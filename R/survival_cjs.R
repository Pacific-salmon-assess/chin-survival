### Initial survival models
## June 20, 2020
## Fit CJS models using individual occurrence data initially only focusing 
# on spatially explicit models (i.e. ind. w/ known stock ID belonging to stocks 
# migrating through study area)
# Information on preliminary model structures in example-models/BPA/Ch07
# Excludes immature individuals and severe injuries
# Linkages between condition, hook location and survival are currently in
# survival logistic


library(tidyverse)
library(rstan)
library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Import occurence matrix (generated in prep_detection_histories.R)
dat_tbl <- readRDS(here::here("data", "generated_data", "det_history_tbl.RDS")) 

# ID tag codes to retain
# 1) no extreme injuries
kept_tags <- dat_tbl %>%
  # unnest and filter out stage 3
  unnest(cols = bio_dat) %>%
  filter(redeploy == "no") %>%
  pull(vemco_code)

dat_tbl_trim <- dat_tbl %>% 
  filter(!agg == "ECVI") %>% 
  mutate(
    bio_dat = purrr::map(bio_dat, function (x) {
      x %>%
        filter(vemco_code %in% kept_tags)
    }),
    wide_array_dat = purrr::map(wide_array_dat, function (x) {
      x %>%
        filter(vemco_code %in% kept_tags)
    })
  ) %>%
  select(agg, bio_dat, wide_array_dat)

# dat_tbl_trim$input_mat <- purrr::map(dat_tbl_trim$wide_array_dat, function(x) {
#   out_mat <- x %>% 
#       select(-vemco_code) %>% 
#       as.matrix()
#   attr(out_mat, "dimnames") <- NULL
#   return(out_mat)
# })


# Prior predictive check -------------------------------------------------------

# simulate distributions for beta_phi and beta_p, representing stage
# specific survival parameters in link space; since 3 beta pars
mu_phi <- rnorm(1000, mean = 0.3, sd = 1.5)
hist(boot::inv.logit(mu_phi),  col=rgb(1,0,0,0.4))
beta_phi <- gamma_phi <- rnorm(1000, mean = 0, sd = 0.5)
beta_p <- rnorm(1000, mean = 0, sd = 1.6)

# transform into survival and detection probabilities and evalute
total_phi <- mu_phi + beta_phi + gamma_phi
hist(boot::inv.logit(total_phi),  col=rgb(1,0,0,0.4))
hist(boot::inv.logit(beta_p), col=rgb(0,0,1,0.4), add = T)


set.seed(123)
n <- 1000
df <- 2
t_samples <- rt(n, df)
sigma <- 0 + 1 * t_samples
# replicates model structure in stan model
surv <- mu_phi + (beta_phi * sigma) + gamma_phi
hist(boot::inv.logit(surv))



# Average survival by segment and year -----------------------------------------

dat_tbl %>% 
  select(agg, long_array_dat) %>% 
  unnest() %>% 
  group_by(agg, year)

mean_det <- purrr::map2(
  dat_tbl_trim$agg,
  dat_tbl_trim$wide_array_dat, 
  ~ .y %>% 
    group_by(year) %>% 
    summarise(
      n = n(),
      across(-c(vemco_code, n), function(x) {sum(x) / n})
    ) %>% 
    ungroup() %>% 
    pivot_longer(
      cols = -c(year, n),
      names_to = "array_number",
      values_to = "ppn_detected"
    ) %>% 
    mutate(
      stock_group = .x
    )
) %>%
  bind_rows()

mean_det_pt <- ggplot(mean_det) +
  geom_point(
    aes(x = array_number, y = ppn_detected, fill = as.factor(year)),
    position = position_dodge(width = 0.5),
    shape = 21
  ) +
  scale_fill_discrete(name = "Year") +
  facet_wrap(~stock_group, scales = "free_x") +
  ggsidekick::theme_sleek() +
  theme(
    legend.position = "top"
  )

png(here::here("figs", "survival", "average_detections.png"), 
    height = 6, width = 9, units = "in", res = 250)
mean_det_pt
dev.off()


# Fit model --------------------------------------------------------------------

## Fixed temporal effects model (survival and probability vary among stages) 
# with each aggregate modeled independently (due to unique migration rates)
# Interaction between year and stage for p because array structure varies, but
# preliminary model comparisons suggest additive (year + stage) has better post
# preds

# Function to convert wide DF to list input for Stan models
prep_cjs_dat <- function(dat, fixp = NULL, grouping_var = NULL) {
  y <- dat %>% 
    dplyr::select(.dots = -c(all_of(grouping_var), "vemco_code"))
  n_occasions <- ncol(y)
  nind <- nrow(y)
  
  out_list <- list(y = y, n_occasions = n_occasions, nind = nind)
  
  # add grouping variables and fixed p values as needed
  if (!is.null(grouping_var)) {
    group <- as.numeric(as.factor(dat[[grouping_var]]))
    g <- length(unique(group))
    out_list <- c(out_list, list(g = g, group = group))
  }
  if (!is.null(fixp) & !is.na(fixp)) {
    out_list <- c(out_list, list(final_fix_p = fixp)) 
  }
  
  return(out_list)
}


# fix terminal value for p for Fraser/Col where detection probability v. high
dat_tbl_trim$fixp <- ifelse(
  grepl("Fraser", dat_tbl_trim$agg) | grepl("Col", dat_tbl_trim$agg), 0.99, NA
  )
dat_tbl_trim$years <- purrr::map(dat_tbl_trim$wide_array_dat, function (x) {
  x$year %>% as.factor() %>% levels()
})
dat_tbl_trim$dat_in <- pmap(list(dat = dat_tbl_trim$wide_array_dat,  
                                 fixp = dat_tbl_trim$fixp),
                            .f = prep_cjs_dat, 
                            grouping_var = "year") 



# Call Stan from R and fit to each aggregate separately 
fixed_mod_sims <- stan_model(
  here::here("R", "stan_models", "cjs_add_fixed_eff.stan")
)
fixed_mod_sims_fixp <- stan_model(
  here::here("R", "stan_models", "cjs_add_fixed_eff_fixp.stan")
)
hier_mod_sims_fixp <- stan_model(
  here::here("R", "stan_models", "cjs_add_hier_eff_fixp.stan")
)



# MCMC settings
n_chains = 4
n_iter = 2000
n_warmup = n_iter / 2
params <- c("beta_phi", "beta_p", #"mu_phi", 
            "gamma_phi", "phi_g", "p_g",
            "sigma_beta_phi",
            "y_hat")

## FOR DEBUGGING
dd <-  dat_tbl_trim$dat_in[[2]]

inits <- lapply(1:n_chains, function (i) {
  list(
    #mu_phi = rnorm(1, 0, 0.5),
    beta_phi = rnorm(dd$g, 0, 0.5),
    gamma_phi = rnorm(dd$n_occasions - 1, 0, 0.5),
    beta_p = matrix(rnorm(dd$g * (dd$n_occasions - 2), 0, 0.5), nrow = dd$g),
    sigma_beta_phi = exp(rnorm(1, 0, 0.5))
  )
})
fit <- sampling(
  hier_mod_sims_fixp, data = dd, pars = params,
  init = inits, chains = n_chains, iter = n_iter, warmup = n_warmup,
  open_progress = FALSE,
  control = list(adapt_delta = 0.95)
)


params_fix <- c("beta_phi", "beta_p", "mu_phi", "gamma_phi", "phi_g", "p_g",
            "y_hat")
inits_fix <- lapply(1:n_chains, function (i) {
  list(
    mu_phi = rnorm(1, 0, 0.5),
    beta_phi = rnorm(dd$g, 0, 0.5),
    gamma_phi = rnorm(dd$n_occasions - 1, 0, 0.5),
    beta_p = matrix(rnorm(dd$g * (dd$n_occasions - 2), 0, 0.5), nrow = dd$g)  )
})
fit2 <- sampling(fixed_mod_sims_fixp, data = dd, pars = params_fix,
                init = inits_fix, chains = n_chains, iter = n_iter, warmup = n_warmup,
                open_progress = FALSE)

fit_list <- list(fit, fit2)

purrr::map(
  fit_list, 
  function(x) {
    dd <- summary(x)$summary
    dd2 <- dd[grepl("gamma_phi", rownames(dd)), ]
    print(dd2)
    sd(dd2[ , "mean"])
  }
)


cjs_fixed_sims <- pmap(
  list(x = dat_tbl_trim$dat_in, fixp = dat_tbl_trim$fixp), 
  .f = function(x, fixp) {
    
    # used fixed p model if fixed p value is provided in tbl and adjust inits
    # for beta p accordingly
    if (!is.na(fixp)) {
      mod <- fixed_mod_sims_fixp
      beta_p_dim <- x$g * (x$n_occasions - 2)
    } else {
      mod <- fixed_mod_sims
      beta_p_dim <- x$g * (x$n_occasions - 1)
    }
    
    # matrix of inits with same dims as estimated parameter matrices
    inits <- lapply(1:n_chains, function (i) {
      list(
        mu_phi = rnorm(1, 0, 0.5),
        beta_phi = rnorm(x$g, 0, 0.5),
        gamma_phi = rnorm(x$n_occasions - 1, 0, 0.5),
        beta_p = matrix(rnorm(beta_p_dim, 0, 0.5), nrow = x$g),
        sigma_beta_phi = exp(rnorm(1, 0, 0.5))
      )
    })
 
    sampling(mod, data = x, pars = params,
             init = inits, chains = n_chains, iter = n_iter, warmup = n_warmup,
             open_progress = FALSE)
})


# add to tibble
dat_tbl_trim$cjs_fixed <- cjs_fixed_sims


## Model checks ----------------------------------------------------------------

# neff 
purrr::map(dat_tbl_trim$cjs_fixed, function (x) {
  tt <- summary(x)
  tt2 <- tt$summary[ , "n_eff"]
  tt2[which(tt2 < 1000)]
})

# rhat
purrr::map(dat_tbl_trim$cjs_fixed, function (x) {
  tt <- bayesplot::rhat(x)
  tt[which(tt > 1.1)]
})


# phi_fit <- cjs_fixed_sims[[2]]
# 
# dum <- summary(phi_fit)$summary 
# phi_ids <- grepl("phi", rownames(dum))
# dum[phi_ids, ]

# chain plots of parameters
pdf(here::here("figs", "survival", "cjs", "diagnostics", 
               "mcmc_trace_fe_add.pdf"),
    height = 5, width = 8)
map(dat_tbl_trim$cjs_fixed, traceplot, pars = c("beta_phi", "beta_p", 
                                                "gamma_phi"))
dev.off()


# posterior predictions check
# compare predicted to observed proportions for each year and stage by aggregate
pp_list <- pmap(list(dat_tbl_trim$dat_in, dat_tbl_trim$cjs_fixed, 
                     dat_tbl_trim$years), 
     function (dat_in, mod, years) {
       obs_dat <- dat_in$y
       n_iters_out <- 500 #dim(post_dat)[1]
       post_dat <- rstan::extract(mod)$y_hat[1:n_iters_out, , ]
       group_ids <- dat_in$group
       n_groups <- length(unique(group_ids))
       n_stages <- dim(post_dat)[3]
       
       # calculate ppns for each the observed data and each iteration/stage of 
       # posterior
       post_ppns_list <- vector(mode = "list", length = n_groups * n_iters_out)
       list_counter <- 1
       for (g in 1:n_groups) {
         g_group <- which(group_ids == g)
         #calc obs ppns
         obs_ppns <- apply(obs_dat[g_group, ], 2, 
                           function(x) sum(x) / nrow(obs_dat[g_group, ]))
         
         for (i in 1:n_iters_out) {
           dum <- post_dat[i, g_group, ] %>% 
             # coerce to matrix to ensure that it's functional even with single
             # obs
             matrix(ncol = n_stages)
           post_ppns_vec <- apply(dum, 2, function(x) sum(x) / nrow(dum))
           post_ppns_list[[list_counter]] <- data.frame(
             array_num = seq(1, n_stages, 1),
             year = years[g],
             post_ppns = post_ppns_vec,
             obs_ppns = obs_ppns
           )
           list_counter <- list_counter + 1
         }
       }
       post_ppns <- post_ppns_list %>% 
         bind_rows() %>% 
         mutate(array_num = as.factor(array_num))
       obs_only <- post_ppns %>% dplyr::select(-post_ppns) %>% distinct()
       
       list(pos = post_ppns, obs = obs_only)
     })

# box plots showing posterior distribution of proportions relative to observed 
pp_plot_list <- map2(pp_list, dat_tbl_trim$agg, function (xx, title) {
  ggplot(xx$pos) +
    geom_boxplot(aes(x = array_num, y = post_ppns)) +
    geom_point(data = xx$obs, aes(x = array_num, y = obs_ppns), 
               color = "red") +
    facet_wrap(~year) +
    labs(title = title) +
    ggsidekick::theme_sleek()
})

pdf(here::here("figs", "survival", "cjs", "diagnostics", 
               "posterior_checks_fe_add.pdf"),
    height = 5, width = 8)
pp_plot_list
dev.off()



## Post-hoc calculations -------------------------------------------------------

## Import survival segment key for labelling plots 
seg_key <- read.csv(here::here("data", "generated_data", 
                               "surv_segment_key_2023.csv")) %>%
  mutate(segment = array_num - 1,
         segment_name = str_replace(segment_name, " ", "\n")) %>% 
  dplyr::select(agg, segment, segment_name) %>% 
  distinct()


# function to adjust final phi estimates to account for imperfect detection
source(here::here("R", "functions", "adjust_phi.R"))

# adjust phi for all stocks except Columbia and Fraser (100% detection rate)
phi_mat <- pmap(
  list(dat_tbl_trim$cjs_fixed, dat_tbl_trim$agg, dat_tbl_trim$years), 
  function(x, y, z) {
  if (grepl("Fraser", y) | grepl("Col", y)) {
    extract(x)[["phi_g"]]
  } else {
    adj_phi(x, group_names = z)
  }
})


# calculate cumulative survival across segments
cum_surv_list <- pmap(
  list(phi_mat, dat_tbl_trim$agg, dat_tbl_trim$years), 
  function (x, yy, group_names) {
    if(dim(x)[2] != length(group_names)) 
      stop("Array dimensions do not match group levels.")
    
    dims <- list(iter = seq(1, dim(x)[1], by = 1),
                 segment = seq(1, dim(x)[3], by = 1))
    
    dumm <- expand.grid(iter = dims$iter,
                        segment = 0, 
                        Freq = 1, 
                        group = group_names)
    
    # calculate the cumulative product across segments for each group and 
    # iteration, then convert to a dataframe
    cumprod_list <- vector(length(group_names), mode = "list") 
    for (i in seq_along(group_names)) {
      cumprod_mat <- t(apply(x[ , i, ], 1, cumprod))
      dimnames(cumprod_mat) = dims[1:2]
      cumprod_list[[i]] <- cumprod_mat %>% 
        as.table() %>% 
        as.data.frame() %>% 
        mutate(group = group_names[i])
    }
    
    cumprod_list %>% 
      bind_rows() %>% 
      rbind(dumm, .) %>%
      mutate(segment = as.integer(as.character(segment)),
             agg = yy) %>%
      rename(est = Freq) %>% 
      left_join(., seg_key, by = c("agg", "segment")) %>% #add segment key
      mutate(segment_name = fct_reorder(as.factor(segment_name), segment),
             hab = case_when(
               segment == max(segment) ~ "river",
               TRUE ~ "marine"
             )) %>% 
      group_by(segment, group) %>% 
      mutate(median = median(est),
             low = quantile(est, 0.05),
             up = quantile(est, 0.95)) %>% 
      ungroup() %>% 
      arrange(segment, group, iter) 
  }
)
dat_tbl_trim$cum_survival <- cum_surv_list

# export
saveRDS(dat_tbl_trim,
        here::here("data", "model_outputs", "fixed_cjs_posterior_tbl.RDS"))


## Visualize posterior ---------------------------------------------------------

source(here::here("R", "functions", "plot_survival.R"))

# import data
dat_tbl_trim <- readRDS(
  here::here("data", "model_outputs", "fixed_cjs_posterior_tbl.RDS")
)

surv_plot_trials <- purrr::map(dat_tbl_trim$cum_survival, function (x) {
  x %>% 
    mutate(agg_name_f = as.factor(agg)) %>% 
    plot_surv(., show_mcmc = T) +
    facet_wrap(~group)
})
surv_plot_clean  <- purrr::map(dat_tbl_trim$cum_survival, function (x) {
  x %>% 
    mutate(agg_name_f = as.factor(agg)) %>%
    plot_surv(., show_mcmc = F) +
    facet_wrap(~group)
})

pdf(here::here("figs", "survival", "cjs", "cum_surv_ind_fixed_trials.pdf"), 
    height = 6, width = 7.5)
surv_plot_trials
dev.off()

pdf(here::here("figs", "survival", "cjs", "cum_surv_ind_fixed_clean.pdf"), 
    height = 6, width = 7.5)
surv_plot_clean
dev.off()


## visualize just Fraser and Upper Columbia river stocks
surv_plot_clean_trim <- dat_tbl_trim %>%
  filter(agg %in% c("Fraser", "Col")) %>%
  pull(cum_survival) %>%
  purrr::map(., function (x) {
  plot_surv(x, show_mcmc = FALSE) +
    facet_wrap(~group) +
    ylim(0.25, 1)
})

surv_fig_trim <- ggpubr::ggarrange(surv_plot_clean_trim[[1]], 
                                   surv_plot_clean_trim[[2]], 
                                   common.legend = TRUE, 
                                   ncol = 1)
cum_surv_trim <- ggpubr::annotate_figure(surv_fig_trim, 
                                         bottom = "Migration Segment",
                                     left = "Cumulative Survival")


png(here::here("figs", "survival", "cjs", "cum_surv_ind_fixed_clean_trim_new.png"), 
    height = 6, width = 9, units = "in", res = 250)
cum_surv_trim
dev.off()


## pared down model for presentation
dat_list <- dat_tbl_trim %>%
  filter(agg %in% c("Fraser", "Col")) %>%
  pull(cum_survival)

plot_list <- purrr::map(
  dat_list, 
  function (x) {
    set.seed(123)
    trials <- sample(1:length(unique(x$iter)), 50, replace = F)
    dat <- x %>% 
      filter(iter %in% trials)
    mu_dat <- x %>% 
      dplyr::select(-iter, -est) %>% 
      distinct() %>% 
      filter(group == "2019")
    
    col_pt <- ifelse(mu_dat$agg == "Col", "#a6cee3", "#1f78b4")
    
   ggplot() +
      geom_pointrange(data = mu_dat, 
                      aes(x = fct_reorder(segment_name, segment), 
                          y = median, ymin = low, ymax = up),
                      fill = col_pt, shape = 21) +
      ggsidekick::theme_sleek() + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank(),
            legend.text=element_text(size = 9),
            legend.title = element_blank()) +
      guides(fill = guide_legend(override.aes = list(size = 1.1))) +
      lims(y = c(0, 1))
  }
)

surv_fig_trim2 <- ggpubr::ggarrange(plot_list[[1]], 
                                    plot_list[[2]], 
                                   common.legend = TRUE, 
                                   ncol = 2)
cum_surv_trim2 <- ggpubr::annotate_figure(surv_fig_trim2, 
                                         bottom = "Migration Segment",
                                         left = "Cumulative Survival")

png(here::here("figs", "survival", "cjs", "cum_surv_ind_fixed_clean_19only.png"), 
    height = 4, width = 9, units = "in", res = 250)
cum_surv_trim2
dev.off()



#-------------------------------------------------------------------------------
## Model with hierarchical, stage-specific effects 

# inputs (include ECVI since stock is treated as RE)
det_dat <- dat_tbl %>%
  select(agg, agg_year, agg_det) %>% 
  unnest(cols = agg_det) %>% 
  mutate(mark = 1) 
det_mat <- det_dat %>% 
  select(mark, int_det, final_det, river_det) %>% 
  as.matrix()
attr(det_mat, "dimnames") <- NULL

colSums(det_mat) / nrow(det_mat) 

#format input list for stan model - exclude individual covariates for now
# fl_z <- as.numeric(scale(det_dat$fl))
stk_id <- as.numeric(as.factor(det_dat$agg_year))
nstks <- length(unique(stk_id))
nstages <- ncol(det_mat)

dat_in <- list(y = det_mat, 
               nind = nrow(det_mat), 
               n_occasions = nstages, 
               group = stk_id, 
               g = nstks)

# mcmc settings
ni <- 1250
nt <- 2
nb <- 400
nc <- 4

# Input data same as fixed effects model
params <- c("mu_phi", "mu_p", "beta_phi", "beta_p", "gamma_phi", "gamma_pp",
            "sigma_phi", "sigma_p", 
            "phi_g", "p_g", "phi_t", "p_t")
inits <- lapply(1:nc, function(i) {
  list(beta_phi = matrix(rnorm(nstks * (nstages - 1)), nrow = nstks),
       beta_p = matrix(rnorm(nstks * (nstages - 1)), nrow = nstks))
})

# fit model (non-centered parameterization necessary to stop div. iters)
hier_mod_nc <- stan_model(here::here("R", "stan_models", 
                                    "cjs_hier_eff_norm_nc.stan"))

fit_hier_nc <- sampling(hier_mod_nc, data = dat_in, init = inits,
                       chains = nc, iter = ni, warmup = nb, thin = 2,
                       control = list(adapt_delta = 0.9,
                                      max_treedepth = 15), 
                       pars = params, seed = 1)

print(fit_hier_nc, digits = 3)

saveRDS(fit_hier_nc, here::here("R", "survival_ex", "hier_nc_cjs_est.RDS"))

shinystan::launch_shinystan(fit_hier_nc)

# code to check divergent iterations (if necessary)
# library(bayesplot)
# posterior_hier <- as.array(fit_hier_nc)
# np_hier <- nuts_params(fit_hier_nc)
# mcmc_parcoord(posterior_hier, np = np_hier,
#               pars = c("sigma_phi[1]", "sigma_phi[2]", "sigma_phi[3]",
#                        "sigma_p[1]", "sigma_p[2]", "sigma_p[3]"))
# mcmc_pairs(posterior_hier, np = np_hier,
#            pars = c("sigma_p[1]", "sigma_p[2]", "sigma_p[3]",
#                     "beta_p[1,1]", "beta_p[1,2]", "beta_p[1,3]"))
# 
# mcmc_scatter(
#   posterior_hier,
#   pars = c("beta_p[1,2]", "sigma_p[2]"),
#   np = np_hier,
#   size = 1
# )

# fit_hier_nc <- readRDS(here::here("R", "survival_ex", "hier_nc_cjs_est.RDS"))


source(here::here("R", "functions", "plot_survival.R"))


# Calculate time varying, but non-stock-specific survival
adj_phi_mat <- adj_phi(mod_obj = fit_hier_nc, phi = "phi_t", p = "p_t")

dims <- list(iter = seq(1, nrow(adj_phi_mat), by = 1),
             segment = seq(1, ncol(adj_phi_mat), by = 1))

#dummy data frame representing release
dumm <- data.frame(iter = dims[1],
                   segment = 0, 
                   est = 1)

avg_surv <- matrix(t(apply(adj_phi_mat, 1, function(x) cumprod(x))),
                   ncol = ncol(adj_phi_mat),
                   dimnames = dims) %>%
  as.table() %>% 
  as.data.frame() %>%
  rename(est = Freq) %>%
  #add dummy row for release
  rbind(dumm, .) %>% 
  mutate(segment = as.integer(as.character(segment))) %>% 
  group_by(segment) %>% 
  mutate(median = median(est),
         low = quantile(est, 0.05),
         up = quantile(est, 0.95)) %>% 
  ungroup() %>% 
  mutate(segment_name = fct_recode(as.factor(segment), 
                                   "Release" = "0", "Intermediate" = "1", 
                                   "Final Marine" = "2", "River" = "3"))

cum_surv_avg <- plot_surv(avg_surv, hier = TRUE)


# Calculate time varying, stock-specific survival via REs
stk_name <- unique(det_dat$agg_year)
cum_dat_hier <- adj_phi(fit_hier_nc, phi = "phi_g", p = "p_g", 
                        stk_names = stk_name) %>% 
  purrr::map(., calc_cum) %>% 
  bind_rows() %>% 
  mutate(segment = fct_recode(as.factor(segment), 
                              "Int" = "1", "Term" = "2", "River" = "3"))

cum_surv_hier <- ggplot(data = cum_dat_hier) +
  geom_pointrange(aes(x = segment, y = mean, ymin = low, ymax = up)) +
  ggsidekick::theme_sleek() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  facet_wrap(~agg)

png(here::here("figs", "survival", "cum_surv_hier_avg.png"), height = 5.5, 
    width = 7, units = "in", res = 300)
cum_surv_avg
dev.off()

png(here::here("figs", "survival", "cum_surv_hier_stock-specific.png"), 
    height = 7, width = 7, units = "in", res = 300)
cum_surv_hier
dev.off()


#-------------------------------------------------------------------------------

## Fixed temporal effects model (survival and probability vary among stages) 
# with aggregates modeled simultaneously and covariates included. 
# Note to ensure unbiased estimates of survival, det probability must also vary
# among stocks because of different migration routes.
# Stock effects currently do not vary among detection stages.
# Generally favor first approach or hierarchical model

det_mat <- det_dat %>% 
  select(mark, int_det, final_det, river_det) %>% 
  as.matrix()
attr(det_mat, "dimnames") <- NULL

colSums(det_mat) / nrow(det_mat) 

#format input list for stan model - exclude individual covariates for now
# fl_z <- as.numeric(scale(det_dat$fl))
stk_id <- as.numeric(as.factor(det_dat$Agg_Name))
nstks <- length(unique(stk_id))
nstages <- ncol(det_mat)

dat_in <- list(y = det_mat, nind = nrow(det_mat), n_occasions = nstages, 
               group = stk_id, g = nstks)

# mcmc settings
ni <- 1250
nt <- 2
nb <- 400
nc <- 4

# parameters to monitor
params <- c("beta_phi", "beta_p", "phi_g", "p_g")

# fit model
fixed_mod2 <- stan_model(here::here("R", "stan_models", 
                                   "cjs_fixed_eff.stan"))
fit_fix <- sampling(fixed_mod2, data = dat_in,
                    chains = nc, iter = ni, warmup = nb, thin = nt,
                    control = list(adapt_delta = 0.85), pars = params)
saveRDS(fit_fix, here::here("R", "survival_ex", "fixed_cjs_est.RDS"))

shinystan::launch_shinystan(fit_fix)


# Calculate cumulative survival for each stock
# adjust phi
stk_name <- unique(det_dat$Agg_Name)


fit_fix <- readRDS(here::here("R", "survival_ex", "fixed_cjs_est.RDS"))

cum_dat_fix <- adj_phi(fit_fix, phi = "phi_g", p = "p_g", 
                       stk_names = stk_name) %>% 
  map(., calc_cum) %>% 
  bind_rows() %>% 
  mutate(segment = fct_recode(as.factor(segment), 
                              "Int" = "1", "Term" = "2", "River" = "3"))

cum_surv_fix <- ggplot(data = cum_dat_fix) +
  geom_pointrange(aes(x = segment, y = mean, ymin = low, ymax = up)) +
  ggsidekick::theme_sleek() + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank()) +
  facet_wrap(~agg)

png(here::here("figs", "survival", "cum_surv_fixed.png"), height = 5.5, 
    width = 7, units = "in", res = 300)
cum_surv_fix
dev.off()


#-------------------------------------------------------------------------------


## Look at linkages between hook location/injuries and whether a fish was 
# detected at least once and whether it succsefully migrated through study area

det_dat <- readRDS(here::here("data", "generatedData", "surv_df.RDS"))
acoustic <- readRDS(here::here("data", "taggingData", 
                               "acousticOnly_GSI.RDS")) %>% 
  mutate(transmitter_id = paste("A69", "9006", acoustic, sep = "-")) %>% 
  select(transmitter_id, Agg_Name, hookLoc, injury, scaleLoss, finDam, meanLogE,
         fl) %>% 
  left_join(., det_dat, by = "transmitter_id") %>% 
  filter(!is.na(det))

det_prop = function(group = "injury") {
  acoustic %>%
    group_by_(group) %>% 
    summarize(total = length(unique(transmitter_id)),
              prop_det = sum(det) / total,
              prop_det_final = sum(final_det) / total)
}
det_prop(group = "injury")
det_prop(group = "scaleLoss")
det_prop(group = "finDam")
