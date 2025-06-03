### CJS survival models
## June 20, 2020
## Fit CJS models using individual occurrence data initially only focusing 
# on spatially explicit models (i.e. ind. w/ known stock ID belonging to stocks 
# migrating through study area)
# Information on preliminary model structures in example-models/BPA/Ch07
# Excludes immature individuals and severe injuries


library(tidyverse)
library(rstan)
library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Import occurence matrix (generated in chinTagging/prep_detection_histories.R 
# and cleaned in data_clean.R)
dat_tbl_trim <- readRDS(here::here("data", "surv_cjs_data.rds")) 

mature_tags <- dat_tbl_trim %>% 
  unnest(bio_dat) %>%
  filter(stage_1 == "mature") %>% 
  pull(vemco_code)

dat_tbl_trim$bio_dat <- purrr::map(
  dat_tbl_trim$bio_dat,
  ~ .x %>% 
    filter(vemco_code %in% mature_tags)
)

dat_tbl_trim$wide_array_dat <- purrr::map(
  dat_tbl_trim$wide_array_dat,
  ~ .x %>% 
    filter(vemco_code %in% mature_tags)
)
dat_tbl_trim$stock_group <- fct_relevel(
  as.factor(dat_tbl_trim$stock_group),
  "Cali", "Low Col.", "Up Col.", "Fraser", "South Puget"
)


# Import duration and distance estimates for scaling survival
array_dist <- read.csv(
  here::here("data", "mean_array_locs_measured.csv")
) %>% 
  select(
    stock_group, array_num, dist
  )


## Import survival segment key for labelling plots 
seg_key <- read.csv(
  here::here("data", "surv_segment_key_2023.csv")
) %>%
  mutate(segment = array_num - 1,
         segment_name = ifelse(
           stock_group == "Up Col." & segment == 6,
           "Above\nBonneville",
           str_replace(segment_name, " ", "\n")
         ),
         array_num_key = paste(segment, segment + 1, sep = "_")) %>% 
  dplyr::select(stock_group, segment, segment_name, array_num, array_num_key,
                max_array_num, terminal) %>% 
  distinct() %>% 
  left_join(
    ., array_dist, by = c("stock_group", "array_num")
  ) 


# Average survival by segment and year -----------------------------------------

mean_det <- purrr::map2(
  dat_tbl_trim$stock_group,
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
      names_to = "array_num",
      values_to = "ppn_detected"
    ) %>% 
    mutate(
      stock_group = .x,
      array_num = as.numeric(array_num)
    )
) %>%
  bind_rows() %>% 
  left_join(
    .,
    seg_key %>% 
      select(stock_group, array_num, segment_name) %>% 
      distinct(),
    by = c("array_num", "stock_group")
  ) 

mean_det_pt <- ggplot(mean_det) +
  geom_point(
    aes(x = fct_reorder(segment_name, array_num), y = ppn_detected, 
        fill = as.factor(year)),
    position = position_dodge(width = 0.5),
    shape = 21
  ) +
  scale_fill_discrete(name = "Year") +
  facet_wrap(~stock_group, scales = "free_x") +
  ggsidekick::theme_sleek() +
  labs(y = "Proportion Tags Detected") +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  )

png(here::here("figs", "average_detections.png"), 
    height = 6, width = 9, units = "in", res = 250)
mean_det_pt
dev.off()

#export for Rmd
saveRDS(mean_det_pt, here::here("figs", "average_detections.rds" ))



# check fraser
left_join(dat_tbl_trim$bio_dat[[2]] %>% 
            select(vemco_code, agg),
          dat_tbl_trim$wide_array_dat[[2]]) %>% 
  group_by(year, agg) %>%
  summarise(
    n = n(),
    across(-c(vemco_code, n), function(x) {sum(x) / n})
  ) %>% 
  pivot_longer(
    cols = -c(year, agg, n),
    names_to = "array_num",
    values_to = "ppn_detected"
  )  %>% 
  ggplot(.) +
  geom_point(
    aes(x = array_num, y = ppn_detected, 
        fill = as.factor(year)),
    position = position_dodge(width = 0.5),
    shape = 21
  ) +
  scale_fill_discrete(name = "Year") +
  facet_wrap(~stock_group, scales = "free_x") +
  ggsidekick::theme_sleek() +
  labs(y = "Proportion Tags Detected") +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  ) +
  facet_wrap(~agg)


# Fit model --------------------------------------------------------------------

## Fixed temporal effects model (survival and probability vary among stages) 
# with each aggregate modeled independently (due to unique migration rates)
# Interaction between year and stage for p because array structure varies, but
# preliminary model comparisons suggest additive (year + stage) has better post
# preds

# Function to convert wide DF to list input for Stan models
prep_cjs_dat <- function(dat, grouping_vars = NULL) {
  y <- dat %>% 
    dplyr::select(matches("[1-9]"))
  n_occasions <- ncol(y)
  nind <- nrow(y)
  
  out_list <- list(y = y, n_occasions = n_occasions, nind = nind)
  
  # add grouping variables and fixed p values as needed
  if (!is.null(grouping_vars)) {
    for(i in seq_along(grouping_vars)) {
      # vector of factor IDs
      out_list[[grouping_vars[i]]] <- as.numeric(
        as.factor(dat[[grouping_vars[i]]])
        )
      # number of factor IDs
      out_list[[paste0("n", grouping_vars[i])]] <- 
        length(unique(dat[[grouping_vars[i]]]))
    }
  }
  
  return(out_list)
}


dat_tbl_trim$years <- purrr::map(dat_tbl_trim$wide_array_dat, function (x) {
  x$year %>% as.factor() %>% levels()
})
dat_tbl_trim$grouping_vars <- purrr::map(
  dat_tbl_trim$wide_array_dat, 
  function(x) {
    dd <- colnames(x)
    # remove columns that include vemco code or array numbers
    dd[!(grepl("^([0-9]+)$", dd) | dd == "vemco_code")]
})
dat_tbl_trim$dat_in <- pmap(
  list(dat = dat_tbl_trim$wide_array_dat,  
       grouping_vars = dat_tbl_trim$grouping_vars),
  .f = prep_cjs_dat
) 



# Call Stan from R and fit to each aggregate separately 
hier_mod_sims <- stan_model(
  here::here("R", "stan_models", "cjs_add_hier_eff_adaptive_v3.stan")
)
hier_mod_sims_stk <- stan_model(
  here::here("R", "stan_models", "cjs_add_hier_eff_adaptive_v3_stock.stan")
)


## TEST FIT --------------------------------------------------------------------

#MCMC settings
# n_chains = 1
# n_iter = 2000
# n_warmup = n_iter / 2
# params <- c(
#   "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
#   "L_Rho_yr", "alpha_p", "alpha_yr_p",
#   # transformed pars or estimated quantities
#   "Rho_yr", "phi_yr", "p_yr", "beta_yr", "y_hat"
# )
# params2 <- c(
#   "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
#   "L_Rho_yr", "alpha_p", "alpha_yr_p",
#   # transformed pars or estimated quantities
#   "Rho_yr", "phi_yr", "p_yr", "y_hat"
# )
# params3 <- c(
#   "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
#   "L_Rho_yr", "alpha_p", "alpha_yr_p",
#   # transformed pars or estimated quantities
#   "Rho_yr", "phi_yr", "p_yr", "beta_yr", "y_hat",
#   # stock specific pars and quants
#   "alpha_stk_phi", "sigma_alpha_stk_phi", "phi_stk",
#   "alpha_yr_phi", "alpha_stk_phi_z"
# )

## FOR DEBUGGING
# dd <-  dat_tbl_trim$dat_in[[2]]
# dd$nstock <- dat_tbl_trim$bio_dat[[2]]$agg %>% unique() %>% length()
# dd$stock <- dat_tbl_trim$bio_dat[[2]]$agg %>%
#   as.factor() %>%
#   droplevels() %>%
#   as.numeric()
# # saveRDS(dd, here::here("data", "model_outputs", "sample_cjs_dat.rds"))
# 
# inits <- lapply(1:n_chains, function (i) {
#   list(
#     alpha_phi = rnorm(1, 0, 0.5),
#     # note z transformed so inverted compared to beta_phi or beta_p
#     alpha_yr_phi_z = matrix(
#       rnorm(dd$nyear * (dd$n_occasions - 1), 0, 0.5),
#       nrow = (dd$n_occasions - 1)
#     ),
#     alpha_t_phi = rnorm(dd$n_occasions - 1, 0, 0.5),
#     sigma_alpha_yr_phi = rexp((dd$n_occasions - 1), 2),
#     # replace with rlkjcorr(XXX, K = 2, eta = 2) from rethinking package
#     L_Rho_yr = matrix(
#       runif((dd$n_occasions - 1)^2, -0.5, 0.5),
#       nrow = (dd$n_occasions - 1)
#     ),
#     alpha_p = rnorm(1, 0, 0.5),
#     alpha_yr_p = matrix(rnorm(dd$nyear * (dd$n_occasions), 0, 0.5), nrow = dd$nyear)
#   )
# })
# inits2 <- lapply(1:n_chains, function (i) {
#   list(
#     alpha_phi = rnorm(1, 0, 0.5),
#     # note z transformed so inverted compared to beta_phi or beta_p
#     alpha_yr_phi_z = matrix(
#       rnorm(dd$nyear * (dd$n_occasions - 1), 0, 0.5), nrow = (dd$n_occasions - 1)
#     ),
#     alpha_t_phi = rnorm(dd$n_occasions - 1, 0, 0.5),
#     sigma_alpha_yr_phi = rexp((dd$n_occasions - 1), 2),
#     # replace with rlkjcorr(XXX, K = 2, eta = 2) from rethinking package
#     L_Rho_yr = matrix(
#       runif((dd$n_occasions - 1)^2, -0.5, 0.5), nrow = (dd$n_occasions - 1)
#     ),
#     alpha_p = rnorm(1, 0, 0.5),
#     alpha_yr_p = matrix(
#       rnorm(dd$nyear * (dd$n_occasions - 1), 0, 0.5), nrow = dd$nyear
#     )
#   )
# })
# inits3 <- lapply(1:n_chains, function (i) {
#   list(
#     alpha_phi = rnorm(1, 0, 0.5),
#     # note z transformed so inverted compared to beta_phi or beta_p
#     alpha_yr_phi_z = matrix(
#       rnorm(dd$nyear * (dd$n_occasions - 1), 0, 0.5), 
#       nrow = (dd$n_occasions - 1)
#     ),
#     alpha_stk_phi = rnorm(dd$nstock, 0, 0.5),
#     alpha_t_phi = rnorm(dd$n_occasions - 1, 0, 0.5),
#     sigma_alpha_yr_phi = rexp((dd$n_occasions - 1), 2),
#     sigma_alpha_stk_phi = rexp(1, 2),
#     # replace with rlkjcorr(XXX, K = 2, eta = 2) from rethinking package
#     L_Rho_yr = matrix(
#       runif((dd$n_occasions - 1)^2, -0.5, 0.5), 
#       nrow = (dd$n_occasions - 1)
#     ),
#     alpha_p = rnorm(1, 0, 0.5),
#     alpha_yr_p = matrix(
#       rnorm(dd$nyear * (dd$n_occasions), 0, 0.5), nrow = dd$nyear
#     )
#   )
# })


# fit <- sampling(
#   hier_mod_sims, data = dd, pars = params,
#   init = inits, chains = n_chains, iter = n_iter, warmup = n_warmup,
#   open_progress = FALSE,
#   control = list(adapt_delta = 0.95)
# )
# fit2 <- sampling(
#   hier_mod_sims_fixp, data = dd, pars = params2,
#   init = inits2, chains = n_chains, iter = n_iter, warmup = n_warmup,
#   open_progress = FALSE,
#   control = list(adapt_delta = 0.95)
# )
# fit3 <- sampling(
#   hier_mod_sims_stk, data = dd, pars = params3,
#   init = inits3, chains = n_chains, iter = n_iter, warmup = n_warmup,
#   open_progress = FALSE,
#   control = list(adapt_delta = 0.95)
# )
# 
# 
# phi_pattern <- "phi_yr"#"phi_yr\\[\\d+,5\\]$"
# p_pattern <- "p_yr"#\\[\\d+,6\\]$"
# 
# fit_post <- summary(fit)$summary
# fit_post[grepl("beta", rownames(fit_post)), ]
# fit_post[grepl(p_pattern, rownames(fit_post)), ]
# fit_post[grepl(phi_pattern, rownames(fit_post)), ]
# fit_post[grepl("sigma_alpha_yr_phi", rownames(fit_post)), ]
# fit_post[grepl("Rho_yr", rownames(fit_post)), ]
# 
# fit_post2 <- summary(fit2)$summary
# fit_post2[grepl(p_pattern, rownames(fit_post2)), ]
# fit_post2[grepl("alpha_yr_phi", rownames(fit_post2)), ]
# 
# fit_post3 <- summary(fit3)$summary
# fit_post3[grepl("alpha_stk_phi", rownames(fit_post3)), ]
# fit_post3[grepl("alpha_yr_phi", rownames(fit_post3)), ]



## REAL FIT --------------------------------------------------------------------

# MCMC settings
n_chains = 4
n_iter = 2000
n_warmup = n_iter / 2
pars_in <- c(
  "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
  "L_Rho_yr", "alpha_p", "alpha_yr_p",
  # transformed pars or estimated quantities
  "Rho_yr", "alpha_yr_phi", "phi_yr", "p_yr", "beta_yr", "y_hat"
)

cjs_hier_sims <- pmap(
  list(x = dat_tbl_trim$dat_in, stock_group = dat_tbl_trim$stock_group,
       bio_dat = dat_tbl_trim$bio_dat),
  .f = function(x, stock_group, bio_dat) {

    if (stock_group == "Fraser") {
      mod <- hier_mod_sims_stk
      x$nstock <- bio_dat$agg %>% unique() %>% length()
      x$stock <- bio_dat$agg %>%
        as.factor() %>%
        droplevels() %>%
        as.numeric()

      pars <- c(pars_in,
                # stock specific pars and quants
                "alpha_stk_phi", "sigma_alpha_stk_phi", "phi_stk",
                "alpha_stk_phi_z")

      inits <- lapply(1:n_chains, function (i) {
        list(
          alpha_phi = rnorm(1, 0, 0.5),
          # note z transformed so inverted compared to beta_phi or beta_p
          alpha_yr_phi_z = matrix(
            rnorm(x$nyear * (x$n_occasions - 1), 0, 0.5),
            nrow = (x$n_occasions - 1)
          ),
          alpha_stk_phi_z = rnorm(x$nstock, 0, 0.5),
          alpha_t_phi = rnorm(x$n_occasions - 1, 0, 0.5),
          sigma_alpha_yr_phi = rexp((x$n_occasions - 1), 2),
          sigma_alpha_stk_phi = rexp(1, 2),
          # replace with rlkjcorr(XXX, K = 2, eta = 2) from rethinking package
          L_Rho_yr = matrix(
            runif((x$n_occasions - 1)^2, -0.5, 0.5), nrow = (x$n_occasions - 1)
          ),
          alpha_p = rnorm(1, 0, 0.5),
          alpha_yr_p = matrix(
            rnorm(x$nyear * (x$n_occasions), 0, 0.5), nrow = x$nyear
          )
        )
      })

    } else{
      mod <- hier_mod_sims

      pars <- pars_in

      # matrix of inits with same dims as estimated parameter matrices
      inits <- lapply(1:n_chains, function (i) {
        list(
          alpha_phi = rnorm(1, 0, 0.5),
          # note z transformed so inverted compared to beta_phi or beta_p
          alpha_yr_phi_z = matrix(
            rnorm(x$nyear * (x$n_occasions - 1), 0, 0.5),
            nrow = (x$n_occasions - 1)
          ),
          alpha_t_phi = rnorm(x$n_occasions - 1, 0, 0.5),
          sigma_alpha_yr_phi = rexp((x$n_occasions - 1), 2),
          # replace with rlkjcorr(XXX, K = 2, eta = 2) from rethinking package
          L_Rho_yr = matrix(
            runif((x$n_occasions - 1)^2, -0.5, 0.5),
            nrow = (x$n_occasions - 1)
          ),
          alpha_p = rnorm(1, 0, 0.5),
          alpha_yr_p = matrix(rnorm(x$nyear * (x$n_occasions), 0, 0.5),
                              nrow = x$nyear)
        )
      })
    }
    
    sampling(
      mod, data = x, pars = pars,
      init = inits, chains = n_chains, iter = n_iter, warmup = n_warmup,
      open_progress = FALSE,
      control = list(adapt_delta = 0.96)
    )
})


# EXPLORE DIVERGENT TRANSITIONS
# sampler_params <- get_sampler_params(dd)
# divergent_indices <- which(sampler_params[[1]][, "divergent__"] == 1)
# posterior_samples <- as.array(dd)
# # Extract samples corresponding to divergent transitions
# divergent_samples <- posterior_samples[divergent_indices, , ]
#
# parameter_samples <- posterior_samples[, , 38]
# divergence_region <- (452 - 5):(452 + 5)
#
# df <- data.frame(
#   Iteration = seq_along(parameter_samples),
#   Parameter = parameter_samples
# ) %>%
#   pivot_longer(
#     starts_with("Parameter"),
#     names_to = "Chain",
#     values_to = "Est"
#   )
#
# ggplot(df[df$Iteration %in% divergence_region, ],
#        aes(x = Iteration, y = Est, Colour = Chain)) +
#   geom_line() +
#   geom_point(data = df %>% filter(Iteration == "452"), color = "red")
#

saveRDS(cjs_hier_sims,
        here::here("data", "model_outputs", "hier_cjs_fit_tbl.RDS"))

cjs_hier_sims <- readRDS(
  here::here("data", "model_outputs", "hier_cjs_fit_tbl.RDS")
  ) 

# add to tibble
dat_tbl_trim$cjs_hier <- cjs_hier_sims



## Model checks ----------------------------------------------------------------

# neff 
purrr::map(dat_tbl_trim$cjs_hier, function (x) {
  tt <- summary(x)
  tt2 <- tt$summary[ , "n_eff"]
  tt2[which(tt2 < 1000)]
})

# rhat
purrr::map(dat_tbl_trim$cjs_hier, function (x) {
  tt <- bayesplot::rhat(x)
  tt[which(tt > 1.05)]
})


# chain plots of parameters
par_list <- list(
  "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
  "L_Rho_yr", "alpha_p", "alpha_yr_p"
)
file_names <- paste(
  c( "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
     "L_Rho_yr", "alpha_p", "alpha_yr_p"),
  "hier_mcmc_trace.pdf",
  sep = "_"
)

for (i in seq_along(par_list)) {
  pdf(here::here("figs", "diagnostics", 
                 file_names[i]),
      height = 5, width = 8)
  for (j in seq_along(dat_tbl_trim$cjs_hier)) {
    print(traceplot(dat_tbl_trim$cjs_hier[[j]], pars = par_list[i]))
  }
  dev.off()  
}


# posterior predictions check
# compare predicted to observed proportions for each year and stage by aggregate
# TODO: modify hier stock model to include alpha_stk in posterior preds
pp_list <- pmap(list(dat_tbl_trim$dat_in, dat_tbl_trim$cjs_hier, 
                     dat_tbl_trim$years), 
     function (dat_in, mod, years) {
       obs_dat <- dat_in$y
       n_iters_out <- 500 #dim(post_dat)[1]
       post_dat <- rstan::extract(mod)$y_hat[1:n_iters_out, , ]
       year_ids <- dat_in$year
       n_years <- length(unique(year_ids))
       n_stages <- dim(post_dat)[3]
       
       # calculate ppns for each the observed data and each iteration/stage of 
       # posterior
       post_ppns_list <- vector(mode = "list", length = n_years * n_iters_out)
       list_counter <- 1
       for (g in 1:n_years) {
         g_year <- which(year_ids == g)
         #calc obs ppns
         obs_ppns <- apply(obs_dat[g_year, ], 2, 
                           function(x) sum(x) / nrow(obs_dat[g_year, ]))
         
         for (i in 1:n_iters_out) {
           dum <- post_dat[i, g_year, ] %>% 
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
pp_plot_list <- map2(pp_list, dat_tbl_trim$stock_group, function (xx, title) {
  ggplot(xx$pos) +
    geom_boxplot(aes(x = array_num, y = post_ppns)) +
    geom_point(data = xx$obs, aes(x = array_num, y = obs_ppns), 
               color = "red") +
    facet_wrap(~year) +
    labs(title = title) +
    ggsidekick::theme_sleek()
})

pdf(here::here("figs", "diagnostics", "posterior_checks_hier.pdf"),
    height = 5, width = 8)
pp_plot_list
dev.off()



## Prior-Posterior Comparisons -------------------------------------------------

# average survival
alpha_phi_prior_df <- data.frame(est = rnorm(4000, 0.8, 1), parameter = "Prior")
purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group, 
  function(x , y) {
    dum1 <- extract(x)[["alpha_phi"]]  
    dum <- data.frame("est" = as.numeric(dum1))
    
    p <- ggplot() +
      geom_density(data = dum, aes(x = est), 
                   fill = "red", colour = "red", alpha = 0.4) +
      geom_density(data = alpha_phi_prior_df, aes(x = est), 
                   fill = "blue", colour = "blue", alpha = 0.4) +
      ggsidekick::theme_sleek() +
      labs(x = "Gamma_Phi Estimate", y = "Kernel Density", title = y) 
    
    file_name <- paste("gamma_phi_", y, ".png", sep = "")
    
    png(here::here("figs", "cjs", "posterior_prior_comp", file_name), 
        height = 4.5, width = 6, units = "in", res = 250)
    print(p)
    dev.off()
  }
)


# average stage specific survival
alpha_phi_t_prior_df <- data.frame(est = rnorm(4000, 0, 0.5), parameter = "Prior")
purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group, 
  function(x , y) {
    dum <- extract(x)[["alpha_t_phi"]] %>% 
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "segment", values_to = "est", 
                   names_prefix = "V")
    
    p <- ggplot() +
      geom_density(data = dum, aes(x = est), 
                   fill = "red", colour = "red", alpha = 0.4) +
      geom_density(data = alpha_phi_t_prior_df, aes(x = est), 
                   fill = "blue", colour = "blue", alpha = 0.4) +
      facet_wrap(~segment) + 
      ggsidekick::theme_sleek() +
      labs(x = "Gamma Phi T Estimate", y = "Kernel Density", title = y) 
    
    file_name <- paste("gamma_phi_t_", y, ".png", sep = "")
    
    png(here::here("figs", "cjs", "posterior_prior_comp", file_name), 
        height = 4.5, width = 6, units = "in", res = 250)
    print(p)
    dev.off()
  }
)

# average detection prob
alpha_p_prior_df <- data.frame(est = rnorm(4000, 0.5, 1.2), parameter = "Prior")
purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group, 
  function(x , y) {
    dum1 <- extract(x)[["alpha_p"]]  
    dum <- data.frame("est" = as.numeric(dum1))
    
    p <- ggplot() +
      geom_density(data = dum, aes(x = est), 
                   fill = "red", colour = "red", alpha = 0.4) +
      geom_density(data = alpha_phi_prior_df, aes(x = est), 
                   fill = "blue", colour = "blue", alpha = 0.4) +
      ggsidekick::theme_sleek() +
      labs(x = "Gamma p Estimate", y = "Kernel Density", title = y) 
    
    file_name <- paste("gamma_p_", y, ".png", sep = "")
    
    png(here::here("figs", "cjs", "posterior_prior_comp", file_name), 
        height = 4.5, width = 6, units = "in", res = 250)
    print(p)
    dev.off()
  }
)


# stage and year specific detection prob
alpha_p_yr_prior_df <- data.frame(est = rnorm(4000, 0, 0.5), parameter = "Prior")
purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group, 
  function(x , y) {
    dum <- extract(x)[["alpha_yr_p"]] %>% 
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "segment", values_to = "est", 
                   names_prefix = "V")
    
    yr_p_mat <- extract(x)[["alpha_yr_p"]] 
    p_sum <- yr_p_mat %>% 
      as.data.frame.table() %>% 
      rename(year = Var2, segment = Var3) %>% 
      mutate(est = Freq,
             year = as.numeric(as.factor(year)) + 2018,
             array_num = as.numeric(as.factor(segment)))
    
    p <- ggplot() +
      geom_density(data = p_sum, aes(x = est), 
                   fill = "red", colour = "red", alpha = 0.4) +
      geom_density(data = alpha_p_yr_prior_df, aes(x = est), 
                   fill = "blue", colour = "blue", alpha = 0.4) +
      facet_wrap(~segment) + 
      ggsidekick::theme_sleek() +
      labs(x = "Gamma p_jt Estimate", y = "Kernel Density", title = y) 
    
    file_name <- paste("gamma_p_jt_", y, ".png", sep = "")
    
    png(here::here("figs", "cjs", "posterior_prior_comp", file_name), 
        height = 4.5, width = 6, units = "in", res = 250)
    print(p)
    dev.off()
  }
)


# among year variability in survival
sigma_yr_prior_df <- data.frame(est = rexp(4000, rate = 2), parameter = "Prior")
purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group, 
  function(x , y) {
    dum <- extract(x)[["sigma_alpha_yr_phi"]] %>% 
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "segment", values_to = "est", 
                   names_prefix = "V")
    
    p <- ggplot() +
      geom_density(data = dum, aes(x = est), 
                   fill = "red", colour = "red", alpha = 0.4) +
      geom_density(data = sigma_yr_prior_df, aes(x = est), 
                   fill = "blue", colour = "blue", alpha = 0.4) +
      facet_wrap(~segment) + 
      ggsidekick::theme_sleek() +
      labs(x = "Sigma Year Estimate", y = "Kernel Density", title = y) 
    
    file_name <- paste("sigma_year_", y, ".png", sep = "")
    
    png(here::here("figs", "cjs", "posterior_prior_comp", file_name), 
        height = 4.5, width = 6, units = "in", res = 250)
    print(p)
    dev.off()
  }
)


## Post-hoc calculations -------------------------------------------------------


# extract phi matrix and swap last col with beta estimates (i.e. combined p and 
# phi) except for fix p models 
# stocks
phi_mat <- purrr::map(
  dat_tbl_trim$cjs_hier, 
  function (x) {
  phi_adj <- extract(x)[["phi_yr"]]
  # replace array corresponding to last stage-specific survival est, w/ beta
  phi_adj[ , , dim(phi_adj)[3]] <- extract(x)[["beta_yr"]]
  return(phi_adj)
})


# calculate cumulative survival across segments
cum_surv_list <- pmap(
  list(phi_mat, dat_tbl_trim$stock_group, dat_tbl_trim$years), 
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
             stock_group = yy) %>%
      rename(est = Freq) %>% 
      #add segment key
      left_join(., seg_key, by = c("stock_group", "segment")) %>% 
      mutate(par = case_when(
               segment == max(segment) ~ "beta",
               TRUE ~ "phi"
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


## calculate average among years
phi_mat_mean <- purrr::map(
  dat_tbl_trim$cjs_hier,
  function(x) {
    alpha <- extract(x)[["alpha_phi"]]
    alpha_t <- extract(x)[["alpha_t_phi"]]
    mean_phi <- apply(
      alpha_t, 2, function (alpha_t) boot::inv.logit(alpha + alpha_t)
      )
    # calculate average beta per iteration among years, then replace mean
    # phi for last stage when p isn't fixed
    mean_beta <- apply(
      extract(x)[["beta_yr"]], 1, mean
    )
    # adjust final stage-specific estimate by detection probability
    mean_phi[, ncol(mean_phi)] <- mean_beta
    return(mean_phi)
  }
)

cum_surv_list_mean <- pmap(
  list(phi_mat_mean, dat_tbl_trim$stock_group), 
  function (x, yy) {
    dims <- list(iter = seq(1, dim(x)[1], by = 1),
                 segment = seq(1, dim(x)[2], by = 1))
    
    dumm <- expand.grid(iter = dims$iter,
                        segment = 0, 
                        Freq = 1)
    
    # calculate the cumulative product across segments for each iteration,
    # then convert to a dataframe
    cumprod_mat <- t(apply(x, 1, cumprod))
    dimnames(cumprod_mat) = dims[1:2]
    cumprod_dat <-  cumprod_mat %>% 
      as.table() %>% 
      as.data.frame()
    
    rbind(dumm, cumprod_dat) %>%
      mutate(segment = as.integer(as.character(segment)),
             stock_group = yy) %>%
      rename(est = Freq) %>% 
      #add segment key
      left_join(., seg_key, by = c("stock_group", "segment")) %>% 
      mutate(segment_name = fct_reorder(as.factor(segment_name), segment),
             par = case_when(
               segment == max(segment) ~ "beta",
               TRUE ~ "phi"
             )) %>% 
      group_by(segment) %>% 
      mutate(median = median(est),
             low = quantile(est, 0.05),
             up = quantile(est, 0.95)) %>% 
      ungroup() %>% 
      arrange(segment, iter) 
  }
)
dat_tbl_trim$phi_mat_mean <- phi_mat_mean
dat_tbl_trim$cum_survival_mean <- cum_surv_list_mean


# export
saveRDS(dat_tbl_trim,
        here::here("data", "model_outputs", "hier_cjs_posterior_tbl.RDS"))


## Parameter estimates ---------------------------------------------------------

dat_tbl_trim <- readRDS(
  here::here("data", "model_outputs", "hier_cjs_posterior_tbl.RDS")
  )


# plots of estimated correlations among stages in year-specific survival 
rho_plot_list <- purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group, 
  function(x , y) {
    rho_mat <- extract(x)[["Rho_yr"]]
    
    rho_list <- vector(mode = "list", length = dim(rho_mat)[3])
    for (i in 1:dim(rho_mat)[3]) {
      dum <- as.data.frame(rho_mat[ , , i])
      colnames(dum) <- paste(
        i, seq(1, ncol(rho_mat[ , , i])), sep = "_"
      )
      rho_list[[i]] <- dum %>% 
        pivot_longer(everything(), names_to = "among_segment_corr",
                     values_to = "est")
    }
    
    rho_list %>% 
      bind_rows() %>% 
      group_by(among_segment_corr) %>% 
      summarize(
        med = median(est),
        low = rethinking::HPDI(est, prob = 0.9)[1],
        up = rethinking::HPDI(est, prob = 0.9)[2]
      ) %>%
      ungroup() %>% 
      ggplot(.) +
      geom_pointrange(aes(x = among_segment_corr,
                          y = med, ymin = low, ymax = up),
                      fill = "white", shape = 21) +
      ggsidekick::theme_sleek() +
      labs(x = "Among Segment Correlation in Beta",
           y = "Rho Estimate", title = y)
  }
)


# plots of estimated among year variability in stage-specific survival 
sigma_plot_list <- purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group, 
  function(x , y) {
    extract(x)[["sigma_alpha_yr_phi"]] %>% 
      as.data.frame() %>%
      pivot_longer(everything(), names_to = "segment", values_to = "est", 
                   names_prefix = "V") %>%
      group_by(segment) %>% 
      summarize(
        med = median(est),
        low = rethinking::HPDI(est, prob = 0.9)[1],
        up = rethinking::HPDI(est, prob = 0.9)[2]
      ) %>% 
      ungroup() %>% 
      mutate(
        segment = as.numeric(segment),
        segment_name = paste("stage", segment, sep = "_")
      ) %>% 
      ggplot(.) +
      geom_pointrange(aes(x = segment_name,
                          y = med, ymin = low, ymax = up),
                      fill = "white", shape = 21) +
      ggsidekick::theme_sleek() +
      labs(x = "Migration Stage", y = "Sigma Estimate", title = y)
  }
)

pdf(here::here("figs", "cjs", "estimated_rho.pdf"), 
    height = 4.5, width = 6)
rho_plot_list
dev.off()

pdf(here::here("figs", "cjs", "estimated_sigma_beta.pdf"), 
             height = 4.5, width = 6)
sigma_plot_list
dev.off()



# posterior estimates of det probability for Upper Col
p_mat_uc <- extract(dat_tbl_trim$cjs_hier[[5]])[["p_yr"]]
hist(p_mat_uc[ , 1:5, 5], main = "Lower River p") 
hist(p_mat_uc[ , 1:5, 6], main = "Bonneville p") 
hist(p_mat_uc[ , 1:5, 7], main = "Upriver p") 

# posterior estimates of det probability for Low Col
p_mat_lc <- extract(dat_tbl_trim$cjs_hier[[3]])[["p_yr"]]
hist(p_mat_lc[ , 1:5, 5], main = "Lower River p") 
hist(p_mat_lc[ , 1:5, 6], main = "Bonneville p")

# posterior estimates of det probability for Fraser
p_mat_f <- extract(dat_tbl_trim$cjs_hier[[2]])[["p_yr"]]
hist(p_mat_f[ , 1:5, 5], main = "Below Mission p") 
hist(p_mat_f[ , 1:5, 6], main = "Above Mission p")


# estimates of stage specific mean survival rates
med_seg_surv <- purrr::map2(
  dat_tbl_trim$phi_mat_mean, dat_tbl_trim$stock_group,
  ~ .x %>% 
    as.table() %>% 
    as.data.frame() %>% 
    mutate(segment = as.integer(Var2),
           iteration = rep(1:nrow(.x), 
                           times = length(unique(Var2))),
           stock_group = .y
    ) %>% 
    group_by(
      segment, stock_group
    ) %>% 
    summarize(
      med = median(Freq),
      lo = rethinking::HPDI(Freq, prob = 0.9)[1],
      up = rethinking::HPDI(Freq, prob = 0.9)[2]
    ) %>% 
    left_join(., seg_key, by = c("stock_group", "segment")) 
) %>% 
  bind_rows() %>% 
  mutate(
    par = ifelse(
      array_num == max_array_num,
      "beta",
      "phi"
    ),
    stock_group = factor(stock_group, levels = levels(dat_tbl_trim$stock_group))
  ) 
  

stage_spec_surv <- ggplot(med_seg_surv %>% filter(!par == "beta")) +
  geom_pointrange(aes(x = fct_reorder(as.factor(segment_name), segment),
                      y = med, ymin = lo, ymax = up)) +
  facet_wrap(~stock_group, scales = "free_x") +
  labs(x = "Migration Segment", y = "Stage-Specific Survival Rate") +
  ggsidekick::theme_sleek() 

png(here::here("figs", "cjs", "phi_ests.png"), 
    height = 5, width = 7.5, units = "in", res = 200)
stage_spec_surv
dev.off()


# year-specific phi estimates
yr_phi_dat <- purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group,
  function(x, y) {
    yr_phi_mat <- extract(x)[["phi_yr"]] 
    p_sum <- yr_phi_mat %>% 
      as.data.frame.table() %>% 
      rename(year = Var2, segment = Var3) %>% 
      mutate(est = Freq,
             year = as.numeric(as.factor(year)) + 2018,
             array_num = as.numeric(as.factor(segment))) %>% 
      group_by(
        array_num, year
      ) %>% 
      reframe(
        med = median(est),
        lo = rethinking::HPDI(est, prob = 0.9)[1],
        up = rethinking::HPDI(est, prob = 0.9)[2],
        stock_group = y
      ) %>% 
      left_join(., seg_key, by = c("stock_group", "array_num")) %>% 
      mutate(
        segment_name = fct_reorder(as.factor(segment_name), segment)
      ) 
  }
) %>% 
  bind_rows() %>%
  filter(!segment_name == "Release") %>% 
  mutate(
    year = as.factor(year),
    segment_name = factor(
      segment_name,
      levels = c(
        "Release", "WCVI/\nSalish\nSea", "Marine", "NW\nWA", "SW\nWA",
        "Central\nCA",
        "Outside\nShelf", "Juan\nde Fuca", "Strait\nof Georgia",
        "Puget\nSound", "Lower\nCol.", "Bonneville", "Above\nBonneville",
        "In\nRiver", "Downstream\nMission", "Upstream\nMission"
      )),
    stock_group = factor(stock_group, levels = levels(dat_tbl_trim$stock_group))
  )

png(here::here("figs", "cjs", "estimated_yearly_phi.png"), 
    height = 5.5, width = 7.5, units = "in", res = 250)
ggplot(yr_phi_dat) +
  geom_pointrange(
    aes(x = segment_name, y = med, ymin = lo, ymax = up, fill = year),
    shape = 21,
    position = position_dodge(width = 0.3)
  ) +
  facet_wrap(~stock_group, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "PuOr", name = "") +
  ggsidekick::theme_sleek() +
  labs(x = "Segment Name", y = "Detection Probability") 
dev.off()

# calculate mean CV among years, then mean among CV
yr_phi_dat %>% 
  group_by(stock_group, segment) %>% 
  summarize(
    cv = sd(med) / mean(med) 
  ) %>% 
  pull(cv) %>% 
  mean()


# year- and stage-specific detection parameter estimates
yr_p_dat <- purrr::map2(
  dat_tbl_trim$cjs_hier, dat_tbl_trim$stock_group,
  function(x, y) {
    yr_p_mat <- extract(x)[["p_yr"]] 
    p_sum <- yr_p_mat %>% 
      as.data.frame.table() %>% 
      rename(year = Var2, segment = Var3) %>% 
      mutate(est = Freq,
             year = as.numeric(as.factor(year)) + 2018,
             array_num = as.numeric(as.factor(segment))) %>% 
      group_by(
        array_num, year
      ) %>% 
      reframe(
        med = median(est),
        lo = rethinking::HPDI(est, prob = 0.9)[1],
        up = rethinking::HPDI(est, prob = 0.9)[2],
        stock_group = y
      ) %>% 
      left_join(., seg_key, by = c("stock_group", "array_num")) %>% 
      mutate(
        segment_name = fct_reorder(as.factor(segment_name), segment)
      )
  }
) %>% 
  bind_rows() %>%
  filter(!segment_name == "Release") %>% 
  mutate(
    year = as.factor(year),
    segment_name = factor(
      segment_name,
      levels = c(
        "Release", "WCVI/\nSalish\nSea", "Marine", "NW\nWA", "SW\nWA",
        "Central\nCA",
        "Outside\nShelf", "Juan\nde Fuca", "Strait\nof Georgia",
        "Puget\nSound", "Lower\nCol.", "Bonneville", "Above\nBonneville",
        "In\nRiver", "Downstream\nMission", "Upstream\nMission"
      )),
    stock_group = factor(stock_group, levels = levels(dat_tbl_trim$stock_group))
  )


png(here::here("figs", "cjs", "estimated_yearly_p.png"), 
    height = 4.5, width = 9, units = "in", res = 250)
ggplot(yr_p_dat) +
  geom_pointrange(
    aes(x = segment_name, y = med, ymin = lo, ymax = up, fill = year),
    shape = 21,
    position = position_dodge(width = 0.3)
  ) +
  facet_wrap(~stock_group, scales = "free_x", ncol = 2) +
  scale_fill_brewer(palette = "PuOr", name = "") +
  ggsidekick::theme_sleek() +
  labs(x = "Segment Name", y = "Detection Probability") 
dev.off()


## Visualize posterior ---------------------------------------------------------

source(here::here("R", "functions", "plot_survival.R"))

# import data
dat_tbl_trim <- readRDS(
  here::here("data", "model_outputs", "hier_cjs_posterior_tbl.RDS")
)

surv_plot_trials <- purrr::map(dat_tbl_trim$cum_survival, function (x) {
  x %>% 
    mutate(agg_name_f = as.factor(stock_group),
           segment_name = fct_reorder(as.factor(segment_name), segment)) %>% 
    plot_surv(., show_mcmc = T) +
    facet_wrap(~group)
})
surv_plot_clean  <- purrr::map(dat_tbl_trim$cum_survival, function (x) {
  x %>% 
    mutate(agg_name_f = as.factor(stock_group),
           segment_name = fct_reorder(as.factor(segment_name), segment)) %>%
    plot_surv(., show_mcmc = F) +
    facet_wrap(~group)
})

mean_surv_dat <- dat_tbl_trim$cum_survival_mean %>% 
  bind_rows() %>% 
  mutate(agg_name_f = NA,
         # segment_name_f = as.factor(segment_name)
         segment_name = factor(
           segment_name,
           levels = c(
             "Release", "WCVI/\nSalish\nSea", "Marine", "NW\nWA", "SW\nWA",
             "Central\nCA",
             "Outside\nShelf", "Juan\nde Fuca", "Strait\nof Georgia",
             "Puget\nSound", "Lower\nCol.", "Bonneville", "Above\nBonneville",
             "In\nRiver", "Downstream\nMission", "Upstream\nMission"
           )),
         stock_group = factor(
           stock_group,
           levels = c(
             "Cali", "Low Col.", "Up Col.", "Fraser", "South Puget"
           )
         )
  ) %>% 
  select(-c(iter, est)) %>% 
  distinct() 

fill_pal <- c("white", "black")
names(fill_pal) <- c("no", "yes")
shape_pal <- c(21, 24)
names(shape_pal) <- c("phi", "beta")

surv_plot_mean <- ggplot(data = mean_surv_dat) +
  geom_pointrange(
    aes(x = fct_reorder(segment_name, segment), 
        y = median, ymin = low, ymax = up, fill = terminal, shape = par)
  ) +
  ggsidekick::theme_sleek() +
  scale_shape_manual(values = shape_pal, name = "Parameter") +
  scale_fill_manual(values = fill_pal, label = c("Marine", "Terminal"),
                    name = "Terminal") +
  theme(legend.text=element_text(size = 9),
        axis.text.x = element_text(size = rel(.8))) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(~stock_group, scales = "free_x", ncol = 2) +
  labs(x = "Segment Name", y = "Cumulative Survival")


# as above but remove beta estimates 
surv_plot_mean_trim <- mean_surv_dat %>% 
  filter(
    par == "phi"
  ) %>% 
  ggplot(.) +
  geom_pointrange(
    aes(x = fct_reorder(segment_name, segment), 
        y = median, ymin = low, ymax = up, fill = terminal),
    shape = 21
  ) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = fill_pal, label = c("Non-terminal", "Terminal")) +
  theme(legend.text=element_text(size = 9),
        legend.title = element_blank(),
        axis.text.x = element_text(size = rel(.8))) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(~stock_group, scales = "free_x", ncol = 2) +
  labs(x = "Segment Name", y = "Cumulative Survival") +
  ylim(c(0.5, 1)) 


size_pal <- c(0.4, 1)
names(size_pal) <- c("no", "yes")

surv_plot_mean_trim_dist <- mean_surv_dat %>% 
  left_join(
    ., array_dist, by = c("stock_group", "array_num")
  ) %>%
  filter(
    par == "phi",
    !(segment == "5" & stock_group == "Up Col.")
  ) %>% 
  group_by(stock_group) %>% 
  mutate(
    cum_dist = cumsum(dist)
  ) %>% 
  ggplot(.) +
  geom_pointrange(
    aes(x = cum_dist, y = median, ymin = low, ymax = up, fill = stock_group,
        size = terminal),
    shape = 21
  ) +
  geom_line(
    aes(x = cum_dist, y = median, colour = stock_group),
    alpha = 1
  ) +
  ggsidekick::theme_sleek() +
  scale_size_manual(values = size_pal) +
  scale_colour_viridis_d() +
  scale_fill_viridis_d() +
  theme(legend.text=element_text(size = 9),
        legend.title = element_blank(),
        axis.text.x = element_text(size = rel(.8))) +
  guides(size = "none") +
  labs(x = "Migration Distance", y = "Cumulative Survival") +
  ylim(c(0.5, 1)) 


pdf(here::here("figs", "cjs", "cum_surv_ind_hier_trials.pdf"), 
    height = 6, width = 7.5)
surv_plot_trials
dev.off()

pdf(here::here("figs", "cjs", "cum_surv_ind_hier_clean.pdf"), 
    height = 6, width = 7.5)
surv_plot_clean
dev.off()

png(here::here("figs", "cjs", "cum_surv_mean_hier_clean.png"), 
    height = 5, width =7.5, units = "in", res = 200)
surv_plot_mean
dev.off()

png(here::here("figs", "cjs", "cum_surv_mean_hier_clean_nobeta.png"), 
    height = 5, width = 7.5, units = "in", res = 200)
surv_plot_mean_trim
dev.off()

png(here::here("figs", "cjs", "cum_surv_mean_hier_clean_nobeta_dist.png"), 
    height = 5, width = 7.5, units = "in", res = 200)
surv_plot_mean_trim_dist
dev.off()



# plot terminal survival rate (absolute and scaled by migration distance) of 
# high detection probability stocks 
term_surv_dat <- dat_tbl_trim$cum_survival_mean %>% 
  bind_rows() %>% 
  filter(
    stock_group == "Fraser" & segment_name == "Downstream\nMission" |
      # stock_group == "Up Col." & segment_name == "Bonneville" |
      stock_group == "Up Col." & segment_name == "Lower\nCol." |
      stock_group == "Low Col." & segment_name == "Lower\nCol." |
      stock_group == "South Puget" & segment_name == "Puget\nSound" |
      stock_group == "Cali" & segment_name == "Central\nCA"
    ) %>% 
  mutate(
    stock_group = factor(stock_group, levels = levels(dat_tbl_trim$stock_group))
  )

p_total <- term_surv_dat %>% 
  select(stock_group, median, low, up) %>% 
  distinct() %>% 
  ggplot(.) +
  geom_pointrange(aes(x = stock_group, y = median, ymin = low, ymax = up)) +
  labs(y = "Cumulative Terminal\nSurvival Rate") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank())

# import terminal migration distance (to avoid in river harvest focus on 
# Lower Columbia for both UC and LW)
term_dist <- read.csv(here::here("data", "terminal_locations.csv")) %>% 
  filter(!location == "bonneville")

term_surv_dat_scaled <- term_surv_dat %>% 
  select(iter, est, stock_group) %>%
  left_join(., term_dist, by = "stock_group") %>%
  mutate(
    scaled_surv = est^(1 / (distance_km / 100)),
    stock_group = factor(stock_group, levels = levels(dat_tbl_trim$stock_group))
  ) %>% 
  group_by(
    stock_group
  ) %>% 
  summarize(
    median = median(scaled_surv),
    low = rethinking::HPDI(scaled_surv, prob = 0.9)[1],
    up = rethinking::HPDI(scaled_surv, prob = 0.9)[2]
  )

# export for comparison with sockeye 
# rbind(
#   term_surv_dat %>% 
#     select(colnames(term_surv_dat_scaled)) %>% 
#     distinct() %>%
#     mutate(
#       data = "unscaled"
#     ),
#   term_surv_dat_scaled %>% 
#     mutate(
#       data= "scaled"
#     )) %>% 
#   saveRDS(.,
#         here::here("data", "term_surv_rate.rds"))


p_dist <- ggplot(term_surv_dat_scaled) +
  geom_pointrange(aes(x = stock_group, y = median, ymin = lo, ymax = up)) +
  labs(y = "Cumulative Terminal\nSurvival Rate per 100 km") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank())


png(here::here("figs", "cjs", "terminal_surv_comp.png"), 
    height = 7.5, width = 5, units = "in", res = 200)
cowplot::plot_grid(p_total, p_dist, ncol = 1)
dev.off()


# plot stage-specific mean survival but scaled by either migration distance
# or duration (calculated in chinTagging repo segment_distance_duration.R)
# scaled_dat_list <- purrr::map2(
#   dat_tbl_trim$phi_mat_mean, dat_tbl_trim$stock_group,
#   ~ .x %>% 
#     as.table() %>% 
#     as.data.frame() %>% 
#     mutate(segment = as.integer(Var2),
#            iter = rep(1:nrow(.x), 
#                       times = length(unique(Var2))),
#            key_stock_group = .y
#     ) %>% 
#     group_by(
#       segment, key_stock_group
#     ) %>%
#     select(key_stock_group, iter, segment, est = Freq) %>% 
#     left_join(., seg_key %>% rename(key_stock_group = stock_group), 
#               by = c("key_stock_group", "segment")) %>% 
#     full_join(., array_dat, 
#               by = c("key_stock_group", "iter", "array_num_key")) %>% 
#     filter(iter %in% array_dat$iter,
#            !segment_name == "Release") %>% 
#     mutate(
#       scaled_surv_dist = est^(1 / (dist_km / 100)),
#       scaled_surv_dur = est^(1 / (time_dat)),
#       par = ifelse(
#         key_stock_group %in% c("Cali", "Low Col.", "WA_OR", "WCVI") & 
#           segment_name == "In\nRiver",
#         "beta",
#         "phi"
#       )
#       ) 
# ) %>% 
#   bind_rows() %>% 
#   select(-dist_km, -time_dat) %>% 
#   pivot_longer(
#     cols = c("scaled_surv_dist", "scaled_surv_dur"),
#     names_to = "metric",
#     values_to = "surv",
#     names_prefix = "scaled_surv_"
#   ) %>% 
#   group_by(
#     key_stock_group, segment, segment_name, par, metric
#   ) %>% 
#   summarize(
#     med = median(surv),
#     lo = rethinking::HPDI(surv, prob = 0.9)[1],
#     up = rethinking::HPDI(surv, prob = 0.9)[2]
#   ) 
# 
# 
# png(here::here("figs", "cjs", "stage_surv_duration.png"), 
#     width = 7.5, height = 5, units = "in", res = 200)
# scaled_dat_list %>% filter(!par == "beta", metric == "dur") %>% 
#   ggplot() +
#   geom_pointrange(aes(x = fct_reorder(as.factor(segment_name), segment),
#                       y = med, ymin = lo, ymax = up)) +
#   labs(y = "Posterior Cumulative Survival Rate per Day") +
#   ggsidekick::theme_sleek() +
#   theme(axis.title.x = element_blank()) +
#   facet_wrap(~key_stock_group, scales = "free")
# dev.off()
#   
# 
# png(here::here("figs", "cjs", "stage_surv_distance.png"), 
#     width = 7.5, height = 5, units = "in", res = 200)
# scaled_dat_list %>% filter(!par == "beta", metric == "dist") %>% 
#   ggplot() +
#   geom_pointrange(aes(x = fct_reorder(as.factor(segment_name), segment),
#                       y = med, ymin = lo, ymax = up)) +
#   labs(y = "Posterior Cumulative Survival Rate per 100 km") +
#   ggsidekick::theme_sleek() +
#   theme(axis.title.x = element_blank()) +
#   facet_wrap(~key_stock_group, scales = "free")
# dev.off()


## Calculate cumulative survival for Fraser by stock ---------------------------

# stock name key
stk_key <- data.frame(stock = dat_tbl_trim$bio_dat[[2]]$agg) %>%
  mutate(
    agg_n = as.factor(stock) %>% 
      droplevels() %>% 
      as.numeric()
  ) %>% 
  distinct() %>% 
  arrange(agg_n)

phi_stk <- extract(dat_tbl_trim$cjs_hier[[2]])[["phi_stk"]]

# calculate cumulative product for each stk and iteration then convert to DF
dims <- list(iter = seq(1, dim(phi_stk)[1], by = 1),
             segment = seq(1, dim(phi_stk)[3], by = 1))

dumm <- expand.grid(iter = dims$iter,
                    segment = 0, 
                    Freq = 1, 
                    agg_n = stk_key$agg_n)

# calculate the cumulative product across segments for each group and 
# iteration, then convert to a dataframe
cumprod_list <- vector(length(stk_key$agg_n), mode = "list") 
for (i in seq_along(stk_key$agg_n)) {
  cumprod_mat <- t(apply(phi_stk[ , i, ], 1, cumprod))
  dimnames(cumprod_mat) = dims[1:2]
  cumprod_list[[i]] <- cumprod_mat %>% 
    as.table() %>% 
    as.data.frame() %>% 
    mutate(agg_n = stk_key$agg_n[i])
}

stk_cumprod_plot <- cumprod_list %>% 
  bind_rows() %>% 
  rbind(dumm, .) %>%
  mutate(segment = as.integer(as.character(segment)),
         stock_group = "Fraser") %>%
  rename(est = Freq) %>% 
  #add segment key
  left_join(., seg_key, by = c("stock_group", "segment")) %>% 
  left_join(., stk_key, by = "agg_n") %>% 
  mutate(segment_name = fct_reorder(as.factor(segment_name), segment),
         hab = case_when(
           segment == max(segment) ~ "river",
           TRUE ~ "marine"
         )) %>% 
  group_by(segment, segment_name, stock) %>% 
  summarize(median = median(est),
         low = quantile(est, 0.05),
         up = quantile(est, 0.95)) %>% 
  ungroup() %>% 
  mutate(agg_name_f = NA) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = fct_reorder(segment_name, segment), 
                      y = median, ymin = low, ymax = up, fill = stock),
                  shape = 21, position = position_dodge(width = 0.5)) +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.text=element_text(size = 9),
        legend.title = element_blank(),
        axis.text.x = element_text(size = rel(.8)),
        legend.position = "top") +
  lims(y = c(0, 1))  

png(here::here("figs", "cjs", "fraser_stock_surv.png"), 
    height = 4.5, width = 5.5, units = "in", res = 200)
stk_cumprod_plot
dev.off()


stk_effect <- extract(dat_tbl_trim$cjs_hier[[2]])[["alpha_stk_phi"]]
colnames(stk_effect) <- stk_key$stock


png(here::here("figs", "cjs", "fraser_stk_effect.png"), 
    height = 4.5, width = 5.5, units = "in", res = 200)
stk_effect %>% 
  as.data.frame() %>% 
  pivot_longer(
    cols = everything(),
    names_to = "stk",
    values_to = "est"
  ) %>% 
  group_by(
    stk
  ) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>% 
  mutate(
    stk = factor(
      stk, 
      levels = c("Fraser Spr. Yr.", "Fraser Sum. Yr.", "Fraser Sum. 4.1", 
                 "Fraser Fall")),
    stk = fct_recode(stk, 
                     "Spring 1.x" = "Fraser Spr. Yr.",
                     "Summer 1.3" = "Fraser Sum. Yr.", 
                     "Summer 0.3" = "Fraser Sum. 4.1", 
                     "Fall 0.3" = "Fraser Fall")
  ) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = stk, y = med, ymin = lo, ymax = up, fill = stk),
                  shape = 21) +
  scale_fill_brewer(palette = "YlGnBu") +
  labs(x = "Fraser River Stock Group", y = "Stock-Specific Survival Effect") +
  geom_hline(yintercept = 0, lty = 2) +
  ggsidekick::theme_sleek() +
  theme(legend.title = element_blank(),
        legend.position = "none")
dev.off()


prior_df <- data.frame(est = rexp(4000, rate = 1), parameter = "Prior")

stk_sigma <- extract(dat_tbl_trim$cjs_hier[[2]])[["sigma_alpha_stk_phi"]]  %>% 
  as.data.frame() %>%
  pivot_longer(everything(), names_to = "segment", values_to = "est", 
               names_prefix = "V")

png(here::here("figs", "cjs", "fraser_stk_sigma_post_prior.png"), 
    height = 4.5, width = 5.5, units = "in", res = 200)
ggplot() +
  geom_density(data = stk_sigma, aes(x = est), 
               fill = "red", colour = "red", alpha = 0.4) +
  geom_density(data = prior_df, aes(x = est), 
               fill = "blue", colour = "blue", alpha = 0.4) +
  facet_wrap(~segment) + 
  ggsidekick::theme_sleek() +
  labs(x = "Sigma Year Estimate", y = "Kernel Density") 
dev.off()


## CALCULATE POSTERIOR MEANS OF TERMINAL DET -----------------------------------

# posterior estimates of det probability for Upper Col
p_mat_uc <- extract(dat_tbl_trim$cjs_hier[[5]])[["p_yr"]]
hist(p_mat_uc[ , 1:5, 5], main = "Lower River p") 
hist(p_mat_uc[ , 1:5, 6], main = "Bonneville p") 
hist(p_mat_uc[ , 1:5, 7], main = "Upriver p") 

# posterior estimates of det probability for Low Col
p_mat_lc <- extract(dat_tbl_trim$cjs_hier[[3]])[["p_yr"]]
phi_mat_lc <- extract(dat_tbl_trim$cjs_hier[[3]])[["phi_yr"]]

p5 <- p_mat_lc[ , 1, 5]
p6 <- p_mat_lc[ , 1, 6]
phi5 <- phi_mat_lc[ , 1, 5]
p_final <- p5 + ((1 - p5) * phi5 * p6)

# calculate posterior detection probability for all stocks except Columbia 
# (more detection arrays)
p_final_list <- purrr::map(
  dat_tbl_trim$cjs_hier[1:4],
  function (x) {
    p_mat <- extract(x)[["p_yr"]]
    phi_mat <- extract(x)[["phi_yr"]]
    
    if (ncol(p_mat) == ncol(phi_mat)) {
      p_final <- matrix(NA, ncol = ncol(p_mat), nrow = nrow(p_mat))
      for (i in 1:ncol(p_mat)) {
        p5 <- p_mat[ , i, 5]
        p6 <- p_mat[ , i, 6]
        phi5 <- phi_mat[ , i, 5]
        p_final[ , i] <- p5 + ((1 - p5) * phi5 * p6)
      }
    }
    return(p_final)
  }
  )
  

# calculate upper columbia separately since has more terminal areas
p_mat_uc <- extract(dat_tbl_trim$cjs_hier[[5]])[["p_yr"]]
phi_mat_uc <- extract(dat_tbl_trim$cjs_hier[[5]])[["phi_yr"]]

p_final_uc <- matrix(NA, ncol = ncol(p_mat_uc), nrow = nrow(p_mat_uc))
for (i in 1:ncol(p_mat_uc)) {
  p5 <- p_mat_uc[ , i, 5]
  p6 <- p_mat_uc[ , i, 6]
  p7 <- p_mat_uc[ , i, 7]
  phi5 <- phi_mat_uc[ , i, 5]
  phi6 <- phi_mat_uc[ , i, 6]
  p_final_uc[ , i] <-  p5 + 
    (1 - p5) * phi5 * p6 + 
    (1 - p5) * phi5 * (1 - p6) * phi6 * p7 - 
    (1 - p5) * phi5 * p6 * phi6 * p7
}

# make column names years
p_final_list2 <- lapply(append(p_final_list, list(p_final_uc)),
                        function(mat) {
  colnames(mat) <- seq(2019, 2023, by = 1)
  return(mat)
})

# add to tbl
post_det_p_tbl <- tibble(
  stock_group = dat_tbl_trim$stock_group,
  final_p_mat = p_final_list2
) %>% 
  mutate(
    p_out = purrr::map(
      final_p_mat, 
      ~ .x %>% 
        rethinking::logit() %>% 
        as.data.frame() %>% 
        pivot_longer(cols = everything(), names_to = "year", values_to = "logit_p")
    )
  )

posterior_dat_out <- post_det_p_tbl %>% 
  select(
    stock_group, p_out
  ) %>% 
  unnest(cols = c(p_out)) %>% 
  group_by(stock_group, year) %>% 
  summarize(
    mean_logit_p = mean(logit_p),
    sd_logit_p = sd(logit_p)
  ) %>% 
  ungroup()

saveRDS(
  posterior_dat_out,
  here::here("data", "posterior_p.rds")
)
