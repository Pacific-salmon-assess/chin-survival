### CJS Sensitivity Analyses
## Compare phi and cumulative survival estimates after a) fixing alternative
# maturation schedules or b) removing individuals with presumed tagging 
# mortality


library(tidyverse)
library(rstan)
library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)

## Import occurence matrix (generated in chinTagging/prep_detection_histories.R 
# and cleaned in data_clean.R)
dat_tbl_trim_in <- readRDS(here::here("data", "surv_cjs_data.rds")) 


# Select alternative maturation stage
mature_tags <- dat_tbl_trim_in %>% 
  unnest(bio_dat) %>%
  filter(stage_2 == "mature") %>% 
  pull(vemco_code)

dat_tbl_trim1 <- dat_tbl_trim_in %>% 
  mutate(
    bio_dat = purrr::map(
      bio_dat, ~.x %>% filter(vemco_code %in% mature_tags)
    ),
    wide_array_dat = purrr::map(
      wide_array_dat, ~.x %>% filter(vemco_code %in% mature_tags)
    ),
    data = "mature_2"
  )


# Remove individuals not observed within five days of tagging
mature_tags1 <- dat_tbl_trim_in %>% 
  unnest(bio_dat) %>%
  filter(stage_1 == "mature") %>% 
  pull(vemco_code)

five_d_tags <- readRDS(
  here::here("data", "detections_all.RDS")) %>% 
  group_by(vemco_code) %>% 
  mutate(
    release_time = min(date_time),
    time_diff_days =  as.numeric(
      difftime(date_time, release_time, units = "days")
      )
  ) %>% 
  filter(
    time_diff_days > 5 & vemco_code %in% mature_tags1
  ) %>% 
  pull(vemco_code) %>% 
  unique()
  
dat_tbl_trim2 <- dat_tbl_trim_in %>% 
  mutate(
    bio_dat = purrr::map(
      bio_dat, ~.x %>% filter(vemco_code %in% five_d_tags)
    ),
    wide_array_dat = purrr::map(
      wide_array_dat, ~.x %>% filter(vemco_code %in% five_d_tags)
    ),
    data = "five_day"
  )


dat_tbl_trim <- rbind(dat_tbl_trim1, dat_tbl_trim2)


dat_tbl_trim$stock_group <- fct_relevel(
  as.factor(dat_tbl_trim$stock_group),
  "Cali", "Low Col.", "Up Col.", "Fraser", "South Puget"
)

## Import survival segment key for labelling plots 
seg_key <- read.csv(here::here("data", 
                               "surv_segment_key_2023.csv")) %>%
  mutate(segment = array_num - 1,
         segment_name = ifelse(
           stock_group == "Up Col." & segment == 6,
           "Above\nBonneville",
           str_replace(segment_name, " ", "\n")
         ),
         array_num_key = paste(segment, segment + 1, sep = "_")) %>% 
  dplyr::select(stock_group, segment, segment_name, array_num, array_num_key,
                max_array_num, terminal) %>% 
  distinct() 


# Import duration and distance estimates for scaling survival
array_dat <- readRDS(here::here("data", "distance_duration_array.rds"))
array_dat$iter <- as.numeric(as.factor(array_dat$iteration))


# Fit model --------------------------------------------------------------------

## Same model structure as survival_cjs_hier.R

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


## FIT -------------------------------------------------------------------------

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

saveRDS(cjs_hier_sims,
        here::here("data", "model_outputs", "hier_cjs_fit_tbl_sens.RDS"))

cjs_hier_sims <- readRDS(
  here::here("data", "model_outputs", "hier_cjs_fit_tbl_sens.RDS")
  )


## CHECKS ----------------------------------------------------------------------

# neff 
purrr::map(cjs_hier_sims, function (x) {
  tt <- summary(x)
  tt2 <- tt$summary[ , "n_eff"]
  tt2[which(tt2 < 1000)]
})

# rhat
purrr::map(cjs_hier_sims, function (x) {
  tt <- bayesplot::rhat(x)
  tt[which(tt > 1.05)]
})



## Post-hoc calculations -------------------------------------------------------

# extract phi matrix and swap last col with beta estimates (i.e. combined p and 
# phi) except for fix p models 
# stocks
phi_mat <- map(
  cjs_hier_sims, 
  function(x) {
    phi_adj <- extract(x)[["phi_yr"]]
    # replace array corresponding to last stage-specific survival est, w/ beta
    phi_adj[ , , dim(phi_adj)[3]] <- extract(x)[["beta_yr"]]
    return(phi_adj)
  })


# calculate cumulative survival across segments for each year
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
phi_mat_mean <- map(
  cjs_hier_sims,
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

# integrate out year effects
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
        here::here("data", "model_outputs", "hier_cjs_posterior_tbl_sens.RDS"))

dat_tbl_trim <- readRDS(
  here::here("data", "model_outputs", "hier_cjs_posterior_tbl_sens.RDS")
  )


# combine new and old tbls
standard_tbl <- readRDS(
  here::here("data", "model_outputs", "hier_cjs_posterior_tbl.RDS")
) %>% 
  mutate(
    data = "standard"
  ) %>% 
  select(
    colnames(dat_tbl_trim)
  )


dat_tbl <- rbind(standard_tbl, dat_tbl_trim) %>% 
  mutate(
    data = factor(data, levels = c("standard", "mature_2", "five_day"),
                  labels = c("standard", "maturity", "tag effect"))
  )


## VISUALIZE POSTERIOR ---------------------------------------------------------


# estimates of stage specific mean survival rates
med_seg_surv <- purrr::pmap(
  list(dat_tbl$phi_mat_mean, dat_tbl$stock_group, dat_tbl$data),
  function(x, y, z) {
    x %>% 
      as.table() %>% 
      as.data.frame() %>% 
      mutate(segment = as.integer(Var2),
             iteration = rep(1:nrow(x), 
                             times = length(unique(Var2))),
             stock_group = y,
             data = z
      ) %>% 
      group_by(
        segment, stock_group, data
      ) %>% 
      summarize(
        med = median(Freq),
        lo = rethinking::HPDI(Freq, prob = 0.9)[1],
        up = rethinking::HPDI(Freq, prob = 0.9)[2]
      ) %>% 
      left_join(., seg_key, by = c("stock_group", "segment")) 
  }
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
  geom_pointrange(
    aes(x = fct_reorder(as.factor(segment_name), segment),
        y = med, ymin = lo, ymax = up, fill = data), 
    shape = 21,
    position = position_dodge(width = 0.4)
  ) +
  facet_wrap(~stock_group, scales = "free_x") +
  labs(x = "Migration Segment", y = "Stage-Specific Survival Rate") +
  ggsidekick::theme_sleek() +
  theme(
    legend.position = "top"
  )

png(here::here("figs", "sens", "phi_ests.png"), 
    height = 5, width = 7.5, units = "in", res = 200)
stage_spec_surv
dev.off()


# cumulative survival
mean_surv_dat <- dat_tbl %>% 
  select(data, cum_survival_mean) %>% 
  unnest(cols= cum_survival_mean) %>% 
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
           ))
  ) %>% 
  select(-c(iter, est)) %>% 
  distinct() %>% 
  filter(
    !par == "beta"
  )

surv_plot_mean <- ggplot(data = mean_surv_dat) +
  geom_pointrange(
    aes(x = fct_reorder(segment_name, segment), 
        y = median, ymin = low, ymax = up, fill = data),
    shape = 21,
    position = position_dodge(width = 0.4)
  ) +
  ggsidekick::theme_sleek() +
  theme(legend.text=element_text(size = 9),
        axis.text.x = element_text(size = rel(.8))) +
  # guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_wrap(~stock_group, scales = "free_x", ncol = 2) +
  labs(x = "Segment Name", y = "Cumulative Survival") +
  theme(
    legend.position = "top"
  )

png(here::here("figs", "sens", "cumulative_surv.png"), 
    height = 5, width = 7.5, units = "in", res = 200)
surv_plot_mean
dev.off()


# difference in cumulative survival
cs_list <- split(dat_tbl, dat_tbl$data) %>% 
  purrr::map(
    ., 
    ~ .x %>% 
      select(cum_survival_mean) %>% 
      unnest(cols = c(cum_survival_mean)) %>% 
      select(iter:array_num, par) %>% 
      filter(!par == "beta") %>% 
      group_by(stock_group) %>% 
      mutate(
        new_max = max(array_num)
      ) %>% 
      filter(
        array_num == new_max
      ) %>% 
      ungroup()
  )

mat_stage_effect <- cs_list[[1]] %>% 
  mutate(
    sens_est = cs_list[[2]]$est,
    diff = sens_est - est,
    data = "Maturation Stage"
  )

tag_effect <- cs_list[[1]] %>% 
  mutate(
    sens_est = cs_list[[3]]$est,
    diff = sens_est - est,
    data = "Tag Effect"
  )


png(here::here("figs", "sens", "cumulative_surv_diff.png"), 
    height = 5, width = 5, units = "in", res = 200)
ggplot(rbind(mat_stage_effect, tag_effect)) +
  geom_histogram(aes(x = diff), bins = 50, fill = "#3690c0") +
  geom_vline(aes(xintercept = 0), lty = 2 , colour = "black", linewidth = 1.25) +
  ggsidekick::theme_sleek() +
  labs(x = "Difference in Cumulative Survival\nRelative to Standard Model") +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank()
  ) +
  facet_wrap(~data, ncol = 1)
dev.off()
