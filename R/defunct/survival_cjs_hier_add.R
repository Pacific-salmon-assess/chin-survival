### CJS survival models
## June 20, 2020
## Same as survival_cjs_hier.R, but uses simpler CJS models without interaction
# between year and stage


library(tidyverse)
library(rstan)
library(shinystan)
library(loo)

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
  here::here("R", "stan_models", "cjs_add_hier_eff_adaptive_v2.stan")
)
hier_mod_sims_stk <- stan_model(
  here::here("R", "stan_models", "cjs_add_hier_eff_adaptive_v2_stock.stan")
)


## FIT --------------------------------------------------------------------

# MCMC settings
n_chains = 4
n_iter = 2000
n_warmup = n_iter / 2
pars_in <- c(
  "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
  "alpha_p", "alpha_yr_p",
  # transformed pars or estimated quantities
  "alpha_yr_phi", "phi_yr", "p_yr", "beta_yr", "y_hat", "log_lik"
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
          alpha_yr_phi_z = rnorm(x$nyear, 0, 0.5),
          alpha_stk_phi_z = rnorm(x$nstock, 0, 0.5),
          alpha_t_phi = rnorm(x$n_occasions - 1, 0, 0.5),
          sigma_alpha_yr_phi = rexp(1, 2),
          sigma_alpha_stk_phi = rexp(1, 2),
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
          alpha_yr_phi_z = rnorm(x$nyear, 0, 0.5),
          alpha_t_phi = rnorm(x$n_occasions - 1, 0, 0.5),
          sigma_alpha_yr_phi = rexp(1, 2),
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
      control = list(adapt_delta = 0.94)
    )
  })


saveRDS(cjs_hier_sims,
        here::here("data", "model_outputs", "hier_cjs_fit_tbl_add.RDS"))

cjs_hier_sims <- readRDS(
  here::here("data", "model_outputs", "hier_cjs_fit_tbl_add.RDS")
) 

# add to tibble
dat_tbl_trim$cjs_hier <- cjs_hier_sims

# extract and save looic for comparison
loo_list <- purrr::map(
  cjs_hier_sims, ~ loo(.x)
)
saveRDS(loo_list,
        here::here("data", "model_outputs", "hier_cjs_add_looic.RDS"))


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


# posterior predictions check
# compare predicted to observed proportions for each year and stage by aggregate
pp_list <- pmap(
  list(dat_tbl_trim$dat_in, dat_tbl_trim$cjs_hier, 
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
  }
)

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

pdf(here::here("figs", "diagnostics", "posterior_checks_hier_add.pdf"),
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
    
    png(here::here("figs", "cjs", "posterior_prior_comp_add", file_name), 
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
    
    png(here::here("figs", "cjs", "posterior_prior_comp_add", file_name), 
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
    
    png(here::here("figs", "cjs", "posterior_prior_comp_add", file_name), 
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
    
    png(here::here("figs", "cjs", "posterior_prior_comp_add", file_name), 
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
    
    png(here::here("figs", "cjs", "posterior_prior_comp_add", file_name), 
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



## Parameter estimates ---------------------------------------------------------

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

phi_stk <- extract(dd[[1]])[["phi_stk"]]

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


stk_effect <- extract(dd[[1]])[["alpha_stk_phi"]]
colnames(stk_effect) <- stk_key$stock

library(ggridges)

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


summ <- rstan::summary(dd[[1]])[[1]]
rbind(summ[which(grepl("alpha_phi", row.names(summ))), ],
      summ[which(grepl("alpha_t_phi", row.names(summ))), ]
)
summ[which(grepl("alpha_stk_phi", row.names(summ))), ]
summ[which(grepl("sigma_alpha_yr_phi", row.names(summ))), ]
sigma_alpha_yr_phi
     


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
