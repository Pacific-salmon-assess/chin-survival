### CJS survival models
## June 5, 2025
## Same as survival_cjs_hier.R, but uses simpler CJS models without interaction
# between year and stage as well as covariate for tagging date


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

dat_tbl_trim$wide_array_dat <- purrr::map2(
  dat_tbl_trim$wide_array_dat, dat_tbl_trim$bio_dat,
  ~ .x %>% 
    filter(vemco_code %in% mature_tags) %>% 
    left_join(
      ., 
      # add z-scored tagging date
      .y %>% 
        select(vemco_code, tag_date_z), 
      by = "vemco_code")
)
dat_tbl_trim$stock_group <- fct_relevel(
  as.factor(dat_tbl_trim$stock_group),
  "Cali", "Low Col.", "Up Col.", "Fraser", "South Puget"
)


# Fit model --------------------------------------------------------------------


# Function to convert wide DF to list input for Stan models
prep_cjs_dat <- function(dat, grouping_vars = NULL) {
  y <- dat %>% 
    dplyr::select(matches("[1-9]"))
  n_occasions <- ncol(y)
  nind <- nrow(y)
  
  tag_date_z <- dat$tag_date_z
  
  out_list <- list(y = y, n_occasions = n_occasions, nind = nind, 
                   tag_date_z = tag_date_z)
  
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
    dd[!(grepl("^([0-9]+)$", dd) | dd == "vemco_code" | dd == "tag_date_z")]
  })
dat_tbl_trim$dat_in <- pmap(
  list(dat = dat_tbl_trim$wide_array_dat,  
       grouping_vars = dat_tbl_trim$grouping_vars),
  .f = prep_cjs_dat
) 

# Call Stan from R and fit to each aggregate separately 
hier_mod_sims <- stan_model(
  here::here("R", "stan_models", "cjs_int_hier_eff_adaptive_date.stan")
)
hier_mod_sims_stk <- stan_model(
  here::here("R", "stan_models", "cjs_int_hier_eff_adaptive_date_stock.stan")
)


## FIT --------------------------------------------------------------------

# MCMC settings
n_chains = 4
n_iter = 2000
n_warmup = n_iter / 2
pars_in <- c(
  "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
  "beta_date_phi", "L_Rho_yr", "alpha_p", "alpha_yr_p",
  # transformed pars or estimated quantities
  "Rho_yr", "alpha_yr_phi", "phi_yr", "p_yr", "beta_yr", "y_hat", "log_lik"
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
          ),
          beta_date_phi = rnorm(1, 0, 0.5)
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
                              nrow = x$nyear),
          beta_date_phi = rnorm(1, 0, 0.5)
        )
      })
    }
    
    sampling(
      mod, data = x, pars = pars,
      init = inits, chains = n_chains, iter = n_iter, warmup = n_warmup,
      open_progress = FALSE,
      control = list(adapt_delta = 0.94)
    )
  }
  )

saveRDS(cjs_hier_sims,
        here::here("data", "model_outputs", "hier_cjs_fit_tbl_int_date.RDS"))

cjs_hier_sims <- readRDS(
  here::here("data", "model_outputs", "hier_cjs_fit_tbl_int_date.RDS")
) 

# extract and save looic for comparison
loo_list <- purrr::map(
  cjs_hier_sims, ~ loo(.x)
)
saveRDS(loo_list,
        here::here("data", "model_outputs", "hier_cjs_int_date_looic.RDS"))

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
  "hier_mcmc_trace_date_int.pdf",
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

pdf(here::here("figs", "diagnostics", "posterior_checks_hier_date_int.pdf"),
    height = 5, width = 8)
pp_plot_list
dev.off()
