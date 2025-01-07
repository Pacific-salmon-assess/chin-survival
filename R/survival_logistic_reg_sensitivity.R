## Logistic Survival Sensitivity Analysis
## Compare parameter estimates after a) fixing alternative maturation schedules
# or b) removing individuals with presumed tagging mortality
# Jan. 6, 2025

library(tidyverse)
library(rethinking)
library(grid)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

det_dat_in <- readRDS(here::here("data", "surv_log_reg_data.rds")) %>% 
  #scale cov
  mutate(
    lipid_z = scale(lipid) %>% as.numeric(),
    fl_z = scale(fl) %>% as.numeric(),
    wt_z = scale(wt) %>% as.numeric(),
    log_wt_z = scale(log(wt)) %>% as.numeric(),
    day_z = scale(year_day) %>% as.numeric(),
    cyer_z = scale(isbm_cyer) %>% as.numeric(),
    terminal_p = as.factor(terminal_p),
    inj = as.integer(as.factor(adj_inj)),
    term_p = as.integer(terminal_p)
  )


# switch maturity stage
det_dat1 <- det_dat_in %>% 
  filter(
    !is.na(isbm_cyer),
    stage_2 == "mature"
  ) %>% 
  mutate(
    year = as.factor(year),
    stock_group = factor(
      stock_group, 
      levels = c(
        "Low Col.", "Up Col.", "WA_OR", "WCVI", "ECVI", 
        "Fraser Spr. Yr.", "Fraser Sum. Yr.", "Fraser Sum. 4.1", "Fraser Fall", 
        "North Puget", "South Puget"
      )) %>% 
      droplevels(),
    yr = as.integer(year),
    stk = as.integer(stock_group)
  ) %>% 
  droplevels()


# five days post-release
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
    time_diff_days > 5
  ) %>% 
  pull(vemco_code) %>% 
  unique()
det_dat2 <- det_dat_in %>% 
  filter(
    !is.na(isbm_cyer),
    vemco_code %in% five_d_tags
  ) %>% 
  mutate(
    year = as.factor(year),
    stock_group = factor(
      stock_group, 
      levels = c(
        "Low Col.", "Up Col.", "WA_OR", "WCVI", "ECVI", 
        "Fraser Spr. Yr.", "Fraser Sum. Yr.", "Fraser Sum. 4.1", "Fraser Fall", 
        "North Puget", "South Puget"
      )) %>% 
      droplevels(),
    yr = as.integer(year),
    stk = as.integer(stock_group)
  ) 


dat_list_in <- list(det_dat1, det_dat2)


## FIT MODELS ------------------------------------------------------------------

path_name <- paste("hier_binomial_cyer", 
                   c("sens_mature.rds", "sens_tag_effect.rds"),
                   sep = "_")


fit_list <- purrr::map2(
  dat_list_in, path_name,
  function (x, y) {
    dat_list <- list(
      surv = as.integer(x$term_det),
      fl_z = x$fl_z,
      lipid_z = x$lipid_z,
      day_z = x$day_z,
      cyer_z = x$cyer_z,
      inj = x$inj,
      term_p = x$term_p,
      yr = x$yr,
      stk = x$stk,
      alpha = rep(2, length(unique(x$inj)) - 1)
    )
    
    fit <- ulam(
      alist(
        # stock-specific sampling dates
        day_z ~ normal(mu_day, sigma_day),
        mu_day <- beta_stk[stk, 4],
        # unobserved latent condition
        c(fl_z, lipid_z) ~ multi_normal(c(mu_fl, mu_lipid), Rho, Sigma),
        mu_fl <- alpha_fl + beta_df * day_z + beta_stk[stk, 1] + beta_yr[yr, 1],
        mu_lipid <- alpha_lipid + beta_dl * day_z + beta_yr[yr, 2] +
          beta_stk[stk, 2],
        # survival
        surv ~ dbinom(1 , p) ,
        logit(p) <- alpha_s[term_p] +
          beta_yr[yr, 3] +
          beta_stk[stk, 3] +
          beta_ds * day_z +
          beta_fs * fl_z +
          beta_ls * lipid_z +
          beta_cs * cyer_z + 
          (beta_ds_cs * day_z * cyer_z) + 
          beta_is * sum(delta_inj[1:inj])
        ,
        # adaptive priors
        transpars> matrix[yr, 3]:beta_yr <-
          compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
        matrix[3, yr]:z_yr ~ normal(0, 0.5),
        transpars> matrix[stk, 4]:beta_stk <-
          compose_noncentered(sigma_stk , L_Rho_stk , z_stk),
        matrix[4, stk]:z_stk ~ normal(0, 0.5),
        # fixed priors
        c(alpha_fl, alpha_lipid) ~ normal(0, 0.5),
        c(beta_df, beta_dl) ~ normal(0.2, 0.5),
        Rho ~ lkj_corr(2),
        Sigma ~ exponential(1),
        sigma_day ~ exponential(1),
        
        cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
        vector[3]:sigma_yr ~ exponential(1),
        cholesky_factor_corr[4]:L_Rho_stk ~ lkj_corr_cholesky(2),
        vector[4]:sigma_stk ~ exponential(1),
        
        # strong rationale for eg negative effect of exploitation and pos effect of 
        # interaction; assign weakly informative priors
        alpha_s[term_p] ~ dnorm(0, 1),
        c(beta_ds, beta_fs, beta_ls, beta_ds_cs) ~ normal(0.2, 0.5),
        c(beta_is, beta_cs) ~ normal(-0.2, 0.5),
        
        # constraints on ordinal effects of injury
        vector[4]: delta_inj <<- append_row(0, delta),
        simplex[3]: delta ~ dirichlet(alpha),
        
        # compute ordinary correlation matrixes from Cholesky factors
        gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr),
        gq> matrix[4, 4]:Rho_stk <<- Chol_to_Corr(L_Rho_stk)
      ),
      data = dat_list, chains = 4, cores = 4, iter = 2000,
      control = list(adapt_delta = 0.95)
    )
    
    saveRDS(fit, 
            here::here("data", "model_outputs", y))
  }
)



fit_list <- readRDS(
  here::here("data", "model_outputs", "hier_binomial_cyer_sens.rds")
  )


## import fit from main text 
original_fit <- readRDS(
  here::here("data", "model_outputs", "hier_binomial_cyer.rds")
  )

fit_list_full <- c(original_fit, fit_list)


## EXAMINE POSTERIOR PARAMETER ESTIMATES ---------------------------------------

dd <- fit_list_full[[1]]

par_names <- "alpha_s", "beta_ds", "beta_fs", "beta_ls", "beta_cs", "beta_ds_cs", "beta_is")
post <- extract.samples(dd)
