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
mature_tags1 <- det_dat_in %>% 
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
    time_diff_days > 5
  ) %>% 
  pull(vemco_code) %>% 
  unique()
det_dat2 <- det_dat_in %>% 
  filter(
    !is.na(isbm_cyer),
    vemco_code %in% five_d_tags & vemco_code %in% mature_tags1
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

post_list <- sens_list <- vector(mode = "list", length = 2)

for (i in seq_along(dat_list_in)) {
  x <- dat_list_in[[i]]
  
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
  
  sens_list[[i]] <- ulam(
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
  post_list[[i]] <- extract.samples(sens_list[[i]])
}

saveRDS(
  sens_list, 
  here::here("data", "model_outputs", "hier_binomial_sens_list.rds")
)
saveRDS(
  post_list, 
  here::here("data", "model_outputs", "hier_binomial_sens_list_samps.rds")
)


## import fit from main text 
post_list <- readRDS(
  here::here("data", "model_outputs", "hier_binomial_sens_list_samps.rds")
)
original_post <- readRDS(
  here::here("data", "model_outputs", "hier_binomial_cyer_samps.rds")
  )

post_list_full <- list(original_post, post_list[[1]], post_list[[2]])


## EXAMINE POSTERIOR PARAMETER ESTIMATES ---------------------------------------

dat_names <- c("standard", "maturity", "tag effect")
par_names <- c("alpha_s", "beta_ds", "beta_fs", "beta_ls", "beta_cs",
               "beta_ds_cs", "beta_is")
par_list <- vector(mode = "list", length = length(post_list_full))

for(i in seq_along(par_list)) {
  par_list[[i]] <- purrr::map(
    par_names,
    function (x) {
      dum <- post_list_full[[i]][[x]]
      # for matrix parameters pivot longer
      if (ncol(dum) > 1) {
        colnames(dum) <- paste(x, seq(1, ncol(dum), by = 1), sep = "_")
        out <- dum %>% 
          as.data.frame() %>%
          pivot_longer(
            cols = everything(), names_to = "parameter", values_to = "est"
          ) %>% 
          group_by(parameter) %>% 
          summarize(
            med = median(est),
            lo = rethinking::HPDI(est, prob = 0.9)[1],
            up = rethinking::HPDI(est, prob = 0.9)[2]
          ) 
      } else {
        data.frame(
          parameter = x,
          med = median(dum),
          lo = rethinking::HPDI(dum, prob = 0.9)[1],
          up = rethinking::HPDI(dum, prob = 0.9)[2]
        )
      }
    }
  ) %>% 
    bind_rows() %>% 
    mutate(
      dataset = dat_names[[i]]
    )
}

all_pars <- par_list %>% 
  bind_rows() %>% 
  mutate(
    dataset = factor(dataset,
                     levels =  c("standard", "maturity", "tag effect")),
    parameter = factor(
      parameter, 
      levels = c("alpha_s_1", "alpha_s_2", "beta_ds",  "beta_fs", "beta_ls",
                 "beta_cs", "beta_ds_cs", "beta_is"),
      labels = c("Low Det. Prob.", "High Det. Prob.", "Date", "Fork Length", 
                 "Lipid", "Exploitation", "Date:Exploitation", "Injury")
      )
  ) 



png(here::here("figs", "sens", "int_ests.png"), 
    height = 3, width = 4.5, units = "in", res = 200)
ggplot(all_pars %>% filter(grepl("Prob.", parameter))) +
  geom_pointrange(
    aes(x = dataset, y = med, ymin = lo, ymax = up, fill = dataset),
    shape = 21
  ) +
  facet_wrap(~parameter, ncol = 2, scales = "free") +
  labs(x = "Dataset", y = "Intercept Estimate") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top",
        legend.title = element_blank())
dev.off()


png(here::here("figs", "sens", "slope_ests.png"), 
    height = 7, width = 4.5, units = "in", res = 200)
ggplot(all_pars %>% filter(!grepl("Prob.", parameter))) +
  geom_pointrange(
    aes(x = dataset, y = med, ymin = lo, ymax = up, fill = dataset),
    shape = 21
  ) +
  facet_wrap(~parameter, ncol = 2, scales = "free") +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Dataset", y = "Parameter Estimate") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top",
        legend.title = element_blank())
dev.off()
