## Logistic Survival Sensitivity Analysis
## Compare parameter estimates after a) fixing alternative maturation schedules
# or b) removing individuals with presumed tagging mortality
# Jan. 6, 2025

library(tidyverse)
library(rstan)
library(tidybayes)
library(grid)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

det_dat_in <- readRDS(here::here("data", "surv_log_reg_data.rds")) %>% 
  filter(
    !is.na(focal_er_adj)
  ) %>% 
  droplevels()


# import posterior estimates of detection probability based on CJS models
# add high detection probability for WCVI in 2021 and 2022 because a) all 
# tagged fish Robertson Creek and b) 100% of in-river detections observed on
# multiple arrays
post_p_wcvi <- data.frame(
  stock_group = "WCVI",
  year = c("2021", "2022"),
  mean_logit_p = 5.64,
  sd_logit_p = 1.09
)
post_p <- readRDS(here::here("data", "posterior_p.rds")) %>% 
  rbind(., post_p_wcvi) %>% 
  mutate(det_group_id = paste(as.character(stock_group), year, sep = "_")) %>% 
  select(-c(stock_group, year))
key <- det_dat_in %>% 
  mutate(
    det_group_id = ifelse(
      grepl("Fraser", stock_group), "Fraser", as.character(stock_group)
    ) %>% 
      paste(., year, sep = "_")
  ) %>% 
  select(stock_group, year, det_group_id) %>% 
  distinct() %>% 
  arrange(stock_group, year, det_group_id) %>% 
  left_join(., post_p, by = "det_group_id") %>% 
  #replace NAs with placeholders (not used in fitting, replaced by priors)
  mutate(
    post = ifelse(is.na(mean_logit_p), 0, 1),
    mean_logit_p = ifelse(is.na(mean_logit_p), 0.001, mean_logit_p),
    sd_logit_p = ifelse(is.na(sd_logit_p), 1, sd_logit_p)
  ) %>% 
  # sort so that unobserved stocks follow observed
  arrange(desc(post), stock_group, year, det_group_id) %>% 
  mutate(
    det_group_id_n = seq(1, nrow(.), by = 1)
  )



det_dat_in2 <- det_dat_in %>% 
  left_join(
    .,
    key %>% select(stock_group, year, det_group_id_n), 
    by = c("stock_group", "year")
  ) %>% 
  #scale cov
  mutate(
    lipid_z = scale(lipid) %>% as.numeric(),
    fl_z = scale(fl) %>% as.numeric(),
    wt_z = scale(wt) %>% as.numeric(),
    log_wt_z = scale(log(wt)) %>% as.numeric(),
    day_z = scale(year_day) %>% as.numeric(),
    cyer_z = scale(focal_er_adj) %>% as.numeric(),
    inj = as.integer(as.factor(adj_inj))
  ) %>% 
  mutate(
    year = as.factor(year),
    stock_group = factor(
      stock_group, 
      levels = c(
        "Low Col.", "Up Col.", "WA_OR", "WCVI", "ECVI", 
        "Fraser Spr. 1.x", "Fraser Sum. 1.2", "Fraser Sum. 0.3", "Fraser Fall", 
        "North Puget", "South Puget"
      )) %>% 
      droplevels(),
    yr = as.integer(year),
    stk = as.integer(stock_group)
  ) 


# switch maturity stage
det_dat1 <- det_dat_in2 %>% 
  filter(
    stage_2 == "mature"
  ) 


# five days post-release
mature_tags1 <- det_dat_in2 %>% 
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
det_dat2 <- det_dat_in2 %>% 
  filter(
    vemco_code %in% five_d_tags & vemco_code %in% mature_tags1
  ) 


dat_list_in <- list(det_dat1, det_dat2)



## FIT MODELS ------------------------------------------------------------------

mod1 <- stan_model(
  here::here("R", "stan_models", "obs_surv_jll_cov2_uninformative.stan")
  )


fit_list <- vector(mode = "list", length = 2)

for (i in seq_along(dat_list_in)) {
  x <- dat_list_in[[i]]
  
  dat_list <- list(
    N = nrow(x),
    N_year = length(unique(x$year)),
    N_stock = length(unique(x$stk)),
    N_det_id = length(unique(key$det_group_id_n)),
    N_det_id_obs = sum(key$post),
    
    s_obs = as.integer(x$term_det),
    
    logit_det_sd = key$sd_logit_p,
    logit_det_mu = key$mean_logit_p,
    use_posterior = key$post,
    det_group_id = x$det_group_id_n,
    
    yr = x$yr,
    stk_n = x$stk,
    inj = x$inj,
    fl_z = x$fl_z,
    lipid_z = x$lipid_z,
    day_z = x$day_z,
    cyer_z = x$cyer_z,
    alpha_i = rep(2, length(unique(x$inj)) - 1)
  )
  
  
  fit_list[[i]] <- sampling(
    mod1, data = dat_list,
    chains = 4, iter = 2000, warmup = 1000,
    control = list(adapt_delta = 0.97)
  )
}


saveRDS(
  fit_list, 
  here::here("data", "model_outputs", "hier_binomial_sens_list_samps.rds")
)


## import fit from main text 
fit_list <- readRDS(
  here::here("data", "model_outputs", "hier_binomial_sens_list_samps.rds")
)
original_fit <- readRDS(
  here::here(
    "data", "model_outputs", "hier_binomial_cyer_stan_uninformative_adj.rds"
    )
  )

post_list_full <- list(original_fit, fit_list[[1]], fit_list[[2]])


## EXAMINE POSTERIOR PARAMETER ESTIMATES ---------------------------------------

dat_names <- c("standard", "maturity", "tag effect")
par_list <- vector(mode = "list", length = length(post_list_full))

par_list <- purrr::map2(
  post_list_full, dat_names,
  function (x, y) {
    spread_draws(
      x, 
      beta_cs, beta_ls, beta_fs, beta_ds, beta_ds_cs, beta_is, alpha_s
    ) %>% 
      select(-starts_with(".")) %>% 
      pivot_longer(
        cols = everything(), names_to = "parameter", values_to = "est"
      ) %>% 
      group_by(parameter) %>% 
      summarize(
        med = median(est),
        lo = rethinking::HPDI(est, prob = 0.9)[1],
        up = rethinking::HPDI(est, prob = 0.9)[2]
      ) %>% 
      ungroup() %>% 
      mutate(
        dataset = y
      )
  }
)

all_pars <- par_list %>% 
  bind_rows() %>% 
  mutate(
    dataset = factor(dataset,
                     levels =  c("standard", "maturity", "tag effect")),
    parameter = factor(
      parameter, 
      levels = c("alpha_s", "beta_ds",  "beta_fs", "beta_ls",
                 "beta_cs", "beta_ds_cs", "beta_is"),
      labels = c("Mean Survival", "Date", "Fork Length", 
                 "Lipid", "Exploitation", "Date:Exploitation", "Injury")
      )
  ) 


png(here::here("figs", "sens", "jll_par_ests.png"), 
    height = 7, width = 4.5, units = "in", res = 200)
ggplot(all_pars) +
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
