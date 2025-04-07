## Logistic Survival Sensitivity Analysis
## Compare parameter estimates after a) fixing alternative maturation schedules
# or b) removing individuals with presumed tagging mortality
# Jan. 6, 2025

library(tidyverse)
library(rethinking)
library(rstan)
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
    cyer_z = scale(focal_er) %>% as.numeric(),
    terminal_p = as.factor(terminal_p),
    inj = as.integer(as.factor(adj_inj))
    )


# switch maturity stage
det_dat1 <- det_dat_in %>% 
  filter(
    !is.na(focal_er),
    stage_2 == "mature"
  )


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
    !is.na(focal_er),
    vemco_code %in% five_d_tags & vemco_code %in% mature_tags1
  ) 


dat_list_in <- list(det_dat1, det_dat2)


# import posterior predictions
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


## FIT MODELS ------------------------------------------------------------------

mod1 <- stan_model(here::here("R", "stan_models", "obs_surv_jll_cov2.stan"))
post_list <- sens_list <- vector(mode = "list", length = 2)

for (i in seq_along(dat_list_in)) {
  x <- dat_list_in[[i]]
  
  # make key within for loop since datasets have different stock-year 
  # combinations which impacts indexing
  key <- x %>% 
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
    #replace NAs with placeholders
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
  
  x2 <- x %>% 
    left_join(
      .,
      key %>% select(stock_group, year, det_group_id_n), 
      by = c("stock_group", "year")
    ) %>% 
    mutate(
      lipid_z = scale(lipid) %>% as.numeric(),
      fl_z = scale(fl) %>% as.numeric(),
      wt_z = scale(wt) %>% as.numeric(),
      log_wt_z = scale(log(wt)) %>% as.numeric(),
      day_z = scale(year_day) %>% as.numeric(),
      cyer_z = scale(focal_er) %>% as.numeric(),
      year = as.factor(year),
      stock_group = factor(
        stock_group, 
        levels = c(
          "Low Col.", "Up Col.", "WA_OR", "WCVI", "ECVI", 
          "Fraser Spr. Yr.", "Fraser Sum. Yr.", "Fraser Sum. 4.1", "Fraser Fall", 
          "North Puget", "South Puget"
        )) %>% 
        droplevels(),
      inj = as.integer(as.factor(adj_inj)),
      term_p = as.integer(terminal_p),
      yr = as.integer(year),
      stk = as.integer(stock_group)
    ) 
  
  dat_list <- list(
    N = nrow(x2),
    N_year = length(unique(x2$year)),
    N_stock = length(unique(x2$stk)),
    N_det_id = length(unique(key$det_group_id_n)),
    N_det_id_obs = sum(key$post),
    
    s_obs = as.integer(x2$term_det),
    
    logit_det_sd = key$sd_logit_p,
    logit_det_mu = key$mean_logit_p,
    use_posterior = key$post,
    det_group_id = x2$det_group_id_n,
    
    yr = x2$yr,
    stk_n = x2$stk,
    inj = x2$inj,
    fl_z = x2$fl_z,
    lipid_z = x2$lipid_z,
    day_z = x2$day_z,
    cyer_z = x2$cyer_z,
    alpha_i = rep(2, length(unique(x2$inj)) - 1)
  )
  
  sens_list[[i]] <- sampling(
    mod1, data = dat_list,
    chains = 4, iter = 2000, warmup = 1000,
    control = list(adapt_delta = 0.96)
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
  here::here("data", "model_outputs", "hier_binomial_cyer_stan.rds")
  ) %>% 
  extract.samples()

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
      dims <- dim(dum)
      ndims <- if (is.null(dims)) 1 else length(dims)
      
      if (ndims == 1) {
        # Scalar or vector
        out <- data.frame(
          parameter = x,
          med = median(dum),
          lo = rethinking::HPDI(dum, prob = 0.9)[1],
          up = rethinking::HPDI(dum, prob = 0.9)[2]
        )
        
      } else if (ndims == 2) {
        # Matrix: iterations × parameter columns
        colnames(dum) <- paste(x, seq_len(ncol(dum)), sep = "_")
        out <- as.data.frame(dum) %>%
          pivot_longer(
            cols = everything(), names_to = "parameter", values_to = "est"
          ) %>%
          group_by(parameter) %>%
          summarize(
            med = median(est),
            lo = rethinking::HPDI(est, prob = 0.9)[1],
            up = rethinking::HPDI(est, prob = 0.9)[2],
            .groups = "drop"
          )
        
      } else {
        warning(paste("Unsupported shape for parameter:", x))
        out <- data.frame(parameter = x, med = NA, lo = NA, up = NA)
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
