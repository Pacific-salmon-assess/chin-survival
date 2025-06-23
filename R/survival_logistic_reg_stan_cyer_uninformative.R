## Equivalent to survival_logistic_reg_stan_cyer.R but with uninformative priors
# April 3, 2025

library(tidyverse)
library(rethinking)
library(grid)
library(rstan)
library(tidybayes)


rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



det_dat1 <- readRDS(here::here("data", "surv_log_reg_data.rds")) %>% 
  filter(
    !is.na(focal_er),
    stage_1 == "mature"
  )


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
key <- det_dat1 %>% 
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


# scale continuous covariates
det_dat <- det_dat1 %>% 
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
    cyer2_z = scale(focal_er_no_ps) %>% as.numeric(),
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


# FIT LOGISTIC REG MODELS ------------------------------------------------------


dat_list <- list(
  N = nrow(det_dat),
  N_year = length(unique(det_dat$year)),
  N_stock = length(unique(det_dat$stk)),
  N_det_id = length(unique(key$det_group_id_n)),
  N_det_id_obs = sum(key$post),
  
  s_obs = as.integer(det_dat$term_det),
  
  logit_det_sd = key$sd_logit_p,
  logit_det_mu = key$mean_logit_p,
  use_posterior = key$post,
  det_group_id = det_dat$det_group_id_n,
  
  yr = det_dat$yr,
  stk_n = det_dat$stk,
  inj = det_dat$inj,
  fl_z = det_dat$fl_z,
  lipid_z = det_dat$lipid_z,
  day_z = det_dat$day_z,
  cyer_z = det_dat$cyer_z,
  alpha_i = rep(2, length(unique(det_dat$inj)) - 1)
)


# mod1 <- stan_model(here::here("R", "stan_models", "obs_surv_jll_cov2_uninformative.stan"))
# m1_stan <- sampling(mod1, data = dat_list,
#                     chains = 4, iter = 2000, warmup = 1000,
#                     control = list(adapt_delta = 0.96))
# saveRDS(m1_stan,
#         here::here("data", "model_outputs", "hier_binomial_cyer_stan_uninformative.rds"))
# # as above but with alternative CYER index
# dat_list$cyer_z <-  det_dat$cyer2_z
# m1_stan_no_ps <- sampling(mod1, data = dat_list,
#                     chains = 4, iter = 2000, warmup = 1000,
#                     control = list(adapt_delta = 0.96))
# saveRDS(
#   m1_stan_no_ps,
#   here::here("data", "model_outputs", "hier_binomial_cyer_stan_uninformative_no_ps.rds")
# )


m1_stan <- readRDS(
  here::here("data", "model_outputs", "hier_binomial_cyer_stan_uninformative.rds"))


# check problematic params
summary_df <- summary(m1_stan)$summary %>% 
  as.data.frame()
summary_df$parameter <- rownames(summary_df)

# Filter based on Rhat and n_eff
subset(summary_df, Rhat > 1.05 | n_eff < 1000)

summary_df %>% 
  filter(grepl("beta", parameter))

# similar effects regardless of whether Puget Sound ISBM fisheries included
loo1 <- loo(m1_stan)
loo2 <- loo(m1_stan_no_ps)


# DETECTION PROBABILITY --------------------------------------------------------

det_p_draws <- as_draws_df(m1_stan) %>%
  spread_draws(p[i]) %>% 
  rename(det_group_id_n = i) %>% 
  left_join(., key, by = "det_group_id_n") %>% 
  ungroup()

terminal_det_p <- ggplot() +
  geom_boxplot(data = det_p_draws, aes(x = as.factor(year), y = p)) +
  geom_point(data = key %>% filter(post == "1"),
             aes(x = as.factor(year), y = inv_logit(mean_logit_p)), 
             colour = "red") +
  labs(y = "Terminal Detection Probability") +
  facet_wrap(~stock_group) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# POSTERIOR INFERENCE  ---------------------------------------------------------

post <- extract.samples(m1_stan)

yday_seq <- c(135, 182, 227) 
day_label <- c("May 15", "Jul 1", "Aug 15")
day_seq <- (yday_seq - mean(det_dat$year_day)) / sd(det_dat$year_day)
cyer_seq <- seq(min(det_dat$cyer_z), max(det_dat$cyer_z), length = 40)
preds <- vector(mode = "list", length = length(day_seq))

for(i in seq_along(day_seq)) {
  pred_x <- sapply(
    cyer_seq, 
    function (x) {
      inv_logit(
        post$alpha_s + post$beta_ds * day_seq[i] +
          post$beta_cs * x + 
          (post$beta_ds_cs * day_seq[i] * x)
      )
    }
  )
  preds[[i]] <- pred_x %>% 
    as.data.frame() %>% 
    set_names(cyer_seq) %>% 
    mutate(iter = seq(1, nrow(.), by = 1)) %>% 
    pivot_longer(
      cols = -iter, names_to = "cyer_z", values_to = "est"
    ) %>% 
    mutate(
      cyer_z = as.numeric(cyer_z),
      cyer = (cyer_z * sd(det_dat$focal_er)) + mean(det_dat$focal_er),
      day = day_label[i]
    )
}

pred_day_cyer_dat <- bind_rows(preds) %>% 
  mutate(
    day = factor(day, levels = c("May 15", "Jul 1", "Aug 15"))
  ) %>% 
  group_by(day, cyer) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  )

pred_day_cyer <- ggplot(
  pred_day_cyer_dat, aes(x = cyer, y = med)
  ) +
  geom_line(
    colour = "#7570b3"
  ) +
  geom_ribbon(
    aes(ymin = lo, ymax = up), fill = "#7570b3", alpha = 0.3
  ) +
  labs(y = "Predicted Survival Rate", x = "Exploitation Rate Index") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.0), expand = c(0, 0)) +
  facet_wrap(~day, ncol = 1) 



## Random Intercepts
# representing among year and among stock variability in fork length,
# lipid content, survival, and date

# sigma_stock 
sigma_stk_mat <- post$sigma_stk 
colnames(sigma_stk_mat) <- c("date", "fl", "lipid", "survival")
sigma_stk <- sigma_stk_mat %>% 
  as.data.frame() %>%
  pivot_longer(
    cols = everything(), names_to = "parameter", values_to = "est"
  ) %>% 
  group_by(parameter) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>% 
  mutate(
    sigma = "stock"
  )

# sigma_year
sigma_yr_mat <- post$sigma_yr 
colnames(sigma_yr_mat) <- c("fl", "lipid", "survival")
sigma_yr <- sigma_yr_mat %>% 
  as.data.frame() %>%
  pivot_longer(
    cols = everything(), names_to = "parameter", values_to = "est"
  ) %>% 
  group_by(parameter) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>% 
  mutate(
    sigma = "year"
  )

sigma_stk_pt <- rbind(sigma_yr, sigma_stk) %>% 
  mutate(
    parameter = factor(parameter, levels = c("date", "fl", "lipid", "survival"))
  ) %>% 
  ggplot() +
  geom_pointrange(aes(x = parameter, y = med, ymin = lo, ymax = up),
                  shape = 21, fill = "#d95f02") +
  facet_wrap(~sigma, ncol = 2) +
  labs(x = "Parameter", y = "Variance Estimate") +
  ggsidekick::theme_sleek() 


# stock and year rho (correlations among random intercepts)
mean_rho_stock <- colMeans(post$Rho_stk)
mean_rho_yr <- colMeans(post$Rho_yr)
colnames(mean_rho_stock) <- rownames(mean_rho_stock) <- c("date", "fl", "lipid",
                                                          "surv")
colnames(mean_rho_yr) <- rownames(mean_rho_yr) <- c("fl", "lipid", "surv")

cor_plot_list <- purrr::map2(
  list(mean_rho_stock, mean_rho_yr),
  c("Stock Intercepts", "Year Intercepts"),
  function (x, y) {
    cor_data <- as.data.frame(as.table(x)) %>% 
      mutate(
        Var1 = fct_relevel(as.factor(Var1), "surv", after = Inf),
        Var2 = fct_relevel(as.factor(Var2), "surv", after = Inf)
      )
    
    # Filter for upper and lower triangles
    cor_data_upper <- cor_data %>% 
      filter(as.numeric(Var1) <= as.numeric(Var2))
    cor_data_lower <- cor_data %>% 
      filter(as.numeric(Var1) >= as.numeric(Var2))
    
    # Plot
    ggplot() +
      geom_tile(
        data = cor_data_upper,
        aes(x = Var1, y = Var2, fill = Freq),
        color = "white"
      ) +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = 0, limit = c(-0.99, .99), space = "Lab",
        name = "Correlation"
      ) +
      geom_text(
        data = cor_data_lower %>% filter(!Var1 == Var2),
        aes(x = Var1, y = Var2, label = round(Freq, 2)),
        color = "black", size = 4
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.title = element_blank()
      ) +
      coord_fixed() +
      labs(title = y)
  }
)


## Submodel effects 

## rho representing correlation between lipid and fork length
# symmetrical matrix so extract single values
rho_hist <- ggplot() +
  geom_histogram(aes(x = post$Rho_sl[ , 2, 1]), bins = 50, fill = "#1b9e77") +
  geom_vline(aes(xintercept = 0), lty = 2 , colour = "black", linewidth = 1.25) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank()
  )

## date effects on lipid content and fork length
samp_date_seq <- seq(from = -2.9, to = 2.7, length.out = 30)

pred_mu_fl <- purrr::map(
  samp_date_seq, ~ post$alpha_f + post$beta_df * .x
) %>% 
  do.call(cbind, .) 
pred_mu_lipid <- purrr::map(
  samp_date_seq, ~ post$alpha_l + post$beta_dl * .x
) %>% 
  do.call(cbind, .) 
colnames(pred_mu_fl) <- colnames(pred_mu_lipid) <- samp_date_seq

# combine posterior preds
pred_mu <- rbind(
  pred_mu_fl %>% 
    as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "yday", values_to = "est"
    ) %>% 
    mutate(
      var = "fl"
    ),
  pred_mu_lipid %>% 
    as.data.frame() %>% 
    pivot_longer(
      cols = everything(), names_to = "yday", values_to = "est"
    ) %>% 
    mutate(
      var = "lipid"
    )
) %>%
  mutate(
    yday_z = as.numeric(yday),
    year_day = (yday_z * sd(det_dat$year_day)) + mean(det_dat$year_day)
  ) %>% 
  group_by(year_day, var) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) 

pred_mu_ribbon <- ggplot(pred_mu, aes(x = year_day, y = med)) +
  geom_line(
    colour = "#1b9e77"
  ) +
  geom_ribbon(
    aes(ymin = lo, ymax = up), fill = "#1b9e77", alpha = 0.3
  ) +
  facet_wrap(~var) +
  scale_x_continuous(breaks = c(135, 182, 226), 
                     labels = c("May 15", "Jul 1", "Aug 15"), 
                     expand = c(0, 0)) +
  labs(y = "Submodel Prediction", x = "Tagging Date") +
  ggsidekick::theme_sleek()



## Survival effects 

# sampling day, including indirect effects on size/lipid
# generate posterior covariance matrix for each draw, combine with pred mu to 
# sample from mvrnorm, then iterate over exp var
sigma <- post$Sigma_fl
rho <- post$Rho_sl
S <- vector(mode = "list", length = nrow(sigma))
sim_surv <- sim_surv_d <- sim_fl <- sim_lipid <- matrix(
  NA, 
  nrow = nrow(sigma), 
  ncol = length(samp_date_seq)
)
for (j in seq_along(samp_date_seq)) {
  for (i in 1:nrow(sigma)) {
    S <- diag(sigma[i,]) %*% rho[i,,] %*% diag(sigma[i,])
    mu <- MASS::mvrnorm(
      n = 1, mu = c(pred_mu_fl[i, j], pred_mu_lipid[i, j]),
      Sigma = S
    )
    sim_fl[i, j] <- mu[1]
    sim_lipid[i, j] <- mu[2]
  }
  # excludes hierarchical intercept; assumes low injury and mean harvest
  sim_eta <- as.numeric(
    post$alpha_s + 
      post$beta_ds * samp_date_seq[j] +
      post$beta_fs * sim_fl[ , j] +
      post$beta_ls * sim_lipid[ , j]
  )
  # as above but sets fl and lipid to zero (i.e. removes them)
  sim_eta_direct <- as.numeric(
    post$alpha_s + 
      post$beta_ds * samp_date_seq[j]
  )
  sim_surv[ , j] <- boot::inv.logit(sim_eta) # Probability of survival
  sim_surv_d[ , j] <- boot::inv.logit(sim_eta_direct) # Probability of survival
}

pred_day_surv_total <- sim_surv %>% 
  as.data.frame() %>% 
  set_names(samp_date_seq) %>% 
  pivot_longer(
    cols = everything(), names_to = "yday_z", values_to = "est"
  ) %>% 
  mutate(
    effect = "total"
  )
pred_day_surv_direct <- sim_surv_d %>% 
  as.data.frame() %>% 
  set_names(samp_date_seq) %>% 
  pivot_longer(
    cols = everything(), names_to = "yday_z", values_to = "est"
  ) %>% 
  mutate(
    effect = "direct"
  )

day_effect_pal <- c(1, 2)
names(day_effect_pal) <- c("direct", "total")

# shared y axis minimum value 
ymin_val <- 0.1

pred_day_ribbon <- rbind(pred_day_surv_total, pred_day_surv_direct) %>% 
  mutate(
    yday_z = as.numeric(yday_z),
    year_day = (yday_z * sd(det_dat$year_day)) + mean(det_dat$year_day)
  ) %>% 
  group_by(year_day, effect) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>% 
  ggplot(
    ., aes(x = year_day, y = med, lty = effect)
  ) +
  geom_line(
    colour = "#7570b3"
  ) +
  geom_ribbon(
    aes(ymin = lo, ymax = up), fill = "#7570b3", alpha = 0.3
  ) +
  scale_linetype_manual(values = day_effect_pal) +
  labs(y = "Predicted Survival Rate", x = "Tagging Date") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(breaks = c(135, 182, 226), 
                     labels = c("May 15", "Jul 1", "Aug 15"), 
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(ymin_val, 1.0), expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
  )


## direct effects of lipid content
lipid_seq <- seq(-3, 4.2, length = 40)
pred_l <- sapply(
  lipid_seq, 
  function (x) {
    inv_logit(post$alpha_s + post$beta_ls * x)
  }
)
pred_lipid_ribbon <- pred_l %>% 
  as.data.frame() %>% 
  set_names(lipid_seq) %>% 
  pivot_longer(
    cols = everything(), names_to = "lipid_z", values_to = "est"
  ) %>% 
  mutate(
    lipid_z = as.numeric(lipid_z),
    lipid = (lipid_z * sd(det_dat$lipid)) + mean(det_dat$lipid)
  ) %>% 
  group_by(lipid) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>% 
  ggplot(
    ., aes(x = lipid, y = med)
  ) +
  geom_line(
    colour = "#7570b3"
  ) +
  geom_ribbon(
    aes(ymin = lo, ymax = up), fill = "#7570b3", alpha = 0.3
  ) +
  labs(y = "Predicted Survival Rate", x = "Lipid Content (% wet weight)") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0), breaks = c(5, 10, 15)) +
  scale_y_continuous(limits = c(ymin_val, 1.0), expand = c(0, 0)) +
  theme(
    axis.title.y = element_blank()
  )


## direct effects of fork length
fl_seq <- seq(-2, 3.4, length = 40)
pred_fl <- sapply(
  fl_seq, 
  function (x) {
    inv_logit(post$alpha_s + post$beta_fs * x)
  }
)
pred_fl_ribbon <- pred_fl %>% 
  as.data.frame() %>% 
  set_names(fl_seq) %>% 
  pivot_longer(
    cols = everything(), names_to = "fl_z", values_to = "est"
  ) %>% 
  mutate(
    fl_z = as.numeric(fl_z),
    fl = (fl_z * sd(det_dat$fl)) + mean(det_dat$fl)
  ) %>% 
  group_by(fl) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>% 
  ggplot(
    ., aes(x = fl, y = med)
  ) +
  geom_line(
    colour = "#7570b3"
  ) +
  geom_ribbon(
    aes(ymin = lo, ymax = up), fill = "#7570b3", alpha = 0.3
  ) +
  labs(y = "Predicted Survival Rate", x = "Fork Length (cm)") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(ymin_val, 1.0), expand = c(0, 0)) +
  theme(
    axis.title.y = element_blank()
  )


pp <- cowplot::plot_grid(
  pred_day_ribbon, pred_lipid_ribbon, pred_fl_ribbon,
  ncol = 1
) 


## Injury effects on survival
# Note scales so that cumsum(delta) * beta = total effect first stage = 0 
# absorbed in intercept
pred_inj <- matrix(
  NA, nrow = nrow(post$delta), ncol = (ncol(post$delta) + 1)
)
for (i in 1:(ncol(post$delta) + 1)) {
  if (i == 1) {
    delta_inj <- 0 
  }
  if (i > 1) {
    delta_inj <- apply(as.matrix(post$delta[ , 1:(i - 1)]), 1, sum) 
  }
  pred_inj[ , i] <- inv_logit(
    post$alpha_s + post$beta_is * delta_inj
  )
}

pred_injury_dat <- pred_inj %>% 
  as.data.frame() %>%
  set_names(
    seq(0, 3, by = 1)
  ) %>% 
  pivot_longer(
    cols = everything(), names_to = "score", values_to = "est"
  )

injury_point <- pred_injury_dat %>% 
  group_by(score) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>% 
  ggplot() +
  geom_pointrange(aes(x = score, y = med, ymin = lo, ymax = up),
                  shape = 21, fill = "#7570b3") +
  labs(x = "Injury Score", y = "Predicted Survival Rate") +
  ggsidekick::theme_sleek() 


diff_dat <- data.frame(
  diff = pred_inj[ , 4] - pred_inj[ , 1]
) 
median(diff_dat$diff)
diff_hist <- ggplot() +
  geom_histogram(data = diff_dat, aes(x = diff), 
                 bins = 50, fill = "#7570b3") +
  geom_vline(aes(xintercept = 0), lty = 2 , colour = "black", linewidth = 1) +
  ggsidekick::theme_sleek() +
  labs(x = "Difference in Survival Rate\nBetween Max/Min Injury Scores") +
  theme(
    axis.title.y = element_blank()
  ) 

n_neg <- diff_dat %>% 
  filter(diff < 0) %>% 
  nrow() 
n_neg / nrow(diff_dat) #80%


# stock, including indirect effects on date/size/lipid, and indirect exploitation
# rate 
# mean date by stock
pred_mu_date <- post$alpha_stk[ , , 1]

stk_seq <- unique(dat_list$stk) %>% sort()

# calculate mean cyer by stk
cyer_seq <- det_dat %>% 
  arrange(stk) %>% 
  group_by(stk) %>% 
  summarize(
    mean_cyer = mean(cyer_z)
  ) %>% 
  pull(mean_cyer)


sigma <- post$Sigma_fl
rho <- post$Rho_sl
S <- vector(mode = "list", length = nrow(sigma))
pred_mu_fl <- pred_mu_lipid <- matrix(
  NA, nrow = nrow(pred_mu_date), ncol = ncol(pred_mu_date)
)
sim_surv <- sim_surv_d <- sim_surv_cyer <- sim_fl <- sim_lipid <- matrix(
  NA, 
  nrow = nrow(sigma), 
  ncol = length(stk_seq)
)
# define low cyer for comparison
low_cyer <- 0
low_cyer_z <- (low_cyer - mean(det_dat$focal_er)) / sd(det_dat$focal_er)


for (j in seq_along(stk_seq)) {
  pred_mu_fl[ , j] <- post$alpha_f + post$beta_df * pred_mu_date[ , j] + 
    post$alpha_stk[ , j, 1]
  pred_mu_lipid[ , j] <- post$alpha_l + post$beta_dl * pred_mu_date[ , j] + 
    post$alpha_stk[ , j, 2]
  
  for (i in 1:nrow(sigma)) {
    S <- diag(sigma[i,]) %*% rho[i,,] %*% diag(sigma[i,])
    mu <- MASS::mvrnorm(
      n = 1, mu = c(pred_mu_fl[i, j], pred_mu_lipid[i, j]),
      Sigma = S
    )
    sim_fl[i, j] <- mu[1]
    sim_lipid[i, j] <- mu[2]
  }
  # excludes hierarchical intercept; assumes low injury and low harvest rate
  sim_eta <- as.numeric(
    post$alpha_s + 
      post$beta_ds * pred_mu_date[ , j] +
      post$beta_fs * sim_fl[ , j] +
      post$beta_ls * sim_lipid[ , j] +
      post$beta_cs * low_cyer_z +
      post$alpha_stk[ , j, 4]
  )
  # as above but sets fl and lipid to zero (i.e. removes them)
  sim_eta_direct <- as.numeric(
    post$alpha_s + post$alpha_stk[ , j, 4] +
      post$beta_cs * low_cyer_z
  )
  # as above but includes stock-specific cyer
  sim_eta_cyer <- as.numeric(
    post$alpha_s + 
      post$beta_ds * pred_mu_date[ , j] +
      post$beta_fs * sim_fl[ , j] +
      post$beta_ls * sim_lipid[ , j] +
      post$alpha_stk[ , j, 4] +
      post$beta_cs * cyer_seq[j] + 
      (post$beta_ds_cs * pred_mu_date[ , j] * cyer_seq[j])
  )
  
  sim_surv[ , j] <- boot::inv.logit(sim_eta) 
  sim_surv_d[ , j] <- boot::inv.logit(sim_eta_direct) 
  sim_surv_cyer[ , j] <- boot::inv.logit(sim_eta_cyer) 
}

stk_key <- det_dat %>% 
  select(stock_group, stk) %>% 
  mutate(stk = as.character(stk)) %>% 
  distinct()

pred_stk_surv_total <- sim_surv %>% 
  as.data.frame() %>% 
  set_names(stk_seq) %>% 
  pivot_longer(
    cols = everything(), names_to = "stk", values_to = "est"
  ) %>% 
  mutate(
    effect = "absent"
  )
# pred_stk_surv_direct <- sim_surv_d %>% 
#   as.data.frame() %>% 
#   set_names(stk_seq) %>% 
#   pivot_longer(
#     cols = everything(), names_to = "stk", values_to = "est"
#   ) %>% 
#   mutate(
#     effect = "direct"
#   )
pred_stk_surv_cyer <- sim_surv_cyer %>% 
  as.data.frame() %>% 
  set_names(stk_seq) %>% 
  pivot_longer(
    cols = everything(), names_to = "stk", values_to = "est"
  ) %>% 
  mutate(
    effect = "present"
  )

alpha_pal <- c("white", #"#7570b3", 
               "black")
names(alpha_pal) <- c("absent", "present")

pred_stk_dat <- rbind(
  #pred_stk_surv_direct, 
  pred_stk_surv_total, pred_stk_surv_cyer
) %>% 
  group_by(stk, effect) %>%
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>%
  left_join(., stk_key, by = "stk") %>% 
  mutate(
    stock_group = fct_recode(stock_group, "Fraser Sum. 0.3" = "Fraser Sum. 4.1",
                             "Fraser Spr. 1.2" = "Fraser Spr. Yr.")
  )

pred_stk_dat %>% 
  group_by(effect) %>% 
  summarize(sd(med) / mean(med))

pred_stk_comb <- ggplot(pred_stk_dat) +
  geom_pointrange(aes(x = stock_group, y = med, ymin = lo, ymax = up, 
                      fill = effect),
                  position = position_dodge(width = 0.4),
                  shape = 21) +
  scale_fill_manual(values = alpha_pal, name = "Exploitation\nRate Index") +
  labs(y = "Predicted Survival Rate") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


## calculate difference in survival when adding in ecological covariates
diff_surv_list <- vector(mode = "list", length = length(stk_seq))
for (i in seq_along(stk_seq)) {
  diff_surv_list[[i]] <- data.frame(
    diff = sim_surv[ , i] - sim_surv_d[ , i],
    stk = as.character(i)
  )
}   

diff_surv_dat <- diff_surv_list %>% 
  bind_rows() %>% 
  left_join(., stk_key, by = "stk") 


stock_pal <- c(#"#b30000",
  "#6baed6", "#08306b", "#fec44f", "#ccece6", "#238b45",
               "#bcbddc", #"#807dba",
  "#54278f", "#3f007d", "#fde0dd", "#f768a1")
names(stock_pal) <- levels(stk_key$stock_group)

diff_survival_hist <- ggplot(data = diff_surv_dat) +
  geom_density(aes(x = diff, fill = stock_group), colour = "black") +
  geom_vline(aes(xintercept = 0), lty = 2 , colour = "black", linewidth = 1) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = stock_pal, guide = "none") +
  labs(x = "Difference in Survival Rate Due to Individual Covariates") +
  theme(
    axis.title.y = element_blank()
  ) +
  facet_wrap(~stock_group)



## export figs
png(here::here("figs", "binomial-glm-cyer-uninformative", "det_prob.png"), 
    units = "in", 
    res = 250, height = 4.5, width = 6.5)
terminal_det_p
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "surv_cyer.png"), 
    units = "in", res = 250, height = 5.25, width = 2.5)
pred_day_cyer
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "ri_sigmas.png"), 
    units = "in", res = 250, height = 3.5, width = 6)
sigma_stk_pt
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "fl_lipid_corr.png"), 
    units = "in", res = 250, height = 3, width = 3.5)
rho_hist
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "stock_ri_corr.png"),
    units = "in", res = 250, height = 5, width = 5.5)
cor_plot_list[[1]]
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "yr_ri_corr.png"), 
    units = "in", res = 250, height = 4.5, width = 5)
cor_plot_list[[2]]
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "fl_lipid_pred.png"),
    units = "in", res = 250, height = 3.5, width = 6)
pred_mu_ribbon
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "surv_pred.png"), 
    units = "in", res = 250, height = 5.25, width = 2.5)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    pp, 
    left = textGrob(
      "Predicted Survival Rate", rot = 90, 
      gp = gpar(fontsize = 11))
  )
)
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "inj_pred.png"), 
    units = "in", res = 250, height = 3, width = 3.5)
injury_point
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "inj_delta.png"),
    units = "in", res = 250, height = 3, width = 3.5)
diff_hist
dev.off()

png(here::here("figs", "binomial-glm-cyer-uninformative", "stock_surv.png"),
    units = "in", res = 250, height = 3.5, width = 6)
pred_stk_comb
dev.off()
png(here::here("figs", "binomial-glm-cyer-uninformative",
               "diff_survival_hist.png"),
    units = "in", res = 250, height = 4.5, width = 6)
diff_survival_hist
dev.off()
