## Logistic survival
# Survival to detection and to terminal arrays using logistic regression (i.e.
# not stage-specific, does not account for detection probability)
# Excludes immature tags and fish with unknown stock ID, but does not exclude
# redeployed tags or injured fish
# Feb. 22, 2021
# Updated July 4, 2022

library(tidyverse)
library(rethinking)
library(grid)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



det_dat1 <- readRDS(here::here("data", "surv_log_reg_data.rds"))


# scale continuous covariates
det_dat <- det_dat1 %>% 
  # filter(
  #   !redeploy == "yes"
  # ) %>% 
  mutate(
    lipid_z = scale(lipid) %>% as.numeric(),
    fl_z = scale(fl) %>% as.numeric(),
    day_z = scale(year_day) %>% as.numeric(),
    cyer_z = scale(isbm_cyer) %>% as.numeric(),
    terminal_p = as.factor(terminal_p),
    year = as.factor(year),
    stock_group = factor(
      stock_group, 
      levels = c("Cali", "Low Col.", "Up Col.", "WA_OR", "WCVI", "ECVI", 
                 "Fraser Sum. 4.1", "Fraser Fall", "Fraser Spr. Yr.", 
                 "Fraser Sum. Yr.", "North Puget", "South Puget")),
    inj = as.integer(as.factor(adj_inj)),
    term_p = as.integer(terminal_p),
    yr = as.integer(year),
    stk = as.integer(stock_group)
  )


# FIT LOGISTIC REG MODELS ------------------------------------------------------

dat_list <- list(
  surv = as.integer(det_dat$term_det),
  fl_z = det_dat$fl_z,
  lipid_z = det_dat$lipid_z,
  day_z = det_dat$day_z,
  # cyer_z = det_dat$cyer_z,
  inj = det_dat$inj,
  term_p = det_dat$terminal_p,
  yr = det_dat$yr,
  stk = det_dat$stk,
  alpha = rep(2, length(unique(det_dat$inj)) - 1)
)


m4 <- ulam(
  alist(
    # stock-specific sampling dates
    day_z ~ normal(mu_day, sigma_day),
    mu_day <- beta_stk[stk, 4],
    # unobserved latent condition
    c(fl_z, lipid_z) ~ multi_normal(c(mu_fl, mu_lipid), Rho, Sigma),
    mu_fl <- alpha_fl + beta_df * day_z + beta_yr[yr, 1] + beta_stk[stk, 1],
    mu_lipid <- alpha_lipid + beta_dl * day_z + beta_yr[yr, 2] + 
      beta_stk[stk, 2],
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_term[term_p] +
      beta_yr[yr, 3] +
      beta_stk[stk, 3] +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z +
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
    c(beta_df, beta_dl) ~ normal(0.1, 0.5),
    Rho ~ lkj_corr(2),
    Sigma ~ exponential(1),
    sigma_day ~ exponential(1),
    
    cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[3]:sigma_yr ~ exponential(1),
    cholesky_factor_corr[4]:L_Rho_stk ~ lkj_corr_cholesky(2),
    vector[4]:sigma_stk ~ exponential(1),
    
    surv_bar ~ normal(0, 1.25),
    beta_term[term_p] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls, beta_is, beta_ss) ~ normal(0, 0.5),
    
    # constraints on ordinal effects of injury
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr),
    gq> matrix[4, 4]:Rho_stk <<- Chol_to_Corr(L_Rho_stk)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)
# as in simpler model, injury effects are very modest (overall effect broadly
# spans zero and scales semi-linearly with injury)

saveRDS(m4, here::here("data", "model_outputs", "hier_binomial.rds"))
m4 <- readRDS(here::here("data", "model_outputs", "hier_binomial.rds"))

post <- extract.samples(m4)


# POSTERIOR PREDICTIONS --------------------------------------------------------


post_pred_foo <- function(dat_in) {
  delta_inj_eff <- cbind(0, post$delta)
  #as matrix necessary for cases where only first column selected
  delta_inj2 <- apply(as.matrix(delta_inj_eff[ , 1:dat_in$inj]), 1, sum)

  logodds <- with(
    post,
    surv_bar + 
      beta_term[ , dat_in$term_p] +
      beta_yr[ , dat_in$yr, 3] +
      beta_stk[ , dat_in$stk, 3] +
      beta_ds * dat_in$day_z +
      beta_fs * dat_in$fl_z +
      beta_ls * dat_in$lipid_z +
      beta_is * delta_inj2
  )
  return( inv_logit(logodds) )
}


sim_mat <- matrix(NA, nrow = nrow(det_dat), ncol = dim(post[[1]])[[1]])
for (i in 1:nrow(det_dat)) {
  sim_mat[i , ] <- post_pred_foo(det_dat[i, ]) 
}
sim_mat_binom <- apply(sim_mat, c(1, 2), function (p) rbinom(1, 1, p))


# Define a function to calculate PIT residuals
calc_pit <- function(y, posterior_pred) {
  # Get the proportion of posterior samples that are less than or equal to the observed value
  n_obs <- length(y)
  pit_residuals <- numeric(n_obs)
  
  for (i in 1:n_obs) {
    # Calculate pmin and pmax for each observation
    y_prime <- posterior_pred[i, ]  # posterior predictions for i-th observation
    pmin_i <- mean(y_prime < y[i])
    pmax_i <- mean(y_prime <= y[i])
    
    # Generate the PIT residuals as a random draw from uniform(pmin, pmax)
    pit_residuals[i] <- runif(1, pmin_i, pmax_i)
  }
  
  return(pit_residuals)
}

# Calculate PIT residuals
pit_residuals <- calc_pit(y = det_dat$term_det, posterior_pred = sim_mat_binom)
qqplot(qunif(ppoints(length(pit_residuals))), pit_residuals,
       main = "QQ-plot of PIT Residuals")



# POSTERIOR INFERENCE  ---------------------------------------------------------

## Random Intercepts
# representing among year and among stock variability in fork length,
# lipid content, survival, and date

# sigma_stock 
sigma_stk_mat <- post$sigma_stk 
colnames(sigma_stk_mat) <- c("fl", "lipid", "survival", "date")
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
  labs(x = "Parameter", y = "Posterior Variance Estimate") +
  ggsidekick::theme_sleek() 



## Submodel effects 

## rho representing correlation between lipid and fork length
# symmetrical matrix so extract single values
rho <- post$Rho[ , 2, 1]
rho_hist <- ggplot() +
  geom_histogram(aes(x = rho), bins = 50, fill = "#1b9e77") +
  geom_vline(aes(xintercept = 0), lty = 2 , colour = "black", linewidth = 1.25) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank()
  )

## date effects on lipid content and fork length
samp_date_seq <- seq(from = -2.9, to = 2.7, length.out = 30)

pred_mu_fl <- purrr::map(
  samp_date_seq, ~ post$alpha_fl + post$beta_df * .x
) %>% 
  do.call(cbind, .) 
pred_mu_lipid <- purrr::map(
  samp_date_seq, ~ post$alpha_lipid + post$beta_dl * .x
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

# TODO: consider back converting pred_mu into real space
pred_mu_ribbon <- ggplot(pred_mu, aes(x = year_day, y = med)) +
  geom_line(
    colour = "#1b9e77"
  ) +
  geom_ribbon(
    aes(ymin = lo, ymax = up), fill = "#1b9e77", alpha = 0.3
  ) +
  facet_wrap(~var) +
  labs(y = "Submodel Prediction", x = "Year Day") +
  ggsidekick::theme_sleek()



## Survival effects 

# sampling day, including indirect effects on size/lipid
# generate posterior covariance matrix for each draw, combine with pred mu to 
# sample from mvrnorm, then iterate over exp var
sigma <- post$Sigma
rho <- post$Rho
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
  # excludes hierarchical intercept; assumes high terminal det prob and low
  # injury
  sim_eta <- as.numeric(
    post$surv_bar + 
      post$beta_ds * samp_date_seq[j] +
      post$beta_fs * sim_fl[ , j] +
      post$beta_ls * sim_lipid[ , j] +
      post$beta_term[ , 2]
  )
  # as above but sets fl and lipid to zero (i.e. removes them)
  sim_eta_direct <- as.numeric(
    post$surv_bar + 
      post$beta_ds * samp_date_seq[j] +
      post$beta_term[ , 2]
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
  labs(y = "Predicted Survival", x = "Year Day") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.2, 1.0), expand = c(0, 0)) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank()
  )


## direct effects of lipid content
lipid_seq <- seq(-3, 4, length = 40)
pred_l <- sapply(
  lipid_seq, 
  function (x) {
    inv_logit(post$surv_bar + post$beta_term[ , 2] + post$beta_ls * x)
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
  labs(y = "Predicted Survival", x = "Lipid Content") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.2, 1.0), expand = c(0, 0)) +
  theme(
    axis.title.y = element_blank()
  )


## direct effects of fork length
fl_seq <- seq(-2.2, 3.3, length = 40)
pred_fl <- sapply(
  fl_seq, 
  function (x) {
    inv_logit(post$surv_bar + post$beta_term[ , 2] + post$beta_fs * x)
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
  labs(y = "Predicted Survival", x = "Fork Length") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.2, 1.0), expand = c(0, 0)) +
  theme(
    axis.title.y = element_blank()
  )


pp <- cowplot::plot_grid(
  pred_day_ribbon, pred_lipid_ribbon, pred_fl_ribbon,
  ncol = 1
) 


## Injury and scale loss effects on survival
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
      post$surv_bar + post$beta_is * delta_inj + post$beta_term[ , 2]
      )
  }

# combine data
# pred_capture_dat <- purrr::map2(
#   list(pred_inj), c("injury"),
#   function (x, y) {
#     x %>% 
#       as.data.frame() %>% 
#       set_names(
#         seq(0, 3, by = 1)
#       ) %>% 
#       pivot_longer(
#         cols = everything(), names_to = "score", values_to = "est"
#       ) %>% 
#       mutate(
#         metric = y
#       )
#   }
# )  %>% 
#   bind_rows()
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
  labs(x = "Injury Score", y = "Survival Probability") +
  ggsidekick::theme_sleek() 


#calculate difference in survival between lowest and highest injury score
# diff_pt <- data.frame(
#   diff_inj = pred_inj[ , 1] - pred_inj[ , 4],
#   diff_scale = pred_scale[ , 1] - pred_scale[ , 4]
# ) %>% 
#   pivot_longer(
#     cols = everything(), names_to = "metric", values_to = "diff", 
#     names_prefix = "diff_"
#   ) %>% 
  # group_by(metric) %>%
  # summarize(
  #   med = median(diff),
  #   lo = rethinking::HPDI(diff, prob = 0.9)[1],
  #   up = rethinking::HPDI(diff, prob = 0.9)[2]
  # ) %>%
  # ggplot() +
  # geom_pointrange(aes(x = metric, y = med, ymin = lo, ymax = up),
  #                 shape = 21, fill = "#7570b3") +
  # labs(x = "Metric",
  #      y = "Difference in Survival Probability Between Max/Min Scores") +
  # ggsidekick::theme_sleek() +
  # geom_hline(aes(yintercept = 0), lty = 2) +
  # theme(
  #   axis.title.y = element_blank()
  # )

diff_hist <- data.frame(
  diff = pred_inj[ , 1] - pred_inj[ , 4]
) %>% 
  ggplot() +
  geom_histogram(aes(x = diff), 
                 bins = 50, fill = "#7570b3") +
  geom_vline(aes(xintercept = 0), lty = 2 , colour = "black", linewidth = 1) +
  ggsidekick::theme_sleek() +
  labs(x = "Difference in Survival Probability Between Max/Min Scores") +
  theme(
    axis.title.y = element_blank()
  ) 



# stock, including indirect effects on date/size/lipid

# mean date by stock
pred_mu_date <- post$beta_stk[ , , 4]

stk_seq <- unique(dat_list$stk) %>% sort()

sigma <- post$Sigma
rho <- post$Rho
S <- vector(mode = "list", length = nrow(sigma))
pred_mu_fl <- pred_mu_lipid <- matrix(
  NA, nrow = nrow(pred_mu_date), ncol = ncol(pred_mu_date)
)
sim_surv <- sim_surv_d <- sim_fl <- sim_lipid <- matrix(
  NA, 
  nrow = nrow(sigma), 
  ncol = length(stk_seq)
)
for (j in seq_along(stk_seq)) {
  pred_mu_fl[ , j] <- post$alpha_fl + post$beta_df * pred_mu_date[ , j] + 
    post$beta_stk[ , j, 1]
  pred_mu_lipid[ , j] <- post$alpha_lipid + post$beta_dl * pred_mu_date[ , j] + 
    post$beta_stk[ , j, 2]
  
  for (i in 1:nrow(sigma)) {
    S <- diag(sigma[i,]) %*% rho[i,,] %*% diag(sigma[i,])
    mu <- MASS::mvrnorm(
      n = 1, mu = c(pred_mu_fl[i, j], pred_mu_lipid[i, j]),
      Sigma = S
    )
    sim_fl[i, j] <- mu[1]
    sim_lipid[i, j] <- mu[2]
  }
  # excludes hierarchical intercept; assumes high terminal det prob and low
  # injury
  sim_eta <- as.numeric(
    post$surv_bar + 
      post$beta_ds * pred_mu_date[ , j] +
      post$beta_fs * sim_fl[ , j] +
      post$beta_ls * sim_lipid[ , j] +
      post$beta_term[ , 2] +
      post$beta_stk[ , j, 3]
  )
  # as above but sets fl and lipid to zero (i.e. removes them)
  sim_eta_direct <- as.numeric(
    post$surv_bar + 
      post$beta_stk[ , j, 3] +
      post$beta_term[ , 2]
  )
  sim_surv[ , j] <- boot::inv.logit(sim_eta) # Probability of survival
  sim_surv_d[ , j] <- boot::inv.logit(sim_eta_direct) # Probability of survival
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
    effect = "total"
  )
pred_stk_surv_direct <- sim_surv_d %>% 
  as.data.frame() %>% 
  set_names(stk_seq) %>% 
  pivot_longer(
    cols = everything(), names_to = "stk", values_to = "est"
  ) %>% 
  mutate(
    effect = "direct"
  )

alpha_pal <- c(0.2, 1)
names(alpha_pal) <- c("direct", "total")

pred_stk_comb <- rbind(pred_stk_surv_total, pred_stk_surv_direct) %>% 
  group_by(stk, effect) %>%
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>%
  left_join(., stk_key, by = "stk") %>% 
  ggplot() +
  geom_pointrange(aes(x = stock_group, y = med, ymin = lo, ymax = up, 
                      alpha = effect),
                  position = position_dodge(width = 0.4),
                  shape = 21, fill = "#7570b3") +
  scale_alpha_manual(values = alpha_pal) +
  labs(y = "Predicted Survival Probability") +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


## export figs
png(here::here("figs", "binomial-glm", "ri_sigmas.png"), units = "in", 
    res = 250, height = 3.5, width = 6)
sigma_stk_pt
dev.off()

png(here::here("figs", "binomial-glm", "fl_lipid_corr.png"), units = "in", 
    res = 250, height = 3, width = 3.5)
rho_hist
dev.off()

png(here::here("figs", "binomial-glm", "fl_lipid_pred.png"), units = "in", 
    res = 250, height = 3.5, width = 6)
pred_mu_ribbon
dev.off()

png(here::here("figs", "binomial-glm", "surv_pred.png"), units = "in", 
    res = 250, height = 6.5, width = 3.5)
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    pp, 
    left = textGrob(
      "Survival Probability", rot = 90, 
      gp = gpar(fontsize = 11))
  )
)
dev.off()

png(here::here("figs", "binomial-glm", "inj_pred.png"), units = "in", 
    res = 250, height = 3, width = 3.5)
injury_point
dev.off()

png(here::here("figs", "binomial-glm", "inj_delta.png"), units = "in", 
    res = 250, height = 3, width = 3.5)
diff_hist
dev.off()

png(here::here("figs", "binomial-glm", "stock_surv.png"), units = "in", 
    res = 250, height = 3.5, width = 6)
pred_stk_comb
dev.off()


## SIMPLIFIED MODEL VERSIONS ---------------------------------------------------

# multiple regression hierarchical model with injury effects 
m1 <- ulam(
  alist(
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_term[term_p] +
      beta_yr[yr] * sigma_yr +
      beta_stk[stk] * sigma_stk +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z +
      beta_is * sum(delta_inj[1:inj])
    ,
    surv_bar ~ normal(0, 1.25),
    beta_yr[yr] ~ normal(0, 0.5),
    beta_stk[stk] ~ normal(0, 0.5),
    beta_term[term_p] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls, beta_is) ~ normal(0, 0.5),
    c(sigma_yr, sigma_stk) ~ exponential(1),
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)

# check priors (look good, note following code doesn't work with dirichlet)
# set.seed(1999)
# prior <- extract.prior( m1 , n=1e3 )
# prior_p_dat <- data.frame(
#   fl_z = seq(-2, 2, length = 30),
#   lipid_z = 0,
#   day_z = 0,
#   term_p = 1,
#   yr = 1,
#   stk = 1
# )
# prior_p <- link(m1, post = prior, data = prior_p_dat)
# hist(prior_p)
# plot(prior_p[1, ] ~ prior_p_dat$fl_z, type = "line", ylim = c(0, 1))
# for (i in 2:50) {
#   lines(prior_p_dat$fl_z, prior_p[i, ])
# }


# as model 1 but assumes date influences size and lipid content, which also 
# covary with one another
m2 <- ulam(
  alist(
    # unobserved latent condition
    c(fl_z, lipid_z) ~ multi_normal(c(mu_fl, mu_lipid), Rho, Sigma),
    mu_fl <- alpha_fl + beta_df * day_z,
    mu_lipid <- alpha_lipid + beta_dl * day_z,
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_term[term_p] +
      beta_yr[yr] * sigma_yr +
      beta_stk[stk] * sigma_stk +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z
    ,
    c(alpha_fl, alpha_lipid) ~ normal(0, 0.5),
    c(beta_df, beta_dl) ~ normal(0.1, 0.5),
    Rho ~ lkj_corr(2),
    Sigma ~ exponential(1),
    surv_bar ~ normal(0, 1.25),
    beta_term[term_p] ~ normal(0, 0.5),
    beta_yr[yr] ~ normal(0, 0.5),
    beta_stk[stk] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls) ~ normal(0, 0.5),
    c(sigma_yr, sigma_stk) ~ exponential(1)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)
# modest covariance among fork length and lipid (~0.3); similar sigmas (0.9);
# moderate day of year effects (0.3); other values unchanged


# as model 2 but assume a) year- and stock-specific intercepts on size and 
# lipid and b) stock-specific effects on tagging date;
# assume year and stock effects in each submodel are drawn from multivariate
m3 <- ulam(
  alist(
    # stock-specific sampling dates
    day_z ~ normal(mu_day, sigma_day),
    mu_day <- beta_stk[stk, 4],
    # unobserved latent condition
    c(fl_z, lipid_z) ~ multi_normal(c(mu_fl, mu_lipid), Rho, Sigma),
    mu_fl <- alpha_fl + beta_df * day_z + beta_yr[yr, 1] + beta_stk[stk, 1],
    mu_lipid <- alpha_lipid + beta_dl * day_z + beta_yr[yr, 2] + 
      beta_stk[stk, 2],
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_term[term_p] +
      beta_yr[yr, 3] +
      beta_stk[stk, 3] +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z,
    # adaptive priors
    transpars> matrix[yr, 3]:beta_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[3, yr]:z_yr ~ normal(0, 0.5),
    transpars> matrix[stk, 4]:beta_stk <-
      compose_noncentered(sigma_stk , L_Rho_stk , z_stk),
    matrix[4, stk]:z_stk ~ normal(0, 0.5),
    # fixed priors
    c(alpha_fl, alpha_lipid) ~ normal(0, 0.5),
    c(beta_df, beta_dl) ~ normal(0.1, 0.5),
    Rho ~ lkj_corr(2),
    Sigma ~ exponential(1),
    sigma_day ~ exponential(1),
    
    cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[3]:sigma_yr ~ exponential(1),
    cholesky_factor_corr[4]:L_Rho_stk ~ lkj_corr_cholesky(2),
    vector[4]:sigma_stk ~ exponential(1),
    
    surv_bar ~ normal(0, 1.25),
    beta_term[term_p] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls) ~ normal(0, 0.5),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr),
    gq> matrix[4, 4]:Rho_stk <<- Chol_to_Corr(L_Rho_stk)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)
# covariance between year effects on survival and lipid/fl much weaker (~0.05)
# covariance between stock effects on survival and lipid/fl much weaker (~0.05)
# greater interstock and interannual variability in lipid content than fork
# length or survival



## FREQUENTIST MODELS ----------------------------------------------------------

dum_fit <- glmmTMB(
  lipid ~ day_z + (1 | stock_group) + (1 | year), 
  data = det_dat
)


### EVALUATE HOW TO DEAL WITH INJURIES -----------------------------------------

# injury vs redeployed tag effects

det_dat_no_redep <- det_dat %>% 
  filter(redeploy == "no")
det_dat_no_3 <- det_dat %>% 
  filter(!injury == "3")

# assume year and stock group are random intercepts
# compare three models for terminal detections: injury (0-3) & redeploy tags
dat_list <- list(
  surv = det_dat$term_det,
  inj = det_dat$injury,
  redeploy = ifelse(det_dat$redeploy == "no", 1, 2),
  yr = as.integer(as.factor(det_dat$year)),
  stock = as.integer(as.factor(det_dat$stock_group)),
  alpha = rep(2, length(unique(det_dat$injury)) - 1)
)
dat_list2 <- list(
  surv = det_dat_no_redep$term_det,
  inj = det_dat_no_redep$injury,
  redeploy = ifelse(det_dat_no_redep$redeploy == "no", 1, 2),
  yr = as.integer(as.factor(det_dat_no_redep$year)),
  stock = as.integer(as.factor(det_dat_no_redep$stock_group)),
  alpha = rep(2, length(unique(det_dat_no_redep$injury)) - 1)
)
dat_list3 <- list(
  surv = det_dat_no_3$term_det,
  inj = det_dat_no_3$injury,
  redeploy = ifelse(det_dat_no_3$redeploy == "no", 1, 2),
  yr = as.integer(as.factor(det_dat_no_3$year)),
  stock = as.integer(as.factor(det_dat_no_3$stock_group)),
  alpha = rep(2, length(unique(det_dat_no_3$injury)) - 1)
)


# null model
m0 <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_yr[yr]*sigma_yr +
      gamma_stock[stock] * sigma_stock,
    # priors
    gamma ~ normal(0, 2.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)


m1 <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_yr[yr] * sigma_yr +
      gamma_stock[stock] * sigma_stock +
      gamma_inj * sum(delta_inj[1:inj]),
    # priors
    gamma ~ normal(0, 2.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    gamma_inj ~ normal(0, 0.5),
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list, chains=4 , cores = 4,#log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)
m1b <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_yr[yr] * sigma_yr +
      gamma_stock[stock] * sigma_stock +
      gamma_inj * sum(delta_inj[1:inj]),
    # priors
    gamma ~ normal(0, 2.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    gamma_inj ~ normal(0, 0.5),
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list2, chains=4 , cores = 4,#log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)
m1c <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_yr[yr] * sigma_yr +
      gamma_stock[stock] * sigma_stock +
      gamma_inj * sum(delta_inj[1:inj]),
    # priors
    gamma ~ normal(0, 2.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    gamma_inj ~ normal(0, 0.5),
    vector[3]: delta_inj <<- append_row(0, delta),
    simplex[2]: delta ~ dirichlet(alpha),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list3, chains=4 , cores = 4,#log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)

post <- extract.samples(m1b)
p_link_inj <- function(injury) {
  if (injury == 0) {
    logodds <- with(post, gamma)
  } else (
    logodds <- with(
      post, gamma + (gamma_inj * rowSums(delta[, 1:injury, drop = FALSE])) 
    )
  )
  return( inv_logit(logodds) )
}
p_raw <- sapply(0:3 , function(i) p_link_inj( i ) )
inj_dat <- purrr::map(
  seq(0, 3, by = 1), function (i) {
    data.frame(
      injury = as.factor(i),
      est = p_raw[ , i + 1])
  }
) %>% 
  bind_rows()
ggplot(inj_dat, aes(x = injury, y = est)) +
  geom_boxplot()


# redeploy model
m2 <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_redeploy[redeploy] +
      gamma_yr[yr]*sigma_yr +
      gamma_stock[stock] * sigma_stock,
    # priors
    gamma ~ normal(0, 2.5),
    gamma_redeploy[redeploy] ~ dnorm(0, 0.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)
m2b <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_redeploy[redeploy] +
      gamma_yr[yr]*sigma_yr +
      gamma_stock[stock] * sigma_stock,
    # priors
    gamma ~ normal(0, 2.5),
    gamma_redeploy[redeploy] ~ dnorm(0, 0.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list3, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)

post2 <- extract.samples(m2b)
p_link_redeploy <- function(redeploy) {
  logodds <- with(
      post2, gamma + gamma_redeploy[ , redeploy] 
    )
  return( inv_logit(logodds) )
}
p_raw <- sapply(1:2, function(i) p_link_redeploy( i ) )
redep_dat <- purrr::map(
  c(1, 2), function (i) {
    data.frame(
      redeploy = as.factor(i),
      est = p_raw[ , i])
  }
) %>% 
  bind_rows()
ggplot(redep_dat, aes(x = redeploy, y = est)) +
  geom_boxplot()

# greatest stabilizing effect from removing redeployed tags
