## Logistic survival
# Survival to detection and to terminal arrays using logistic regression (i.e.
# not stage-specific, does not account for detection probability)
# Excludes immature tags and fish with unknown stock ID, but does not exclude
# redeployed tags or injured fish
# Identical to survival_logistic_reg.R but excludes stocks without exploitation
# rate estimates and models interaction between exploitation and capture date
# Feb. 22, 2021
# Updated July 4, 2022

library(tidyverse)
library(rethinking)
library(grid)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



det_dat1 <- readRDS(here::here("data", "surv_log_reg_data.rds")) %>% 
  filter(
    !is.na(isbm_cyer),
    stage_1 == "mature"
  )


# scale continuous covariates
det_dat <- det_dat1 %>% 
  mutate(
    lipid_z = scale(lipid) %>% as.numeric(),
    fl_z = scale(fl) %>% as.numeric(),
    wt_z = scale(wt) %>% as.numeric(),
    log_wt_z = scale(log(wt)) %>% as.numeric(),
    day_z = scale(year_day) %>% as.numeric(),
    cyer_z = scale(isbm_cyer) %>% as.numeric(),
    terminal_p = as.factor(terminal_p),
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

det_dat <- det_dat %>% filter()

dat_list <- list(
  surv = as.integer(det_dat$term_det),
  fl_z = det_dat$fl_z,
  lipid_z = det_dat$lipid_z,
  day_z = det_dat$day_z,
  cyer_z = det_dat$cyer_z,
  inj = det_dat$inj,
  term_p = det_dat$term_p,
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
    mu_fl <- alpha_fl + beta_df * day_z + beta_stk[stk, 1] + beta_yr[yr, 1] ,
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
  data=dat_list, chains = 4, cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)

saveRDS(m4, here::here("data", "model_outputs", "hier_binomial_cyer.rds"))
m4 <- readRDS(here::here("data", "model_outputs", "hier_binomial_cyer.rds"))


# PRIOR PREDICTIONS ------------------------------------------------------------

day_seq <- c(-2, 0, 2)
cyer_seq <- seq(-1.5, 4.5, length = 40)
preds <- vector(mode = "list", length = length(day_seq))
prior_sbar <- rnorm(1000, 0, 0.5)
# prior_term <- rnorm(1000, 0, 0.5)
# prior_ds <- rnorm(1000, 0, 0.5)
prior_cs <- rnorm(1000, -0.5, 0.2)
prior_ds_cs <- rnorm(1000, 0.2, 0.2)


for(i in seq_along(day_seq)) {
  pred_x <- sapply(
    cyer_seq, 
    function (x) {
      inv_logit(
        prior_sbar +
          # prior_term + prior_ds * day_seq[i] +
          prior_cs * x + 
          (prior_ds_cs * day_seq[i] * x)
      )
    }
  )
  preds[[i]] <- pred_x[1:20, ] %>% 
    as.data.frame() %>% 
    set_names(cyer_seq) %>% 
    mutate(iter = seq(1, nrow(.), by = 1)) %>% 
    pivot_longer(
      cols = -iter, names_to = "cyer", values_to = "est"
    ) %>% 
    mutate(
      cyer = as.numeric(cyer),
      day_z = day_seq[i]
    )
}

bind_rows(preds) %>% 
  ggplot(
    ., aes(x = cyer, y = est, group = iter)
  ) +
  geom_line(
  ) +
  ggsidekick::theme_sleek() +
  facet_wrap(~day_z)


# POSTERIOR PREDICTIONS --------------------------------------------------------

post <- extract.samples(m4)

post_pred_foo <- function(dat_in) {
  delta_inj_eff <- cbind(0, post$delta)
  #as matrix necessary for cases where only first column selected
  delta_inj2 <- apply(as.matrix(delta_inj_eff[ , 1:dat_in$inj]), 1, sum)

  logodds <- with(
    post,
    alpha_s[ , dat_in$term_p] +
      beta_yr[ , dat_in$yr, 3] +
      beta_stk[ , dat_in$stk, 3] +
      beta_ds * dat_in$day_z +
      beta_fs * dat_in$fl_z +
      beta_ls * dat_in$lipid_z +
      beta_cs * dat_in$cyer_z + 
      (beta_ds_cs * dat_in$day_z * dat_in$cyer_z) +
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
abline(0, 1)


png(here::here("figs", "binomial-glm-cyer", "qq_plot.png"), units = "in", 
    res = 250, height = 3.5, width = 3.5)
qqplot(qunif(ppoints(length(pit_residuals))), pit_residuals,
       main = "QQ-plot of PIT Residuals")
abline(0, 1)
dev.off()


# POSTERIOR INFERENCE  ---------------------------------------------------------

## Calendar year exploitation rate (all other plots are shared w/ 
# survival_logistic_reg.R)

yday_seq <- c(120, 181, 243) 
day_label <- c("May 1", "July 1", "Sep 1")
day_seq <- (yday_seq - mean(det_dat$year_day)) / sd(det_dat$year_day)
cyer_seq <- seq(-1.5, 4.5, length = 40)
preds <- vector(mode = "list", length = length(day_seq))

for(i in seq_along(day_seq)) {
  pred_x <- sapply(
    cyer_seq, 
    function (x) {
      inv_logit(
        post$alpha_s[ , 2] + post$beta_ds * day_seq[i] +
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
      cyer = (cyer_z * sd(det_dat$isbm_cyer)) + mean(det_dat$isbm_cyer),
      day = day_label[i]
    )
}

pred_day_cyer <- bind_rows(preds) %>% 
  mutate(
    day = factor(day, levels = c("May 1", "July 1", "Sep 1"))
  ) %>% 
  group_by(day, cyer) %>% 
  summarize(
    med = median(est),
    lo = rethinking::HPDI(est, prob = 0.9)[1],
    up = rethinking::HPDI(est, prob = 0.9)[2]
  ) %>% 
  ggplot(
    ., aes(x = cyer, y = med)
  ) +
  geom_line(
    colour = "#7570b3"
  ) +
  geom_ribbon(
    aes(ymin = lo, ymax = up), fill = "#7570b3", alpha = 0.3
  ) +
  labs(y = "Predicted Survival", x = "ISBM CYER") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.0), expand = c(0, 0)) +
  facet_wrap(~day, ncol = 1) 



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


# stock and year rho (correlations among random intercepts)
mean_rho_stock <- colMeans(post$Rho_stk)
mean_rho_yr <- colMeans(post$Rho_yr)
colnames(mean_rho_stock) <- rownames(mean_rho_stock) <- c("fl", "lipid", "surv", "date")
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
  geom_histogram(aes(x = post$Rho[ , 2, 1]), bins = 50, fill = "#1b9e77") +
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
  # injury and mean harvest
  sim_eta <- as.numeric(
    post$alpha_s[ , 2] + 
      post$beta_ds * samp_date_seq[j] +
      post$beta_fs * sim_fl[ , j] +
      post$beta_ls * sim_lipid[ , j]
  )
  # as above but sets fl and lipid to zero (i.e. removes them)
  sim_eta_direct <- as.numeric(
    post$alpha_s[ , 2] + 
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
  labs(y = "Predicted Survival", x = "Tagging Date") +
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
lipid_seq <- seq(-3, 4, length = 40)
pred_l <- sapply(
  lipid_seq, 
  function (x) {
    inv_logit(post$alpha_s[ , 2] + post$beta_ls * x)
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
  scale_y_continuous(limits = c(ymin_val, 1.0), expand = c(0, 0)) +
  theme(
    axis.title.y = element_blank()
  )


## direct effects of fork length
fl_seq <- seq(-2.2, 3.3, length = 40)
pred_fl <- sapply(
  fl_seq, 
  function (x) {
    inv_logit(post$alpha_s[ , 2] + post$beta_fs * x)
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
      post$alpha_s[ , 2] + post$beta_is * delta_inj
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
  labs(x = "Injury Score", y = "Survival Probability") +
  ggsidekick::theme_sleek() 


diff_hist <- data.frame(
  diff = pred_inj[ , 1] - pred_inj[ , 4]
) %>% 
  ggplot() +
  geom_histogram(aes(x = diff), 
                 bins = 50, fill = "#7570b3") +
  geom_vline(aes(xintercept = 0), lty = 2 , colour = "black", linewidth = 1) +
  ggsidekick::theme_sleek() +
  labs(x = "Difference in Survival Prob.\nBetween Max/Min Scores") +
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
    post$alpha_s[ , 2] + 
      post$beta_ds * pred_mu_date[ , j] +
      post$beta_fs * sim_fl[ , j] +
      post$beta_ls * sim_lipid[ , j] +
      post$beta_stk[ , j, 3]
  )
  # as above but sets fl and lipid to zero (i.e. removes them)
  sim_eta_direct <- as.numeric(
    post$alpha_s[ , 2] + post$beta_stk[ , j, 3]
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
png(here::here("figs", "binomial-glm-cyer", "surv_cyer.png"), units = "in", 
    res = 250, height = 6.5, width = 3.5)
pred_day_cyer
dev.off()

png(here::here("figs", "binomial-glm-cyer", "ri_sigmas.png"), units = "in", 
    res = 250, height = 3.5, width = 6)
sigma_stk_pt
dev.off()

png(here::here("figs", "binomial-glm-cyer", "fl_lipid_corr.png"), units = "in", 
    res = 250, height = 3, width = 3.5)
rho_hist
dev.off()

png(here::here("figs", "binomial-glm-cyer", "stock_ri_corr.png"), units = "in", 
    res = 250, height = 5, width = 5.5)
cor_plot_list[[1]]
dev.off()

png(here::here("figs", "binomial-glm-cyer", "yr_ri_corr.png"), units = "in", 
    res = 250, height = 4.5, width = 5)
cor_plot_list[[2]]
dev.off()

png(here::here("figs", "binomial-glm-cyer", "fl_lipid_pred.png"), units = "in", 
    res = 250, height = 3.5, width = 6)
pred_mu_ribbon
dev.off()

png(here::here("figs", "binomial-glm-cyer", "surv_pred.png"), units = "in", 
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

png(here::here("figs", "binomial-glm-cyer", "inj_pred.png"), units = "in", 
    res = 250, height = 3, width = 3.5)
injury_point
dev.off()

png(here::here("figs", "binomial-glm-cyer", "inj_delta.png"), units = "in", 
    res = 250, height = 3, width = 3.5)
diff_hist
dev.off()

png(here::here("figs", "binomial-glm-cyer", "stock_surv.png"), units = "in", 
    res = 250, height = 3.5, width = 6)
pred_stk_comb
dev.off()

