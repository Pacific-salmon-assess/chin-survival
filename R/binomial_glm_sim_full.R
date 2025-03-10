## Data Structure Simulation
# Adapts binomial_glm_sim.R to account for complexities in data structure
# Use a mix of glmmTMB and rethinking models to evaluate how model accuracy
# degrades across increasingly realistic data structures

library(tidyverse)
library(glmmTMB)
library(rethinking)


rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# Set seed for reproducibility
set.seed(123)

n <- 1000  # Number of observations
n_years <- 5
n_pops <- 5
yrs <- seq(1, 5, by = 1)
pops <- c("a", "b", "c", "d", "e")


# CYER from MVN to represent covariance among stocks with different means
cyer_yr_rho <- 0.5 
cyer_mu_yr <- runif(n_pops, 0.05, 0.6) #range of mean stock ERs
cyer_yr_sigmas <- rep(0.3, n_pops)
cyer_Rho <- diag(1, 5) + matrix(cyer_yr_rho, 5, 5) - diag(cyer_yr_rho, 5)
cyer_Sigma <- diag(cyer_yr_sigmas) %*% cyer_Rho %*% diag(cyer_yr_sigmas) 
cyer_mat <- MASS::mvrnorm(n_years, log(cyer_mu_yr), cyer_Sigma) %>% exp()
colnames(cyer_mat) <- pops
cyer_dat <- as.data.frame(cyer_mat) %>%
  mutate(yr = yrs) %>% 
  pivot_longer(cols = a:e, names_to = "pop", values_to = "cyer") 


# detection probability varies among years and pops
yr_pop_key <- expand.grid(
  yr = seq(1, 5, by = 1),
  pop = c("a", "b", "c", "d", "e")
) %>% 
  mutate(
    pop_n = as.numeric(as.factor(pop)),
    det = case_when(
      pop_n < 3 ~ 0,
      pop_n < 3 & yr %in% c(1, 5) ~ 1,
      TRUE ~ 1
    )
  ) %>% 
  left_join(., cyer_dat, by = c("pop", "yr"))


## yr intercepts (condition, survival)
# correlated versions
yr_rho <- 0.9
yr_mu <- c(0, 0) #average effects both centered on zero
yr_sigmas <- c(0.4, 0.2)
yr_Rho <- matrix(c(1, yr_rho, yr_rho, 1), nrow = 2)
yr_Sigma <- diag(yr_sigmas) %*% yr_Rho %*% diag(yr_sigmas)
yr_ints <- MASS::mvrnorm(n_years, yr_mu, yr_Sigma)
alpha_yr_cond <- yr_ints[ , 1]
alpha_yr_surv <- yr_ints[ , 2]
# uncorrelated version
alpha_yr_cond2 <- rnorm(5, 0, 0.4)
alpha_yr_surv2 <- rnorm(5, 0, 0.2)


# pop intercepts (condition, date, survival)
# correlated versions
pop_rho <- 0.9
pop_mu <- c(0, 0, 0) #average effects both centered on zero
pop_sigmas <- c(0.4, 0.5, 0.5)
pop_Rho <- diag(1, 3) + matrix(pop_rho, 3, 3) - diag(pop_rho, 3)
pop_Sigma <- diag(pop_sigmas) %*% pop_Rho %*% diag(pop_sigmas)
pop_ints <- MASS::mvrnorm(n_pops, pop_mu, pop_Sigma)
alpha_pop_cond <- pop_ints[ , 1]
alpha_pop_date <- pop_ints[ , 2]
alpha_pop_surv <- pop_ints[ , 3]
# uncorrelated version
alpha_pop_cond2 <- rnorm(5, 0, 0.4)
alpha_pop_date2 <- rnorm(5, 0, 0.5)
alpha_pop_surv2 <- rnorm(5, 0, 0.5)



## Make dataframe with full structure of correlations
dat <- data.frame(
  yr = sample(yrs, n, replace = TRUE), 
  pop = sample(pops, n, replace =  TRUE)
) %>% 
  left_join(
    ., yr_pop_key, by = c("yr", "pop")
  ) %>% 
  mutate(
    #scale ER data
    cyer = scale(cyer) %>% as.numeric(),
    samp_date = NA,
    cond = NA,
    size = NA,
    lipid = NA,
    yr_f = as.factor(yr),
    pop_f = as.factor(pop)
  ) 

# simulate continuous covariates
beta_dc <- 0.7
beta_cl <- 0.9 #results in a mean correlation of ~0.6
beta_cf <- 0.3

for (i in 1:nrow(dat)) { 
  dat$samp_date[i] <- rnorm(1, alpha_pop_date[dat$pop_n[i]], 1)
  dat$cond[i] <- rnorm(
    1, 
    alpha_pop_cond[dat$pop_n[i]] + alpha_yr_cond[dat$yr[i]] + 
      beta_dc * dat$samp_date[i], 
    1
  )
  # generate observations from latent state 
  dat$size[i] <- rnorm(1, beta_cf * dat$cond[i], 0.5)
  dat$lipid[i] <- rnorm(1, beta_cl * dat$cond[i], 0.5)
  dat$alpha_yr[i] <- alpha_yr_surv[dat$yr[i]]
  dat$alpha_pop[i] <- alpha_pop_surv[dat$pop_n[i]]
}


## make opposite, simplified dataframe where no correlation structure
dat_simple <- dat %>% 
  mutate(
    det = sample(c(0, 1), n, replace =  TRUE),
    cyer = runif(n, 0.2, 0.8) %>% 
      scale() %>% 
      as.numeric()
  )
dat_simple$samp_date <- rnorm(n, 0, 1)
dat_simple$cond <- rnorm(n, 0, 1)
# generate observations from latent state 
dat_simple$size <- rnorm(n, beta_cf * dat_simple$cond, 0.5)
dat_simple$lipid <- rnorm(n, beta_cl * dat_simple$cond, 0.5)
for (i in 1:nrow(dat_simple)) { 
  dat_simple$alpha_yr[i] <- alpha_yr_surv2[dat_simple$yr[i]]
  dat_simple$alpha_pop[i] <- alpha_pop_surv2[dat_simple$pop_n[i]]
}


## fixed effect parameters
surv_bar <- -1  # intercept
beta_cs <- 0.5 # condition effect
beta_ds <- 0.4 # date effect
beta_cyer <- -1 # er effect
beta_d_cyer <- 0.5 # date ER interaction
beta_ps <- 1.5 # detection probability effect

true_pars <- data.frame(
  par = c("(Intercept)", "samp_date", "cyer", "cond", "size", "lipid", "det",
          "samp_date:cyer"),
  true = c(surv_bar, beta_ds, beta_cyer, beta_cs, 0.3 * beta_cs, 0.9 * beta_cs,
           beta_ps, beta_d_cyer)
)


## FREQUENTIST -----------------------------------------------------------------

## DATASET 1
# uncorrelated variables, only fixed effects, no latent variable
eta <- surv_bar + 
  beta_cs * dat_simple$cond +
  beta_ds * dat_simple$samp_date +
  beta_cyer * dat_simple$cyer +
  (beta_d_cyer * dat_simple$samp_date * dat_simple$cyer) +
  beta_ps * dat_simple$det

p <- 1 / (1 + exp(-(eta)))  # Probability of survival
y <- rbinom(n, 1, p)  # Binary outcome

d1 <- dat_simple %>% 
  mutate(surv = y)

fit1 <- glmmTMB(
  surv ~ samp_date * cyer + cond + det,  
  data = d1,
  family = binomial()
)


## DATASET 2 
# as above but replace condition (latent variable) with correlated observed 
# variables
fit2 <- glmmTMB(
  surv ~ samp_date * cyer + size + lipid + det,  
  data = d1,
  family = binomial()
)


## DATASET 3
# as above but with RIs that are uncorrelated with another, detection 
# probability and CYER
eta3 <- surv_bar + 
  dat_simple$alpha_yr + dat_simple$alpha_pop +
  beta_cs * dat_simple$cond +
  beta_ds * dat_simple$samp_date +
  beta_cyer * dat_simple$cyer +
  (beta_d_cyer * dat_simple$samp_date * dat_simple$cyer) +
  beta_ps * dat_simple$det

p3 <- 1 / (1 + exp(-(eta3)))  # Probability of survival
y3 <- rbinom(n, 1, p3)  # Binary outcome

d3 <- dat_simple %>% 
  mutate(surv = y3)

fit3 <- glmmTMB(
  surv ~ samp_date * cyer + size + lipid + det + (1 | yr_f) + (1 | pop_f),  
  data = d3,
  family = binomial()
)


## DATASET 4
# as above but with RIs that covary with another but are independent of 
# detection probability / CYER
eta4 <- surv_bar + 
  dat$alpha_yr + dat$alpha_pop +
  beta_cs * dat_simple$cond +
  beta_ds * dat_simple$samp_date +
  beta_cyer * dat_simple$cyer +
  (beta_d_cyer * dat_simple$samp_date * dat_simple$cyer) +
  beta_ps * dat_simple$det

p4 <- 1 / (1 + exp(-(eta4)))  # Probability of survival
y4 <- rbinom(n, 1, p4)  # Binary outcome

d4 <- dat_simple %>% 
  mutate(
    alpha_yr = dat$alpha_yr,
    alpha_pop = dat$alpha_pop,
    surv = y4
    )

fit4 <- glmmTMB(
  surv ~ samp_date * cyer + size + lipid + det + (1 | yr_f) + (1 | pop_f),  
  data = d4,
  family = binomial()
)


## DATASET 5
# as above but with cascading effects of pop and year on sample date/condition
eta5 <- surv_bar + 
  dat$alpha_yr + dat$alpha_pop +
  beta_cs * dat$cond +
  beta_ds * dat$samp_date +
  beta_cyer * dat_simple$cyer +
  (beta_d_cyer * dat$samp_date * dat_simple$cyer) +
  beta_ps * dat_simple$det

p5 <- 1 / (1 + exp(-(eta5)))  # Probability of survival
y5 <- rbinom(n, 1, p5)  # Binary outcome

d5 <- dat %>% 
  mutate(
    cyer = dat_simple$cyer,
    det = dat_simple$det,
    surv = y5
    )

fit5 <- glmmTMB(
  surv ~ samp_date * cyer + size + lipid + det + (1 | yr_f) + (1 | pop_f),  
  data = d5,
  family = binomial()
)


## DATASET 6
# as above but with RIs that that are confounded with P (but not CYER); similar 
# effect if either switched out
eta6 <- surv_bar + 
  dat$alpha_yr + dat$alpha_pop +
  beta_cs * dat$cond +
  beta_ds * dat$samp_date +
  beta_cyer * dat_simple$cyer +
  (beta_d_cyer * dat$samp_date * dat_simple$cyer) +
  beta_ps * dat$det

p6 <- 1 / (1 + exp(-(eta6)))  # Probability of survival
y6 <- rbinom(n, 1, p6)  # Binary outcome

d6 <- dat %>% 
  mutate(
    cyer = dat_simple$cyer,
    surv = y6
    )

fit6 <- glmmTMB(
  surv ~ samp_date * cyer + size + lipid + det + (1 | yr_f) + (1 | pop_f),  
  data = d6,
  family = binomial()
)


## DATASET 7
# as above but with RIs that that are confounded with CYER and P
eta7 <- surv_bar + 
  dat$alpha_yr + dat$alpha_pop +
  beta_cs * dat$cond +
  beta_ds * dat$samp_date +
  beta_cyer * dat$cyer +
  (beta_d_cyer * dat$samp_date * dat$cyer) +
  beta_ps * dat$det

p7 <- 1 / (1 + exp(-(eta7)))  # Probability of survival
y7 <- rbinom(n, 1, p7)  # Binary outcome

d7 <- dat %>% 
  mutate(surv = y7)

fit7 <- glmmTMB(
  surv ~ samp_date * cyer + size + lipid + det + (1 | yr_f) + (1 | pop_f),  
  data = d7,
  family = binomial()
)


## DATASET 8
# as above but ignoring variability in P
fit8 <- glmmTMB(
  surv ~ samp_date * cyer + size + lipid + (1 | yr_f) + (1 | pop_f),  
  data = d7,
  family = binomial()
)



fit_list <- list(fit1, fit2, fit3, fit4, fit5, fit6, fit7, fit8)

plot_list <- purrr::map(
  fit_list, 
  function (x) {
    ests <- summary(x)$coefficients$cond
    est_df1 <- data.frame(
      par = rownames(ests),
      est = ests[, "Estimate"] %>% as.numeric(),
      se = ests[, "Std. Error"] %>% as.numeric()
    ) %>% 
      left_join(., true_pars, by = "par")
    
    ggplot(est_df1, aes(x = par)) +
      geom_point(aes(y = est), size = 3) +
      geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.2) +
      geom_point(aes(y = true), color = "red", size = 3, alpha = 0.5) +
      theme_minimal() 
  }
)

plot_list[[1]]
plot_list[[4]]
plot_list[[5]]
plot_list[[6]]
plot_list[[7]]
plot_list[[8]]


## BAYESIAN 1 ------------------------------------------------------------------

# Fit equivalent Bayesian models to 1 and 7

# model 1
dat_list1 <- list(
  surv = d1$surv,
  pop_n = d1$pop_n,
  yr = d1$yr,
  cyer = d1$cyer %>% scale() %>% as.numeric(),
  samp_date = d1$samp_date %>% scale() %>% as.numeric(),
  lipid = d1$lipid %>% scale() %>% as.numeric(),
  size = d1$size %>% scale() %>% as.numeric(),
  det = d1$det + 1,
  cond = d1$cond %>% scale() %>% as.numeric()
)
m1 <- ulam(
  alist(
    surv ~ dbinom(1 , p) ,
    logit(p) <- alpha[det] + beta_ds * samp_date +
      beta_cyer * cyer + (beta_d_cyer * samp_date * cyer) +
      beta_cs * cond,
    
    beta_cyer ~ normal(-0.5, 0.5),
    c(beta_ds, beta_cs, beta_d_cyer) ~ normal(0, 0.5),
    alpha[det] ~ dnorm(0, 1)
  ),
  data = dat_list1, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)

precis(m1)

post <- extract.samples(m1)

post_list <- purrr::map2(
  post[1:4], 
  c("cyer", "samp_date:cyer", "cond", "samp_date"),
  function (x, y) {
    data.frame(
      par = y,
      est = mean(x),
      se = sd(x)
    ) 
  }
  ) %>% 
  bind_rows()
post_int <- data.frame(
  par = "(Intercept)",
  est = mean(post$alpha[ , 1]),
  se = sd(post$alpha[ , 1])
) 
post_det <- data.frame(
  par = "det",
  est = mean(post$alpha[ , 2] - post$alpha[ , 1]),
  se = sd(post$alpha[ , 2] - post$alpha[ , 1])
) 

list(post_list, post_int, post_det) %>% 
  bind_rows() %>% 
  left_join(., true_pars, by = "par") %>% 
  ggplot(., aes(x = par)) +
  geom_point(aes(y = est), size = 3) +
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.2) +
  geom_point(aes(y = true), color = "red", size = 3, alpha = 0.5) +
  theme_minimal() 
plot_list[[1]]
# effectively identical


# model 7
dat_list7 <- list(
  surv = d7$surv,
  pop_n = d7$pop_n,
  yr = d7$yr,
  cyer = d7$cyer %>% scale() %>% as.numeric(),
  samp_date = d7$samp_date %>% scale() %>% as.numeric(),
  lipid = d7$lipid %>% scale() %>% as.numeric(),
  size = d7$size %>% scale() %>% as.numeric(),
  det = d7$det + 1,
  det_n = d7$det,
  cond = d7$cond %>% scale() %>% as.numeric()
)
m7 <- ulam(
  alist(
    surv ~ dbinom(1, p),
    logit(p) <- alpha[det] + beta_ds * samp_date +
      beta_cyer * cyer + (beta_d_cyer * samp_date * cyer) +
      beta_ls * lipid +
      beta_fs * size +
      alpha_pop[pop_n] * sigma_pop +
      alpha_yr[yr] * sigma_yr,
    
    beta_cyer ~ normal(-0.5, 0.5),
    c(beta_ds, beta_ls, beta_fs, beta_d_cyer) ~ normal(0, 0.5),
    alpha[det] ~ dnorm(0, 1),
    alpha_yr[yr] ~ dnorm(0, 0.5),
    alpha_pop[pop_n] ~ dnorm(0, 0.5),
    c(sigma_yr, sigma_pop) ~ exponential(1)
  ),
  data = dat_list7, chains = 4 , cores = 4,
  control = list(adapt_delta = 0.95)
)

precis(m7, depth = 2)
post7 <- extract.samples(m7)

## check priors (looks reasonable)
# prior7 <- extract.prior(m7)
# eta_prior <- prior7$surv_bar + 
#   prior7$beta_fs * d7$size + 
#   prior7$beta_ls * d7$lipid +
#   prior7$beta_ds * d7$samp_date +
#   prior7$beta_cyer * d7$cyer +
#   (prior7$beta_d_cyer * d7$samp_date * d7$cyer) +
#   prior7$beta_ps * d7$det

post_list7 <- purrr::map2(
  post7[1:5], 
  c("cyer", "samp_date:cyer", "size", "lipid", "samp_date"),
  function (x, y) {
    data.frame(
      par = y,
      est = mean(x),
      se = sd(x)
    ) 
  }
) %>% 
  bind_rows()
post_int7 <- data.frame(
  par = "(Intercept)",
  est = mean(post7$alpha[ , 1]),
  se = sd(post7$alpha[ , 1])
) 
post_det7 <- data.frame(
  par = "det",
  est = mean(post7$alpha[ , 2] - post7$alpha[ , 1]),
  se = sd(post7$alpha[ , 2] - post7$alpha[ , 1])
) 

list(post_list7, post_int7, post_det7) %>% 
  bind_rows() %>% 
  left_join(., true_pars, by = "par") %>% 
  ggplot(., aes(x = par)) +
  geom_point(aes(y = est), size = 3) +
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.2) +
  geom_point(aes(y = true), color = "red", size = 3, alpha = 0.5) +
  theme_minimal() 
plot_list[[7]]
# similar directionally, but fixed effects biased weak; unclear what the issue 
# is since an m3 version matches plot_list[[3]] well


## BAYESIAN 2 ------------------------------------------------------------------

# fit full version accounting for covariance parameters and see if performance
# improves relative to frequentist

m8 <- ulam(
  alist(
    # detection probability
    det_n ~ dbinom(1, p_det),
    logit(p_det) <- alpha_det_pop[pop_n] * sigma_det_pop + 
      alpha_det_yr[yr] * sigma_det_yr,
    
    # cyer
    cyer ~ dnorm(mu_cyer, sigma_cyer),
    mu_cyer <- alpha_cyer_pop[pop_n] * sigma_cyer_pop + 
      alpha_cyer_yr[yr] * sigma_cyer_yr,
    
    # date
    samp_date ~ dnorm(mu_date, sigma_date),
    mu_date <- alpha_pop[pop_n, 1],
    
    # covariance among size and lipid
    c(size, lipid) ~ multi_normal(c(mu_size, mu_lipid), Rho_sl, Sigma_sl),
    mu_size <- alpha_yr[yr, 1] + alpha_pop[pop_n, 2] + beta_df * samp_date,
    mu_lipid <- alpha_yr[yr, 2] + alpha_pop[pop_n, 3] + beta_dl * samp_date,
    
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- alpha[det] + alpha_pop[pop_n, 4] + alpha_yr[yr, 3] +
      beta_ds * samp_date +
      beta_fs * size +
      beta_ls * lipid +
      beta_cyer * cyer + (beta_d_cyer * samp_date * cyer) ,
    
    # adaptive priors
    transpars> matrix[yr, 3]:alpha_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[3, yr]:z_yr ~ normal(0 , 1),
    transpars> matrix[pop_n, 4]:alpha_pop <-
      compose_noncentered(sigma_pop , L_Rho_pop , z_pop),
    matrix[4, pop_n]:z_pop ~ normal(0 , 1),
    
    # priors
    c(beta_df, beta_dl) ~ normal(0, 1),
    beta_cyer ~ normal(-0.5, 0.5),
    c(beta_ds, beta_fs, beta_ls, beta_d_cyer) ~ normal(0, 0.5),
    alpha[det] ~ dnorm(0, 1),
    Rho_sl ~ lkj_corr( 2 ),
    cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[3]:sigma_yr ~ exponential(1),
    cholesky_factor_corr[4]:L_Rho_pop ~ lkj_corr_cholesky(2),
    vector[4]:sigma_pop ~ exponential(1),
    
    alpha_det_pop[pop_n] ~ dnorm(0, 1),
    alpha_det_yr[yr] ~ dnorm(0, 1),
    alpha_cyer_pop[pop_n] ~ dnorm(0, 1),
    alpha_cyer_yr[yr] ~ dnorm(0, 1),
    
    c(Sigma_sl, sigma_date, sigma_cyer, sigma_cyer_pop, sigma_cyer_yr, 
      sigma_det_pop, sigma_det_yr) ~ exponential(1),
    # compute ordinary correlation matrices from Cholesky factors
    gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr),
    gq> matrix[4, 4]:Rho_pop <<- Chol_to_Corr(L_Rho_pop)
  ),
  data = dat_list7, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)

m8b <- ulam(
  alist(
    # date
    samp_date ~ dnorm(mu_date, sigma_date),
    mu_date <- alpha_pop[pop_n, 1],
    
    # covariance among size and lipid
    c(size, lipid) ~ multi_normal(c(mu_size, mu_lipid), Rho_sl, Sigma_sl),
    mu_size <- alpha_yr[yr, 1] + alpha_pop[pop_n, 2] + beta_df * samp_date,
    mu_lipid <- alpha_yr[yr, 2] + alpha_pop[pop_n, 3] + beta_dl * samp_date,
    
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- alpha[det] + alpha_pop[pop_n, 4] + alpha_yr[yr, 3] +
      beta_ds * samp_date +
      beta_fs * size +
      beta_ls * lipid +
      beta_cyer * cyer + (beta_d_cyer * samp_date * cyer),
    
    # adaptive priors
    transpars> matrix[yr, 3]:alpha_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[3, yr]:z_yr ~ normal(0 , 1),
    transpars> matrix[pop_n, 4]:alpha_pop <-
      compose_noncentered(sigma_pop , L_Rho_pop , z_pop),
    matrix[4, pop_n]:z_pop ~ normal(0 , 1),
    
    # priors
    c(beta_df, beta_dl) ~ normal(0, 1),
    beta_cyer ~ normal(-0.5, 0.5),
    c(beta_ds, beta_fs, beta_ls, beta_d_cyer) ~ normal(0, 0.5),
    alpha[det] ~ dnorm(0, 1),
    Rho_sl ~ lkj_corr( 2 ),
    cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[3]:sigma_yr ~ exponential(1),
    cholesky_factor_corr[4]:L_Rho_pop ~ lkj_corr_cholesky(2),
    vector[4]:sigma_pop ~ exponential(1),
    c(Sigma_sl, sigma_date) ~ exponential(1),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr),
    gq> matrix[4, 4]:Rho_pop <<- Chol_to_Corr(L_Rho_pop)
  ),
  data = dat_list7, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)



post8 <- extract.samples(m8)
post8b <- extract.samples(m8b)

post_list8 <- purrr::map2(
  list(
    post8$beta_cyer, post8$beta_d_cyer, post8$beta_fs, post8$beta_ls, 
    post8$beta_ds
  ), 
  c("cyer", "samp_date:cyer", "size", "lipid", "samp_date"),
  function (x, y) {
    data.frame(
      par = y,
      est = mean(x),
      se = sd(x)
    ) 
  }
) %>% 
  bind_rows()
post_int8 <- data.frame(
  par = "(Intercept)",
  est = mean(post8$alpha[ , 1]),
  se = sd(post8$alpha[ , 1])
) 
post_det8 <- data.frame(
  par = "det",
  est = mean(post8$alpha[ , 2] - post8$alpha[ , 1]),
  se = sd(post8$alpha[ , 2] - post8$alpha[ , 1])
) 

list(post_list8, post_int8, post_det8) %>% 
  bind_rows() %>% 
  left_join(., true_pars, by = "par") %>% 
  ggplot(., aes(x = par)) +
  geom_point(aes(y = est), size = 3) +
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.2) +
  geom_point(aes(y = true), color = "red", size = 3, alpha = 0.5) +
  theme_minimal() 
plot_list[[7]]



# examine RIs specifically
fit_list_bayes <- list(post7, post8b, post8)

ri_pop_dat <- purrr::map2(
  fit_list_bayes, 
  c("simp", "nest1", "nest2"),
  function (x, x2) {
    data.frame(
      par = paste("alpha_pop", seq(1, 5, by = 1), sep = "_"),
      est = apply(post8$alpha_pop[ , , 4], 2, mean),
      se = apply(post8$alpha_pop[ , , 4], 2, sd),
      true = alpha_pop_surv,
      model = x2
    ) 
  }
) %>% 
  bind_rows()

ggplot(ri_pop_dat, aes(x = model)) +
  geom_point(aes(y = est), size = 3) +
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.2) +
  geom_point(aes(y = true), color = "red", size = 3, alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~par)


ri_yr_dat <- purrr::map2(
  fit_list_bayes, 
  c("simp", "nest1", "nest2"),
  function (x, x2) {
    data.frame(
      par = paste("alpha_yr", seq(1, 5, by = 1), sep = "_"),
      est = apply(post8$alpha_yr[ , , 3], 2, mean),
      se = apply(post8$alpha_yr[ , , 3], 2, sd),
      true = alpha_yr_surv,
      model = x2
    ) 
  }
) %>% 
  bind_rows()

ggplot(ri_yr_dat, aes(x = model)) +
  geom_point(aes(y = est), size = 3) +
  geom_errorbar(aes(ymin = est - se, ymax = est + se), width = 0.2) +
  geom_point(aes(y = true), color = "red", size = 3, alpha = 0.5) +
  theme_minimal() +
  facet_wrap(~par)
# minimal differences in means for RIs; though positive correlation among RIs 
# is estimated, it's biased low (perhaps due to small number of levels)

# focus on estimating parameters of interest only, including detection 
# probability



## ORDERED CATEGORICAL RESPONSE ------------------------------------------------

## Test model estimating size effects on injury for potential incorporation as 
# submodel

library(rethinking)

# Simulate example data
set.seed(123)
N <- 200  # Sample size
X <- rnorm(N, mean = 0, sd = 1)  # Continuous predictor

# Generate latent variable and assign ordered categories
latent_Y <- 1.5 * X + rnorm(N, mean = 0, sd = 0.4)  # Linear effect with noise
thresholds <- quantile(latent_Y, probs = c(0.25, 0.50, 0.75))  # Define cutpoints

# Convert latent Y into ordered categorical response
Y <- cut(latent_Y, breaks = c(-Inf, thresholds, Inf), labels = 1:4,
         ordered_result = TRUE)
Y <- as.integer(Y)  # Convert to numeric for modeling (1, 2, 3, 4)

ggplot(data.frame(X, Y = factor(Y)), aes(x = X, fill = Y)) +
  geom_histogram(bins = 30, position = "stack") +
  scale_fill_brewer(palette = "Set2") +
  labs(x = "X (Continuous Predictor)", y = "Count", fill = "Category",
       title = "Category Distribution Across X") +
  theme_minimal()

data.frame(X = X, Y = Y) %>%
  mutate(
    X_quartile = cut(X, breaks = quantile(X, probs = seq(0, 1, by = 0.25)), 
                     include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4"))
  ) %>% 
  ggplot(., aes(x = X_quartile, fill = as.factor(Y))) +
  geom_bar(position = "stack") +
  labs(x = "X Quartile", y = "Count of Y", fill = "Y") +
  ggtitle("Stacked Bar Plot of Y Across Quartiles of X") +
  theme_minimal()


# Create the data list for Stan
data_list <- list(
  N = N,
  X = X,
  Y = Y
)

# Define the ordinal regression model using the log-cumulative-odds function
m_ord <- ulam(
  alist(
    Y ~ dordlogit(phi, cutpoints),  # Ordered logistic likelihood
    phi <- a + b * X,  # Linear predictor
    a ~ normal(0, 1),  # Intercept prior
    b ~ normal(0, 1),  # Slope prior
    cutpoints ~ normal(0, 1)  # Priors for category thresholds
  ), data = data_list, chains = 4, cores = 4
)

# Summarize the model
precis(m_ord, depth = 2)


# Define new X values at which we want predictions
new_X <- c(-1, 0, 1)

# Extract posterior samples of model parameters
post <- extract.samples(m_ord)

# Compute the latent variable phi for new values of X
phi_samples <- sapply(new_X, function(x) {
  post$a + post$b * x
})

# Convert phi into category probabilities using the ordinal logistic function
pred_probs <- lapply(1:ncol(phi_samples), function(i) {  # Loop over X values
  t(sapply(1:nrow(phi_samples), function(j) {  # Loop over posterior samples
    pordlogit(1:4, phi_samples[j, i], post$cutpoints[j, ])  # Use correct cutpoints
  }))
})

cat_probs <- lapply(pred_probs, function(p) {
  t(sapply(1:nrow(p), function(j) {  # Loop over posterior samples
    c(p[j, 1], diff(p[j, ]))
  }))
})


# Compute median probabilities across posterior samples
median_probs <- sapply(cat_probs, function(p) apply(p, 2, median))

# Convert to long format for ggplot
df_plot <- as.data.frame(median_probs) %>%
  mutate(Category = factor(1:4)) %>%
  pivot_longer(starts_with("V"), names_to = "X", values_to = "Probability") %>% 
  mutate(
    X = factor(X, levels = c("V1", "V2", "V3"), labels = c("-1", "0", "1"))
  )
ggplot(df_plot, aes(x = X, y = Probability, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "X Value", y = "Median Predicted Proportion", fill = "Y Category") +
  ggtitle("Stacked Bar Chart of Median Predicted Proportions") +
  theme_minimal()


## EXPLORE OBSERVATION ERROR ---------------------------------------------------

## true survival is a function of one covariate and observed survival includes
# error
surv_bar <- -1  # intercept
beta_cs <- 0.5 # condition effect

eta <- surv_bar #+ beta_cs * dat$cond 
p <- 1 / (1 + exp(-(eta)))  # Probability of survival
set.seed(123)
s_true <- rbinom(n, 1, p)  # Binary outcome

# detection probability varies across five groups with the first three having
# posterior estimates from a previous model that can be passed as data
dat_obs <- dat %>% 
  mutate(
    group_id = sample(c(1, 5), size = nrow(dat), replace = TRUE),
    det_p = ifelse(group_id == "1", 0.9, 0.3),
    s_obs = ifelse(s_true == 1 & runif(n) < det_p, 1, 0)
  )

dat_list <- list(
  s_obs = dat_obs$s_obs,
  group_id = dat_obs$group_id
)

# m1b <- ulam(
#   alist(
#     # obs model
#     s_obs | s_obs == 1 ~ custom( 
#       log_p + log(det_p) 
#       ),
#     s_obs | s_obs == 0 ~ custom(
#       log_sum_exp(log1m_exp(log_p), log_p + log1m(det_p)) 
#       ),
#     # Prior for detection probability
#     det_p ~ beta(19, 2),  # Adjusted to favor ~0.95 detection probability
#     # det_p ~ beta(3, 2),  # Adjusted to favor ~0.6 detection probability
#     
#     # Process model
#     log_p <- log_inv_logit(surv_bar + beta_cs * cond),
#     
#     # Priors
#     beta_cs ~ normal(0, 0.5),
#     surv_bar ~ dnorm(0, 2)
#   ),
#   data = dat_list, chains=4 , cores=4,
#   control = list(adapt_delta = 0.95)
# )

## not possible to estimate separate detection probabilities with ulam,
# write in stan instead
# also modify detection probability so that it varies across five groups with 
# the first three having posterior estimates from a previous model that can be 
# passed as data
# NOTE remove condition effect since it reduces precision of estimates

# posterior matrix
set.seed(123)  # For reproducibility

M <- 4000  # Number of posterior samples
det_p_posterior <- matrix(0, nrow = M, ncol = 5)
det_p_posterior[, 1] <- rbeta(M, shape1 = 40, shape2 = 7)   # Median ~0.85
det_p_posterior[, 2] <- rbeta(M, shape1 = 50, shape2 = 4)   # Median ~0.92
det_p_posterior[, 3] <- rbeta(M, shape1 = 75, shape2 = 1.5) # Median ~0.98


# true det_p
det_p_true <- data.frame(
  group_id = seq(1, 5, 1),
  det_p = c(.85, .92, .98, 0.3, 0.2), #median
  use_posterior = ifelse(det_p > 0.5, 1, 0)
)

dat_obs <- dat %>% 
  mutate(
    group_id = sample(seq(1, 5, 1), size = nrow(dat), replace = TRUE)
  ) %>%  
  left_join(., det_p_true, by = "group_id") %>%
  mutate(
    s_obs = ifelse(s_true == 1 & runif(n) < det_p, 1, 0)
  )

ggplot(dat_obs) +
  geom_bar(aes(x = s_obs)) +
  facet_wrap(~group_id)

dat_list <- list(
  N = nrow(dat_obs),
  s_obs = dat_obs$s_obs,
  group_id = dat_obs$group_id,
  G = length(unique(dat_obs$group_id)),
  use_posterior = det_p_true$use_posterior,
  # P = sum(det_p_true$use_posterior), # number of groups using posterior
  M = nrow(det_p_posterior),
  det_p_posterior = det_p_posterior
)

library(rstan)

samp_mod <- stan_model(here::here("R", "stan_models", "obs-error-example.stan"))

m1_stan <- sampling(samp_mod, data = dat_list, 
                    # chains = 1, iter = 2000, warmup = 1000, 
                    chains = 4, iter = 4000, warmup = 1000,
                    control = list(adapt_delta = 0.95))

print(m1_stan, pars = c(#"beta_cs", 
                        "surv_bar", "det_p_out"))
