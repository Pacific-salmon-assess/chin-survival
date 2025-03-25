## Data Structure Simulation
# Adapts binomial_glm_sim_full.R to account for imperfect detection efficiency
# Requires adjusting model to be estimated in stan not rethinking
# March 22, 2025


library(tidyverse)
library(glmmTMB)
library(rethinking)
library(rstan)
library(tidybayes)


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
# simplifying assumption that posterior is used as true mean
det_key1 <- data.frame(
  pop = c("a", "b", "c", "d", "e"),
  mu_det_prob = logit(c(.85, .92, .98, 0.2, 0.3))
) 
det_key2 <- vector(mode = "list", length = length(unique(det_key1$pop)))
for (i in 1:nrow(det_key1)) {
  det_key2[[i]] <- expand.grid(
    pop = det_key1$pop[i],
    yr = yrs
  ) %>% 
    mutate(
      logit_det_prob = rnorm(length(unique(yr)), det_key1$mu_det_prob[i], sd = 0.4),
      det_prob = inv_logit(logit_det_prob)
    )
}
det_key <- bind_rows(det_key2) %>% 
  mutate(
    sd_logit_det_prob = 0.4,
    # define groups that have posterior samples
    post = ifelse(pop %in% c("d", "e"), 0, 1),
    det_group_id = paste(pop, yr, sep = "_") %>% 
      as.factor() %>% 
      as.numeric()
  ) 


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


## Make dataframe with full structure of correlations
dat <- data.frame(
  yr = sample(yrs, n, replace = TRUE), 
  pop = sample(pops, n, replace =  TRUE)
) %>% 
  left_join(
    ., det_key, by = c("yr", "pop")
  ) %>% 
  left_join(
    ., cyer_dat, by = c("yr", "pop")
  ) %>%
  mutate(
    #scale ER data
    cyer = scale(cyer) %>% as.numeric(),
    samp_date = NA,
    cond = NA,
    size = NA,
    lipid = NA,
    yr_f = as.factor(yr),
    pop_f = as.factor(pop),
    pop_n = as.numeric(pop_f)
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


## fixed effect parameters
surv_bar <- -0.5  # intercept
beta_cs <- 0.5 # condition effect
beta_ds <- 0.4 # date effect
beta_cyer <- -1 # er effect
beta_d_cyer <- 0.5 # date ER interaction

true_pars <- data.frame(
  par = c("(Intercept)", "beta_ds", "beta_cyer", "beta_cs", "beta_fs", "beta_ls",
          "beta_d_cyer"),
  true = c(surv_bar, beta_ds, beta_cyer, beta_cs, 0.3 * beta_cs, 0.9 * beta_cs,
           beta_d_cyer)
)


## simulate true survival
eta <- surv_bar + 
  dat$alpha_yr + dat$alpha_pop +
  beta_cs * dat$cond +
  beta_ds * dat$samp_date +
  beta_cyer * dat$cyer +
  (beta_d_cyer * dat$samp_date * dat$cyer) 

p_true <- 1 / (1 + exp(-(eta)))  # Probability of survival
dat$s_true <- rbinom(n, 1, p_true)  # Binary outcome


# adjust observed survival to reflect detection probability
dat$s_obs <- ifelse(
  dat$s_true == 1 & runif(n) < dat$det_prob, 1, 0
)


dat_list <- list(
  N = length(dat$s_true),
  N_year = length(unique(dat$yr)),
  N_stock = length(unique(dat$pop_n)),
  N_det_id = length(unique(det_key$det_group_id)),
  N_det_id_obs = det_key %>% 
    filter(post == 1) %>% 
    pull(det_group_id) %>% 
    unique() %>% 
    length(),
  use_posterior = det_key$post,
  
  s_true = dat$s_true,
  s_obs = dat$s_obs,
  
  logit_det_mu = dat$logit_det_prob,
  logit_det_sd = dat$sd_logit_det_prob,
  
  pop_n = dat$pop_n,
  yr = dat$yr,
  det_group_id = dat$det_group_id,

  cyer = dat$cyer %>% scale() %>% as.numeric(),
  samp_date = dat$samp_date %>% scale() %>% as.numeric(),
  lipid = dat$lipid %>% scale() %>% as.numeric(),
  size = dat$size %>% scale() %>% as.numeric()
  )



## MODEL FIT -------------------------------------------------------------------

## rethinking version with no observation process
fit <- ulam(
  alist(
    # date
    samp_date ~ dnorm(mu_date, sigma_date),
    mu_date <- alpha_pop[pop_n, 1],
    
    # covariance among size and lipid
    c(size, lipid) ~ multi_normal(c(mu_size, mu_lipid), Rho_sl, Sigma_sl),
    mu_size <- alpha_yr[yr, 1] + alpha_pop[pop_n, 2] + beta_df * samp_date,
    mu_lipid <- alpha_yr[yr, 2] + alpha_pop[pop_n, 3] + beta_dl * samp_date,
    
    # survival
    s_true ~ dbinom(1 , p) ,
    logit(p) <- alpha_bar + alpha_pop[pop_n, 4] + alpha_yr[yr, 3] +
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
    alpha_bar ~ dnorm(0, 1),
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
  # set iterations and chains low because just want to extract stan code
  data = dat_list, chains = 4, cores = 4, iter = 2000,
  control = list(adapt_delta = 0.95)
)

precis(fit, depth = 2)

stancode(fit)


# true survival model
mod1 <- stan_model(here::here("R", "stan_models", "true_surv_jll.stan"))

m1_stan <- sampling(mod1, data = dat_list, 
                    chains = 4, iter = 2000, warmup = 1000,
                    control = list(adapt_delta = 0.95))

# true survival model
mod2 <- stan_model(here::here("R", "stan_models", "obs_surv_jll.stan"))

m2_stan <- sampling(mod2, data = dat_list, 
                    chains = 4, iter = 2000, warmup = 1000,
                    control = list(adapt_delta = 0.95))




alpha_pars <- names(m1_stan)[grepl("alpha", names(m1_stan))]
alpha_pop_pars <- names(m1_stan)[grepl("alpha_pop", names(m1_stan))]
sigma_pars <- names(m1_stan)[grepl("sigma", names(m1_stan))]
beta_pars <- names(m1_stan)[grepl("beta", names(m1_stan))]

print(m1_stan, 
      pars = c(alpha_pars, beta_pars, sigma_pars))
print(m1_stan, 
      pars = alpha_pop_pars)


post <- as_draws_df(m1_stan)

# Extract parameters (excluding transformed or latent variables if needed)
draws <- post %>%
  # spread_draws(alpha_bar, alpha_yr[i, j], alpha_pop[k, l]) %>% 
  spread_draws(beta_dl, beta_df, beta_cyer, beta_d_cyer, beta_ls, beta_fs, 
               beta_ds) %>% 
  pivot_longer(starts_with("beta_"), names_to = "par", 
               values_to = "value") %>% 
  left_join(., true_pars, by = "par")

ggplot(draws %>% filter(!is.na(true))) +
  geom_boxplot(aes(x = par, y = value)) +
  geom_point(aes(x = par, y = true), colour = "red") +
  ggsidekick::theme_sleek()

