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
  mu_det_prob = logit(c(.85, .92, .98, .8, .7))
  # mu_det_prob = logit(c(.85, .92, .98, 0.2, 0.3))
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
    post = 1,#ifelse(pop %in% c("d", "e"), 0, 1),
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


# pop intercepts (date, lipid, size, survival)
# correlated versions
pop_rho <- 0.9
pop_mu <- c(0, 0, 0, 0) #average effects both centered on zero
pop_sigmas <- c(0.3, 0.4, 0.4, 0.5)
pop_Rho <- diag(1, 4) + matrix(pop_rho, 4, 4) - diag(pop_rho, 4)
pop_Sigma <- diag(pop_sigmas) %*% pop_Rho %*% diag(pop_sigmas)
pop_ints <- MASS::mvrnorm(n_pops, pop_mu, pop_Sigma)
alpha_pop_date <- pop_ints[ , 1]
alpha_pop_size <- pop_ints[ , 2]
alpha_pop_lipid <- pop_ints[ , 3]
alpha_pop_surv <- pop_ints[ , 4]


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
beta_cl <- 0.9 
beta_cf <- 0.4

for (i in 1:nrow(dat)) { 
  dat$samp_date[i] <- rnorm(1, alpha_pop_date[dat$pop_n[i]], 0.4)
  dat$cond[i] <- rnorm(
    1, 
    alpha_yr_cond[dat$yr[i]] + beta_dc * dat$samp_date[i], 
    0.4
  )
  # generate observations from latent state 
  dat$size[i] <- rnorm(1, beta_cf * dat$cond[i] + 
                         alpha_pop_size[dat$pop_n[i]], 0.4)
  dat$lipid[i] <- rnorm(1, beta_cl * dat$cond[i] + 
                          alpha_pop_lipid[dat$pop_n[i]], 0.4)
  dat$alpha_yr[i] <- alpha_yr_surv[dat$yr[i]]
  dat$alpha_pop[i] <- alpha_pop_surv[dat$pop_n[i]]
}


## fixed effect parameters
surv_bar <- 1  # intercept
# beta_cs <- 0.5 # condition effect
beta_fs <- 0.2 # size effect
beta_ls <- 0.4 # lipid effect
beta_ds <- 0.5 # date effect
beta_cyer <- -1 # er effect
beta_d_cyer <- 0.5 # date ER interaction

true_pars <- data.frame(
  par = c("beta_dc", "beta_cl", "beta_cf", "alpha_bar", "beta_ds", "beta_cyer",
          "beta_cs", "beta_fs", "beta_ls", "beta_d_cyer"),
  true = c(beta_dc, beta_cl, beta_cf, surv_bar, beta_ds, beta_cyer, beta_cs,
           beta_fs, beta_ls, beta_d_cyer)
)

# focus on survival only
true_pars_year_ri <- data.frame(
    year = seq(1, 5, by = 1),
    true = alpha_yr_surv
  )
true_pars_pop_ri <- data.frame(
  pop = seq(1, 5, by = 1),
  true = alpha_pop_surv
)


## simulate true survival
eta <- surv_bar + 
  dat$alpha_yr + dat$alpha_pop +
  # beta_cs * dat$cond +
  beta_fs * dat$size +
  beta_ls * dat$lipid +
  beta_ds * dat$samp_date +
  beta_cyer * dat$cyer +
  (beta_d_cyer * dat$samp_date * dat$cyer)

p_true <- 1 / (1 + exp(-(eta)))  # Probability of survival
dat$s_true <- rbinom(n, 1, p_true)  # Binary outcome


# adjust observed survival to reflect detection probability
dat$s_obs <- ifelse(
  dat$s_true == 1 & runif(n) < dat$det_prob, 1, 0
)

ggplot(dat) +
  geom_bar(aes(x = s_obs)) +
  facet_wrap(~pop)
ggplot(dat) +
  geom_bar(aes(x = s_true)) +
  facet_wrap(~pop)


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
  
  logit_det_mu = det_key$logit_det_prob,
  logit_det_sd = det_key$sd_logit_det_prob,
  
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
# fit <- ulam(
#   alist(
#     # date
#     samp_date ~ dnorm(mu_date, sigma_date),
#     mu_date <- alpha_pop[pop_n, 1],
#     
#     # covariance among size and lipid
#     c(size, lipid) ~ multi_normal(c(mu_size, mu_lipid), Rho_sl, Sigma_sl),
#     mu_size <- alpha_yr[yr, 1] + alpha_pop[pop_n, 2] + beta_df * samp_date,
#     mu_lipid <- alpha_yr[yr, 2] + alpha_pop[pop_n, 3] + beta_dl * samp_date,
#     
#     # survival
#     s_true ~ dbinom(1 , p) ,
#     logit(p) <- alpha_bar + alpha_pop[pop_n, 4] + alpha_yr[yr, 3] +
#       beta_ds * samp_date +
#       beta_fs * size +
#       beta_ls * lipid +
#       beta_cyer * cyer + (beta_d_cyer * samp_date * cyer),
#     
#     # adaptive priors
#     transpars> matrix[yr, 3]:alpha_yr <-
#       compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
#     matrix[3, yr]:z_yr ~ normal(0 , 1),
#     transpars> matrix[pop_n, 4]:alpha_pop <-
#       compose_noncentered(sigma_pop , L_Rho_pop , z_pop),
#     matrix[4, pop_n]:z_pop ~ normal(0 , 1),
#     
#     # priors
#     c(beta_df, beta_dl) ~ normal(0, 1),
#     beta_cyer ~ normal(-0.5, 0.5),
#     c(beta_ds, beta_fs, beta_ls, beta_d_cyer) ~ normal(0, 0.5),
#     alpha_bar ~ dnorm(0, 1),
#     Rho_sl ~ lkj_corr( 2 ),
#     cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
#     vector[3]:sigma_yr ~ exponential(1),
#     cholesky_factor_corr[4]:L_Rho_pop ~ lkj_corr_cholesky(2),
#     vector[4]:sigma_pop ~ exponential(1),
#     c(Sigma_sl, sigma_date) ~ exponential(1),
#     
#     # compute ordinary correlation matrixes from Cholesky factors
#     gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr),
#     gq> matrix[4, 4]:Rho_pop <<- Chol_to_Corr(L_Rho_pop)
#   ),
#   # set iterations and chains low because just want to extract stan code
#   data = dat_list, chains = 4, cores = 4, iter = 2000,
#   control = list(adapt_delta = 0.95)
# )
# 
# precis(fit, depth = 2)
# 
# stancode(fit)


# true survival model
mod1 <- stan_model(here::here("R", "stan_models", "true_surv_jll.stan"))

m1_stan <- sampling(mod1, data = dat_list, 
                    chains = 4, iter = 2000, warmup = 1000,
                    control = list(adapt_delta = 0.95))

# observed survival model
mod2 <- stan_model(here::here("R", "stan_models", "obs_surv_jll.stan"))

m2_stan <- sampling(mod2, data = dat_list, 
                    chains = 4, iter = 2000, warmup = 1000,
                    control = list(adapt_delta = 0.95))



# simplified observed survival model (no RIs)
# mod3 <- stan_model(here::here("R", "stan_models", "obs_surv_jll_FE.stan"))
# 
# m3_stan <- sampling(mod3, data = dat_list, 
#                     chains = 4, iter = 2000, warmup = 1000,
#                     control = list(adapt_delta = 0.95))


alpha_pars <- names(m2_stan)[grepl("alpha", names(m2_stan))]
alpha_pop_pars <- names(m2_stan)[grepl("alpha_pop", names(m2_stan))]
sigma_pars <- names(m2_stan)[grepl("sigma", names(m2_stan))]
beta_pars <- names(m2_stan)[grepl("beta", names(m2_stan))]

print(m1_stan, 
      pars = "alpha_bar")
print(m2_stan, 
      pars = "alpha_bar")

fit_list <- list(m1_stan, m2_stan)
post_list <- purrr::map(fit_list, ~ as_draws_df(.x))


alpha_dat <- purrr::map2(
  post_list[1:2],
  c("m1", "m2"),
  function(x, y) {
    a_yr <- x %>%
      spread_draws(alpha_yr[i, j]) %>%
      filter(j == "3") %>% 
      rename(year = i, value = alpha_yr) %>% 
      left_join(., true_pars_year_ri, by = "year") %>% 
      mutate(
        par = paste("alpha_yr", year, sep = "_")
      ) %>%
      ungroup() %>% 
      select(value, true, par)
    a_pop <- x %>%
      spread_draws(alpha_pop[k, l]) %>%
      filter(l == "4") %>%
      rename(pop = k, value = alpha_pop) %>%
      left_join(., true_pars_pop_ri, by = "pop") %>%
      mutate(
        par = paste("alpha_pop", pop, sep = "_")
      ) %>%
      ungroup() %>%
      select(colnames(a_yr))
    a_bar <- x %>%
      spread_draws(alpha_bar) %>%
      rename(value = alpha_bar) %>%
      mutate(
        par = "alpha_bar"
      ) %>%
      left_join(., true_pars, by = "par") %>%
      ungroup() %>%
      select(value, true, par)
    
    list(a_yr, a_pop, a_bar) %>%
      bind_rows() %>%
      mutate(
        model = y
      )
  }
) %>% 
  bind_rows()

ggplot(alpha_dat) +
  geom_boxplot(aes(x = model, y = value)) +
  geom_point(aes(x = model, y = true), colour = "red") +
  facet_wrap(~par) +
  ggsidekick::theme_sleek()

# beta_dat <- purrr::map2(
#   post_list[[2]],
#   c("m2"),
#   function(x, y) {
#     x %>%
#       spread_draws(beta_dc, beta_cl, beta_cf, beta_ds, beta_cyer,
#                    beta_fs, beta_ls, beta_d_cyer) %>% 
#       pivot_longer(starts_with("beta_"), names_to = "par", 
#                    values_to = "value") %>% 
#       mutate(
#         model = y
#       )
#   }
# ) %>% 
#   bind_rows()

beta_dat <- post_list[[2]] %>%
  spread_draws(beta_dc, beta_cl, beta_cf, beta_ds, beta_cyer,
               beta_fs, beta_ls, beta_d_cyer) %>% 
  pivot_longer(starts_with("beta_"), names_to = "par", 
               values_to = "value") %>% 
  mutate(
    model = "m2"
  )

ggplot() +
  geom_boxplot(data = beta_dat, aes(x = par, y = value)) +
  geom_point(data = true_pars %>% filter(par %in% beta_dat$par),
             aes(x = par, y = true), colour = "red") +
  ggsidekick::theme_sleek() +
  facet_wrap(~model)


det_p_draws <- post_list[[2]] %>%
  spread_draws(p[i]) %>% 
  rename(det_group_id = i) %>% 
  left_join(., det_key, by = "det_group_id") %>% 
  ungroup()

ggplot() +
  geom_boxplot(data = det_p_draws, aes(x = as.factor(yr), y = p)) +
  geom_point(data = det_key, aes(x = as.factor(yr), y = det_prob), 
             colour = "red") +
  facet_wrap(~pop) +
  ggsidekick::theme_sleek()


# posterior detection
dat$ind <- seq(1, nrow(dat), by = 1)
mu_draws <- post_list[[2]] %>%
  spread_draws(s_obs_rep[i]) %>%
  rename(ind = i) %>% 
  ungroup() %>% 
  group_by(ind) %>% 
  summarize(mean_s = mean(s_obs_rep))


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
post_rep <- as.matrix(m2_stan, pars = "s_obs_rep")
pit_residuals <- calc_pit(y = dat$s_obs, posterior_pred = post_rep)
qqplot(qunif(ppoints(length(pit_residuals))), pit_residuals,
       main = "QQ-plot of PIT Residuals")
abline(0, 1)

post_rep1 <- as.matrix(m1_stan, pars = "s_obs_rep")
length(post_rep1[post_rep1 == "0"]) / length(post_rep1)
length(dat$s_obs[dat$s_obs == "0"]) / length(dat$s_obs)

post_rep3 <- as.matrix(m3_stan, pars = "s_obs_rep")
length(post_rep3[post_rep3 == "0"]) / length(post_rep3)
length(dat$s_obs[dat$s_true == "0"]) / length(dat$s_true)


length(post_rep[post_rep == "0"]) / length(post_rep)
length(dat$s_true[dat$s_true == "0"]) / length(dat$s_true)
