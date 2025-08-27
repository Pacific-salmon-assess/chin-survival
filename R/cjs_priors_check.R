### CJS Priors Check
## July 8, 2024
## Conduct prior predictive checks using hierarchical CJS models without data

library(tidyverse)
library(rstan)
library(shinystan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
set.seed(123)


# simple sim 
# simulate distributions for beta_phi and beta_p, representing stage
# specific survival parameters in link space; since 3 beta pars
mu_phi <- rnorm(1000, mean = 0.8, sd = 1)
hist(boot::inv.logit(mu_phi), col=rgb(1,0,0,0.4))
beta_phi <- rnorm(1000, mean = 0, sd = 1)
gamma_phi <- rnorm(1000, mean = 0, sd = 0.5)
beta_p <- rnorm(1000, mean = 0.25, sd = 1)
beta_p_j <- rnorm(1000, mean = 0, sd = 1)
beta_p_j2 <- rnorm(1000, mean = 0, sd = 0.5)

# transform into survival and detection probabilities and evaluate
total_phi <- mu_phi + beta_phi + gamma_phi
hist(boot::inv.logit(mu_phi),  col=rgb(0,0,1,0.4))
hist(boot::inv.logit(total_phi),  col=rgb(1,0,0,0.4), add = T)
hist(boot::inv.logit(beta_p ), col=rgb(0,0,1,0.4), )
hist(boot::inv.logit(beta_p + beta_p_j), col=rgb(1,0,0,0.4), add = T)



#import sample data genereted in survival_cjs_hier.R
# dd <- readRDS(here::here("data", "model_outputs", "sample_cjs_dat.rds"))
dd <- readRDS(here::here("data", "model_outputs", "sample_cjs_dat_date.rds"))

# stan code 
hier_mod_sims_priors_only <- stan_model(
  here::here("R", "stan_models", "cjs_add_hier_eff_adaptive_date_priorsONLY.stan")
)


n_chains = 1
n_iter = 1000
n_warmup = n_iter / 2
params_fixp <- c(
  "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
  "beta_date_phi", "alpha_p", "alpha_yr_p", 
  # transformed pars or estimated quantities
  "alpha_yr_phi", "phi_yr", "p_yr", "beta_yr"
)
alpha_yr_p_dim <- dd$nyear * (dd$n_occasions)

inits <- lapply(1:n_chains, function (i) {
  list(
    alpha_phi = rnorm(1, 0, 0.5),
    # note z transformed so inverted compared to beta_phi or beta_p
    alpha_yr_phi_z = rnorm(dd$nyear, 0, 0.5),
    alpha_t_phi = rnorm(dd$n_occasions - 1, 0, 0.5),
    sigma_alpha_yr_phi = rexp(1, 2),
    beta_date_phi = rnorm(1, 0, 0.5),
    alpha_p = rnorm(1, 0, 0.5),
    alpha_yr_p = matrix(rnorm(dd$nyear * (dd$n_occasions), 0, 0.5),
                        nrow = dd$nyear)
  )
})

fit_priors <- sampling(hier_mod_sims_priors_only, data = dd, pars = params_fixp,
         init = inits, chains = n_chains, iter = n_iter, warmup = n_warmup,
         open_progress = FALSE,
         control = list(adapt_delta = 0.95))



phi_est <- extract(fit_priors)[["phi_yr"]]
# drop first column fixed to one
p_est <- extract(fit_priors)[["p_yr"]][ , , 2:6]
mean_phi_est <- extract(fit_priors)[["alpha_phi"]]
mean_p_est <- extract(fit_priors)[["alpha_p"]]
hist(phi_est)
hist(p_est)
hist(boot::inv.logit(mean_phi_est))
hist(boot::inv.logit(mean_p_est))
