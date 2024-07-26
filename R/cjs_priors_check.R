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
mu_phi <- rnorm(1000, mean = 0.3, sd = 1)
hist(boot::inv.logit(mu_phi),  col=rgb(1,0,0,0.4))
beta_phi <- gamma_phi <- rnorm(1000, mean = 0, sd = 0.5)
beta_p <- rnorm(1000, mean = 0, sd = 1.6)

# transform into survival and detection probabilities and evalute
total_phi <- mu_phi + beta_phi + gamma_phi
hist(boot::inv.logit(total_phi),  col=rgb(1,0,0,0.4))
hist(boot::inv.logit(beta_p), col=rgb(0,0,1,0.4), add = T)



#import sample data genereted in survival_cjs_hier.R
dd <- readRDS(here::here("data", "generated_data", "sample_cjs_dat.rds"))

params <- c(
  "alpha_phi", "alpha_t_phi", "alpha_yr_phi_z", "sigma_alpha_yr_phi",
  "L_Rho_yr", "alpha_p", "alpha_yr_p",
  # transformed pars or estimated quantities
  "Rho_yr", "phi_yr", "p_yr", "beta_yr"
)

## FOR DEBUGGING
inits <- list(
    alpha_phi = rnorm(1, 0, 0.5),
    # note z transformed so inverted compared to beta_phi or beta_p
    alpha_yr_phi_z = matrix(
      rnorm(dd$nyear * (dd$n_occasions - 1), 0, 0.5), nrow = (dd$n_occasions - 1)
    ),
    alpha_t_phi = rnorm(dd$n_occasions - 1, 0, 0.5),
    sigma_alpha_yr_phi = rexp((dd$n_occasions - 1), 2),
    # replace with rlkjcorr(XXX, K = 2, eta = 2) from rethinking package
    L_Rho_yr = matrix(
      runif((dd$n_occasions - 1)^2, -0.5, 0.5), nrow = (dd$n_occasions - 1)
    ),
    alpha_p = rnorm(1, 0, 0.5),
    alpha_yr_p = matrix(
      rnorm(dd$nyear * (dd$n_occasions), 0, 0.5), nrow = dd$nyear
    )
  )

# stan code 
hier_mod_sims_priors_only <- stan_model(
  here::here("R", "stan_models", "cjs_add_hier_eff_adaptive_v2_priorsONLY.stan")
)

fit_priors <- sampling(
  hier_mod_sims_priors_only, data = dd, pars = params,
  init = inits, chains = 1, iter = 4000, warmup = 50,
  open_progress = FALSE
)

phi_est <- extract(fit_priors)[["phi_yr"]]
p_est <- extract(fit_priors)[["p_yr"]]
mean_phi_est <- extract(fit_priors)[["alpha_phi"]]
mean_p_est <- extract(fit_priors)[["alpha_p"]]
hist(phi_est)
hist(p_est)
hist(boot::inv.logit(mean_phi_est))
hist(boot::inv.logit(mean_p_est))

     