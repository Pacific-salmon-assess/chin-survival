## Simulation test beta GLMs
# July 5, 2024

library(tidyverse)
library(rethinking)


# Set seed for reproducibility
set.seed(123)

n <- 1000  # Number of observations
n_years <- 5
yr <- sample(1:n_years, n, replace = TRUE)  # group index


## EXAMPLE 1 -------------------------------------------------------------------

# represents Y -> L & Y + L -> S 
# model years as hierarchical intercepts that covary 

# create empty data frame
dat <- data.frame(
  yr = yr,
  x = rep(NA, n_years), # intercept to predict x as a function of yr
  beta_yr = rep(NA, n_years)
)

# Simulate predictor variable (length) as a function of categorical variable
# (year)
x_sigma <- 0.6 
# make covariance matrix for year effects
yr_rho <- 0.6 
mu_yr <- c(0, 0) #average effects both centered on zero
yr_sigmas <- c(0.4, 0.2)
Rho <- matrix(c(1, yr_rho, yr_rho, 1), nrow = 2)
Sigma <- diag(yr_sigmas) %*% Rho %*% diag(yr_sigmas) 
yr_ints <- MASS::mvrnorm(n_years, mu_yr, Sigma)
alpha_yr <- yr_ints[ , 1]
beta_yr <- yr_ints[ , 2]

# simulate covariate data
for (i in seq_along(dat$yr)) { 
  dat$x[i] <- rnorm(1, alpha_yr[yr[i]], x_sigma)
  dat$beta_yr[i] <- beta_yr[yr[i]]
}

# Simulate binary outcome based on int (beta), x and yr (beta_yr)
beta <- -1  # Intercept
beta_x <- 2   # Slope
eta <- beta + 
  dat$beta_yr +
  beta_x * dat$x #calculate eta
p <- 1 / (1 + exp(-(eta)))  # Probability of survival
y <- rbinom(n, 1, p)  # Binary outcome
dat$y <- y

# Fit binomial GLM
# dat$yr_f <- as.factor(dat$yr)
# model <- glm(y ~ x, data = dat, family = binomial)
# model <- glm(y ~ 0 + x + yr_f, data = dat, family = binomial)
# summary(model)

# Fit rethinking alternative
dat_list <- list(
  y = dat$y,
  x = dat$x,
  yr = dat$yr
)
m1 <- ulam(
  alist(
    # length
    x ~ dnorm(mu, sigma_x),
    mu <- alpha_bar + alpha_yr[yr]*sigma_yr2,
    alpha_yr[yr] ~ dnorm(0, 1),
    sigma_yr2 ~ exponential(1),
    
    # survival
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta_bar + beta_yr[yr]*sigma_yr + beta_x * x,
    
    # priors
    beta_yr[yr] ~ dnorm(0, 0.5),
    alpha_bar ~ normal(0, 2.5),
    beta_bar ~ normal(0, 2.5),
    beta_x ~ normal(0, 1),
    sigma_yr ~ exponential(1),
    sigma_x ~ exponential(1)
    ),
  data=dat_list, chains=4 , log_lik=TRUE,
  control = list(adapt_delta = 0.95)
  )
precis(m1 , depth=2)

# as above but accounting for covariance among year intercepts
m1b <- ulam(
  alist(
    # length
    x ~ dnorm(mu, sigma_x),
    mu <- alpha_bar + alpha_yr[yr, 1],
    # survival
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta_bar + alpha_yr[yr, 2] + beta_x * x,
    # adaptive priors
    transpars> matrix[yr,2]:alpha_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[2, yr]:z_yr ~ normal(0 , 1),
    # priors
    alpha_bar ~ normal(0, 1.5),
    beta_bar ~ normal(0, 1.5),
    beta_x ~ normal(0, 0.5),
    vector[2]:sigma_yr ~ exponential(1),
    cholesky_factor_corr[2]:L_Rho_yr ~ lkj_corr_cholesky(2),
    sigma_x ~ exponential(1),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[2,2]:Rho_yr <<- Chol_to_Corr(L_Rho_yr)
  ),
  data=dat_list, chains=4 , log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)
precis(m1b, depth=3)

# misspecified model that ignores length submodel
m1c <- ulam(
  alist(
    # survival
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta_bar + beta_yr[yr]*sigma_yr + beta_x * x,
    # priors
    beta_bar ~ normal(0, 2.5),
    beta_yr[yr] ~ dnorm(0, 0.5),
    beta_x ~ normal(0, 0.5),
    sigma_yr ~ exponential(1)
  ),
  data=dat_list, chains=4 , log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)
precis(m1c, depth = 2)


## EXAMPLE 2 -------------------------------------------------------------------

# as above but include capture date covariate as well 

# create empty data frame
dat <- data.frame(
  yr = yr,
  x = rep(NA, n_years), # intercept to predict x as a function of yr
  cap_date = rep(NA, n_years),
  beta_yr = rep(NA, n_years)
)

# Simulate predictor variable (length) as a function of categorical variable
# (year)
x_sigma <- 0.6 
# make covariance matrix for year effects
yr_rho <- 0.6 
mu_yr <- c(0, 0, 0) #average effects both centered on zero
yr_sigmas <- c(0.2, 0.2, 0.3)
Rho <- matrix(0.6, 3, 3)
diag(Rho) <- 1
Sigma <- diag(yr_sigmas) %*% Rho %*% diag(yr_sigmas) 
yr_ints <- MASS::mvrnorm(n_years, mu_yr, Sigma)
alpha_yr <- yr_ints[ , 1]
gamma_yr <- yr_ints[ , 2]
beta_yr <- yr_ints[ , 3]

# simulate covariate data
for (i in seq_along(dat$yr)) { 
  dat$x[i] <- rnorm(1, alpha_yr[yr[i]], x_sigma)
  dat$cap_date[i] <- rnorm(1, gamma_yr[yr[i]], x_sigma)
  dat$beta_yr[i] <- beta_yr[yr[i]]
}

# Simulate binary outcome based on int (beta), x and yr (beta_yr)
beta <- -1  # Intercept
beta_x <- 2   # Slope
beta_date <- -1
eta <- beta + 
  dat$beta_yr +
  beta_x * dat$x +
  beta_date * dat$cap_date
p <- 1 / (1 + exp(-(eta)))  # Probability of survival
y <- rbinom(n, 1, p)  # Binary outcome
dat$y <- y


# Fit binomial GLM
dat$yr_f <- as.factor(dat$yr)
model2 <- glm(y ~ 0 + x + cap_date + yr_f, data = dat, family = binomial)
summary(model2)

dat_list <- list(
  y = dat$y,
  x = dat$x,
  yr = dat$yr,
  cap_date = dat$cap_date
)
m2 <- ulam(
  alist(
    # length
    x ~ dnorm(mu_yr, sigma_x),
    mu_yr <- alpha_bar + alpha_yr[yr, 1],
    # date
    cap_date ~ dnorm(mu_date, sigma_date),
    mu_date <- gamma_bar + alpha_yr[yr, 2],
    # survival
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta_bar + alpha_yr[yr, 3] + beta_x * x + beta_date * cap_date,
    # adaptive priors
    transpars> matrix[yr,3]:alpha_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[3, yr]:z_yr ~ normal(0 , 1),
    # priors
    c(alpha_bar, beta_bar, gamma_bar) ~ normal(0, 1.5),
    c(beta_x, beta_date) ~ normal(0, 0.5),
    vector[3]:sigma_yr ~ exponential(1),
    cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
    c(sigma_x, sigma_date) ~ exponential(1),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr)
  ),
  data=dat_list, chains=4 , log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)
precis(m2 , depth=2)
