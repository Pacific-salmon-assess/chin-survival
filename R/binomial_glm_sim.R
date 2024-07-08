## Simulation test beta GLMs
# July 5, 2024

library(tidyverse)
library(rethinking)


## SIMULATE HIERARCHICAL EFFECTS -----------------------------------------------

# represents Y -> L -> S

# Set seed for reproducibility
set.seed(123)

n <- 1000  # Number of observations
n_years <- 5
yr <- sample(1:n_years, n, replace = TRUE)  # group index


# Simulate predictor variable (length) as a function of categorical variable
# (year)
alpha_x_yr <- rnorm(yr, 0, 2.5)
sigma_x <- 1
x <- rnorm(n, alpha_x_yr[yr], sigma_x)

# Combine into a data frame
dat <- data.frame(x = x, yr = yr, beta_yr = NA)

# Simulate binary outcome
beta <- -1  # Intercept
beta_x <- 2   # Slope
beta_yr_sigma <- 1
beta_yr <- rnorm(n_years, 0, beta_yr_sigma) # yr intercept
for (i in seq_along(yr)) { # add to dataframe
  dat$beta_yr[i] <- beta_yr[yr[i]]
}
eta <- beta + dat$beta_yr + beta_x * x #calculate eta
p <- 1 / (1 + exp(-(eta)))  # Probability of survival
y <- rbinom(n, 1, p)  # Binary outcome
dat$y <- y

# Fit binomial GLM
dat$yr_f <- as.factor(dat$yr)
# model <- glm(y ~ x, data = dat, family = binomial)
model <- glm(y ~ 0 + x + yr_f, data = dat, family = binomial)
summary(model)

# Fit rethinking alternative
dat_list <- list(
  y = y,
  x = x,
  yr = yr
)
m1 <- ulam(
  alist(
    # length
    x ~ dnorm(mu, sigma_x), 
    mu <- alpha_x_yr[yr],
    alpha_x_yr[yr] ~ dnorm(0, 2.5),
    sigma_x ~ dexp(1),
    # survival
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta_yr[yr] + beta_x * x,
    beta_yr[yr] ~ dnorm(beta_yr_bar, beta_yr_sigma),
    beta_yr_bar ~ dnorm(0, 1.5),
    beta_yr_sigma ~ dexp(1),
    beta_x ~ dnorm(0, 0.5)
  ),
  data=dat_list, chains=4 , log_lik=TRUE 
  )
precis(m1 , depth=2)

prior <- extract.prior( m1 , n=1e4 )
p <- inv_logit( prior$beta + prior$beta_x )
dens( p , adj=0.1 )


# as above but include capture date covariate as well 
alpha_x_yr <- rnorm(yr, 0, 2.5)
sigma_x <- 1
samp_date <- rnorm(yr, 0, 1) #sampling dates
alpha_x_date <- -0.5
eta_x <- alpha_x_yr[yr] + samp_date * alpha_x_date
x <- rnorm(n, eta_x, sigma_x)

# Combine into a data frame
dat <- data.frame(x = x, yr = yr, samp_date = samp_date, beta_yr = NA)

# Simulate binary outcome
beta <- -1  # Intercept
beta_x <- 2   # Slope
beta_date <- 0.5
beta_yr_sigma <- 1
beta_yr <- rnorm(n_years, 0, beta_yr_sigma) # yr intercept
for (i in seq_along(yr)) { # add to dataframe
  dat$beta_yr[i] <- beta_yr[yr[i]]
}
#calculate eta
eta <- beta + dat$beta_yr + beta_x * x + beta_date * samp_date
p <- 1 / (1 + exp(-(eta)))  # Probability of survival
y <- rbinom(n, 1, p)  # Binary outcome
dat$y <- y


# Fit binomial GLM
dat$yr_f <- as.factor(dat$yr)
model2 <- glm(y ~ 0 + x + samp_date + yr_f, data = dat, family = binomial)
summary(model2)

dat_list <- list(
  y = y,
  x = x,
  yr = yr,
  samp_date = samp_date
)
m2 <- ulam(
  alist(
    # length
    x ~ dnorm(mu, sigma_x), 
    mu <- alpha_x_yr[yr] + alpha_x_date * samp_date,
    alpha_x_yr[yr] ~ dnorm(0, 2.5),
    alpha_x_date ~ dnorm(0, 0.5),
    sigma_x ~ dexp(1),
    # survival
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta_yr[yr] + beta_x * x + beta_date * samp_date,
    beta_yr[yr] ~ dnorm(beta_yr_bar, beta_yr_sigma),
    beta_yr_bar ~ dnorm(0, 1.5),
    beta_yr_sigma ~ dexp(1),
    beta_x ~ dnorm(0, 0.5),
    beta_date ~ dnorm(0, 0.5)
  ),
  data=dat_list, chains=4 , log_lik=TRUE 
)
precis(m2 , depth=2)
