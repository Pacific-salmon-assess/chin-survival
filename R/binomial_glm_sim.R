## Simulation test beta GLMs
# July 5, 2024

library(tidyverse)
library(rethinking)


## SIMULATE FIXED EFFECTS ------------------------------------------------------

# represents Y -> L -> S

# Set seed for reproducibility
set.seed(123)

# Simulate predictor variable
n <- 1000  # Number of observations
x <- rnorm(n)  # Predictor variable

# Simulate binary outcome
beta0 <- -1  # Intercept
beta1 <- 2   # Slope
eta <- beta0 + beta1 * x
p <- 1 / (1 + exp(-(eta)))  # Probability of success
y <- rbinom(n, 1, p)  # Binary outcome

# Combine into a data frame
dat <- data.frame(x = x, y = y)

# Fit binomial GLM
model <- glm(y ~ x, data = dat, family = binomial)
summary(model)

# Fit rethinking alternative
dat_list <- list(
  y = dat$y,
  x = dat$x
)
m1 <- ulam(
  alist(
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta0 + beta1 * x,
    beta0 ~ dnorm( 0 , 2 ),
    beta1 ~ dnorm( 0 , 2)
  ) , data=dat_list , chains=4 , log_lik=TRUE )
precis( m1 , depth=2 )


m1 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list , chains=4 , log_lik=TRUE )
precis( m11.4 , depth=2 )


data(chimpanzees)
d <- chimpanzees
d$treatment <- 1 + d$prosoc_left + 2*d$condition
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment) )
