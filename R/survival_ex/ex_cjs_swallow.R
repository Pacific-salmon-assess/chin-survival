## Example CJS w/ covariates model from:
#Section 14.5 of Bayesian Data Analysis in Ecology Using Linear Models with R, 
#BUGS, and Stan 

library(blmeco)
library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# import survival data
data(survival_swallows)
datax <- survival_swallows


# ------------------------------------------------------------------------------
## Model fitting

swallows_mod <- stan_model(here::here("R", "stan_models", "CJS_swallows.stan"))
# fit_mod_sw <- sampling(swallows_mod, data = datax, chains = 4, iter = 5000) 
# saveRDS(fit_mod_sw, here::here("R", "survival_ex", "cjs_swallow_ex_fits.RDS"))
fit_mod_sw <- readRDS(here::here("R", "survival_ex", "cjs_swallow_ex_fits.RDS"))

print(fit_mod_sw, c("a", "a1", "b0", "b1", "sigmayearphi", "sigmaphi", "sigmap"))


# ------------------------------------------------------------------------------
## Model checks

# example of diagnostic plot
historyplot(fit_mod_sw, "a")

## Posterior predictive check
modsims <- extract(fit_mod_sw) # extract values from posterior distribution

# define arrays for replicated obs. (y) and replicated state var. (z)
nind <- datax$I
nobs <- datax$K
nsim <- dim(modsims$phi)[1] # extract the number of simulations
yrep <- array(dim = c(nsim, nind, nobs)) 
zrep <- array(dim = c(nsim, nind, nobs))

# specify that first capture occasion is always 1 (i.e. marking)
yrep[ , , 1] <- 1
zrep[ , , 1] <- 1

# simulates states based on previous time step and estimated phi, then obs based
# on states and estimated p
for(j in 1:(nobs - 1)) {
  zrep[ , , j + 1] <- rbinom(nsim * datax$I, size = zrep[ , , j], 
                             prob = modsims$phi[ , , j])
  yrep[ , , j + 1] <- rbinom(nsim * datax$I, size = zrep[ , , j + 1], 
                             prob = modsims$p[ , , j])
}

# compare simulated and observed by assess whether between family varaince 
# assumed in model fits data
# obs data
nobspind <- apply(datax$CH, 1, sum) # number of observations per individual
mpfam <- tapply(nobspind, datax$family, mean)
minpfam <- tapply(nobspind, datax$family, min)
maxpfam <- tapply(nobspind, datax$family, max)

# sim data
repnpind <- apply(yrep, c(1, 2), sum)
repmpfam <- matrix(nrow = nsim, ncol = datax$nfam)
repminpfam <- matrix(nrow = nsim, ncol = datax$nfam)
repmaxpfam <- matrix(nrow = nsim, ncol = datax$nfam)
for(f in 1:nsim) {
  repmpfam[f, ] <- tapply(repnpind[f, ], datax$family, mean)
  repminpfam[f, ] <- tapply(repnpind[f, ], datax$family, min)
  repmaxpfam[f, ] <- tapply(repnpind[f, ], datax$family, max)
}

mpfam2 <- sort(mpfam)
repmpfam2 <- apply(repmpfam[ , order(mpfam2)], 2, mean)

plot(mpfam2, type = "l")
lines(repmpfam2, col = "blue")


# ------------------------------------------------------------------------------
## Summarize posterior

# estimates of daily survival
ests <- plogis(apply(modsims$a, 2, mean))
ests.lwr <- plogis(apply(modsims$a, 2, quantile, prob = 0.025))
ests.upr <- plogis(apply(modsims$a, 2, quantile, prob = 0.975))

plot(ests, type = "l", col = "blue", ylim = c(min(ests.lwr), max(ests.upr)))
lines(ests.lwr, col = "blue", lty = 2)
lines(ests.upr, col = "blue", lty = 2)

# effect of parental care
ma <- apply(modsims$a[ , 1:12], 1, mean) # mean over the first 12 days
newdat <- data.frame(carez = seq(-2, 2, length = 100))
b <- c(mean(ma), mean(modsims$a1)) # model coefficients
Xmat <- model.matrix(~carez, data = newdat)
newdat$fit <- plogis(Xmat %*% b)
nsim <- nrow(modsims$a)
fitmat <- matrix(ncol = nsim, nrow = nrow(newdat))
for(i in 1:nsim) {
  fitmat[ , i] <- plogis(Xmat %*% c(ma[i], modsims$a1[1]))
}
newdat$lwr <- apply(fitmat, 1, quantile, prob = 0.025)
newdat$upr <- apply(fitmat, 1, quantile, prob = 0.975)

plot(fit ~ carez, type = "l", col = "blue", data = newdat,
     ylim = c(min(newdat$lwr), max(newdat$upr)))
lines(x = newdat$carez, y = newdat$lwr, col = "blue", lty = 2)
lines(x = newdat$carez, y = newdat$upr, col = "blue", lty = 2)
