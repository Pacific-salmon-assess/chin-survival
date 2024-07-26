### Example CJS survival models in stan
## Dec. 19, 2019

library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# SIM DATA ---------------------------------------------------------------------

# From rjags CJS example: 
# https://sites.google.com/site/uwrandbayesianmodeling/schedule/7-capture-mark-
# recapture-modeling/7-3-open-dynamic-cmr-models/7-3-2-cjs-model-in-jags-with-
# random-effects

# Functions to simulate data

#simulate CJS data for 1 group
simul.cjs <- function(phi, p, marked) {
  n.occasions <- length(p) + 1
  Phi <- matrix(phi, n.occasions - 1, nrow = sum(marked), byrow = T)
  P <- matrix(p, n.occasions - 1, nrow = sum(marked), byrow = T)
  
  CH <- matrix(0, ncol = n.occasions, nrow = sum(marked))
  #define a vector with marking occasion
  mark.occ <- rep(1:length(marked), marked[1:length(marked)])
  #fill in CH
  for (i in 1:sum(marked)) {
    CH[i, mark.occ[i]] <- 1
    if (mark.occ[i] == n.occasions) next
    for(t in (mark.occ[i] + 1):n.occasions) {
      #survive?
      sur <- rbinom(1, 1, Phi[i, t-1])
      if(sur == 0) break #move to next
      #recaptured?
      rp <- rbinom(1, 1, P[i, t - 1])
      if(rp == 1) CH[i, t] <- 1
    } #t
  } #i
  return(CH)
}

#function to get occasion of marking for each animal
get.first <- function(x) min(which(x != 0))

#function to create capture history character strings (only need this for RMark 
# runs)
pasty <- function(x) {
  k <- ncol(x)
  n <- nrow(x)
  out <- array(dim = n)
  for (i in 1:n) {
    out[i] <- paste(x[i, ], collapse = "")
  }
  return(out)
}

#function to add to data known states (where we know z=1)
known.states.cjs <- function(ch) {
  state <- ch
  for (i in 1:dim(ch)[1]) {
    n1 <- min(which(ch[i, ] == 1))
    n2 <- max(which(ch[i, ] == 1))
    state[i, n1:n2] <- 1
    state[i, n1] <- NA
  }
  state[state == 0] <- NA
  return(state)
}

## initialize z states have to be NA where observed
cjs.init.z <- function(ch,f){
  for (i in 1:dim(ch)[1]){
    if (sum(ch[i, ]) == 1) next
    n2 <- max(which(ch[i, ] == 1))
    ch[i, f[i]:n2] <- NA
  }
  for (i in 1:dim(ch)[1]) { 
    ch[i, 1:f[i]] <- NA
  }
  return(ch)
}

## Simulate data
n.occas <- 6 #number of sampling events
marked <- rep(50, n.occas - 1) #number of marked individuals
mean.phi <- 0.65
var.phi <- 1

p.input <- rep(.4, n.occas - 1)
logit.phi <- rnorm(n.occas - 1, qlogis(mean.phi), var.phi^0.5)
#plogis is a built in inverse logit transformation
phi.input <- plogis(logit.phi)

CH <- simul.cjs(phi.input, p.input, marked)
f <- apply(CH, 1, get.first)

#subset to focus on first 50
trim_mat <- CH[1:50, ]

#model list to pass to stan
dat <- list(`T` = ncol(trim_mat),
            I = nrow(trim_mat),
            y = trim_mat)
mod1 <- stan_model(here::here("R", "survival_ex", "cjs_individual.stan"))

fitModel <- sampling(mod,
                     data = dat,
                     iter = 4000, chains = 4, cores = 1,
                     control = list(adapt_delta = 0.95, max_treedepth = 20))

