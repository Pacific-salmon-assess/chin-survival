## Simulation test beta GLMs
# July 5, 2024

library(tidyverse)
library(rethinking)


rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


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


## EXAMPLE 3 -------------------------------------------------------------------

# simulate condition as a latent process that determines size and lipid content,
# condition is influenced by sampling date, and lipid content, size, and date
# all influence survival

set.seed(73)


samp_date <- rnorm(n)
alpha_date <- 0.7
cond <- rnorm(n, alpha_date * samp_date, 0.5)
beta_cond <- 0.6 #assume same effect on length and lipid for simplicity's sake
size <- rnorm(n, beta_cond * cond, 0.5)
lipid <- rnorm(n, beta_cond * cond, 0.5)
gamma <- -1  # Intercept
gamma_date <- 0.5   # Slope
gamma_size <- 0.25
gamma_lipid <- 2
eta <- gamma + 
  gamma_date * samp_date +
  gamma_size * size +
  gamma_lipid * lipid
p <- boot::inv.logit(eta) # Probability of survival
surv <- rbinom(n, 1, p)  # Binary outcome

dat_list <- list(
  samp_date = samp_date,
  size = size,
  lipid = lipid,
  surv = surv
  )

m3 <- ulam(
  alist(
    # survival
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + gamma_date * samp_date +
      gamma_size * size +
      gamma_lipid * lipid,
    # priors
    gamma ~ normal(0, 2.5),
    c(gamma_date, gamma_size, gamma_lipid) ~ normal(0, 0.5)
  ),
  data=dat_list, chains=4 , log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)

m3b <- ulam(
  alist(
    # covariance among size and lipid
    c(size, lipid) ~ multi_normal(c(mu_size, mu_lipid), Rho, Sigma),
    mu_size <- alpha_size + beta_ds * samp_date,
    mu_lipid <- alpha_lipid + beta_dl * samp_date,
    # survival
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + gamma_date * samp_date +
      gamma_size * size +
      gamma_lipid * lipid,
    # priors
    c(alpha_size, alpha_lipid) ~ normal(0, 0.2),
    c(beta_ds, beta_dl) ~ normal(0, 0.5),
    Rho ~ lkj_corr( 2 ),
    Sigma ~ exponential( 1 ),
    gamma ~ normal(0, 2.5),
    c(gamma_date, gamma_size, gamma_lipid) ~ normal(0, 0.5)
  ),
  data=dat_list, chains=4 , #log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)


# sim doesn't work with large number of variables so calc manually accounting 
# for covariance
samp_date_seq <- seq( from=-2 , to=2 , length.out=30 )

post <- extract.samples(m3b)
pred_mu_size <- purrr::map(
  samp_date_seq, ~ post$alpha_size + post$beta_ds * .x
)
pred_mu_lipid <- purrr::map(
  samp_date_seq, ~ post$alpha_lipid + post$beta_dl * .x
)
# generate posterior covariance matrix for each draw, combine with pred mu to 
# sample from mvrnorm, then iterate over exp var
sigma <- post$Sigma
rho <- post$Rho
S <- vector(mode = "list", length = nrow(sigma))
sim_surv <- sim_size <- sim_lipid <-  matrix(NA, nrow = nrow(sigma), 
                                      ncol = length(samp_date_seq))
for (j in seq_along(samp_date_seq)) {
  for (i in 1:nrow(sigma)) {
    S <- diag(sigma[i,]) %*% rho[i,,] %*% diag(sigma[i,])
    mu <- MASS::mvrnorm(
      n = 1, mu=c(pred_mu_size[[j]][i], pred_mu_lipid[[j]][i]),
      Sigma = S
    )
    sim_size[i, j] <- mu[1]
    sim_lipid[i, j] <- mu[2]
  }
  sim_eta <- as.numeric(post$gamma + post$gamma_date * samp_date_seq[j] +
    post$gamma_size * sim_size[ , j] +
    post$gamma_lipid * sim_lipid[ , j]
    )
  p <- boot::inv.logit(sim_eta) # Probability of survival
  sim_surv[ , j] <- rbinom(length(p), 1, p)
}
plot(samp_date_seq, colMeans(sim_surv) , ylim=c(0, 1) , type="l" ,
     xlab="date" , ylab="survival"  )


# parameter estimates are nearly identical but what happens to estimate of
# sampling day effect
post3 <- extract.samples(m3)
sim_surv3 <- matrix(NA, nrow = nrow(sigma), ncol = length(samp_date_seq))
for (j in seq_along(samp_date_seq)) {
  sim_eta1 <- as.numeric(post3$gamma + post3$gamma_date * samp_date_seq[j])
  p1 <- boot::inv.logit(sim_eta1) # Probability of survival
  sim_surv3[ , j] <- rbinom(length(p1), 1, p1)
}
plot(samp_date_seq , colMeans(sim_surv3) , ylim=c(0, 1) , type="l" ,
     xlab="date" , ylab="survival"  )



## EXAMPLE 4 -------------------------------------------------------------------

# as above but include year effects on all covariates

set.seed(73)

## correlated annual intercepts for condition, date, survival
Rho <- diag(3)
Rho[upper.tri(Rho)] <- c(0.2, 0.7, 0.8)
Rho[lower.tri(Rho)] <- t(Rho)[lower.tri(Rho)]
mu_yr <- c(0, 0, 0) #average effects both centered on zero
yr_sigmas <- c(0.4, 0.2, 0.3)
Sigma <- diag(yr_sigmas) %*% Rho %*% diag(yr_sigmas) 
yr_ints <- MASS::mvrnorm(n_years, mu_yr, Sigma)
alpha_yd <- yr_ints[ , 1]
alpha_yc <- yr_ints[ , 2]
alpha_ys <- yr_ints[ , 3]

dat <- data.frame(
  yr = yr
) %>% 
  mutate(
    samp_date = NA,
    cond = NA,
    size = NA,
    lipid = NA
  )

# simulate covariate data
beta_dc <- 0.7
beta_cl <- beta_cf <- 0.6
for (i in seq_along(dat$yr)) { 
  dat$samp_date[i] <- rnorm(1, alpha_yd[yr[i]], 1)
  dat$cond[i] <- rnorm(1, alpha_yc[yr[i]] + beta_dc * dat$samp_date[i], 1)
  dat$size[i] <- rnorm(1, beta_cf * dat$cond[i], 1)
  dat$lipid[i] <- rnorm(1, beta_cl * dat$cond[i], 1)
  dat$alpha_ys[i] <- alpha_ys[yr[i]]
}

surv_bar <- -1  # Intercept
beta_ds <- 0.5   # Slope
beta_fs <- 0.25
beta_ls <- 2
eta <- surv_bar + dat$alpha_ys +
  beta_ds * dat$samp_date +
  beta_fs * dat$size +
  beta_ls * dat$lipid
p <- boot::inv.logit(eta) # Probability of survival
dat$surv <- rbinom(n, 1, p)  # Binary outcome

dat_list <- list(
  samp_date = dat$samp_date,
  yr = dat$yr, 
  size = dat$size,
  lipid = dat$lipid,
  surv = dat$surv
)

m4 <- ulam(
  alist(
    # date
    samp_date ~ dnorm(mu_date, sigma_date),
    mu_date <- alpha_yr[yr, 1],
    # covariance among size and lipid
    c(size, lipid) ~ multi_normal(c(mu_size, mu_lipid), Rho_sl, Sigma_sl),
    mu_size <- alpha_yr[yr, 2] + beta_df * samp_date,
    mu_lipid <- alpha_yr[yr, 3] + beta_dl * samp_date,
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + alpha_yr[yr, 4] +
      beta_ds * samp_date +
      beta_fs * size +
      beta_ls * lipid,
    # adaptive priors
    transpars> matrix[yr, 4]:alpha_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[4, yr]:z_yr ~ normal(0 , 1),
    # priors
    c(beta_df, beta_dl) ~ normal(0, 1),
    c(beta_ds, beta_fs, beta_ls) ~ normal(0, 0.5),
    c(surv_bar) ~ normal(0, 1.5),
    Rho_sl ~ lkj_corr( 2 ),
    cholesky_factor_corr[4]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[4]:sigma_yr ~ exponential(1),
    c(Sigma_sl, sigma_date) ~ exponential(1),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4, 4]:Rho_yr <<- Chol_to_Corr(L_Rho_yr)
  ),
  data=dat_list, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)
precis(m4, depth = 2)

# as above but ignores covariance between lipid and size
m4b <- ulam(
  alist(
    # date
    samp_date ~ dnorm(mu_date, sigma_date),
    mu_date <- alpha_yr[yr, 1],
    # covariance among size and lipid
    size ~ dnorm(mu_size, sigma_size),
    mu_size <- alpha_yr[yr, 2] + beta_df * samp_date,
    lipid ~ dnorm(mu_lipid, sigma_lipid),
    mu_lipid <- alpha_yr[yr, 3] + beta_dl * samp_date,
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + alpha_yr[yr, 4] +
      beta_ds * samp_date +
      beta_fs * size +
      beta_ls * lipid,
    # adaptive priors
    transpars> matrix[yr, 4]:alpha_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[4, yr]:z_yr ~ normal(0 , 1),
    # priors
    c(beta_df, beta_dl) ~ normal(0, 1),
    c(beta_ds, beta_fs, beta_ls) ~ normal(0, 0.5),
    c(surv_bar) ~ normal(0, 1.5),
    cholesky_factor_corr[4]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[4]:sigma_yr ~ exponential(1),
    c(sigma_date, sigma_size, sigma_lipid) ~ exponential(1),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4, 4]:Rho_yr <<- Chol_to_Corr(L_Rho_yr)
  ),
  data=dat_list, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)
precis(m4b, depth = 2)

# sim doesn't work with large number of variables so calc manually accounting 
# for covariance (excludes year effects)
samp_date_seq <- seq(-2, 2, length.out = 30)
post <- extract.samples(m4)
pred_mu_size <- purrr::map(
  samp_date_seq, ~ post$size_bar + post$beta_ds * .x
)
pred_mu_lipid <- purrr::map(
  samp_date_seq, ~ post$lipid_bar + post$beta_dl * .x
)
# generate posterior covariance matrix for each draw, combine with pred mu to 
# sample from mvrnorm, then iterate over exp var
sigma <- post$Sigma
rho <- post$Rho_sl
S <- vector(mode = "list", length = nrow(sigma))
sim_surv <- sim_size <- sim_lipid <- matrix(NA, nrow = nrow(sigma), 
                                            ncol = length(samp_date_seq))
for (j in seq_along(samp_date_seq)) {
  for (i in 1:nrow(sigma)) {
    S <- diag(sigma[i,]) %*% rho[i,,] %*% diag(sigma[i,])
    mu <- MASS::mvrnorm(
      n = 1, mu=c(pred_mu_size[[j]][i], pred_mu_lipid[[j]][i]),
      Sigma = S
    )
    sim_size[i, j] <- mu[1]
    sim_lipid[i, j] <- mu[2]
  }
  sim_eta <- as.numeric(post$surv_bar + post$beta_ds * samp_date_seq[j] +
                          post$beta_fs * sim_size[ , j] +
                          post$beta_ls * sim_lipid[ , j]
  )
  p <- boot::inv.logit(sim_eta) # Probability of survival
  sim_surv[ , j] <- rbinom(length(p), 1, p)
}
plot(samp_date_seq, colMeans(sim_surv) , ylim=c(0, 1) , type="l" ,
     xlab="date" , ylab="survival"  )


## EXAMPLE 5 -------------------------------------------------------------------


# as example 3 but includes tagging injuries following Dirichlet distribution

set.seed(73)

# injuries range from 1 to 4
n_inj <- 4
inj <- sample(0:3, n, replace = TRUE)  # group index
gamma_inj <- -1.5
delta_inj <- c(0.1, 0.2, 0.7)
# injury = sum of each level * gamma (i.e. level 4 = full effect)
inj_eff <- rep(NA, times = n)
for (i in 1:length(inj)) {
  inj_eff[i] <- ifelse(
    inj[i] == 0,
    0,
    gamma_inj * sum(delta_inj[1:inj[i]])
  )
}

samp_date <- rnorm(n)
alpha_date <- 0.7
cond <- rnorm(n, alpha_date * samp_date, 1)
beta_cond <- 0.6 #assume same effect on length and lipid for simplicity's sake
size <- rnorm(n, beta_cond * cond, 1)
lipid <- rnorm(n, beta_cond * cond, 1)
gamma <- 1  # Intercept
gamma_date <- 0.5   # Slope
gamma_size <- 0.25
gamma_lipid <- 2
eta <- gamma + 
  gamma_date * samp_date +
  gamma_size * size +
  gamma_lipid * lipid +
  inj_eff
p <- boot::inv.logit(eta) # Probability of survival
surv <- rbinom(n, 1, p)  # Binary outcome

dat_list <- list(
  samp_date = samp_date,
  size = size,
  lipid = lipid,
  surv = surv,
  inj = inj,
  alpha = rep(2, n_inj - 1)
)


m5 <- ulam(
  alist(
    # covariance among size and lipid
    c(size, lipid) ~ multi_normal(c(mu_size, mu_lipid), Rho, Sigma),
    mu_size <- alpha_size + beta_ds * samp_date,
    mu_lipid <- alpha_lipid + beta_dl * samp_date,
    # survival
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + gamma_date * samp_date +
      gamma_size * size +
      gamma_lipid * lipid +
      gamma_inj * sum(delta_inj[1:inj])
    ,
    # priors
    c(alpha_size, alpha_lipid) ~ normal(0, 0.2),
    c(beta_ds, beta_dl) ~ normal(0, 0.5),
    Rho ~ lkj_corr( 2 ),
    Sigma ~ exponential( 1 ),
    gamma ~ normal(0, 2.5),
    gamma_inj ~ normal(0, 0.5),
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha),
    c(gamma_date, gamma_size, gamma_lipid) ~ normal(0, 0.5)
  ),
  data=dat_list, chains=4 , cores = 4,#log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)

# note order of estimates depends on order in which inj index appears
precis(m5, depth = 2)


## EXAMPLE 6 -------------------------------------------------------------------

# Model interaction between capture date, exploitation rate and survival

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

# choose 10 plausible exploitation rate vals then sample to length dat
cyer_vals <- runif(10, 0.3, 0.8)
cyer_vec <- sample(cyer_vals, size = nrow(dat), replace = TRUE)
dat$cyer <- scale(cyer_vec) %>% as.numeric()


# Simulate binary outcome based on int (beta), x and yr (beta_yr)
# define cyer to have negative imapct on survival and negative interaction
# with date (i.e. effect of CYER weaker late in year)
beta <- -1  # Intercept
beta_x <- 0.5 # Slope
beta_cyer <- -1
beta_x_cyer <- 0.5

eta <- beta + 
  dat$beta_yr +
  beta_x * dat$x +
  beta_cyer * dat$cyer +
  (beta_x_cyer * dat$x * dat$cyer)
  
p <- 1 / (1 + exp(-(eta)))  # Probability of survival
y <- rbinom(n, 1, p)  # Binary outcome
dat$y <- y

dat_list <- list(
  y = dat$y,
  x = dat$x,
  yr = dat$yr,
  cyer = dat$cyer
)

m6 <- ulam(
  alist(
    # length
    x ~ dnorm(mu, sigma_x),
    mu <- alpha_bar + alpha_yr[yr]*sigma_yr2,
    alpha_yr[yr] ~ dnorm(0, 1),
    sigma_yr2 ~ exponential(1),
    
    # survival
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta_bar + beta_yr[yr]*sigma_yr + beta_x * x + 
      beta_cyer * cyer + (beta_x_cyer * x * cyer),
    
    # priors
    beta_yr[yr] ~ dnorm(0, 0.5),
    alpha_bar ~ normal(0, 1.25),
    beta_bar ~ normal(0, 1.25),
    c(beta_x, beta_cyer, beta_x_cyer) ~ normal(0, 0.5),
    sigma_yr ~ exponential(1),
    sigma_x ~ exponential(1)
  ),
  data=dat_list, chains=4 , log_lik=TRUE, cores = 4,
  control = list(adapt_delta = 0.95)
)
precis(m6)

# equivalent model with informative prior
m6b <- ulam(
  alist(
    # length
    x ~ dnorm(mu, sigma_x),
    mu <- alpha_bar + alpha_yr[yr]*sigma_yr2,
    alpha_yr[yr] ~ dnorm(0, 1),
    sigma_yr2 ~ exponential(1),
    
    # survival
    y ~ dbinom( 1 , p ) ,
    logit(p) <- beta_bar + beta_yr[yr]*sigma_yr + beta_x * x + 
      beta_cyer * cyer + (beta_x_cyer * x * cyer),
    
    # priors
    beta_yr[yr] ~ dnorm(0, 0.5),
    alpha_bar ~ normal(0, 1.25),
    beta_bar ~ normal(0, 1.25),
    c(beta_x, beta_x_cyer) ~ normal(0, 0.5),
    beta_cyer ~ normal(-0.5, 0.5),
    sigma_yr ~ exponential(1),
    sigma_x ~ exponential(1)
  ),
  data=dat_list, chains=4 , log_lik=TRUE, cores = 4,
  control = list(adapt_delta = 0.95)
)
precis(m6b , depth=2)

post <- extract.samples(m6)
prior <- extract.prior(m6)

x_seq <- c(-2, 0, 2)
cyer_seq <- seq(-1.3, 1.3, length = 40)
preds <- vector(mode = "list", length = length(x_seq))

for(i in seq_along(x_seq)) {
  pred_x <- sapply(
    cyer_seq, 
    function (x) {
      inv_logit(
        post$beta_bar + post$beta_x * x_seq[i] + post$beta_cyer * x +
                  (post$beta_x_cyer * x_seq[i] * x)
                )
    }
  )
  preds[[i]] <- pred_x[1:20, ] %>% 
    as.data.frame() %>% 
    set_names(cyer_seq) %>% 
    mutate(iter = seq(1, nrow(.), by = 1)) %>% 
    pivot_longer(
      cols = -iter, names_to = "cyer", values_to = "est"
    ) %>% 
    mutate(
      cyer = as.numeric(cyer),
      x = x_seq[i]
    )
}

bind_rows(preds) %>% 
  ggplot(
    ., aes(x = cyer, y = est, group = iter)
  ) +
  geom_line(
  ) +
  ggsidekick::theme_sleek() +
  facet_wrap(~x)


## EXAMPLE 7 -------------------------------------------------------------------

library(glmmTMB)

set.seed(73)

pop_int <- data.frame(
  pop = c("a", "b", "c", "d", "e")
) %>% 
  mutate(
    alpha = rnorm(length(unique(pop)), 0, 0.5)
  )
surv_bar <- -1  # Intercept
beta_qs <- 0.7 # detection probability effect


# detection probability qq
det_key <- expand.grid(
  # yr = seq(1, 5, by = 1),
  pop = c("a", "b", "c", "d", "e")
) %>% 
  mutate(
    pop_n = as.numeric(as.factor(pop)),
    qq = case_when(
      pop_n < 3 ~ 0,
      TRUE ~ 1
    ), 
    qq_alpha = ifelse(qq == "0", 0, beta_qs)
  ) 

dat <- data.frame(
  pop = sample(c("a", "b", "c", "d", "e"), size = 1000, replace = T)
) %>% 
  left_join(., det_key, by = "pop") %>% 
  left_join(., pop_int, by = "pop") %>% 
  mutate(
    #alternative random detection probability
    qq = as.factor(qq),
    qq2 = sample(c(0, 1), size = 1000, replace = T) %>% 
      as.factor(),
    qq2_alpha = ifelse(qq2 == "0", 0, beta_qs)
  )

nsim <- 100
sim_list <- vector(mode = "list", length = 100)


for (i in 1:30) {
  eta <- surv_bar + dat$qq_alpha + dat$alpha
  p <- boot::inv.logit(eta) # Probability of survival
  dat$surv <- rbinom(nrow(dat), 1, p) 
  
  # as above but for random det prob effect
  eta2 <-  surv_bar + dat$qq2_alpha + dat$alpha
  p2 <- boot::inv.logit(eta2) # Probability of survival
  dat$surv2 <- rbinom(nrow(dat), 1, p2) 
  
  
  fit1 <- glmmTMB(
    surv ~ qq + (1 | pop),
    data = dat,
    family = binomial()
  )
  fes <- fixef(fit1)$cond 
  ris <- ranef(fit1)$cond %>% as.data.frame()
  d1 <- data.frame(
    par = c("int", "qq", paste("ri_", seq(1, 5, by = 1), sep = "")),
    est = c(fes, ris[ , 1])
  ) %>% 
    mutate(
      model = "1"
    )
  fit1b <- glmmTMB(
    surv ~ qq, #+ (1 | pop),
    data = dat,
    family = binomial()
  )
  fes1b <- fixef(fit1b)$cond 
  d1b <- data.frame(
    par = c("int", "qq"),
    est = c(fes1b)
  ) %>% 
    mutate(
      model = "1b"
    )
  
  
  fit2 <- glmmTMB(
    surv2 ~ qq2 + (1 | pop),
    data = dat,
    family = binomial()
  )
  fes2 <- fixef(fit2)$cond 
  ris2 <- ranef(fit2)$cond %>% as.data.frame()
  d2 <- data.frame(
    par = c("int", "qq", paste("ri_", seq(1, 5, by = 1), sep = "")),
    est = c(fes2, ris2[ , 1])
  ) %>% 
    mutate(
      model = "2"
    )
  fit2b <- glmmTMB(
    surv2 ~ qq2, #+ (1 | pop),
    data = dat,
    family = binomial()
  )
  fes2b <- fixef(fit2b)$cond 
  d2b <- data.frame(
    par = c("int", "qq"),
    est = c(fes2b)
  ) %>% 
    mutate(
      model = "2b"
    )
  
  sim_list[[i]] <- list(d1, d1b, d2, d2b) %>% bind_rows()
}

sim_dat <- sim_list[1:20] %>% bind_rows()
true_dat <- data.frame(
  par = c("int", "qq", paste("ri_", seq(1, 5, by = 1), sep = "")),
  est = c(-1, 0.7, pop_int$alpha)
) 

ggplot() +
  geom_boxplot(
    data = sim_dat, aes(x = par, y = est)
  ) +
  geom_point(
    data = true_dat, aes(x = par, y = est), color = "red",
  ) +
  facet_wrap(
    ~ model
  )
# parameter estimates are identical between b and a, estimates are marginally
# more accurate in 2 than 1, but effects are modest


## Bayesian version of above to answer whether including submodel changes
# parameter estimates or not
dat_list <- list(
  surv = dat$surv,
  pop_n = dat$pop_n,
  qq = as.numeric(dat$qq) - 1,
  qq_n = as.numeric(dat$qq)
)

m8a <- ulam(
  alist(
    # detection probability
    qq ~ dbinom(1, p_qq),
    logit(p_qq) <- alpha_qq_pop[pop_n],
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- alpha_pop[pop_n] * sigma_pop + alpha_qq[qq_n],
    # adaptive priors
    # priors
    alpha_pop[pop_n] ~ dnorm(0, 0.5),
    alpha_qq_pop[pop_n] ~ dnorm(0, 0.5),
    alpha_qq[qq_n] ~ normal(0, 1),
    sigma_pop ~ exponential(1)
  ),
  data=dat_list, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)
m8b <- ulam(
  alist(
    surv ~ dbinom(1 , p) ,
    logit(p) <- alpha_pop[pop_n] * sigma_pop + alpha_qq[qq_n],
    # priors
    alpha_pop[pop_n] ~ dnorm(0, 0.5),
    alpha_qq[qq_n] ~ normal(0, 1),
    sigma_pop ~ exponential(1)
  ),
  data=dat_list, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)

precis(m8a, depth = 2)
precis(m8b, depth = 2)
# inclusion of submodel doesn't impact estimates



# model where stock is also present as random intercept, detection probability
# influences apparent survival,  and is confounded with stock/year; otherwise 
# same as model 4

set.seed(73)

pop <- sample(c("a", "b", "c", "d", "e"), size = 1000, replace = T)
pop_int <- data.frame(
  pop = unique(pop) %>% sort()
) %>% 
  mutate(
    alpha = rnorm(length(unique(pop)), 0, 0.5)
  )


## correlated annual intercepts for condition, date, survival
Rho <- diag(3)
Rho[upper.tri(Rho)] <- c(0.2, 0.7, 0.8)
Rho[lower.tri(Rho)] <- t(Rho)[lower.tri(Rho)]
mu_yr <- c(0, 0, 0) #average effects both centered on zero
yr_sigmas <- c(0.4, 0.2, 0.3)
Sigma <- diag(yr_sigmas) %*% Rho %*% diag(yr_sigmas)
yr_ints <- MASS::mvrnorm(n_years, mu_yr, Sigma)
alpha_yd <- yr_ints[ , 1]
alpha_yc <- yr_ints[ , 2]
alpha_ys <- yr_ints[ , 3]

# detection probability qq
det_key <- expand.grid(
  # yr = seq(1, 5, by = 1),
  pop = unique(pop)
) %>% 
  mutate(
    pop_n = as.numeric(as.factor(pop)),
    qq = case_when(
      yr < 3 & pop_n < 3 ~ 0,
      yr == "4" & pop_n == "3" ~ 0,
      TRUE ~ 1
    )
  ) 

dat <- data.frame(
  yr = yr,
  pop = pop
) %>%
  mutate(
    samp_date = NA,
    cond = NA,
    size = NA,
    lipid = NA
  ) %>%
  left_join(
    ., det_key, by = c("yr", "pop")
  )

# simulate covariate data
beta_dc <- 0.7
beta_cl <- beta_cf <- 0.6
for (i in seq_along(dat$yr)) { 
  dat$samp_date[i] <- rnorm(1, alpha_yd[yr[i]], 1)
  dat$cond[i] <- rnorm(1, alpha_yc[yr[i]] + beta_dc * dat$samp_date[i], 1)
  dat$size[i] <- rnorm(1, beta_cf * dat$cond[i], 1)
  dat$lipid[i] <- rnorm(1, beta_cl * dat$cond[i], 1)
  dat$alpha_ys[i] <- alpha_ys[yr[i]]
}

# random intercepts for pop survival 
alpha_ps <- rnorm(length(unique(pop)), 0, 1)
for (i in seq_along(dat$pop)) {
  dat$alpha_ps[i] <- alpha_ps[dat$pop_n[i]]
}


surv_bar <- -1  # Intercept
beta_ds <- 0.5   # date effect
beta_fs <- 0.25 # fork length
beta_ls <- 1 # lipid
beta_qs <- 0.5 # detection probability effect
eta <- surv_bar + dat$alpha_ys + dat$alpha_ps +
  beta_ds * dat$samp_date +
  beta_fs * dat$size +
  beta_ls * dat$lipid + 
  beta_qs * dat$qq
p <- boot::inv.logit(eta) # Probability of survival
dat$surv <- rbinom(n, 1, p)  # Binary outcome

dat_list <- list(
  samp_date = dat$samp_date,
  yr = dat$yr, 
  size = dat$size,
  lipid = dat$lipid,
  surv = dat$surv,
  pop_n = dat$pop_n,
  qq = dat$qq
)

m8 <- ulam(
  alist(
    # date
    samp_date ~ dnorm(mu_date, sigma_date),
    mu_date <- alpha_yr[yr, 1],
    # covariance among size and lipid
    c(size, lipid) ~ multi_normal(c(mu_size, mu_lipid), Rho_sl, Sigma_sl),
    mu_size <- alpha_yr[yr, 2] + beta_df * samp_date,
    mu_lipid <- alpha_yr[yr, 3] + beta_dl * samp_date,
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + alpha_yr[yr, 4] + alpha_pop[pop_n] * sigma_pop +
      beta_ds * samp_date +
      beta_fs * size +
      beta_ls * lipid +
      beta_qs * qq,
    # adaptive priors
    transpars> matrix[yr, 4]:alpha_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[4, yr]:z_yr ~ normal(0 , 1),
    # priors
    alpha_pop[pop_n] ~ dnorm(0, 1),
    c(beta_df, beta_dl) ~ normal(0, 1),
    c(beta_ds, beta_fs, beta_ls, beta_qs) ~ normal(0, 0.5),
    c(surv_bar) ~ normal(0, 1.5),
    Rho_sl ~ lkj_corr( 2 ),
    cholesky_factor_corr[4]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[4]:sigma_yr ~ exponential(1),
    c(Sigma_sl, sigma_date, sigma_pop) ~ exponential(1),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4, 4]:Rho_yr <<- Chol_to_Corr(L_Rho_yr)
  ),
  data=dat_list, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)
precis(m8, depth = 2)


mu <- alpha_bar + alpha_yr[yr]*sigma_yr2
alpha_yr[yr] ~ dnorm(0, 1)
sigma_yr2 ~ exponential(1)

