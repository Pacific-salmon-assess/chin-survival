## Logistic survival
# Survival to detection and to terminal arrays using logistic regression (i.e.
# not stage-specific, does not account for detection probability)
# Excludes immature tags and fish with unknown stock ID, but does not exclude
# redeployed tags or injured fish
# Feb. 22, 2021
# Updated July 4, 2022

library(tidyverse)
library(rethinking)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



det_dat1 <- readRDS(here::here("data", "surv_log_reg_data.rds"))


# remove redeployed tags and scale continuous covariates
det_dat <- det_dat1 %>% 
  filter(
    !redeploy == "yes"
  ) %>% 
  mutate(
    lipid_z = scale(lipid) %>% as.numeric(),
    fl_z = scale(fl) %>% as.numeric(),
    day_z = scale(year_day) %>% as.numeric(),
    cyer_z = scale(isbm_cyer) %>% as.numeric(),
    terminal_p = as.factor(det_dat$terminal_p),
    year = as.factor(det_dat$year),
    stock_group = as.factor(det_dat$stock_group)
  )


# PLOTS OF RAW DATA ------------------------------------------------------------

ppn_foo <- function(group, response) {
  group_exp <- c("vemco_code", group)
  labs <- det_dat %>% 
    dplyr::select_at(group_exp) %>% 
    distinct() %>% 
    group_by_at(group) %>% 
    tally() %>% 
    mutate(year = as.factor(year))
  
  dat_out <- det_dat %>% 
    group_by_at(group) %>% 
    summarize(n = length(unique(vemco_code)),
              ppn = sum(.data[[response]] / n),
              se = sqrt((ppn * (1 - ppn)) / n),
              up = pmin(1, ppn + (1.96 * se)),
              lo = pmax(0, ppn - (1.96 * se)),
              .groups = "drop") %>% 
    mutate(year = as.factor(year))
  
  list("labs" = labs, "dat" = dat_out)
}


det_ppns <- ppn_foo(group = c("stock_group", "year"), response = "term_det") 
ggplot(data = det_ppns$dat, aes(y = ppn)) +
  geom_pointrange(aes(x = year, ymin = lo, ymax = up)) +
  facet_wrap(~stock_group) +
  ggsidekick::theme_sleek() +
  geom_text(data = det_ppns$lab, aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

inj_ppns <- ppn_foo(group = c("injury", "year"), response = "term_det")
ggplot(data = inj_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~injury) +
  ggsidekick::theme_sleek()+
  geom_text(data = inj_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

inj_ppns2 <- ppn_foo(group = c("comp_inj", "year"), response = "term_det")
ggplot(data = inj_ppns2$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~comp_inj) +
  ggsidekick::theme_sleek()+
  geom_text(data = inj_ppns2$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)


loc_ppns <- ppn_foo(group = c("hook_loc", "year"), response = "term_det")
ggplot(data = loc_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~hook_loc) +
  ggsidekick::theme_sleek()+
  geom_text(data = loc_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

fin_ppns <- ppn_foo(group = c("fin_dam", "year"), response = "term_det")
ggplot(data = fin_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~fin_dam) +
  ggsidekick::theme_sleek()+
  geom_text(data = fin_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

scale_ppns <- ppn_foo(group = c("scale_loss", "year"), response = "term_det")
ggplot(data = scale_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~scale_loss) +
  ggsidekick::theme_sleek()+
  geom_text(data = scale_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

redep_ppns <- ppn_foo(group = c("redeploy", "year"), response = "det")
ggplot(data = redep_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~redeploy) +
  ggsidekick::theme_sleek()+
  geom_text(data = redep_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)


ggplot(det_dat, aes(x = lipid_z, y = final_det)) +
  geom_point() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()
ggplot(det_dat, aes(x = trough_time, y = det)) +
  geom_point() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()
ggplot(det_dat, aes(x = mean_log_e, y = final_det)) +
  geom_point() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()
ggplot(det_dat, aes(x = cyer_z, y = final_det)) +
  geom_point() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()



# FIT LOGISTIC REG MODELS ------------------------------------------------------

dat_list <- list(
  surv = as.integer(det_dat$term_det),
  fl_z = det_dat$fl_z,
  lipid_z = det_dat$lipid_z,
  day_z = det_dat$day_z,
  # cyer_z = det_dat$cyer_z,
  inj = as.integer(det_dat$injury),
  term_p = as.integer(det_dat$terminal_p),
  yr = as.integer(det_dat$year),
  stk = as.integer(det_dat$stock_group),
  alpha = rep(2, length(unique(det_dat$injury)) - 1)
)


# multiple regression hierarchical model with injury effects 
m1 <- ulam(
  alist(
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_term[term_p] +
      beta_yr[yr] * sigma_yr +
      beta_stk[stk] * sigma_stk +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z +
      beta_is * sum(delta_inj[1:inj])
      ,
    surv_bar ~ normal(0, 1.25),
    beta_yr[yr] ~ normal(0, 0.5),
    beta_stk[stk] ~ normal(0, 0.5),
    beta_term[term_p] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls, beta_is) ~ normal(0, 0.5),
    c(sigma_yr, sigma_stk) ~ exponential(1),
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)

# check priors (look good, note following code doesn't work with dirichlet)
# set.seed(1999)
# prior <- extract.prior( m1 , n=1e3 )
# prior_p_dat <- data.frame(
#   fl_z = seq(-2, 2, length = 30),
#   lipid_z = 0,
#   day_z = 0,
#   term_p = 1,
#   yr = 1,
#   stk = 1
# )
# prior_p <- link(m1, post = prior, data = prior_p_dat)
# hist(prior_p)
# plot(prior_p[1, ] ~ prior_p_dat$fl_z, type = "line", ylim = c(0, 1))
# for (i in 2:50) {
#   lines(prior_p_dat$fl_z, prior_p[i, ])
# }


# as model 1 but assumes date influences size and lipid content, which also 
# covary with one another
m2 <- ulam(
  alist(
    # unobserved latent condition
    c(fl_z, lipid_z) ~ multi_normal(c(mu_fl, mu_lipid), Rho, Sigma),
    mu_fl <- alpha_fl + beta_df * day_z,
    mu_lipid <- alpha_lipid + beta_dl * day_z,
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_term[term_p] +
      beta_yr[yr] * sigma_yr +
      beta_stk[stk] * sigma_stk +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z
    ,
    c(alpha_fl, alpha_lipid) ~ normal(0, 0.5),
    c(beta_df, beta_dl) ~ normal(0.1, 0.5),
    Rho ~ lkj_corr(2),
    Sigma ~ exponential(1),
    surv_bar ~ normal(0, 1.25),
    beta_term[term_p] ~ normal(0, 0.5),
    beta_yr[yr] ~ normal(0, 0.5),
    beta_stk[stk] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls) ~ normal(0, 0.5),
    c(sigma_yr, sigma_stk) ~ exponential(1)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)
# modest covariance among fork length and lipid (~0.3); similar sigmas (0.9);
# moderate day of year effects (0.3); other values unchanged


# as model 2 but assume a) year- and stock-specific intercepts on size and 
# lipid and b) stock-specific effects on tagging date;
# assume year and stock effects in each submodel are drawn from multivariate
m3 <- ulam(
  alist(
    # stock-specific sampling dates
    day_z ~ normal(mu_day, sigma_day),
    mu_day <- beta_stk[stk, 4],
    # unobserved latent condition
    c(fl_z, lipid_z) ~ multi_normal(c(mu_fl, mu_lipid), Rho, Sigma),
    mu_fl <- alpha_fl + beta_df * day_z + beta_yr[yr, 1] + beta_stk[stk, 1],
    mu_lipid <- alpha_lipid + beta_dl * day_z + beta_yr[yr, 2] + 
      beta_stk[stk, 2],
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_term[term_p] +
      beta_yr[yr, 3] +
      beta_stk[stk, 3] +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z,
    # adaptive priors
    transpars> matrix[yr, 3]:beta_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[3, yr]:z_yr ~ normal(0, 0.5),
    transpars> matrix[stk, 4]:beta_stk <-
      compose_noncentered(sigma_stk , L_Rho_stk , z_stk),
    matrix[4, stk]:z_stk ~ normal(0, 0.5),
    # fixed priors
    c(alpha_fl, alpha_lipid) ~ normal(0, 0.5),
    c(beta_df, beta_dl) ~ normal(0.1, 0.5),
    Rho ~ lkj_corr(2),
    Sigma ~ exponential(1),
    sigma_day ~ exponential(1),
    
    cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[3]:sigma_yr ~ exponential(1),
    cholesky_factor_corr[4]:L_Rho_stk ~ lkj_corr_cholesky(2),
    vector[4]:sigma_stk ~ exponential(1),
    
    surv_bar ~ normal(0, 1.25),
    beta_term[term_p] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls) ~ normal(0, 0.5),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr),
    gq> matrix[4, 4]:Rho_stk <<- Chol_to_Corr(L_Rho_stk)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)
# covariance between year effects on survival and lipid/fl much weaker (~0.05)
# covariance between stock effects on survival and lipid/fl much weaker (~0.05)
# greater interstock and interannual variability in lipid content than fork
# length or survival


# as model 3 but includes injury effects
m4 <- ulam(
  alist(
    # stock-specific sampling dates
    day_z ~ normal(mu_day, sigma_day),
    mu_day <- beta_stk[stk, 4],
    # unobserved latent condition
    c(fl_z, lipid_z) ~ multi_normal(c(mu_fl, mu_lipid), Rho, Sigma),
    mu_fl <- alpha_fl + beta_df * day_z + beta_yr[yr, 1] + beta_stk[stk, 1],
    mu_lipid <- alpha_lipid + beta_dl * day_z + beta_yr[yr, 2] + 
      beta_stk[stk, 2],
    # survival
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_term[term_p] +
      beta_yr[yr, 3] +
      beta_stk[stk, 3] +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z +
      beta_is * sum(delta_inj[1:inj])
    ,
    # adaptive priors
    transpars> matrix[yr, 3]:beta_yr <-
      compose_noncentered(sigma_yr , L_Rho_yr , z_yr),
    matrix[3, yr]:z_yr ~ normal(0, 0.5),
    transpars> matrix[stk, 4]:beta_stk <-
      compose_noncentered(sigma_stk , L_Rho_stk , z_stk),
    matrix[4, stk]:z_stk ~ normal(0, 0.5),
    # fixed priors
    c(alpha_fl, alpha_lipid) ~ normal(0, 0.5),
    c(beta_df, beta_dl) ~ normal(0.1, 0.5),
    Rho ~ lkj_corr(2),
    Sigma ~ exponential(1),
    sigma_day ~ exponential(1),
    
    cholesky_factor_corr[3]:L_Rho_yr ~ lkj_corr_cholesky(2),
    vector[3]:sigma_yr ~ exponential(1),
    cholesky_factor_corr[4]:L_Rho_stk ~ lkj_corr_cholesky(2),
    vector[4]:sigma_stk ~ exponential(1),
    
    surv_bar ~ normal(0, 1.25),
    beta_term[term_p] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls, beta_is) ~ normal(0, 0.5),
    
    # constraints on ordinal effects of injury
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha),
    
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[3, 3]:Rho_yr <<- Chol_to_Corr(L_Rho_yr),
    gq> matrix[4, 4]:Rho_stk <<- Chol_to_Corr(L_Rho_stk)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)
# as in simpler model, injury effects are very modest (overall effect broadly
# spans zero and scales semi-linearly with injury)


## posterior predictions
# sampling day, including indirect effects on size/lipid
samp_date_seq <- seq(from = -2, to = 2, length.out = 30)

post <- extract.samples(m3)
pred_mu_fl <- purrr::map(
  samp_date_seq, ~ post$alpha_fl + post$beta_df * .x
)
pred_mu_lipid <- purrr::map(
  samp_date_seq, ~ post$alpha_lipid + post$beta_dl * .x
)

# generate posterior covariance matrix for each draw, combine with pred mu to 
# sample from mvrnorm, then iterate over exp var
sigma <- post$Sigma
rho <- post$Rho
S <- vector(mode = "list", length = nrow(sigma))
sim_surv_d <- sim_fl <- sim_lipid <- matrix(
  NA, 
  nrow = nrow(sigma), 
  ncol = length(samp_date_seq)
)
for (j in seq_along(samp_date_seq)) {
  for (i in 1:nrow(sigma)) {
    S <- diag(sigma[i,]) %*% rho[i,,] %*% diag(sigma[i,])
    mu <- MASS::mvrnorm(
      n = 1, mu=c(pred_mu_fl[[j]][i], pred_mu_lipid[[j]][i]),
      Sigma = S
    )
    sim_fl[i, j] <- mu[1]
    sim_lipid[i, j] <- mu[2]
  }
  sim_eta <- as.numeric(post$surv_bar + post$beta_ds * samp_date_seq[j] +
                          post$beta_fs * sim_fl[ , j] +
                          post$beta_ls * sim_lipid[ , j]
  )
  sim_surv_d[ , j] <- boot::inv.logit(sim_eta) # Probability of survival
  # sim_surv_d[ , j] <- rbinom(length(sim_eta), 1, boot::inv.logit(sim_eta))
}
pred_d_mu <- apply(sim_surv_d, 2, mean)
pred_d_pi <- apply(sim_surv_d, 2, PI)
plot( NULL , xlab="Date Scaled" , ylab="Proportion Terminal Det",
      ylim=c(0,1) , xaxt="n" , xlim=c(-2, 2) )
lines(seq(-2, 2, length = 30) , pred_d_mu )
shade( pred_d_pi , seq(-2, 2, length = 30))


# plot posterior preds
# post <- extract.samples(m1)
# link_foo <- function(pred_dat) {
#   logodds <- with(
#     post,
#     surv_bar + beta_ds * pred_dat$day_z + beta_fs * pred_dat$fl_z + beta_ls *
#       pred_dat$lipid_z 
#   )
#   return( inv_logit(logodds) )
# }


pred_l <- sapply(
  seq(-2, 2, length = 30), 
  function (x) {
    inv_logit(post$surv_bar + post$beta_ls * x)
  }
)
pred_l_mu <- apply(pred_l, 2, mean)
pred_l_pi <- apply(pred_l, 2, PI)
plot( NULL , xlab="Lipid Scaled" , ylab="Proportion Terminal Det",
      ylim=c(0,1) , xaxt="n" , xlim=c(-2, 2) )
lines(seq(-2, 2, length = 30) , pred_l_mu )
shade( pred_l_pi , seq(-2, 2, length = 30))


pred_f <- sapply(
  seq(-2, 2, length = 30), 
  function (x) {
    inv_logit(post$surv_bar + post$beta_fs * x)
  }
)
pred_f_mu <- apply(pred_f, 2, mean)
pred_f_pi <- apply(pred_f, 2, PI)
plot( NULL , xlab="FL Scaled" , ylab="Proportion Terminal Det",
      ylim=c(0,1) , xaxt="n" , xlim=c(-2, 2) )
lines(seq(-2, 2, length = 30) , pred_f_mu )
shade( pred_f_pi , seq(-2, 2, length = 30))


# pred_d <- sapply(
#   seq(-2, 2, length = 30), 
#   function (x) {
#     inv_logit(post$surv_bar + post$beta_ds * x)
#   }
# )
# pred_d_mu <- apply(pred_d, 2, mean)
# pred_d_pi <- apply(pred_d, 2, PI)
# plot( NULL , xlab="Yday Scaled" , ylab="Proportion Terminal Det",
#       ylim=c(0,1) , xaxt="n" , xlim=c(-2, 2) )
# lines(seq(-2, 2, length = 30) , pred_d_mu )
# shade( pred_d_pi , seq(-2, 2, length = 30))



## FREQUENTIST MODELS ----------------------------------------------------------

dum_fit <- glmmTMB(
  lipid ~ day_z + (1 | stock_group) + (1 | year), 
  data = det_dat
)


### EVALUATE HOW TO DEAL WITH INJURIES -----------------------------------------

# injury vs redeployed tag effects

det_dat_no_redep <- det_dat %>% 
  filter(redeploy == "no")
det_dat_no_3 <- det_dat %>% 
  filter(!injury == "3")

# assume year and stock group are random intercepts
# compare three models for terminal detections: injury (0-3) & redeploy tags
dat_list <- list(
  surv = det_dat$term_det,
  inj = det_dat$injury,
  redeploy = ifelse(det_dat$redeploy == "no", 1, 2),
  yr = as.integer(as.factor(det_dat$year)),
  stock = as.integer(as.factor(det_dat$stock_group)),
  alpha = rep(2, length(unique(det_dat$injury)) - 1)
)
dat_list2 <- list(
  surv = det_dat_no_redep$term_det,
  inj = det_dat_no_redep$injury,
  redeploy = ifelse(det_dat_no_redep$redeploy == "no", 1, 2),
  yr = as.integer(as.factor(det_dat_no_redep$year)),
  stock = as.integer(as.factor(det_dat_no_redep$stock_group)),
  alpha = rep(2, length(unique(det_dat_no_redep$injury)) - 1)
)
dat_list3 <- list(
  surv = det_dat_no_3$term_det,
  inj = det_dat_no_3$injury,
  redeploy = ifelse(det_dat_no_3$redeploy == "no", 1, 2),
  yr = as.integer(as.factor(det_dat_no_3$year)),
  stock = as.integer(as.factor(det_dat_no_3$stock_group)),
  alpha = rep(2, length(unique(det_dat_no_3$injury)) - 1)
)


# null model
m0 <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_yr[yr]*sigma_yr +
      gamma_stock[stock] * sigma_stock,
    # priors
    gamma ~ normal(0, 2.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list, chains=4 , cores = 4,
  control = list(adapt_delta = 0.95)
)


m1 <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_yr[yr] * sigma_yr +
      gamma_stock[stock] * sigma_stock +
      gamma_inj * sum(delta_inj[1:inj]),
    # priors
    gamma ~ normal(0, 2.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    gamma_inj ~ normal(0, 0.5),
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list, chains=4 , cores = 4,#log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)
m1b <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_yr[yr] * sigma_yr +
      gamma_stock[stock] * sigma_stock +
      gamma_inj * sum(delta_inj[1:inj]),
    # priors
    gamma ~ normal(0, 2.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    gamma_inj ~ normal(0, 0.5),
    vector[4]: delta_inj <<- append_row(0, delta),
    simplex[3]: delta ~ dirichlet(alpha),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list2, chains=4 , cores = 4,#log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)
m1c <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_yr[yr] * sigma_yr +
      gamma_stock[stock] * sigma_stock +
      gamma_inj * sum(delta_inj[1:inj]),
    # priors
    gamma ~ normal(0, 2.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    gamma_inj ~ normal(0, 0.5),
    vector[3]: delta_inj <<- append_row(0, delta),
    simplex[2]: delta ~ dirichlet(alpha),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list3, chains=4 , cores = 4,#log_lik=TRUE,
  control = list(adapt_delta = 0.95)
)

post <- extract.samples(m1b)
p_link_inj <- function(injury) {
  if (injury == 0) {
    logodds <- with(post, gamma)
  } else (
    logodds <- with(
      post, gamma + (gamma_inj * rowSums(delta[, 1:injury, drop = FALSE])) 
    )
  )
  return( inv_logit(logodds) )
}
p_raw <- sapply(0:3 , function(i) p_link_inj( i ) )
inj_dat <- purrr::map(
  seq(0, 3, by = 1), function (i) {
    data.frame(
      injury = as.factor(i),
      est = p_raw[ , i + 1])
  }
) %>% 
  bind_rows()
ggplot(inj_dat, aes(x = injury, y = est)) +
  geom_boxplot()


# redeploy model
m2 <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_redeploy[redeploy] +
      gamma_yr[yr]*sigma_yr +
      gamma_stock[stock] * sigma_stock,
    # priors
    gamma ~ normal(0, 2.5),
    gamma_redeploy[redeploy] ~ dnorm(0, 0.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)
m2b <- ulam(
  alist(
    surv ~ dbinom( 1 , p ) ,
    logit(p) <- gamma + 
      gamma_redeploy[redeploy] +
      gamma_yr[yr]*sigma_yr +
      gamma_stock[stock] * sigma_stock,
    # priors
    gamma ~ normal(0, 2.5),
    gamma_redeploy[redeploy] ~ dnorm(0, 0.5),
    gamma_yr[yr] ~ dnorm(0, 0.5),
    gamma_stock[stock] ~ dnorm(0, 0.5),
    c(sigma_yr, sigma_stock) ~ exponential(1)
  ),
  data=dat_list3, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95)
)

post2 <- extract.samples(m2b)
p_link_redeploy <- function(redeploy) {
  logodds <- with(
      post2, gamma + gamma_redeploy[ , redeploy] 
    )
  return( inv_logit(logodds) )
}
p_raw <- sapply(1:2, function(i) p_link_redeploy( i ) )
redep_dat <- purrr::map(
  c(1, 2), function (i) {
    data.frame(
      redeploy = as.factor(i),
      est = p_raw[ , i])
  }
) %>% 
  bind_rows()
ggplot(redep_dat, aes(x = redeploy, y = est)) +
  geom_boxplot()

# greatest stabilizing effect from removing redeployed tags
