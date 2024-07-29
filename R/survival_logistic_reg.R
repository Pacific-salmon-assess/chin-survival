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

det_tbl <- readRDS(here::here("data", "det_history_tbl.RDS"))

chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS"))

# pull agg names for fraser
agg_names <- det_tbl$bio_dat[[2]] %>% 
  select(vemco_code, agg) 
  
  
det_dat1 <- det_tbl %>% 
  dplyr::select(stock_group, agg_det) %>% 
  unnest(cols = c(agg_det)) %>%
  left_join(., agg_names, by = "vemco_code") %>% 
  mutate(
    term_det = ifelse(final_det + river_det > 0, 1, 0),
    stock_group = ifelse(stock_group == "Fraser", agg, stock_group) %>% 
      as.factor()
  ) %>% 
  left_join(., 
            chin %>% 
              select(vemco_code = acoustic_year, year, acoustic_type, 
                     fl, lipid, year_day, hook_loc, fin_dam, injury, comment),
            by = "vemco_code") %>%
  mutate(
    redeploy = ifelse(acoustic_type %in% c("V13P", "V13"), "no", "yes")
  ) %>% 
  arrange(
    year, stock_group
  )


# small number of tags (<2% missing lipid data; impute)
bio_dat <- det_dat1 %>% 
  select(vemco_code, fl, lipid, stock_group, year) %>%
  distinct() 
interp_lipid <- bio_dat %>%
  select(-vemco_code) %>% 
  VIM::kNN(., k = 5) %>% 
  select(-ends_with("imp")) 
det_dat1$lipid <- interp_lipid$lipid


# remove redeployed tags and scale continuous covariates
det_dat <- det_dat1 %>% 
  filter(
    !redeploy == "yes"
  ) %>% 
  mutate(
    lipid_z = scale(lipid) %>% as.numeric(),
    fl_z = scale(fl) %>% as.numeric(),
    day_z = scale(year_day) %>% as.numeric()
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


# FIT LOGISTIC REG MODELS ------------------------------------------------------


dat_list <- list(
  surv = as.integer(det_dat$term_det),
  fl_z = det_dat$fl_z,
  lipid_z = det_dat$lipid_z,
  day_z = det_dat$day_z,
  # inj = det_dat$injury,
  yr = as.integer(as.factor(det_dat$year)),
  stk = as.integer(as.factor(det_dat$stock_group))#,
  # alpha = rep(2, length(unique(det_dat$injury)) - 1)
)


m1 <- ulam(
  alist(
    surv ~ dbinom(1 , p) ,
    logit(p) <- surv_bar + 
      beta_yr[yr] * sigma_yr +
      beta_stk[stk] * sigma_stk +
      beta_ds * day_z +
      beta_fs * fl_z +
      beta_ls * lipid_z
      ,
    surv_bar ~ normal(0, 1.25),
    beta_yr[yr] ~ normal(0, 0.5),
    beta_stk[stk] ~ normal(0, 0.5),
    c(beta_ds, beta_fs, beta_ls) ~ normal(0, 0.5),
    c(sigma_yr, sigma_stk) ~ exponential(1)
  ),
  data=dat_list, chains = 4 , cores = 4, iter = 2000,
  control = list(adapt_delta = 0.95)
)


# check priors
# set.seed(1999)
# prior <- extract.prior( m1 , n=1e4 )
# prior_p_dat <- data.frame(
#   fl_z = seq(-2, 2, length = 30),
#   lipid_z = 0,
#   day_z = 0,
#   yr = 1,
#   stk = 1
# )
# prior_p <- link(m1, post = prior, data = prior_p_dat)
# hist(prior_p)
# plot(prior_p[1, ] ~ prior_p_dat$fl_z, type = "line", ylim = c(0, 1))
# for (i in 2:50) {
#   lines(prior_p_dat$fl_z, prior_p[i, ])
# }


# plot posterior preds
post <- extract.samples(m1)
link_foo <- function(pred_dat) {
  logodds <- with(
    post,
    surv_bar + beta_ds * pred_dat$day_z + beta_fs * pred_dat$fl_z + beta_ls *
      pred_dat$lipid_z 
  )
  return( inv_logit(logodds) )
}


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


pred_d <- sapply(
  seq(-2, 2, length = 30), 
  function (x) {
    inv_logit(post$surv_bar + post$beta_ds * x)
  }
)
pred_d_mu <- apply(pred_d, 2, mean)
pred_d_pi <- apply(pred_d, 2, PI)
plot( NULL , xlab="Yday Scaled" , ylab="Proportion Terminal Det",
      ylim=c(0,1) , xaxt="n" , xlim=c(-2, 2) )
lines(seq(-2, 2, length = 30) , pred_d_mu )
shade( pred_d_pi , seq(-2, 2, length = 30))



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
