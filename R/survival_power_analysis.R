## Power analyses for stock-specific survival for 2023 AUP and stage-specific
# survival for sockeye LOI

library(tidyverse)

## Simulate 5 years of data, 30 individuals per year with survival varying 
# randomly among years

fr_mean_surv <- 0.7
col_mean_surv <- 0.8
yr_sig <- 0.15
# rho <- 0.8
# cov <- 0.2 * 0.2 * rho
# cov_mat <- matrix(c(yr_sig, cov, cov, yr_sig), nrow = 2)
n_years <- 5
n_iters <- 500

future::plan(future::multisession, workers = 3)

p_list <- furrr::future_map(col_mean_surv, function (x) {
  
  dum_list <- vector(n_years, mode = "list")
  p_val <- rep(NA, n_iters)
  
  for (j in seq_len(n_iters)) {
    for (i in seq_len(n_years)) {
      # random draw representing interannual variability
      surv_anom <- rnorm(1, 0, yr_sig)
      # calculate annual survival in logit space then back transform
      fraser_surv <- plogis(qlogis(fr_mean_surv) + surv_anom)
      col_surv <- plogis(qlogis(x) + surv_anom)
      dum_list[[i]] <- data.frame(
        year = i,
        surv_fraser = rbinom(30, 1, fraser_surv),
        surv_col = rbinom(30, 1, col_surv)
      )
    }
    dum_dat <- dum_list %>% 
      bind_rows() %>% 
      pivot_longer(cols = starts_with("surv"),
                   names_prefix = "surv_",
                   names_to = "stock",
                   values_to = "surv") %>% 
      mutate(year = as.factor(year))
    fit <- lme4::glmer(
      surv ~ stock + (1 | year),
      data = dum_dat,
      family = binomial
    )  
    p_val[j] <- coef(summary(fit))[2, 4]  
  }
  
  sum(p_val <= 0.05) / n_iters
},
.options = furrr::furrr_options(seed = 123)
)

# note singularity warnings are a function of the relatively small number of 
# random intercept levels (n = 5); using a Bayesian model with weakly informative
# priors can address these issues but was impractical for this application due 
# to longer run times


## Sockeye 

qcs_mean_surv <- 0.7
sog_mean_surv <- c(0.8, 0.85, 0.95)
qcs_n <- 150
yr_sig <- 0.15
n_years <- 2
n_iters <- 500

future::plan(future::multisession, workers = 3)

p_list <- furrr::future_map(sog_mean_surv, function (x) {
  
  dum_list <- vector(n_years, mode = "list")
  p_val <- rep(NA, n_iters)
  
  for (j in seq_len(n_iters)) {
    for (i in seq_len(n_years)) {
      # random draw representing interannual variability
      surv_anom <- rnorm(1, 0, yr_sig)
      # calculate annual survival in logit space then back transform
      qcs_surv <- plogis(qlogis(qcs_mean_surv) + surv_anom)
      sog_surv <- plogis(qlogis(x) + surv_anom)
      surv_qcs1 = rbinom(qcs_n, 1, qcs_surv)
      # since sog occurs after qcs, sample size will be the number of individuals
      # that survive stage 1
      surv_sog1 = rbinom(sum(surv_qcs1), 1, sog_surv)
      dum_list[[i]] <- data.frame(
        year = i,
        surv_qcs = sum(surv_qcs1),
        surv_sog = sum(surv_sog1),
        n = qcs_n
      )
    }
    dum_dat <- dum_list %>% 
      bind_rows() %>% 
      mutate(qcs_n = surv_qcs) %>% 
      pivot_longer(cols = starts_with("surv"),
                   names_prefix = "surv_",
                   names_to = "stage",
                   values_to = "surv") %>% 
      mutate(
        year = as.factor(year),
        # account for sog dead being a function of n - qcs_dead
        dead = ifelse(
          stage == "qcs",
          n - surv,
          qcs_n - surv
        )
      )
    fit <- glm(
      cbind(surv, dead) ~ stage + year,
      data = dum_dat,
      family = binomial
    )  
    p_val[j] <- coef(summary(fit))[2, 4]  
  }
  
  sum(p_val <= 0.05) / n_iters
},
.options = furrr::furrr_options(seed = 123)
)


# Chinook tagging two; goal is simply to estimate survival rates
fr_mean_surv <- 0.65
ps_mean_surv <- 0.85
yr_sig <- 0.0
n_years <- 2
n_iters <- 500


dum_list <- vector(n_years, mode = "list")
ps_surv <- fr_surv <- p_val <- rep(NA, n_iters)

for (j in seq_len(n_iters)) {
  for (i in seq_len(n_years)) {
    # random draw representing interannual variability
    surv_anom <- rnorm(1, 0, yr_sig)
    # calculate annual survival in logit space then back transform
    fraser_surv <- plogis(qlogis(fr_mean_surv) + surv_anom)
    col_surv <- plogis(qlogis(ps_mean_surv) + surv_anom)
    dum_list[[i]] <- data.frame(
      year = i,
      surv_fraser = rbinom(16, 1, fraser_surv),
      surv_col = rbinom(32, 1, col_surv)
    )
  }
  dum_dat <- dum_list %>% 
    bind_rows() %>% 
    pivot_longer(cols = starts_with("surv"),
                 names_prefix = "surv_",
                 names_to = "stock",
                 values_to = "surv") %>% 
    mutate(year = as.factor(year))
  fit <- glm(
    surv ~ stock,
    data = dum_dat,
    family = binomial
  )  
  p_val[j] <- coef(summary(fit))[2, 4]
  ps_surv[j] <- boot::inv.logit(coef(summary(fit))[1, 1])
  fr_surv[j] <- boot::inv.logit(coef(summary(fit))[1, 1] + coef(summary(fit))[2, 1])
}

sum(p_val <= 0.05) / n_iters
sum(ps_surv < ps_mean_surv - 0.1 |
      ps_surv > ps_mean_surv + 0.1) / n_iters
sum(fr_surv < fr_mean_surv - 0.1 |
      fr_surv > fr_mean_surv + 0.1) / n_iters
