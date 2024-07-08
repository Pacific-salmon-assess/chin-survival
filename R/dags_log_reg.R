## Logistic Regression DAGs
# Jun 16, 2024
# Explore conditional dependence of variables influencing Chinook salmon 
# survival

library(dagitty)

dag_simple <-dagitty( 
  "dag {
  L <- D <- A -> S
  A -> F -> S
  A -> L -> S
  D -> F
  D -> S
  S <- Y -> D -> F
  L <- Y -> F
  F -> L 
  }")
adjustmentSets(dag_simple, exposure = "L", outcome = "S")
adjustmentSets(dag_simple, exposure = "F", outcome = "S")
adjustmentSets(dag_simple, exposure = "D", outcome = "S")
adjustmentSets(dag_simple, exposure = "A", outcome = "S")
impliedConditionalIndependencies(dag_simple)

d_model <- bf(D ~ 1 + A)
f_model <- bf(FF ~ 1 + A + Y + D)
l_model <- bf(L ~ 1 + A + Y + D + FF)
s_model <- bf(S ~ 1 + A + Y + D + FF + L)


dag_full <-dagitty( 
  "dag {
  L <- D <- A -> S
  A -> F -> S
  A -> L -> S
  D -> F
  D -> S
  S <- Y -> D -> F
  L <- Y -> F
  F -> L 
  F -> I -> S
  F -> H -> S
  H -> I
  }")
adjustmentSets(dag_full, exposure = "L", outcome = "S")
adjustmentSets(dag_full, exposure = "I", outcome = "S")


# Following 14.4, which parameters should be correlated?
# stock effects on survival, size, and lipid
# year effects on survival, size, and lipid
# size and lipid effects on survival 


################################################################################

## experiment with multimodel inference 

library(rethinking)
library(brms)
library(ggdist)


## Example 1: Chapter 5 H2 homework 
# (https://sr2-solutions.wjakethompson.com/causes-confounds-colliders)

data("WaffleDivorce")

dat <- WaffleDivorce %>%
  select(A = MedianAgeMarriage,
         D = Divorce,
         M = Marriage) %>%
  mutate(across(everything(), standardize))

d_model <- bf(D ~ 1 + A)
a_model <- bf(A ~ 1 + M)

b5h2 <- brm(d_model + a_model + set_rescor(FALSE),
            data = dat, family = gaussian,
            prior = c(prior(normal(0, 0.2), class = Intercept, resp = D),
                      prior(normal(0, 0.5), class = b, resp = D),
                      prior(exponential(1), class = sigma, resp = D),
                      
                      prior(normal(0, 0.2), class = Intercept, resp = A),
                      prior(normal(0, 0.5), class = b, resp = A),
                      prior(exponential(1), class = sigma, resp = A)),
            iter = 4000, warmup = 2000, chains = 4, cores = 4, seed = 1234)

summary(b5h2)


as_draws_df(b5h2) %>%
  as_tibble() %>%
  select(.draw, b_D_Intercept:sigma_A) %>% 
  expand(nesting(.draw, b_D_Intercept, b_A_Intercept, b_D_A, b_A_M,
                 sigma_D, sigma_A),
         m = seq(from = -2, to = 2, length.out = 30)) %>%
  mutate(a_sim = rnorm(n(), mean = b_A_Intercept + b_A_M * m, sd = sigma_A),
         d_sim = rnorm(n(), mean = b_D_Intercept + b_D_A * a_sim, sd = sigma_D)) %>%
  pivot_longer(ends_with("_sim"), names_to = "name", values_to = "value") %>%
  group_by(m, name) %>%
  mean_qi(value, .width = c(0.89)) %>%
  ungroup() %>%
  mutate(name = case_when(name == "a_sim" ~ "Counterfactual M &rarr; A",
                          TRUE ~ "Counterfactual M &rarr; A &rarr; D")) %>%
  ggplot(aes(x = m, y = value, ymin = .lower, ymax = .upper)) +
  facet_wrap(~name, nrow = 1) +
  geom_smooth(stat = "identity") +
  labs(x = "Manipulated M", y = "Counterfactual")


## Chimpanzees adaptive priors example

data(chimpanzees, package = "rethinking")
d <- chimpanzees 
d2 <- d %>% 
  mutate(actor     = factor(actor),
         block     = factor(block),
         treatment = factor(1 + prosoc_left + 2 * condition),
         # this will come in handy, later
         labels    = factor(treatment,
                            levels = 1:4,
                            labels = c("r/n", "l/n", "r/p", "l/p")))

b14.3 <- brm(data = d2, 
      family = binomial,
      pulled_left | trials(1) ~ 0 + treatment + (0 + treatment | actor) #+
        # (0 + treatment | block)
      ,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(1), class = sd, group = actor),
                # prior(exponential(1), class = sd, group = block),
                prior(lkj(2), class = cor, group = actor)#,
                # prior(lkj(2), class = cor, group = block)
                ),
      iter = 2000, warmup = 1000, chains = 4, cores = 4,  
      seed = 4387510)
stancode(b14.3)

d$block_id <- d$block
d$treatment <- 1L + d$prosoc_left + 2L*d$condition
dat <- list(
  L = d$pulled_left,
  tid = d$treatment,
  actor = d$actor,
  block_id = as.integer(d$block_id)
)

m14.3 <- ulam(
  alist(
    L ~ binomial(1,p),
    logit(p) <- g[tid] + alpha[actor,tid] + beta[block_id,tid],
    # adaptive priors - non-centered
    transpars> matrix[actor,4]:alpha <-
      compose_noncentered( sigma_actor , L_Rho_actor , z_actor ),
    transpars> matrix[block_id,4]:beta <-
      compose_noncentered( sigma_block , L_Rho_block , z_block ),
    matrix[4,actor]:z_actor ~ normal( 0 , 1 ),
    matrix[4,block_id]:z_block ~ normal( 0 , 1 ),
    # fixed priors
    g[tid] ~ normal(0,1),
    vector[4]:sigma_actor ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_actor ~ lkj_corr_cholesky( 2 ),
    vector[4]:sigma_block ~ dexp(1),
    cholesky_factor_corr[4]:L_Rho_block ~ lkj_corr_cholesky( 2 ),
    # compute ordinary correlation matrixes from Cholesky factors
    gq> matrix[4,4]:Rho_actor <<- Chol_to_Corr(L_Rho_actor),
    gq> matrix[4,4]:Rho_block <<- Chol_to_Corr(L_Rho_block)
  ) , data=dat , chains=4 , cores=4 , log_lik=TRUE )
rethinking::stancode(m14.3)

precis(m14.3,3,pars=c("alpha"))
precis(m14.3,3,pars=c("Rho_actor"))

