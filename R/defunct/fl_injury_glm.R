## Body size vs. injury effects
# Estimate effect of body size on injury category for FULL dataset (i.e., not
# just tagged individuals), to determine whether relationship should be included
# in final DAG
# March 18, 2025

library(tidyverse)
library(rethinking)


chin_dat <- read.csv(
  here::here("data", "tagChinData.csv")
) %>% 
  janitor::clean_names() %>% 
  mutate(
    # define injury scores as per SJ's analyses
    inj = case_when(
      hook_loc == "eye" ~ 4,
      injury == "3" | fin_dam == "3" | scale_loss == "3" ~ 3,
      injury == "2" | fin_dam == "2" | scale_loss == "2" ~ 2,
      TRUE ~ 1
    ),
    fl_z = scale(fl) %>% as.numeric()
  ) %>% 
  filter(
    # remove fish smaller than tagging size
    !fl < 60
  )


chin_dat %>%
  mutate(
    size_bin = case_when(
      fl <= 65 ~ "xsmall",
      fl <= 70 ~ "small",
      fl <= 80 ~ "med",
      fl > 80 ~ "big"
    )
  ) %>% 
  ggplot(., aes(x = size_bin, fill = as.factor(inj))) +
  geom_bar(position = "stack") +
  ggsidekick::theme_sleek()


dat_list <- list(
  fl_z = chin_dat$fl_z,
  inj = chin_dat$inj
)

m_ord <- ulam(
  alist(
    inj ~ dordlogit(phi, cutpoints),  # Ordered logistic likelihood
    phi <- a + b * fl_z,  # Linear predictor
    a ~ normal(0, 1),  # Intercept prior
    b ~ normal(0, 1),  # Slope prior
    cutpoints ~ normal(0, 1)  # Priors for category thresholds
  ), 
  data = dat_list, chains = 4, cores = 4
)

precis(m_ord, depth = 2)


# Define new X values at which we want predictions
new_fl <- c(-1, 0, 1)

# Extract posterior samples of model parameters
post <- extract.samples(m_ord)

# Compute the latent variable phi for new values of X
phi_samples <- sapply(new_fl, function(x) {
  post$a + post$b * x
})

# Convert phi into category probabilities using the ordinal logistic function
pred_probs <- lapply(1:ncol(phi_samples), function(i) {  # Loop over X values
  t(sapply(1:nrow(phi_samples), function(j) {  # Loop over posterior samples
    pordlogit(1:4, phi_samples[j, i], post$cutpoints[j, ])  # Use correct cutpoints
  }))
})

cat_probs <- lapply(pred_probs, function(p) {
  t(sapply(1:nrow(p), function(j) {  # Loop over posterior samples
    c(p[j, 1], diff(p[j, ]))
  }))
})


# Compute median probabilities across posterior samples
median_probs <- sapply(cat_probs, function(p) apply(p, 2, median))

# Convert to long format for ggplot
df_plot <- as.data.frame(median_probs) %>%
  mutate(Category = factor(1:4)) %>%
  pivot_longer(starts_with("V"), names_to = "X", values_to = "Probability") %>% 
  mutate(
    X = factor(X, levels = c("V1", "V2", "V3"), labels = c("-1", "0", "1"))
  )
ggplot(df_plot, aes(x = X, y = Probability, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Scaled Fork Length", y = "Median Predicted Proportion", 
       fill = "Injury Category") +
  theme_minimal()
