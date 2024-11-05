## Simulate covariance among stage- and year-intercepts

# each sigma represents the among year SD for a given stage
sigma_s <- c(0.2, 0.3, 0.7)

# assume five years per sigma
set.seed(123)
nyr <- 5
alpha_yr_s <- matrix(NA, nrow = nyr, ncol = length(sigma_s))
for(i in seq_along(sigma_s)) {
  alpha_yr_s[ , i] <- rnorm(nyr, 0, sigma_s[i])
}

# as above but add covariance
Rho <- matrix(
  c(
    1, 0.5, -0.9,
    0.5, 1, -0.2,
    -0.9, -0.2, 1
  ),
  nrow = 3
)
eigen_values <- eigen(Rho)$values
print(eigen_values)

sigma_s2 <- diag(sigma_s) %*% Rho %*% diag(sigma_s)

set.seed(123)
alpha_yr_s2 <- MASS::mvrnorm(5, c(0,0,0), sigma_s2)

cor(alpha_yr_s2[, 1], alpha_yr_s2[,3])
