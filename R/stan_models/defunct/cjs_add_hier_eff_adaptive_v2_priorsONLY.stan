// This models is derived from section 12.3 of "Stan Modeling Language
// User's Guide and Reference Manual"
// Based on cjs_add.stan from example-models repo, but includes time-varying p and further adjusted to include t steps for p (as in Stan manual), not t-1 (BPA book)

functions {
  int first_capture(int[] y_i) {
    for (k in 1:size(y_i))
      if (y_i[k])
        return k;
    return 0;
  }

  int last_capture(int[] y_i) {
    for (k_rev in 0:(size(y_i) - 1)) {
      // Compoud declaration was enabled in Stan 2.13
      int k = size(y_i) - k_rev;
      if (y_i[k])
        return k;
    }
    return 0;
  }

  matrix prob_uncaptured(int nind, int n_occasions,
                         matrix p, matrix phi) {
    matrix[nind, n_occasions] chi;

    for (i in 1:nind) {
      chi[i, n_occasions] = 1.0;
      for (t in 1:(n_occasions - 1)) {
        // Compoud declaration was enabled in Stan 2.13
        int t_curr = n_occasions - t;
        int t_next = t_curr + 1;
        
        chi[i, t_curr] = (1 - phi[i, t_curr])
                        + phi[i, t_curr] * (1 - p[i, t_next]) * chi[i, t_next];
      }
    }
    return chi;
  }
}

data {
  int<lower=0> nind;            // Number of individuals
  int<lower=2> n_occasions;     // Number of capture occasions
  int<lower=0,upper=1> y[nind, n_occasions];    // Capture-history
  int<lower=1> nyear;               // Number of years
  int<lower=1,upper=nyear> year[nind];     // Year
}

transformed data {
  int n_occ_minus_1 = n_occasions - 1;
  int<lower=0,upper=n_occasions> first[nind];
  int<lower=0,upper=n_occasions> last[nind];
  
  for (i in 1:nind)
    first[i] = first_capture(y[i]);
  for (i in 1:nind)
    last[i] = last_capture(y[i]);
}

parameters {
  real alpha_phi;                            // Mean phi
  vector[n_occ_minus_1] alpha_t_phi;        // Time effect phi
  matrix[n_occ_minus_1, nyear] alpha_yr_phi_z;    // Year/time random int for phi; note z so inverted before transformation
  vector<lower=0>[n_occ_minus_1] sigma_alpha_yr_phi;   // SD among years
  cholesky_factor_corr[n_occ_minus_1] L_Rho_yr;    // for covariance among year-stage intercepts
  real alpha_p;                              // Mean det prob
  matrix[nyear, n_occasions] alpha_yr_p;        // Year/time intercepts for p
}

transformed parameters {
  // Construct non-centered adaptive priors for hierarchical effects
  matrix[nyear, n_occ_minus_1] alpha_yr_phi;
  alpha_yr_phi = (diag_pre_multiply(sigma_alpha_yr_phi, L_Rho_yr) * alpha_yr_phi_z)';


  matrix<lower=0,upper=1>[nind, n_occ_minus_1] phi;
  matrix<lower=0,upper=1>[nind, n_occasions] p;
  matrix<lower=0,upper=1>[nind, n_occasions] chi;

  // Constraints
  for (i in 1:nind) {
    for (t in 1:(first[i] - 1)) {
      phi[i, t] = 0;
      p[i, t] = 0;
    }
    for (t in first[i]:n_occ_minus_1) {
      phi[i, t] = inv_logit(alpha_phi + alpha_yr_phi[year[i], t] + alpha_t_phi[t]);
    }
    for (tt in first[i]:n_occasions) {
      p[i, tt] = inv_logit(alpha_p + alpha_yr_p[year[i], tt]);
    }
  }

  chi = prob_uncaptured(nind, n_occasions, p, phi);
}

model {
  // Priors
  to_vector(alpha_yr_p) ~ normal(0, 1.5);
  to_vector(alpha_yr_phi_z) ~ normal(0, 0.5);
  alpha_phi ~ normal(0.3, 1);
  alpha_p ~ normal(0, 0.5); 
  alpha_t_phi ~ normal(0, 0.5);
  sigma_alpha_yr_phi ~ exponential(1);
  L_Rho_yr ~ lkj_corr_cholesky(2);

  // Likelihood
  // for (i in 1:nind) {
  //   if (first[i] > 0) {
  //     for (t in (first[i] + 1):last[i]) {
  //       1 ~ bernoulli(phi[i, t - 1]);
  //       y[i, t] ~ bernoulli(p[i, t]);
  //     }
  //     1 ~ bernoulli(chi[i, last[i]]);
  //   }
  // }
}

generated quantities {
  matrix<lower=0,upper=1>[nyear, n_occ_minus_1] phi_yr;
  matrix<lower=0,upper=1>[nyear, n_occasions] p_yr;
  matrix[n_occ_minus_1, n_occ_minus_1] Rho_yr;
  // phi[n_occ_minus_1] and p[n_occasions] can't be identified simultaneously; define their product instead based on year-specific estimates
  vector[nyear] beta_yr; 
  
  // posterior estimates of year_specific phi
  for (yy in 1:nyear) {
    for (t in 1:n_occ_minus_1) {
      phi_yr[yy, t] = inv_logit(alpha_phi + alpha_yr_phi[yy, t] + alpha_t_phi[t]);
    }
    for (tt in 1:n_occasions) {
      p_yr[yy, tt] = inv_logit(alpha_p + alpha_yr_p[yy, tt]);
    }
    beta_yr[yy] = phi_yr[yy, n_occ_minus_1] * p_yr[yy, n_occasions];
  }
  Rho_yr = multiply_lower_tri_self_transpose(L_Rho_yr);
}
