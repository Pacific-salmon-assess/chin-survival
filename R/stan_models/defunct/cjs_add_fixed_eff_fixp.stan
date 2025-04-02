// This models is derived from section 12.3 of "Stan Modeling Language
// User's Guide and Reference Manual"
// Based on cjs_add.stan from example-models repo, but includes time-varying p

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
                        + phi[i, t_curr] * (1 - p[i, t_next - 1]) * chi[i, t_next];
      }
    }
    return chi;
  }
}

data {
  int<lower=0> nind;            // Number of individuals
  int<lower=2> n_occasions;     // Number of capture occasions
  int<lower=0,upper=1> y[nind, n_occasions];    // Capture-history
  int<lower=1> g;               // Number of groups
  int<lower=1,upper=g> group[nind];     // Groups
  real<lower=0> final_fix_p;     // fixed values for final detection prob
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
  real mu_phi;                            // Mean phi
  vector[n_occ_minus_1] gamma_phi;        // Time effect phi
  vector[g] beta_phi;                     // Group effect phi
  matrix[g, n_occ_minus_1 - 1] beta_p;    // Group/time interacting effect p (truncated relative to phi because last observation fixed in transformed parameters block)
}

transformed parameters {
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] phi;
  matrix<lower=0,upper=1>[nind, n_occ_minus_1] p;
  matrix<lower=0,upper=1>[nind, n_occasions] chi;


  // Constraints
  for (i in 1:nind) {
    for (t in 1:(first[i] - 1)) {
      phi[i, t] = 0;
      p[i, t] = 0;
    }
    for (t in first[i]:n_occ_minus_1) {
      phi[i, t] = inv_logit(mu_phi + beta_phi[group[i]] + gamma_phi[t]);
    //  p[i, t] = inv_logit(beta_p[group[i], t]);
    }
    // fix p in last time step
    for (t in first[i]:(n_occ_minus_1 - 1)) {
      p[i, t] = inv_logit(beta_p[group[i], t]);
    }
    p[i, n_occ_minus_1] = final_fix_p;
  }

  chi = prob_uncaptured(nind, n_occasions, p, phi);
}

model {
  // Priors
  // Choose slightly diffuse priors
  for(t in 1:(n_occ_minus_1 - 1)) {
    for(gg in 1:g) {
      beta_p[gg, t] ~ normal(0, 1.5);
    }
  }
  mu_phi ~ normal(0.3, 1.5);
  beta_phi ~ normal(0, 0.5);
  gamma_phi ~ normal(0, 0.5);
     

  // Likelihood
  for (i in 1:nind) {
    if (first[i] > 0) {
      for (t in (first[i] + 1):last[i]) {
        1 ~ bernoulli(phi[i, t - 1]);
        y[i, t] ~ bernoulli(p[i, t - 1]);
      }
      1 ~ bernoulli(chi[i, last[i]]);
    }
  }
}

generated quantities {
  matrix<lower=0,upper=1>[g, n_occ_minus_1] phi_g;
  matrix<lower=0,upper=1>[g, n_occ_minus_1] p_g;

  // matrices to store obs
  int<lower=0,upper=1> y_hat[nind, n_occasions];
  real<lower=0,upper=1> mu_state[nind, n_occasions];
  real<lower=0,upper=1> mu_obs[nind, n_occasions];
  int<lower=0,upper=1> z[nind, n_occasions];

  // posterior estimates of group_specific phi
  for (gg in 1:g) {
    for (t in 1:n_occ_minus_1) {
      phi_g[gg, t] = inv_logit(mu_phi + beta_phi[gg] + gamma_phi[t]);
      //p_g[gg, t] = inv_logit(beta_p[gg, t]);
    }
    for (tt in 1:n_occ_minus_1 - 1) {
      p_g[gg, tt] = inv_logit(beta_p[gg, tt]);
    }
    p_g[gg, n_occ_minus_1] = final_fix_p;
  }

  // posterior predictions of observations
  for (i in 1:nind) {
    y_hat[i, 1] = 1;
    mu_obs[i, 1] = 1;
    mu_state[i, 1] = 1;
    z[i, 1] = 1;

    for (t in 2:n_occasions) {
      // state process
      mu_state[i, t] = phi_g[group[i], t - 1] * z[i, t - 1];
      z[i, t] = bernoulli_rng(mu_state[i, t]);
      
      // obs process
      mu_obs[i, t] = p_g[group[i], t - 1] * z[i, t];
      y_hat[i, t] = bernoulli_rng(mu_obs[i, t]);
    }
  }
}
