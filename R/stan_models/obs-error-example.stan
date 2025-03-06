// Model based on rethinking ulam function

data{
    int<lower=1> N;                  // Number of observations
	array[N] int group_id;
    array[N] int s_obs;
    // vector[N] cond;

    // Number of groups 
	int<lower=1> G;  // Total groups
  	int<lower=0, upper=1> use_posterior[G];  // Flag: 1 if posterior, 0 if estimate

    // Posterior estimates for groups 1-3, placeholders for remainder
  	int<lower=1> M;                   // Number of posterior samples per group
 	matrix[M, G] det_p_posterior;      // Posterior samples for groups 
}
// NOTE: following converges but is not sampling from posterior each iteration, only once when the model is compiled
transformed data {
    real det_p_data[G]; // Precomputed posterior samples for groups with use_posterior == 1
    int post_sample_idx[G]; // Stores sampled indices for posterior groups

    // Sample posterior values before inference begins
	for (g in 1:G) {
        post_sample_idx[g] = categorical_rng(rep_vector(1.0 / M, M));  // Randomly select a sample
        if (use_posterior[g] == 1) {
            det_p_data[g] = det_p_posterior[post_sample_idx[g], g]; // Use sampled posterior
        } else {
            det_p_data[g] = 0;  // Placeholder for estimated groups
        }
    }
}
parameters{
  	real<lower=0, upper=1> det_p_estimate[G];  // Detection probabilities for groups that don't have posterior
  	// real beta_cs;
    real surv_bar;
}
model{
    // priors
    // det_p_estimate ~ beta(2, 10); // Weak prior for estimated groups
  	for (g in 1:G) {
        if (use_posterior[g] == 0) {
            det_p_estimate[g] ~ beta(2, 10); // Weak prior for estimated groups
        }
    }

  	surv_bar ~ normal( 0 , 2 );

    vector[N] log_p;    
	for (i in 1:N) {
    	log_p[i] = log_inv_logit(surv_bar);
    	real det_p_i;  // Select the correct det_p for each group
    	if (use_posterior[group_id[i]] == 1) {
            det_p_i = det_p_data[group_id[i]];  // Use precomputed posterior sample
        } else {
            det_p_i = det_p_estimate[group_id[i]];  // Use estimated value
        }

	    if (s_obs[i] == 0)
    	  target += log_sum_exp(log1m_exp(log_p[i]), log_p[i] + log1m(det_p_i));

	    if (s_obs[i] == 1)
    	  target += log_p[i] + log(det_p_i);
	}
}
generated quantities {
    // Sample posterior values each iteration
	real det_p_data[G]; // Precomputed posterior samples for groups with use_posterior == 1
    int post_sample_idx[G]; // Stores sampled indices for posterior groups

    for (g in 1:G) {
        post_sample_idx[g] = categorical_rng(rep_vector(1.0 / M, M));  // Randomly select a sample
        if (use_posterior[g] == 1) {
            det_p_data[g] = det_p_posterior[post_sample_idx[g], g]; // Use sampled posterior
        } else {
            det_p_data[g] = 0;  // Placeholder for estimated groups
        }
    }

    // Export draws from input data (use_posterior == 1) or from posterior dist (use_posterior == 0)
    real det_p_out[G];  // Store final det_p_data values
    for (g in 1:G) {
        if (use_posterior[g] == 1) {
            det_p_out[g] = det_p_data[g];  // Save the detection probability input data used
        } else {
            det_p_out[g] = det_p_estimate[g];  // Use estimated value for reporting
        }
    }
}

