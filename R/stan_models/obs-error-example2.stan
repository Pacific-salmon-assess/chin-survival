// Model based on rethinking ulam function

data{
    int<lower=1> N;                  // Number of observations
	array[N] int group_id;
    array[N] int s_obs;

    // Number of groups 
	int<lower=1> G;  // Total groups
  	int<lower=1> G_obs;  // Total groups with posterior estimates
    int<lower=0, upper=1> use_posterior[G];  // Flag: 1 if posterior, 0 if estimate

    // Posterior estimates for groups 1-3, placeholders for remainder
  	vector[G_obs] logit_p_obs;     // Mean logit detection prob from posterior
    vector<lower=0>[G_obs] logit_p_obs_sd;     // SD logit detection prob from posterior
}
parameters{
    vector[G] logit_p_true_unobs;  // Detection probabilities from unobserved systems
    vector[G_obs] logit_p_true_obs;  // Detection probabilities for observed systems
    real surv_bar;
}
transformed parameters{
    vector<lower=0, upper = 1>[G] p;  // p in real space
    for (g in 1:G) {
        if (use_posterior[g] == 1) {
            p[g] = inv_logit(logit_p_true_obs[g]);
        }
        if (use_posterior[g] == 0) {
            p[g] = inv_logit(logit_p_true_unobs[g]);
        } 
    }
}
model{
    // logit_p_true_obs ~ normal(0, 1.25); // Weak prior centered on 0.5 after inv logit
  	logit_p_true_unobs ~ normal(-1.7, 1); // prior for systems with low observation efficiency, most mass < 0.4
    
    // replace with mean observation if posterior reliable
    for (g in 1:G_obs) {
        logit_p_obs[g] ~ normal(logit_p_true_obs[g], logit_p_obs_sd[g]);
    }

  	surv_bar ~ normal(0, 1.5);

    vector[N] log_surv;    
	for (i in 1:N) {
    	log_surv[i] = log_inv_logit(surv_bar);

	    if (s_obs[i] == 0)
    	  target += log_sum_exp(log1m_exp(log_surv[i]), log_surv[i] + log1m(p[group_id[i]]));

	    if (s_obs[i] == 1)
    	  target += log_surv[i] + log(p[group_id[i]]);
	}
}

