// Model based on rethinking ulam function

data{
    array[1000] int group_id;
    array[1000] int s_obs;
    vector[1000] cond;
}
parameters{
  	 vector<lower=0,upper=1>[2] det_p; // Detection probability for two groups
     real beta_cs;
     real surv_bar;
}
model{
    vector[1000] log_p;
    surv_bar ~ normal( 0 , 2 );
    beta_cs ~ normal( 0 , 0.5 );
    for ( i in 1:1000 ) {
        log_p[i] = log_inv_logit(surv_bar + beta_cs * cond[i]);
    }
    det_p[1] ~ beta(22, 2);  // Group 1: High detection (~0.9)
  	det_p[2] ~ beta(2, 10);   // Group 2: Lower detection (mostly <0.5)
    for (i in 1:1000) {
    real det_p_i = det_p[group_id[i]];  // Extract det_p value for individual i

    if (s_obs[i] == 0)
        target += log_sum_exp(log1m_exp(log_p[i]), log_p[i] + log1m(det_p_i));

    if (s_obs[i] == 1)
        target += log_p[i] + log(det_p_i);
}
}
