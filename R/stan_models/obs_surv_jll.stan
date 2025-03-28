// Joint likelihood state space model estimating true survival while 
// incorporating posterior estimates of detection probability
data{
    int<lower=1> N;                 // Number of observations
    int<lower=1> N_year;            // Number of years
    int<lower=1> N_stock;           // Number of stocks
    int<lower=1> N_det_id;          // Year-stock specific id for assigning posterior det probs
    int<lower=1> N_det_id_obs;         // Number of stock-years with posterior det prob data 
    
    array[N] int s_obs;             // Observed survival
    
    // Posterior detection probability estimates
    vector[N_det_id] logit_det_sd;
    vector[N_det_id] logit_det_mu;
    vector[N_det_id] use_posterior;  // Flag: 1 if posterior available

    array[N] int det_group_id;      // Year-stock specific id for assigning posterior det probs
    array[N] int yr;
    array[N] int pop_n;

    vector[N] cyer;
    vector[N] lipid;
    vector[N] size;
    vector[N] samp_date;
}
parameters{
    matrix[2, N_year] z_yr;
    matrix[4, N_stock] z_pop;
    real beta_dc;
    real beta_cf;
    real beta_cl;
    real beta_cyer;
    real beta_d_cyer;
    real beta_ls;
    real beta_fs;
    real beta_ds;
    real alpha_bar;
    cholesky_factor_corr[2] L_Rho_yr;
    vector<lower=0>[2] sigma_yr;
    cholesky_factor_corr[4] L_Rho_pop;
    vector<lower=0>[4] sigma_pop;
    real<lower=0> sigma_date;
    real<lower=0> sigma_size;
    real<lower=0> sigma_lipid;
    real<lower=0> sigma_cond;

    vector[N] cond;

    vector[N_det_id] logit_p_true_unobs;  // Unobserved systems
    vector[N_det_id_obs] logit_p_true_obs;  // Observed systems
}
transformed parameters{
    matrix[N_year, 2] alpha_yr;
    matrix[N_stock, 4] alpha_pop;
    alpha_pop = (diag_pre_multiply(sigma_pop, L_Rho_pop) * z_pop)';
    alpha_yr = (diag_pre_multiply(sigma_yr, L_Rho_yr) * z_yr)';

    vector<lower=0, upper = 1>[N_det_id] p;  // detection probability in real space
    for (g in 1:N_det_id) {
        if (use_posterior[g] == 1) {
            p[g] = inv_logit(logit_p_true_obs[g]);
        }
        if (use_posterior[g] == 0) {
            p[g] = inv_logit(logit_p_true_unobs[g]);
        } 
    }
}
model{
    vector[N] mu_date;
    vector[N] mu_cond;
    vector[N] mu_size;
    vector[N] mu_lipid;
    vector[N] logit_phi;
    vector[N] log_phi;
    
    sigma_date ~ exponential( 2 );
    sigma_size ~ exponential( 2 );
    sigma_cond ~ exponential( 2 );
    sigma_lipid ~ exponential( 2 );
    sigma_pop ~ exponential( 2 );
    L_Rho_pop ~ lkj_corr_cholesky( 2 );
    sigma_yr ~ exponential( 2 );
    L_Rho_yr ~ lkj_corr_cholesky( 2 );
    alpha_bar ~ normal( 0.5, 1 );
    beta_ds ~ normal( 0 , 0.5 );
    beta_fs ~ normal( 0 , 0.5 );
    beta_ls ~ normal( 0 , 0.5 );
    beta_d_cyer ~ normal( 0 , 0.5 );
    beta_cyer ~ normal( -0.5 , 0.5 );
    beta_dc ~ normal( 0 , 1 );
    beta_cl ~ normal( 0 , 1 );
    beta_cf ~ normal( 0 , 1 );
    to_vector( z_pop ) ~ normal( 0 , 1 );
    to_vector( z_yr ) ~ normal( 0 , 1 );

    // detection probability submodel
    logit_p_true_unobs ~ normal(-1.2, 1); // prior for systems with low observation efficiency, most mass < 0.4
    
    // replace with mean observation if posterior reliable
    for (g in 1:N_det_id_obs) {
        logit_det_mu[g] ~ normal(logit_p_true_obs[g], logit_det_sd[g]);
    }
    
    for ( i in 1:N ) {
        logit_phi[i] = alpha_bar + alpha_pop[pop_n[i], 4] + alpha_yr[yr[i], 2] + beta_ds * samp_date[i] + beta_fs * size[i] + beta_ls * lipid[i] + beta_cyer * cyer[i] + (beta_d_cyer * samp_date[i] * cyer[i]);

        log_phi[i] = log_inv_logit(logit_phi[i]);

        if (s_obs[i] == 0)
          target += log_sum_exp(log1m_exp(log_phi[i]), log_phi[i] + log1m(p[det_group_id[i]]));

        if (s_obs[i] == 1)
          target += log_phi[i] + log(p[det_group_id[i]]);
    }
        
    // observation model for capture date, size, lipid, and condition
    for ( i in 1:N ) {
        mu_date[i] = alpha_pop[pop_n[i], 1];
        mu_cond[i] = alpha_yr[yr[i], 1] + beta_dc * samp_date[i];      
        mu_size[i] = beta_cf * cond[i] + alpha_pop[pop_n[i], 2];
        mu_lipid[i] = beta_cl * cond[i] + alpha_pop[pop_n[i], 3];
    }

    samp_date ~ normal(mu_date, sigma_date);
    cond ~ normal(mu_cond, sigma_cond);
    size ~ normal(mu_size, sigma_size);
    lipid ~ normal(mu_lipid, sigma_lipid);
}
generated quantities{
    matrix[2,2] Rho_yr;
    matrix[4,4] Rho_pop;
    Rho_pop = multiply_lower_tri_self_transpose(L_Rho_pop);
    Rho_yr = multiply_lower_tri_self_transpose(L_Rho_yr);

    int s_obs_rep[N];

    for (i in 1:N) {
      real phi = inv_logit(alpha_bar + beta_ds * samp_date[i] + 
      beta_fs * size[i] + beta_ls * lipid[i] + beta_cyer * cyer[i] + 
       beta_d_cyer * samp_date[i] * cyer[i]);
      
      real p_det = p[det_group_id[i]];
    
      // Probability of observing survival: survived AND detected
      real prob_s_obs_1 = phi * p_det;
    
      s_obs_rep[i] = bernoulli_rng(prob_s_obs_1);
    }
}
