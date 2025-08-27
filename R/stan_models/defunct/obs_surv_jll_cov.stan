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
    matrix[3, N_year] z_yr;
    matrix[4, N_stock] z_pop;
    real beta_dl;
    real beta_df;
    real beta_cyer;
    real beta_d_cyer;
    real beta_ls;
    real beta_fs;
    real beta_ds;
    real alpha_bar;
    corr_matrix[2] Rho_sl;
    cholesky_factor_corr[3] L_Rho_yr;
    vector<lower=0>[3] sigma_yr;
    cholesky_factor_corr[4] L_Rho_pop;
    vector<lower=0>[4] sigma_pop;
    real<lower=0> sigma_date;
    vector<lower=0>[2] Sigma_sl;

    vector[N_det_id] logit_p_true_unobs;  // Unobserved systems
    vector[N_det_id_obs] logit_p_true_obs;  // Observed systems
}
transformed parameters{
    matrix[N_year, 3] alpha_yr;
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
    vector[N] mu_size;
    vector[N] mu_lipid;
    vector[N] logit_phi;
    vector[N] log_phi;
    
    Sigma_sl ~ exponential( 2 );
    sigma_date ~ exponential( 2 );
    sigma_pop ~ exponential( 2 );
    L_Rho_pop ~ lkj_corr_cholesky( 2 );
    sigma_yr ~ exponential( 2 );
    L_Rho_yr ~ lkj_corr_cholesky( 2 );
    Rho_sl ~ lkj_corr( 2 );
    alpha_bar ~ normal( 0.9, 1.5 );
    beta_ds ~ normal( 0.25 , 0.5 );
    beta_fs ~ normal( 0 , 0.5 );
    beta_ls ~ normal( 0 , 0.5 );
    beta_d_cyer ~ normal( 0 , 0.5 );
    beta_cyer ~ normal( -0.5 , 0.5 );
    beta_df ~ normal( 0.25 , 1 );
    beta_dl ~ normal( 0.25 , 1 );
    to_vector( z_pop ) ~ normal( 0 , 1 );
    to_vector( z_yr ) ~ normal( 0 , 1 );

    // detection probability submodel
    logit_p_true_unobs ~ normal(-1.2, 1); // prior for systems with low observation efficiency, most mass < 0.4
    
    // replace with mean observation if posterior reliable
    for (g in 1:N_det_id_obs) {
        logit_det_mu[g] ~ normal(logit_p_true_obs[g], logit_det_sd[g]);
    }
    
    for ( i in 1:N ) {
        logit_phi[i] = alpha_bar + alpha_pop[pop_n[i], 4] + alpha_yr[yr[i], 3] + beta_ds * samp_date[i] + beta_fs * size[i] + beta_ls * lipid[i] + beta_cyer * cyer[i] + (beta_d_cyer * samp_date[i] * cyer[i]);

        log_phi[i] = log_inv_logit(logit_phi[i]);

        if (s_obs[i] == 0)
          target += log_sum_exp(log1m_exp(log_phi[i]), log_phi[i] + log1m(p[det_group_id[i]]));

        if (s_obs[i] == 1)
          target += log_phi[i] + log(p[det_group_id[i]]);
    }
        
    for ( i in 1:N ) {
        mu_lipid[i] = alpha_yr[yr[i], 2] + alpha_pop[pop_n[i], 3] + beta_dl * samp_date[i];
    }
    
    for ( i in 1:N ) {
        mu_size[i] = alpha_yr[yr[i], 1] + alpha_pop[pop_n[i], 2] + beta_df * samp_date[i];
    }
    
    {
        array[N] vector[2] YY;
        array[N] vector[2] MU;
        for ( j in 1:N ) MU[j] = [ mu_size[j] , mu_lipid[j] ]';
        for ( j in 1:N ) YY[j] = [ size[j] , lipid[j] ]';
        YY ~ multi_normal_cholesky(MU, diag_pre_multiply(Sigma_sl, cholesky_decompose(Rho_sl)));
    }

    for ( i in 1:N ) {
        mu_date[i] = alpha_pop[pop_n[i], 1];
    }

    samp_date ~ normal( mu_date , sigma_date );
}
generated quantities{
    matrix[3,3] Rho_yr;
    matrix[4,4] Rho_pop;
    Rho_pop = multiply_lower_tri_self_transpose(L_Rho_pop);
    Rho_yr = multiply_lower_tri_self_transpose(L_Rho_yr);

    int s_obs_rep[N];
    int s_true_rep[N];

    for (i in 1:N) {
      real phi = inv_logit(alpha_bar + alpha_pop[pop_n[i], 4] + alpha_yr[yr[i], 3] + beta_ds * samp_date[i] + beta_fs * size[i] + beta_ls * lipid[i] + beta_cyer * cyer[i] + (beta_d_cyer * samp_date[i] * cyer[i]));
      
      real p_det = p[det_group_id[i]];
    
      // Probability of observing survival: survived AND detected
      real prob_s_obs_1 = phi * p_det;
    
      s_obs_rep[i] = bernoulli_rng(prob_s_obs_1);
      s_true_rep[i] = bernoulli_rng(phi);
    }
}
