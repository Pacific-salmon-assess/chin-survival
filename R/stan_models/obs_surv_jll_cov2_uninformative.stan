// Joint likelihood state space model estimating true survival while 
// incorporating posterior estimates of detection probability
// Adds injury effects
data{
    int<lower=1> N;                 // Number of observations
    int<lower=1> N_year;            // Number of years
    int<lower=1> N_stock;           // Number of stocks
    int<lower=1> N_det_id;          // Year-stock specific id for assigning posterior det probs
    int<lower=1> N_det_id_obs;         // Number of stock-years with posterior det prob data 
    vector[3] alpha_i;              // Initial values for dirichlet

    array[N] int s_obs;             // Observed survival
    
    // Posterior detection probability estimates
    vector[N_det_id] logit_det_sd;
    vector[N_det_id] logit_det_mu;
    vector[N_det_id] use_posterior;  // Flag: 1 if posterior available

    array[N] int det_group_id;      // Year-stock specific id for assigning posterior det probs
    array[N] int yr;
    array[N] int stk_n;
    array[N] int inj;

    vector[N] cyer_z;
    vector[N] lipid_z;
    vector[N] fl_z;
    vector[N] day_z;
}
parameters{
    matrix[3, N_year] z_yr;
    matrix[4, N_stock] z_stk;
    real beta_dl;
    real beta_df;
    real beta_cs;
    real beta_ds_cs;
    real beta_ls;
    real beta_fs;
    real beta_ds;
    real beta_is;
    real alpha_s;
    real alpha_f;
    real alpha_l;
    corr_matrix[2] Rho_sl;
    cholesky_factor_corr[3] L_Rho_yr;
    vector<lower=0>[3] sigma_yr;
    cholesky_factor_corr[4] L_Rho_stk;
    vector<lower=0>[4] sigma_stk;
    real<lower=0> sigma_day;
    vector<lower=0>[2] Sigma_fl;

    vector[N_det_id] logit_p_true_unobs;  // Unobserved systems
    vector[N_det_id_obs] logit_p_true_obs;  // Observed systems

    simplex[3] delta; // for calculating categorical injury effects
}
transformed parameters{
    matrix[N_year, 3] alpha_yr;
    matrix[N_stock, 4] alpha_stk;
    vector[4] delta_inj;
    alpha_stk = (diag_pre_multiply(sigma_stk, L_Rho_stk) * z_stk)';
    alpha_yr = (diag_pre_multiply(sigma_yr, L_Rho_yr) * z_yr)';

    delta_inj = append_row(0, delta);

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
    vector[N] mu_day;
    vector[N] mu_fl;
    vector[N] mu_lipid;
    vector[N] logit_phi;
    vector[N] log_phi;
    
    Sigma_fl ~ exponential( 2 );
    sigma_day ~ exponential( 2 );
    sigma_stk ~ exponential( 2 );
    L_Rho_stk ~ lkj_corr_cholesky( 2 );
    sigma_yr ~ exponential( 2 );
    L_Rho_yr ~ lkj_corr_cholesky( 2 );
    Rho_sl ~ lkj_corr( 2 );
    alpha_s ~ normal( 0.9, 1.5 );
    alpha_f ~ normal( 0 , 1 );
    alpha_l ~ normal( 0 , 1 );
    beta_ds ~ normal( 0 , 0.5 );
    beta_fs ~ normal( 0 , 0.5 );
    beta_ls ~ normal( 0 , 0.5 );
    beta_is ~ normal( 0 , 0.5 );
    beta_ds_cs ~ normal( 0 , 0.5 );
    beta_cs ~ normal( 0 , 0.5 );
    beta_df ~ normal( 0 , 1 );
    beta_dl ~ normal( 0 , 1 );
    to_vector( z_stk ) ~ normal( 0 , 1 );
    to_vector( z_yr ) ~ normal( 0 , 1 );

    delta ~ dirichlet(alpha_i);

    // detection probability submodel
    logit_p_true_unobs ~ normal(-1.2, 1); // prior for systems with low observation efficiency, most mass < 0.4
    
    // replace with mean observation if posterior reliable
    for (g in 1:N_det_id_obs) {
        logit_det_mu[g] ~ normal(logit_p_true_obs[g], logit_det_sd[g]);
    }
    
    for ( i in 1:N ) {
        logit_phi[i] = alpha_s + alpha_stk[stk_n[i], 4] + alpha_yr[yr[i], 3] + beta_ds * day_z[i] + beta_fs * fl_z[i] + beta_ls * lipid_z[i] + beta_cs * cyer_z[i] + (beta_ds_cs * day_z[i] * cyer_z[i]) + beta_is * sum(delta_inj[1:inj[i]]);

        log_phi[i] = log_inv_logit(logit_phi[i]);

        if (s_obs[i] == 0)
          target += log_sum_exp(log1m_exp(log_phi[i]), log_phi[i] + log1m(p[det_group_id[i]]));

        if (s_obs[i] == 1)
          target += log_phi[i] + log(p[det_group_id[i]]);
    }
        
    for ( i in 1:N ) {
        mu_lipid[i] = alpha_l + alpha_yr[yr[i], 2] + alpha_stk[stk_n[i], 3] + beta_dl * day_z[i];
    }
    
    for ( i in 1:N ) {
        mu_fl[i] = alpha_f + alpha_yr[yr[i], 1] + alpha_stk[stk_n[i], 2] + beta_df * day_z[i];
    }
    
    {
        array[N] vector[2] YY;
        array[N] vector[2] MU;
        for ( j in 1:N ) MU[j] = [ mu_fl[j] , mu_lipid[j] ]';
        for ( j in 1:N ) YY[j] = [ fl_z[j] , lipid_z[j] ]';
        YY ~ multi_normal_cholesky(MU, diag_pre_multiply(Sigma_fl, cholesky_decompose(Rho_sl)));
    }

    for ( i in 1:N ) {
        mu_day[i] = alpha_stk[stk_n[i], 1];
    }

    day_z ~ normal( mu_day , sigma_day );
}
generated quantities{
    matrix[3,3] Rho_yr;
    matrix[4,4] Rho_stk;
    Rho_stk = multiply_lower_tri_self_transpose(L_Rho_stk);
    Rho_yr = multiply_lower_tri_self_transpose(L_Rho_yr);

    int s_obs_rep[N];
    int s_true_rep[N];

    for (i in 1:N) {
      real phi = inv_logit(alpha_s + alpha_stk[stk_n[i], 4] + alpha_yr[yr[i], 3] + beta_ds * day_z[i] + beta_fs * fl_z[i] + beta_ls * lipid_z[i] + beta_cs * cyer_z[i] + (beta_ds_cs * day_z[i] * cyer_z[i]));
      
      real p_det = p[det_group_id[i]];
    
      // Probability of observing survival: survived AND detected
      real prob_s_obs_1 = phi * p_det;
    
      s_obs_rep[i] = bernoulli_rng(prob_s_obs_1);
      s_true_rep[i] = bernoulli_rng(phi);
    }

    // loglikelihood for loo
    vector[N] log_lik;

    for (i in 1:N) {
      real phi = inv_logit(
        alpha_s +
        alpha_stk[stk_n[i], 4] +
        alpha_yr[yr[i], 3] +
        beta_ds * day_z[i] +
        beta_fs * fl_z[i] +
        beta_ls * lipid_z[i] +
        beta_cs * cyer_z[i] +
        (beta_ds_cs * day_z[i] * cyer_z[i]) +
        beta_is * sum(delta_inj[1:inj[i]])
      );
    
      real p_det = p[det_group_id[i]];
    
      if (s_obs[i] == 1)
        log_lik[i] = log(phi) + log(p_det);
      else
        log_lik[i] = log1m(phi * p_det);
    }
}
