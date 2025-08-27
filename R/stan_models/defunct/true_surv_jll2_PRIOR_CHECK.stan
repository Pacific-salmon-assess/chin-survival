// Joint likelihood model estimating TRUE survival only based on 
// rethinking ulam function
data{
    int<lower=1> N;             // Number of observations
    int<lower=1> N_year;             // Number of years
    int<lower=1> N_stock;             // Number of stocks
    
    array[N] int s_true;
    
    vector[N] cyer;
    vector[N] lipid;
    vector[N] size;
    vector[N] samp_date;

    array[N] int yr;
    array[N] int pop_n;
}
parameters{
    matrix[4,N_stock] z_pop;
    real beta_dl;
    real beta_df;
    real beta_cyer;
    real beta_d_cyer;
    real beta_ls;
    real beta_fs;
    real beta_ds;
    real alpha_bar;
    corr_matrix[2] Rho_sl;
    cholesky_factor_corr[4] L_Rho_pop;
    vector<lower=0>[4] sigma_pop;
    real<lower=0> sigma_date;
    vector<lower=0>[2] Sigma_sl;
}
transformed parameters{
    matrix[N_stock,4] alpha_pop;
    alpha_pop = (diag_pre_multiply(sigma_pop, L_Rho_pop) * z_pop)';
}
model{
    vector[N] mu_date;
    vector[N] mu_size;
    vector[N] mu_lipid;
    vector[N] p;
    
    Sigma_sl ~ exponential( 2 );
    sigma_date ~ exponential( 2 );
    sigma_pop ~ exponential( 2 );
    L_Rho_pop ~ lkj_corr_cholesky( 2 );
    Rho_sl ~ lkj_corr( 2 );
    alpha_bar ~ normal( 0.9, 1.5);
    beta_ds ~ normal( 0.25 , 0.5 );
    beta_fs ~ normal( 0 , 0.5 );
    beta_ls ~ normal( 0 , 0.5 );
    beta_d_cyer ~ normal( 0 , 0.5 );
    beta_cyer ~ normal( -0.5 , 0.5 );
    beta_df ~ normal( 0.25 , 1);
    beta_dl ~ normal( 0.25 , 1);
    to_vector( z_pop ) ~ normal( 0 , 1 );
    
    for ( i in 1:N ) {
        p[i] = alpha_bar + alpha_pop[pop_n[i], 4] + beta_ds * samp_date[i] + beta_fs * size[i] + beta_ls * lipid[i] + beta_cyer * cyer[i] + (beta_d_cyer * samp_date[i] * cyer[i]);
        p[i] = inv_logit(p[i]);
    }
    
    // s_true ~ binomial( 1 , p );
    
    for ( i in 1:N ) {
        mu_lipid[i] = alpha_pop[pop_n[i], 3] + beta_dl * samp_date[i];
    }
    
    for ( i in 1:N ) {
        mu_size[i] = alpha_pop[pop_n[i], 2] + beta_df * samp_date[i];
    }
    
    {
        array[N] vector[2] YY;
        array[N] vector[2] MU;
        for ( j in 1:N ) MU[j] = [ mu_size[j] , mu_lipid[j] ]';
        for ( j in 1:N ) YY[j] = [ size[j] , lipid[j] ]';
     //   YY ~ multi_normal_cholesky(MU, diag_pre_multiply(Sigma_sl, cholesky_decompose(Rho_sl)));
    }

    for ( i in 1:N ) {
        mu_date[i] = alpha_pop[pop_n[i], 1];
    }

    // samp_date ~ normal( mu_date , sigma_date );
}
generated quantities{
    matrix[4,4] Rho_pop;
    Rho_pop = multiply_lower_tri_self_transpose(L_Rho_pop);

    int s_obs_rep[N];

    for (i in 1:N) {
      real phi = inv_logit(alpha_bar + beta_ds * samp_date[i] + alpha_pop[pop_n[i], 4] + beta_fs * size[i] + beta_ls * lipid[i] + beta_cyer * cyer[i] + (beta_d_cyer * samp_date[i] * cyer[i]));
      
      s_obs_rep[i] = bernoulli_rng(phi);
    }
}
