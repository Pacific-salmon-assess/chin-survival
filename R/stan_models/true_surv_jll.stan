// Joint likelihood model estimating TRUE survival only based on 
// rethinking ulam function
data{
    int<lower=1> N;             // Number of observations
    int<lower=1> N_year;             // Number of years
    int<lower=1> N_stock;             // Number of stocks
    vector[N] logit_det_sd;
    vector[N] logit_det_mu;
    vector[N] cond;
    array[N] int det_group_id;
    array[N] int s_true;
    vector[N] cyer;
    vector[N] lipid;
    vector[N] size;
    vector[N] samp_date;
    array[N] int yr;
    array[N] int pop_n;
}
parameters{
    matrix[3,5] z_yr;
    matrix[4,5] z_pop;
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
}
transformed parameters{
    matrix[N_year,3] alpha_yr;
    matrix[N_stock,4] alpha_pop;
    alpha_pop = (diag_pre_multiply(sigma_pop, L_Rho_pop) * z_pop)';
    alpha_yr = (diag_pre_multiply(sigma_yr, L_Rho_yr) * z_yr)';
}
model{
    vector[N] mu_date;
    vector[N] mu_size;
    vector[N] mu_lipid;
    vector[N] p;
    
    Sigma_sl ~ exponential( 1 );
    sigma_date ~ exponential( 1 );
    sigma_pop ~ exponential( 1 );
    L_Rho_pop ~ lkj_corr_cholesky( 2 );
    sigma_yr ~ exponential( 1 );
    L_Rho_yr ~ lkj_corr_cholesky( 2 );
    Rho_sl ~ lkj_corr( 2 );
    alpha_bar ~ normal( 0 , 1 );
    beta_ds ~ normal( 0 , 0.5 );
    beta_fs ~ normal( 0 , 0.5 );
    beta_ls ~ normal( 0 , 0.5 );
    beta_d_cyer ~ normal( 0 , 0.5 );
    beta_cyer ~ normal( -0.5 , 0.5 );
    beta_df ~ normal( 0 , 1 );
    beta_dl ~ normal( 0 , 1 );
    to_vector( z_pop ) ~ normal( 0 , 1 );
    to_vector( z_yr ) ~ normal( 0 , 1 );
    
    for ( i in 1:N ) {
        p[i] = alpha_bar + alpha_pop[pop_n[i], 4] + alpha_yr[yr[i], 3] + beta_ds * samp_date[i] + beta_fs * size[i] + beta_ls * lipid[i] + beta_cyer * cyer[i] + (beta_d_cyer * samp_date[i] * cyer[i]);
        p[i] = inv_logit(p[i]);
    }
    
    s_true ~ binomial( 1 , p );
    
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
        YY ~ multi_normal( MU , quad_form_diag(Rho_sl , Sigma_sl) );
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
}
