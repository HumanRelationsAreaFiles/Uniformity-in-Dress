functions{
  /* compute monotonic effects, function lifted from brms package
  * Args:
    *   scale: a simplex parameter
  *   i: index to sum over the simplex
  * Returns:
    *   a scalar between 0 and 1
  */
    real mo(vector scale, int i) {
      if (i > 0) return sum(scale[1:i]);
      else return 0;
    }
}

data{
  int N; // number of societies
  int J; // number of outcomes
  int N_obs; // number of observations (N societies x J outcomes)
  int N_TL; // number of TL observations
  int N_RS; // number of RS observations
  int TL_j[N_TL]; // indicator of which tight-loose variable
  int RS_j[N_RS]; // indicator of which resource stress variable
  int K[J]; // number of ordered categories in each CA outcome
  int factor[J]; // indicator of which factor an outcome belongs to
  int fix[J]; // indicator of whether to fix the loading
  int N_f; // num latent factors
  int resp[N_obs]; // index of which response for each observation
  int id[N_obs]; // index of which society
  int y[N_obs]; // response for each CA outcome
  int y_TL[N_TL]; // TL codes
  int y_RS[N_RS]; // RS codes
  real TL_diff[N_TL]; // TL diff from nearest integer
  int id_TL[N_TL]; // society identifier 
  int id_RS[N_RS]; // same
  int records[N]; // indicator of whether there are written records (1) or not (0)
  int strat[N]; // indicator of whether society is relatively egalitarian (0) or stratified (1)
  int comm_size[N]; // mean size of local community, ordinal variable
  int pol_int[N]; // political integration, ordinal variable
  matrix[N,N] cor_phy; // phylogenetic distances
  matrix[N,N] cor_EP; // temporal distances
  matrix[N,N] cor_geo; // geographic distances
}

parameters{
  ordered[max(K)-1] c[J];  // cutpoints for CA ordered logit, unpruned
  ordered[max(y_RS)-1] c_RS[max(RS_j)]; // cutpoints for RS variables
  ordered[max(y_RS)-1] c_TL[max(TL_j)]; // cutpoints for TL variables

  matrix[(N_f + 1),N] eta_z; // latent factor values, unscaled and uncorrelated
  vector<lower=0>[N_f + 1] sigma_eta; // standard deviation of latent factors, sqrt(total variance)
  simplex[4] S[N_f + 1]; // for the CA latent factors, we'll scale the phylogenetic, geographic, EP, and residual variance using a simplex
  cholesky_factor_corr[N_f + 1] L_eta;  // factor correlation matrix
  vector[J-N_f] lambda; // factor loadings, -1 for each that we fix equal to 1
  
  vector[N] eta_RS_z; // latent RS factor values, unscaled
  real<lower=0> sigma_eta_RS; // standard deviation of RS latent factor

  vector[N_f + 1] alpha_eta; // intercepts for CA latent factors and TL
  
  vector[N_f + 1] b_comm; // overall effect of commnuity size
  simplex[max(comm_size)] scale_comm[N_f + 1]; //  monotonic intervals of commnuity size 
  
  vector[max(RS_j)-1] lambda_RS; // factor loadings for resource stress, leaving out chronic scarcity which we fix equal to 1
  vector[max(TL_j)-1] lambda_TL; // factor loadings for TL, leaving out 1 which we fix

  matrix[N,N_f + 1] phy_z;  // phylogenetic random effects, unscaled and uncorrelated
  matrix[N,N_f + 1] EP_z;  // temporal random effects, unscaled and uncorrelated
  matrix[N,N_f + 1] geo_z;  // temporal random effects, unscaled and uncorrelated
  
  cholesky_factor_corr[N_f + 1] L_EP;  // correlations between EP effects
  cholesky_factor_corr[N_f + 1] L_phy;  // correlations between phy effects
  cholesky_factor_corr[N_f + 1] L_geo;  // correlations between geo effects
}

transformed parameters{
  matrix[N,(N_f + 1)] eta_v;  // residual cor random effects, scaled and correlated (note the transpose here from eta_z)
  vector[N] eta_RS_v; // 
  matrix[J,N_f] lambda_hat; // matrix of factor loadings that we'll fill with unknown variables, 1's and 0's. Only for the CA factors
  
  matrix[N,N_f + 1] phy_v; // phy random effects, scaled and correlated
  matrix[N,N_f + 1] EP_v; // EP random effects, scaled and correlated
  matrix[N,N_f + 1] geo_v; // geographic random effects, scaled and correlated
  
  vector[N_f + 1] sigma_res; // residual sd
  vector[N_f + 1] sigma_phy; // phylogenetic sd
  vector[N_f + 1] sigma_EP; // EP sd
  vector[N_f + 1] sigma_geo; // geo sd
  
  // Scale each sd by its corresponding simplex parameter (sqrt so it's equivalent to scaling the variance)
  for (f in 1:N_f + 1) {
  sigma_phy[f] = ( sqrt(S[f,1])*sigma_eta[f] );
  sigma_EP[f]  = ( sqrt(S[f,2])*sigma_eta[f] );
  sigma_geo[f] = ( sqrt(S[f,3])*sigma_eta[f] );
  sigma_res[f] = ( sqrt(S[f,4])*sigma_eta[f] );
  }
 
 
 // Scale RS factor
 eta_RS_v = eta_RS_z * sigma_eta_RS;
  
  // Phylogenetic covariance function
    for (f in 1:N_f + 1) {
    matrix[N,N] L_phy_cov; // cholesky decomposition matrix

    L_phy_cov = cholesky_decompose(cor_phy);
    phy_v[,f] = L_phy_cov * phy_z[,f];
    }
    /// Also correlate phy effects across factors, matrix normal sampling
    phy_v = phy_v * diag_pre_multiply(sigma_phy, L_phy)';
  
    // Temporal covariance function
    for (f in 1:N_f + 1) {
    matrix[N,N] L_EP_cov; // cholesky decomposition matrix
    
    L_EP_cov = cholesky_decompose(cor_EP);
    EP_v[,f] = L_EP_cov * EP_z[,f];
    }
    /// Also correlate EP effects across factors, matrix normal sampling
    EP_v = EP_v * diag_pre_multiply(sigma_EP, L_EP)';
  
    // Geographic covariance function
    for (f in 1:N_f + 1) {
    matrix[N,N] L_geo_cov; // cholesky decomposition matrix
    
    L_geo_cov = cholesky_decompose(cor_geo);
    geo_v[,f] = L_geo_cov * geo_z[,f];
  }
    /// Also correlate geo effects across factors, matrix normal sampling
    geo_v = geo_v * diag_pre_multiply(sigma_geo, L_geo)';
  
  // Scale and correlate residual variance
  eta_v = (diag_pre_multiply(sigma_res, L_eta) * eta_z)';
  
  
  // Center latent factors based on predictors
  for (n in 1:N) {
  for (f in 1:N_f + 1) {
    eta_v[n,f] = alpha_eta[f] + eta_v[n,f] + phy_v[n,f] + EP_v[n,f] + geo_v[n,f] + b_comm[f]*mo(scale_comm[f],comm_size[n]);
  }
  }
  
  {
  int ticker = 1; // keep track of how we fill the matrix
  for (f in 1:N_f)
  for (j in 1:J) {
    if (factor[j] == f && fix[j] == 0 ) {
      lambda_hat[j,f] = lambda[ticker];
      ticker = ticker + 1;
    }
    
    else if (factor[j] == f && fix[j] == 1) {
      lambda_hat[j,f] = 1;
    }
    
    else lambda_hat[j,f] = 0;
  }
}
  
} // end transformed parameters block

model{
///// priors
// cutpoints for CA variables
for (j in 1:J) {
c[j,] ~ std_normal(); // N(0,1)
}

// cutpoints for resource stress variables
for (j in 1:max(RS_j)) {
  c_RS[j,] ~ std_normal();
}

// cutpoints for TL variables
for (j in 1:max(TL_j)) {
  c_TL[j,] ~ std_normal();
}

alpha_eta ~ std_normal();
b_comm ~ std_normal();

lambda ~ std_normal();
lambda_RS ~ std_normal();
lambda_TL ~ std_normal();
to_vector(eta_z) ~ std_normal();
eta_RS_z ~ std_normal();
L_eta ~ lkj_corr_cholesky(2);
sigma_eta ~ exponential(1);
sigma_eta_RS ~ exponential(1);

// simplex priors
for (f in 1:N_f) {
  S[f] ~ dirichlet(rep_vector(2, 4)); 
}
for (f in 1:(N_f+1)) {
  scale_comm[f] ~ dirichlet(rep_vector(2, max(comm_size))); 
}

to_vector(phy_z) ~ std_normal();
to_vector(EP_z) ~ std_normal();
to_vector(geo_z) ~ std_normal();
L_phy ~ lkj_corr_cholesky(2);
L_EP ~ lkj_corr_cholesky(2);
L_geo ~ lkj_corr_cholesky(2);

  
//// Resource stress models ///////////////
  for (i in 1:N_RS) {
    real mu; // expected value latent scale
    
    if (RS_j[i] == 1) mu = 1*eta_RS_v[id_RS[i]]; // fixed loading for identifiability
    else mu = lambda_RS[RS_j[i]-1]*eta_RS_v[id_RS[i]];
    
    // if RS variable missing
    if (y_RS[i] == -99) {
    vector[max(y_RS)] lp; // log likelihood of each possible state
    vector[max(y_RS)] pk; // prob of each possible state
    vector[max(y_RS)-1] cum_pr = inv_logit( c_RS[RS_j[i],] - mu ); // cumulative probability of each response
    
    // loop over possible responses
    for (k in 1:max(y_RS)) {
      if (k == 1) pk[k] = cum_pr[k]; // cumulative prob for first step is the exact prob
      else if (k < max(y_RS)) pk[k] = cum_pr[k] - cum_pr[k-1];
      else pk[k] = 1 - cum_pr[k-1];
      
      lp[k] = log(pk[k]) + ordered_logistic_lpmf( k | mu, c_RS[RS_j[i],] );
    }
    target += log_sum_exp(lp);
    }
    
    // not missing
    else{
      y_RS[i] ~ ordered_logistic( mu, c_RS[RS_j[i],] );
    }
    }
//////////////////////////////////////////

//// TL models ///////////////////////////
  for (i in 1:N_TL) {
    real mu; // expected value latent scale
    
    if (TL_j[i] == 1) mu = 1*eta_v[id_TL[i],(N_f + 1)]; // fixed loading for identifiability
    else mu = lambda_TL[TL_j[i]-1]*eta_v[id_TL[i],(N_f + 1)];
    
    if (y_TL[i] == -99) {
    vector[max(y_TL)] lp; // log likelihood of each possible state
    vector[max(y_TL)] pk; // prob of each possible state
    vector[max(y_TL)-1] cum_pr = inv_logit( c_TL[TL_j[i],] - mu ); // cumulative probability of each response
    
    // loop over possible responses
    for (k in 1:max(y_TL)) {
      if (k == 1) pk[k] = cum_pr[k]; // cumulative prob for fiTLt step is the exact prob
      else if (k < max(y_TL)) pk[k] = cum_pr[k] - cum_pr[k-1];
      else pk[k] = 1 - cum_pr[k-1];
      
      lp[k] = log(pk[k]) + ordered_logistic_lpmf( k | mu, c_TL[TL_j[i],] );
    }
    target += log_sum_exp(lp);
    }
    
    // If there's uncertainty in the TL code, mix over two possibilities
    else if (TL_diff[i] > 0) {
    real pr; // probabiltiy that the code is actually 1 greater than the rounded value
    target += log_mix( TL_diff[i], ordered_logistic_lpmf( y_TL[i]+1 | mu, c_TL[TL_j[i],] ), ordered_logistic_lpmf( y_TL[i] | mu, c_TL[TL_j[i],] ) );
    }
    
    // If there's uncertainty in the TL code, mix over two possibilities
    else if (TL_diff[i] < 0) {
    real pr; // probabiltiy that the code is actually 1 less than the rounded value
    target += log_mix( (1 + TL_diff[i]), ordered_logistic_lpmf( y_TL[i]-1 | mu, c_TL[TL_j[i],] ), ordered_logistic_lpmf( y_TL[i] | mu, c_TL[TL_j[i],] ) );
    }
    
    // If TL codes observed and is an integer
    else {
      y_TL[i] ~ ordered_logistic( mu, c_TL[TL_j[i],] );
    }
    
    }
//////////////////////////////////////

//////////////////////////////////////
// Run model loop of main CA variables
for (i in 1:N_obs) {
  vector[K[resp[i]]-1] cp; // pruned cutpoints
  
  for (cut in 1:(K[resp[i]]-1)) cp[cut] = c[resp[i],cut];
  
  // If response is missing, need to marginalize over possible states
  if (y[i] == -99) {
    vector[K[resp[i]]] lp;
    vector[K[resp[i]]] pk; // probability of each response
    vector[K[resp[i]]-1] cum_pr = inv_logit( cp ); // cumulative probability of each response
    
    // loop over possible responses
    for (k in 1:K[resp[i]]) {
      if (k == 1) pk[k] = cum_pr[k]; // cumulative prob for first step is the exact prob
      else if (k < K[resp[i]]) pk[k] = cum_pr[k] - cum_pr[k-1];
      else pk[k] = 1 - cum_pr[k-1];
      
      lp[k] = log(pk[k]) + ordered_logistic_lpmf( k | eta_v[id[i],1]*lambda_hat[resp[i],1] + eta_v[id[i],2]*lambda_hat[resp[i],2] + eta_v[id[i],3]*lambda_hat[resp[i],3] + eta_v[id[i],4]*lambda_hat[resp[i],4], cp );
    }
    target += log_sum_exp(lp);
  }
  
  // If response not missing
  else {
    y[i] ~ ordered_logistic( eta_v[id[i],1]*lambda_hat[resp[i],1] + eta_v[id[i],2]*lambda_hat[resp[i],2] + eta_v[id[i],3]*lambda_hat[resp[i],3] + eta_v[id[i],4]*lambda_hat[resp[i],4], cp );
  }
} // end loop over observations
  
} // end model block

generated quantities{
  matrix[(N_f + 1),(N_f + 1)] Rho_eta;
  matrix[(N_f + 1),(N_f + 1)] Rho_phy;
  matrix[(N_f + 1),(N_f + 1)] Rho_EP;
  matrix[(N_f + 1),(N_f + 1)] Rho_geo;
  
  Rho_eta = L_eta * L_eta';
  Rho_phy = L_phy * L_phy';
  Rho_EP = L_EP * L_EP';
  Rho_geo = L_geo * L_geo';
}

