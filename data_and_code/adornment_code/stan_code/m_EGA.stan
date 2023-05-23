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
  int K[J]; // number of ordered categories in each outcome
  int resp[N_obs]; // index of which response for each observation
  int id[N_obs]; // index of which society
  int y[N_obs]; // response for each outcome
  matrix[N,N] cor_phy; // phylogenetic distances
  matrix[N,N] cor_EP; // temporal distances
  matrix[N,N] cor_geo; // geographic distances
}

parameters{
  ordered[max(K)-1] c[J];  // cutpoints for CA ordered logit, unpruned
  matrix[J,N] eta_z; // latent factor values, unscaled and uncorrelated
    vector<lower=0>[J] sigma_eta; // standard deviation of latent factors, sqrt(total variance)
  cholesky_factor_corr[J] L_eta;  // factor correlation matrix
  
    vector[J] alpha_eta; // intercepts for CA latent factors
    simplex[4] S[J]; // for the CA latent factors, we'll scale the phylogenetic, geographic, EP, and residual variance using a simplex
  
  matrix[N,J] phy_z;  // phylogenetic random effects, unscaled and uncorrelated
  matrix[N,J] EP_z;  // temporal random effects, unscaled and uncorrelated
  matrix[N,J] geo_z;  // temporal random effects, unscaled and uncorrelated
  
  cholesky_factor_corr[J] L_EP;  // correlations between EP effects
  cholesky_factor_corr[J] L_phy;  // correlations between phy effects
  cholesky_factor_corr[J] L_geo;  // correlations between geo effects
}

transformed parameters{
  matrix[N,J] eta_v;  // residual cor random effects, scaled and correlated (note the transpose here from eta_z)

  matrix[N,J] phy_v; // phy random effects, scaled and correlated
  matrix[N,J] EP_v; // EP random effects, scaled and correlated
  matrix[N,J] geo_v; // geographic random effects, scaled and correlated
  
  vector[J] sigma_res; // residual sd
  vector[J] sigma_phy; // phylogenetic sd
  vector[J] sigma_EP; // EP sd
  vector[J] sigma_geo; // geo sd
  
  // Scale each sd by its corresponding simplex parameter (sqrt so it's equivalent to scaling the variance)
  for (f in 1:J) {
  sigma_phy[f] = ( sqrt(S[f,1])*sigma_eta[f] );
  sigma_EP[f]  = ( sqrt(S[f,2])*sigma_eta[f] );
  sigma_geo[f] = ( sqrt(S[f,3])*sigma_eta[f] );
  sigma_res[f] = ( sqrt(S[f,4])*sigma_eta[f] );
  }
  
    // Phylogenetic covariance function
    for (f in 1:J) {
    matrix[N,N] L_phy_cov; // cholesky decomposition matrix
    
    L_phy_cov = cholesky_decompose(cor_phy);
    phy_v[,f] = L_phy_cov * phy_z[,f];
    }
    /// Also correlate phy effects across factors, matrix normal sampling
    phy_v = phy_v * diag_pre_multiply(sigma_phy, L_phy)';
    
    // Temporal covariance function
    for (f in 1:J) {
    matrix[N,N] L_EP_cov; // cholesky decomposition matrix
    
    L_EP_cov = cholesky_decompose(cor_EP);
    EP_v[,f] = L_EP_cov * EP_z[,f];
    }
    /// Also correlate EP effects across factors, matrix normal sampling
    EP_v = EP_v * diag_pre_multiply(sigma_EP, L_EP)';
  
    // Geographic covariance function
    for (f in 1:J) {
    matrix[N,N] geo_cov; // phylogenetic/geographic/temporal correlation matrices
    matrix[N,N] L_geo_cov; // cholesky decomposition matrix
    
    for ( q in 1:N )
    geo_cov[q,q] = 1 + 0.001; // adding small constant to diagonal to keep positive semidefinite
    
    L_geo_cov = cholesky_decompose(cor_geo);
    geo_v[,f] = L_geo_cov * geo_z[,f];
  }

  // Scale and correlate residual variance
  eta_v = (diag_pre_multiply(sigma_eta, L_eta) * eta_z)';
  
    // Center latent factors based on predictors
  for (n in 1:N) {
  for (f in 1:J) {
    eta_v[n,f] = alpha_eta[f] + eta_v[n,f] + phy_v[n,f] + EP_v[n,f] + geo_v[n,f];
  }
  
} // end transformed parameters block
}

model{
///// priors
// cutpoints for CA variables
for (j in 1:J) {
c[j,] ~ std_normal(); // N(0,1)
}

to_vector(eta_z) ~ std_normal(); 
L_eta ~ lkj_corr_cholesky(1);
sigma_eta ~ std_normal();

// simplex priors
for (f in 1:J) {
  S[f] ~ dirichlet(rep_vector(2, 4)); 
}

alpha_eta ~ std_normal();
to_vector(phy_z) ~ std_normal();
to_vector(EP_z) ~ std_normal();
to_vector(geo_z) ~ std_normal();
L_phy ~ lkj_corr_cholesky(2);
L_EP ~ lkj_corr_cholesky(2);
L_geo ~ lkj_corr_cholesky(2);

//////////////////////////////////////
// Run model loop of main CA variables
for (i in 1:N_obs) {
  vector[K[resp[i]]-1] cp; // pruned cutpoints
  
  for (cut in 1:(K[resp[i]]-1)) cp[cut] = c[resp[i],cut];
  
  // If response is missing, need to marginalize over possible states
  if (y[i] == -99) {
    vector[K[resp[i]]] lp;
    vector[K[resp[i]]] pk; // probability of each response
    vector[K[resp[i]]-1] cum_pr = inv_logit( cp - eta_v[id[i],resp[i]] ); // cumulative probability of each response
    
    // loop over possible responses
    for (k in 1:K[resp[i]]) {
      if (k == 1) pk[k] = cum_pr[k]; // cumulative prob for first step is the exact prob
      else if (k < K[resp[i]]) pk[k] = cum_pr[k] - cum_pr[k-1];
      else pk[k] = 1 - cum_pr[k-1];
      
      lp[k] = log(pk[k]) + ordered_logistic_lpmf( k | eta_v[id[i],resp[i]], cp );
    }
    target += log_sum_exp(lp);
  }
  
  // If response not missing
  else {
    y[i] ~ ordered_logistic( eta_v[id[i],resp[i]], cp );
  }
} // end loop over observations
  
} // end model block

generated quantities{
  matrix[J,J] Rho_eta;
  matrix[J,J] Rho_phy;
  matrix[J,J] Rho_EP;
  matrix[J,J] Rho_geo;
  
  Rho_eta = L_eta * L_eta';
  Rho_phy = L_phy * L_phy';
  Rho_EP = L_EP * L_EP';
  Rho_geo = L_geo * L_geo';
}

