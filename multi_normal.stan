//Multivariate Normal distribution

data {
  int<lower=0> N; // length of data set
  int<lower=0> n_prot; //number of proteins
  vector[n_prot] y[N];  // data
}

parameters {
  vector[n_prot] mu; 
  cov_matrix[n_prot] Sigma;
}

model {
  
  matrix[n_prot, n_prot] identity;
  // prior distribution
  //mu ~ normal(0, 5000);
  
  // for(k in 1:n_prot){
  // Sigma[k] ~ normal(0,1000);
  // }
  
  //identity = diag_matrix(rep_vector(1.0,n_prot)); 
  //Sigma ~ inv_wishart(n_prot+1, identity);
  // likelihood
  
  for(n in 1:N) {
    
  //y[n] ~ multi_normal(mu, Sigma);
  
  target += multi_normal_cholesky_lpdf(y[n] | mu, Sigma);
  
  }
}
