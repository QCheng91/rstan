data {
  int<lower=0> N; // Sample size
  int<lower=0> num_proteins; // Sample size
  int x[N]; // Value for each sample on each dimension
  int y[N, num_proteins]; // Value for each sample on each dimension
}

parameters {
  real<lower=0> mu;
  real<lower=0> phi;
  vector<lower=0>[num_proteins] a;
  matrix<lower=0>[N, num_proteins] lambda;
}

model {
  a   ~ normal(0, 0.1);
  mu  ~ normal(0, 5000);
  phi ~ normal(0, 5000);
  
  for (n in 1:num_proteins) {
    
    lambda[:,n] ~ normal(0.,5000.);
  }
    
  // for(n in 1:N)
   
  // target += log(poisson_lpmf( y[n,] | lambda[n,]));
  
  for (m in 1:num_proteins)
  
  y[,m] ~ poisson(lambda[,m]);
  
  x   ~ neg_binomial_2(mu + lambda*a, phi);
}

// data {
//   int<lower=0> N; // Sample size
//   int x[N]; // Value for each sample on each dimension
// }
// 
// parameters {
//   real<lower=0, upper=10> mu1;
//   real<lower=20, upper=1000> mu2;
//   
//   real<lower=0> phi1;
//   real<lower=0> phi2;
//   
//   real<lower=0> alpha; //percentage of noise data
// }
// 
// model {
//   
//   alpha0 = 0.1
//   
//   mu1  ~ normal(0, 5000);
//   mu2  ~ normal(0, 5000);
//   
//   phi1 ~ normal(0, 5000);
//   phi2  ~ normal(0, 5000);
//   
//   x   ~ alpha0*neg_binomial_2(mu1, phi1) + (1-alpha0)*neg_binomial_2(mu2, phi2);
// }
// 
