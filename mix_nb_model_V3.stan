//mixture model for noise+data+0spikes

data {
  int<lower=0> N; // Sample size
  int x[N]; // Value for each sample on each dimension
  int<lower=0> d; // Day Index
  int day_index[N];
}

parameters {
  
  real<lower =0, upper = 40>  mu_n;
  real<lower = 2.0*mu_n> mu[d];
  
  real<lower=0> phi_n;
  real<lower=0.01> phi[d];
  
  real<lower=0, upper =1>   theta[d];
  real<lower=0, upper =1>   theta_0[d];
  
  real<lower=0> mu0;
  real<lower=0> sigma;
  
  // real<lower=0, upper=1> alpha[d];
  // real<lower=0> mun0;
  // real<lower=0, upper = 10> sigma0;
  
}

model {
  
  real log_noise;
  //real log_noise; 
  
  mu_n    ~ normal(0, 40);
  phi_n   ~ normal(0, 1000);
  
  // mun0     ~ normal(0, 100);
  // sigma0   ~ normal(0, 10);
  // 
  // mu_n ~ normal(0, 20);
  
  // for(t in 1:d) {
  // 
  // mu_h[t,1]   ~ normal(mun0, sigma0);
  // mu_h[t,2]   ~ normal(mu0, sigma);
  // 
  // }
  
  //phi_n   ~ normal(0, 1000);
  // mu_h    ~ normal(0, 1000);
  // sigma_h ~ normal(0, 100);
  // mu      ~ normal(mu_h, sigma_h);
  // alpha ~ beta(1,10);

  // for(t in 1:d) {
  // 
  // target += log_mix(alpha[t], 
  // beta_lpdf(theta[t] | 2,2),
  // beta_lpdf(theta[t] | 1,1000));
  // 
  // }
   
  theta ~ beta(1, 3);
  
  phi      ~ normal(0, 1000);
  
  mu0      ~ normal(250, 200);
  sigma    ~ normal(0, 10);
  
  mu ~ normal(mu0, sigma);
  
  for(n in 1:N) {
    
    if(x[n] == 0) log_noise = log_sum_exp(bernoulli_lpmf(1 | theta_0[day_index[n]]), bernoulli_lpmf(0 | theta_0[day_index[n]]) + neg_binomial_2_lpmf(0 | mu_n, phi_n));
    
    else log_noise = bernoulli_lpmf(0 | theta_0[day_index[n]]) + neg_binomial_2_lpmf(x[n] | mu_n, phi_n);
    
    target += log_mix(theta[day_index[n]], neg_binomial_2_lpmf(x[n] | mu[day_index[n]], phi[day_index[n]]), log_noise);
    
    // if(x[n] == 0) signal[n] = log_sum_exp(bernoulli_lpmf(1 | theta_0[day_index[n]]), bernoulli_lpmf(0 | theta_0[day_index[n]]) + neg_binomial_2_lpmf(x[n] | mu_n, phi_n));
    // 
    // else signal[n] = bernoulli_lpmf(0 | theta_0[day_index[n]]) + neg_binomial_2_lpmf(x[n] | mu_n, phi_n);
    // 
    // target += log_mix(theta[day_index[n]], signal[n], neg_binomial_2_lpmf(x[n] | mu[day_index[n]], phi[day_index[n]]));

  }
  
  // 
  //   for(n in 1:N) {
  //   
  //   signal[n] = log_mix(theta[day_index[n]], neg_binomial_2_lpmf(x[n] | mu_n, phi_n), neg_binomial_2_lpmf(x[n] | mu[day_index[n]], phi[day_index[n]]));
  //   
  //   if(x[n] == 0) 
  //   
  //   target += log_sum_exp(bernoulli_lpmf(1 | theta_0), bernoulli_lpmf(0 | theta_0) + signal[n]);
  //   
  //   else 
  //   
  //   target += bernoulli_lpmf(0 | theta_0) + signal[n];
  //   
  // }
    
  // log_mix(theta[day_index[n]], neg_binomial_2_lpmf(x[n] | mu_n, phi_n), neg_binomial_2_lpmf(x[n] | mu[day_index[n]], phi[day_index[n]]));
    //target += log_mix(theta, gamma_lpdf(x[n]| a[1], b[1]), gamma_lpdf(x[n]| a[2], b[2]));
}


//log_sum_exp(bernoulli_lpmf(1 | theta), bernoulli_lpmf(0 | theta) + poisson_lpmf(y[n] | lambda));


//mixture model three distributions

// data {
//   int<lower=0> N; // Sample size
//   int x[N]; // Value for each sample on each dimension
//   int<lower=0> d; // Day Index
//   int day_index[N];
//   int n_groups;
// }
// 
// parameters {
//   positive_ordered[n_groups] mu;
//   real<lower=0>  phi[n_groups];
// 
//   //matrix<lower=0, upper=1>[d, n_groups] theta;
//   simplex[n_groups] theta[d];
// }
// 
// model {
//   
//   vector[n_groups] contributions;
//   vector[N] result[n_groups];
//   //matrix[n_groups, N] result;
//   
//   //priors
//   mu[1]   ~   normal(0, 10);
//   mu[2]  ~   normal(10,5000);
//   mu[3]  ~   normal(100,5000);
//   phi  ~   normal(0, 5000);
//   
//   for(m in 1:n_groups) {
//     
//   theta[m] ~ dirichlet(rep_vector(2., n_groups));
//   
//   }
//  
//   // // likelihood
//   
//  // result[1] = neg_binomial_2_lpmf(x | mu[1], phi[1]);
//   // result[2] = neg_binomial_2_lpmf(x | mu[2], phi[2]);
//   // result[3] = neg_binomial_2_lpmf(x | mu[3], phi[3]);
//   // result = neg_binomial_2_lpmf(x | mu[1], phi[1]);
//    // result[1][N] = neg_binomial_2_lpmf(x | mu[1], phi[1]);
//    //  result[2][N] = neg_binomial_2_lpmf(x | mu[2], phi[2]);
//    //  result[3][N] = neg_binomial_2_lpmf(x | mu[3], phi[3]);
// 
//   for(n in 1:N){
//     
//     for(k in 1:n_groups) {
//       
//       contributions[k] = log(theta[day_index[n]][k]) + neg_binomial_2_lpmf(x[n] | mu[k], phi[k]);
//     }
//     
//     target += log_sum_exp(contributions);
//     
//   }
//   // 
//   // for(n in 1:N) {
//   //   
//   //   for(k in 1:n_groups) {
//   //     
//   //     contributions[k] = log(theta[day_index[n]][k]) + neg_binomial_2_lpmf(x[n] | mu[k], phi[k]);
//   //   }
//   //   target += log_sum_exp(contributions);
//   //   
//   // }
//   
//     //log_mix(theta[day_index[n]], neg_binomial_2_lpmf(x[n] | mu_n, phi_n), neg_binomial_2_lpmf(x[n] | mu[day_index[n]], phi[day_index[n]]));
//     //target += log_mix(theta, gamma_lpdf(x[n]| a[1], b[1]), gamma_lpdf(x[n]| a[2], b[2]));
// }

//mixture model for noise+data

// data {
//   int<lower=0> N; // Sample size
//   int x[N]; // Value for each sample on each dimension
//   int<lower=0> d; // Day Index
//   int day_index[N];
// }
// 
// parameters {
//   real<lower=0, upper=20>  mu_n;
//   real<lower=5> mu[d];
//   real<lower=0> phi_n;
//   real<lower=0.01> phi[d];
//   real<lower=0, upper =1>   theta[d];
//   //
//   // positive_ordered[2]  phi;
//   // vector<lower=0.>[2]  mu;
// }
// 
// model {
//   mu_n    ~ normal(0, 10);
//   phi_n   ~ normal(0, 5000);
//   mu      ~ normal(0, 5000);
//   phi     ~ normal(0, 5000);
//   theta   ~ beta(2,2);
// 
//   for(n in 1:N)
//     target += log_mix(theta[day_index[n]], neg_binomial_2_lpmf(x[n] | mu_n, phi_n), neg_binomial_2_lpmf(x[n] | mu[day_index[n]], phi[day_index[n]]));
//     //target += log_mix(theta, gamma_lpdf(x[n]| a[1], b[1]), gamma_lpdf(x[n]| a[2], b[2]));
// }

///////////////////////////////
// data {
//   int<lower=0> N; // Sample size
//   int x[N]; // Value for each sample on each dimension
// }
// 
// parameters {
//   positive_ordered[2]  mu;
//   vector<lower=0.>[2]  phi;
//   real<lower=0, upper =1>   theta;
//   // 
//   // positive_ordered[2]  phi;
//   // vector<lower=0.>[2]  mu;
// }
// 
// model {
//   // mu[1]   ~ normal(0, 10);
//   // mu[2]   ~ normal(20,5000);
//   mu      ~ normal(0,5000);
//   phi     ~ normal(0, 5000);
//   theta   ~ beta(2,2);
//  
//   for(n in 1:N)
//     target += log_mix(theta, neg_binomial_2_lpmf(x[n] | mu[1], phi[1]), neg_binomial_2_lpmf(x[n] | mu[2], phi[2]));
//     target += log_mix(theta, gamma_lpdf(x[n]| a[1], b[1]), gamma_lpdf(x[n]| a[2], b[2]));
// }
