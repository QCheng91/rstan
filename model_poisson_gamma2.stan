data {
  int<lower=0> N; // Sample size
  int<lower=0> num_proteins; // Sample size
  int x[N]; // Value for each sample on each dimension
  matrix[N, num_proteins] y; // Value for each sample on each dimension
  int day_index[N]; 
  int<lower=0> d;
}

parameters {
  real<lower=0> mu[d];
  real<lower=0> phi[d];
  vector<lower=0>[num_proteins] a;
}

model {
  a   ~ normal(0, 0.0025);
  mu  ~ normal(0, 5000);
  phi ~ normal(0, 5000);
  
  for(n in 1:N)
    
    x[n]   ~ neg_binomial_2(mu[day_index[n]] + y[n]*a, phi[day_index[n]]);
}

