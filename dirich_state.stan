//State transition Dirichlet distribution

data {

  int<lower=0> n_state; //number of states
  int<lower=0> n_time; //number of time points
  vector[n_state] x_t_in[n_time]; //data at t
  vector[n_state] y[n_state]; //Hamming distance
}

transformed data {
  vector[n_state] x_t[n_time]; //data at t
  for(t in 1:n_time) {
   x_t[t] = x_t_in[t] + 1e-4;
   x_t[t] = x_t[t]/sum(x_t[t]);
  }
}

parameters {
  vector<lower=0>[n_state-1] lambda;
  // vector[n_state] Sigma;
  real<lower = 0> Sigma;
}

model {

  matrix[n_state, n_state] alpha;
  vector[n_state] x_t_1[n_time];
  
  for(i in 1:n_state){
    for(j in 1:n_state){
      
      alpha[i,j] = 0;
    }
  }
  
  for(k in 1:n_state) {
    
    if(k==1){
      
      alpha[k, k] = -lambda[k];
      
    }
    
    else if(k==n_state) {
      
      alpha[k, k-1] = lambda[k-1];
    }
    
    else {
      
      alpha[k, k-1] = lambda[k-1];
      alpha[k, k] = -lambda[k];
    }
    
  }
  //print(alpha);
  
     // prior distribution
  lambda ~ normal(0,10);
  Sigma ~ normal(0, 1000);
  
 // for(t in 1:n_time) 
  
  for(t in 2:n_time) {
    
    
     x_t_1[t] = matrix_exp((2*t) * alpha) * x_t[1];
     // print(x_t_1[t]* Sigma + 1e-5);
     //Sigma
     x_t[t] ~ dirichlet(x_t_1[t]*Sigma + 1e-5);
   
  }
}

generated quantities{
  vector[n_state] x_hat;
  // vector[n_state] x_t_1[n_time];
  real t =1;
  simplex[n_state] x_t_1[n_time];

  matrix[n_state, n_state] alpha;

  for(i in 1:n_state){
    for(j in 1:n_state){

      alpha[i,j] = 0;
    }
  }

  for(k in 1:n_state) {

    if(k==1){

      alpha[k, k] = -lambda[k];

    }

    else if(k==n_state) {

      alpha[k, k-1] = lambda[k-1];
    }

    else {

      alpha[k, k-1] = lambda[k-1];
      alpha[k, k] = -lambda[k];
    }

  }
  //print(alpha);
  x_t_1[1] = x_t[1];
  ind_loop = 1;
  while(t < n_time){
    
    if(t<5.5) {
      x_hat    = matrix_exp((2*t) * alpha) * x_t[1];
      x_t_1[t] = dirichlet_rng(x_hat*Sigma);
      
    }
    
    if(t>5.5 && t< 7.5) {
       x_hat    = matrix_exp((2*(t-4)) * alpha) * x_t_1[5];
       x_t_1[t] = dirichlet_rng(x_hat*Sigma);
      
    }
    
    if(t > 7.5) {
      
       x_hat    = matrix_exp((2*(t-6)) * alpha) * x_t_1[7];
       x_t_1[ind_loop] = dirichlet_rng(x_hat*Sigma);
      
    }
    
    t = t + 0.5;
    ind_loop = ind_loop+1;
      // for(n in 1:n_state) x_t_1[t, n] = normal_rng(x_hat[n], Sigma[n]);
  }

}


// generated quantities{
//   vector[n_state] x_hat;
//   // vector[n_state] x_t_1[n_time];
//   real t =1;
//   simplex[n_state] x_t_1[n_time];
// 
//   matrix[n_state, n_state] alpha;
// 
//   for(i in 1:n_state){
//     for(j in 1:n_state){
// 
//       alpha[i,j] = 0;
//     }
//   }
// 
//   for(k in 1:n_state) {
// 
//     if(k==1){
// 
//       alpha[k, k] = -lambda[k];
// 
//     }
// 
//     else if(k==n_state) {
// 
//       alpha[k, k-1] = lambda[k-1];
//     }
// 
//     else {
// 
//       alpha[k, k-1] = lambda[k-1];
//       alpha[k, k] = -lambda[k];
//     }
// 
//   }
//   //print(alpha);
//   x_t_1[1] = x_t[1];
//   ind_loop = 1;
//   while(t < n_time){
//     
//     if(t<5.5) {
//       x_hat    = matrix_exp((2*t) * alpha) * x_t[1];
//       x_t_1[t] = dirichlet_rng(x_hat*Sigma);
//       
//     }
//     
//     if(t>5.5 && t< 7.5) {
//        x_hat    = matrix_exp((2*(t-4)) * alpha) * x_t_1[5];
//        x_t_1[t] = dirichlet_rng(x_hat*Sigma);
//       
//     }
//     
//     if(t > 7.5) {
//       
//        x_hat    = matrix_exp((2*(t-6)) * alpha) * x_t_1[7];
//        x_t_1[ind_loop] = dirichlet_rng(x_hat*Sigma);
//       
//     }
//     
//     t = t + 0.5;
//     ind_loop = ind_loop+1;
//       // for(n in 1:n_state) x_t_1[t, n] = normal_rng(x_hat[n], Sigma[n]);
//   }
// 
// }

// //State transition Dirichlet distribution
// 
// data {
// 
//   int<lower=0> n_state; //number of proteins
//   int<lower=0> n_time; //number of time points
//   vector[n_state] x_t[n_time]; //data at t
//   vector[n_state] x_t_1[n_time]; //data at t+1
//   vector[n_state] y[n_state]; //Hamming distance
// }
// 
// parameters {
//   vector<lower=0>[n_state-1] lambda; 
//   vector[n_state] Sigma;
// }
// 
// model {
//   
//   matrix[n_state, n_state] alpha;
//   
//   for(i in 1:n_state){
//     for(j in 1:n_state){
//       
//       alpha[i,j] = 0;
//     }
//   }
//   
//   // prior distribution
//   Sigma ~ normal(0, 10);
//   lambda ~ normal(0,10);
//   
//   for(k in 1:n_state) {
//     
//     if(k==1){
//       
//       alpha[k, k] = -lambda[k];
//       
//     }
//     
//     else if(k==n_state) {
//       
//       alpha[k, k-1] = lambda[k-1];
//     }
//     
//     else {
//       
//       alpha[k, k-1] = lambda[k-1];
//       alpha[k, k] = -lambda[k];
//     }
//     
//   }
// 
//   // likelihood
//   for(n in 1:n_time) {
//   
//   // print(alpha');
//   // print(alpha);
//   // print(alpha'*x_t[n]);
//   // print(alpha*x_t[n]);
//   //   print((alpha)'*x_t[n]);
//     
//   x_t_1[n] ~ normal(alpha*(x_t[n]), Sigma);
//   
//   }
// }
