// stan_mlogit_simulation 2 - using categorical_logit_lpmf instead of categorical_lpmf

data {
  int N;
  int K;
  int D;
  array[N] int y;
  matrix[N, D] x;
}


parameters {
  // matrix[D, K] beta;
  real beta1;
  real beta3;
}

transformed parameters{
  matrix[D, K] beta = rep_matrix(0, D, K);
  beta[1,1] = beta1;
  beta[1,3] = beta3;
}

model {
  matrix[N, K] x_beta = x * beta;
  
  // Fix 2nd and 4th spots to zero
  for (i in 1:N){
    x_beta[i,2] = 0;
    x_beta[i,4] = 0;
  }

  
  for (n in 1:N){
    y[n] ~ categorical_logit(x_beta[n]');
  }
  
}
