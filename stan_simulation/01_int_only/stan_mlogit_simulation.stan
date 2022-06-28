// stan_mlogit_simulation

data {
  array[100] int y;
  int N;
}


parameters {
  real beta1;
  real beta3;
}

transformed parameters {
    // Populate p vector
      vector[5] num_vec;
      // vector[5] p_vec;
      simplex[5] p_vec;

      num_vec[1] = exp(beta1);
      num_vec[3] = exp(beta3);
      // Calculate the denominator
      real denom;
      denom = 1 + num_vec[1] + num_vec[3];

      p_vec = num_vec/denom;
      
      // Assign the loss probability
      p_vec[5] = 1 - p_vec[1] - p_vec[3];
  
  // Add in the 0s
  // Would floating point zero help?
  p_vec[2] = 0.0;
  p_vec[4] = 0.0;
  
  // print("p_vec, final: ", p_vec);
  
  
}

model {

  for (n in 1:N){
    print("y[n]: ", y[n]);
    print("p_vec, final: ", p_vec);
    target += categorical_lpmf(y[n] | p_vec);
    print("target: ", target());
  }
  
  
}
