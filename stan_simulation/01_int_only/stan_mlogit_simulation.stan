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
      // for (k in 1:9){
      //   num_vec[k] = possible_states[states_mat[i,j],k] * exp(b0_matrix[states_mat[i,j],k]);
      // }
      num_vec[1] = exp(beta1);
      num_vec[3] = exp(beta3);
      // Calculate the denominator
      real denom;
      denom = 1 + num_vec[1] + num_vec[3];
      // print("denom: ", denom);
      
      // Populate p_vec
      // for (l in 1:9){
      //   p_vec[l] = num_vec[l]/denom;
      // }
      p_vec = num_vec/denom;
      
      // Assign the loss probability
      p_vec[5] = 1 - p_vec[1] - p_vec[3];
  
  // Add in the 0s
  p_vec[2] = 0;
  p_vec[4] = 0;
  
  // print("p_vec, final: ", p_vec);
  
  
}

model {

  for (n in 1:N){
    // print("y[n]: ", y[n]);
    // print("p_vec, final: ", p_vec);
    target += categorical_lpmf(y[n] | p_vec);
  }
  
  
}
