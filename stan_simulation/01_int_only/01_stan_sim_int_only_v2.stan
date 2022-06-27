// 01_stan_sim_int_only_v2

data {
  // array[10, 41, 1200] int y; // array with dimensions 10 states (rows), 48 columns (maximum number of site visits), 1200 matrix slices (number of fish); has to be int for multinomial logit
  array[1200,41] int y;
  int n_ind; // number of individuals (this is 1200)
  array[1200] int n_obs; // n_obs is a vector that contains the number of site visits for each individual - but needs to be declared as array to take integer values
  vector[10] possible_movements; // a vector containing the number of possible transitions out of each state
  array[1200, 40] int states_mat; // a matrix (array to take integer values) with rows = each fish and columns = the site visits
  array[16, 2] int movements; // a matrix that contains the possible transitions out of each state
  array[74,2] int not_movements; // a matrix containing all non-allowed state transitions
  int nmovements; // an integer value containing the number of possible transitions (same as rows in movements data)
  // array[1200, 48] int dates; // a matrix containing dates (as integers) where rows = number of fish and columns = site visits
  int n_notmovements; // an integer value containing the number of non- possible transitions (same as rows in not_movements data)
  matrix[10, 10] possible_states; // the transition matrix with 1s for allowable transitions and 0s for non-allowable
  array[1200, 3] int cat_X_mat; // a matrix with rows = individual fish and columns = various categorical covariates (rear and origin)
}

// The parameters accepted by the model. Our model has only one matrix of parameters
parameters {
  // matrix[10,10] b0_matrix; // parameters we are monitoring - in this case it's the intercept matrix
  // Let's instead only select certain parameters, then put them in a matrix later
  real b0_matrix_2_1;
  real b0_matrix_1_2;
  real b0_matrix_3_2;
  real b0_matrix_6_2;
  real b0_matrix_7_2;
  real b0_matrix_2_3;
  real b0_matrix_4_3;
  real b0_matrix_5_3;
  real b0_matrix_9_3;
  real b0_matrix_3_4;
  real b0_matrix_3_5;
  real b0_matrix_8_5;
  real b0_matrix_2_6;
  real b0_matrix_2_7;
  real b0_matrix_5_8;
  real b0_matrix_3_9;
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  matrix[10,10] b0_matrix;
  // Set all of the elements of the b0 matrix to -999 (the non-zero elements will be overwritten)
  b0_matrix = rep_matrix(-999, 10, 10);
  
  
  // Set the priors for each of the non-zero elements of the b0 matrix
  b0_matrix_2_1 ~ normal(0,10);
  b0_matrix_1_2 ~ normal(0,10);
  b0_matrix_3_2 ~ normal(0,10);
  b0_matrix_6_2 ~ normal(0,10);
  b0_matrix_7_2 ~ normal(0,10);
  b0_matrix_2_3 ~ normal(0,10);
  b0_matrix_4_3 ~ normal(0,10);
  b0_matrix_5_3 ~ normal(0,10);
  b0_matrix_9_3 ~ normal(0,10);
  b0_matrix_3_4 ~ normal(0,10);
  b0_matrix_3_5 ~ normal(0,10);
  b0_matrix_8_5 ~ normal(0,10);
  b0_matrix_2_6 ~ normal(0,10);
  b0_matrix_2_7 ~ normal(0,10);
  b0_matrix_5_8 ~ normal(0,10);
  b0_matrix_3_9 ~ normal(0,10);
  
  // Fill in all of the non-zero elements into the b0_matrix
  b0_matrix[2,1] = b0_matrix_2_1;
  b0_matrix[1,2] = b0_matrix_1_2;
  b0_matrix[3,2] = b0_matrix_3_2;
  b0_matrix[6,2] = b0_matrix_6_2;
  b0_matrix[7,2] = b0_matrix_7_2;
  b0_matrix[2,3] = b0_matrix_2_3;
  b0_matrix[4,3] = b0_matrix_4_3;
  b0_matrix[5,3] = b0_matrix_5_3;
  b0_matrix[9,3] = b0_matrix_9_3;
  b0_matrix[3,4] = b0_matrix_3_4;
  b0_matrix[3,5] = b0_matrix_3_5;
  b0_matrix[8,5] = b0_matrix_8_5;
  b0_matrix[2,6] = b0_matrix_2_6;
  b0_matrix[2,7] = b0_matrix_2_7;
  b0_matrix[5,8] = b0_matrix_5_8;
  b0_matrix[3,9] = b0_matrix_3_9;
  
  // Loop through the detection matrices for each individual
  // for (i in 1:n_ind){
    for (i in 1:1){
    // Loop through each of the observations, stopping at the loss column
    
    // print("fish # ", i);
    // for (j in 1:n_obs[i]){
      for (j in 1:1){
      
      // I think that stan does let you exponentiate a vector, so this will be more compact
      // I think that the categorical_logit distribution works here
      // Create a vector where it's c(alive transitions, loss)
      // First, declare the vector for numerators and total probabilities
      vector[10] num_vec;
      vector[10] p_vec;
      // Assign all of the non-loss probabilities
      // p_vec[1:9] = possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j]]) / 
      // (1 + sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j]]))); 
      // Use a loop here
      for (k in 1:9){
        num_vec[k] = possible_states[states_mat[i,j],k] * exp(b0_matrix[states_mat[i,j],k]);
      }
      // print("num_vec: ", num_vec);
      
      // Calculate the denominator
      real denom;
      denom = 1 + sum(num_vec[1:9]);
      // print("denom: ", denom);
      
      // Populate p_vec
      // for (l in 1:9){
      //   p_vec[l] = num_vec[l]/denom;
      // }
      
      p_vec = num_vec/denom;
      
      // Assign the loss probability
      p_vec[10] = 1 - sum(p_vec[1:9]);
      
      // print("p_vec: ", p_vec);
      // Evaluate the categorical logit
      // From stan reference manual: As of Stan 2.18, the categorical-logit distribution is not 
      // vectorized for parameter arguments, so the loop is required.
      // Well, we are using the multinomial_logit_lpmf function instead of categorical_logit, so I don't think we have to loop it
      // y[,j+1,i] ~ categorical_logit(p_vec);
      // for (z in 1:10){
        // y[z,j+1,i] ~ categorical_logit(p_vec[z]);
        // y[1:10,j+1,i] ~ multinomial_logit_lpmf(10 | p_vec)
        // Do I actually want multinomial_lpmf? Since we're already doing the logit piece by hand?
        // target += multinomial_lpmf(y[1:10,j+1,i] | p_vec);
        // target += multinomial_logit_lpmf(y[1:10,j+1,i] | p_vec);
        // Or is categorical the most appropriate here, since it's a single trial?
        // target += categorical_logit_lpmf(y[1:10,j+1,i] | p_vec);
        print("p_vec: ", p_vec);
        print("y[i,j]: ", y[i,j+1]);
        target += categorical_lpmf(y[i,j+1] | p_vec);
        // y[i,j+1] ~ categorical(p_vec);
      // }
      
      // y[,j+1,i] ~ categorical_logit(c(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]) / 
      // (1 + sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]))),
      // 1 - sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]) / 
      // (1 + sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]))))));
      
    }
    
  }
  

}

