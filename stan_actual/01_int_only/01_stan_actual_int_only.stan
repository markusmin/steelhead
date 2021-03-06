// 01_stan_sim_int_only_v2

data {
  // array[10, 41, 1200] int y; // array with dimensions 10 states (rows), 48 columns (maximum number of site visits), 1200 matrix slices (number of fish); has to be int for multinomial logit
  int max_visits;
  array[1200,max_visits] int y; // array with dimensions 1200 (number of fish), 41 (maximum number of site visits)
  int n_ind; // number of individuals (this is 1200)
  array[1200] int n_obs; // n_obs is a vector that contains the number of site visits for each individual - but needs to be declared as array to take integer values
  vector[10] possible_movements; // a vector containing the number of possible transitions out of each state
  array[1200, max_visits-1] int states_mat; // a matrix (array to take integer values) with rows = each fish and columns = the site visits
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

transformed parameters {
  // Declare a matrix to store b0 params
  matrix[10,10] b0_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  b0_matrix = rep_matrix(-100000, 10, 10);
  
  // Populate this matrix with betas
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
        

}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

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
  
  // Loop through the detection matrices for each individual
  for (i in 1:n_ind){
    // Loop through each of the observations, stopping at the loss column
    
    // print("fish # ", i);
    for (j in 1:n_obs[i]){

        // vector for logits
        vector[10] logits;
        
        // derived proportions
        // simplex[10] p_vec;
        // So it looks like you're not allowed to use a simplex as a local variable. That's annoying.
        vector[10] p_vec;
        

        
        // Get the index of the current state
        int current;
        current = states_mat[i,j];
        

        // Populate each of the first nine (non-loss)
        for (k in 1:9){
          logits[k] = b0_matrix[current, k];
        }
        
        // loss param
        logits[10] = 0;
        // proper proportion vector
        p_vec = softmax(logits);
            
        // print("p_vec: ", p_vec);
        // print("y[i,j]: ", y[i,j+1]);
        target += categorical_lpmf(y[i,j+1] | p_vec);
      
    }
    
  }
  

}

