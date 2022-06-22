// 01_stan_sim_int_only

data {
  array[10, 48, 1200] real y; // array with dimensions 10 states (rows), 48 columns (maximum number of site visits), 1200 matrix slices (number of fish)
  int n_ind; // number of individuals (this is 1200)
  int vector[1200] n_obs; // n_obs is a vector that contains the number of site visits for each individual
  vector[10] possible_movements; // a vector containing the number of possible transitions out of each state
  matrix[1200, 47] states_mat; // a matrix with rows = each fish and columns = the site visits
  matrix[16, 2] movements; // a matrix that contains the possible transitions out of each state
  matrix[74,2] not_movements; // a matrix containing all non-allowed state transitions
  int nmovements; // an integer value containing the number of possible transitions (same as rows in movements data)
  matrix[1200, 48] dates; // a matrix containing dates (as integers) where rows = number of fish and columns = site visits
  int n_notmovements; // an integer value containing the number of non- possible transitions (same as rows in not_movements data)
  matrix[10, 10] possible_states; // the transition matrix with 1s for allowable transitions and 0s for non-allowable
  matrix[1200, 3] cat_X_mat; // a matrix with rows = individual fish and columns = various categorical covariates (rear and origin)
}

// The parameters accepted by the model. Our model has only one matrix of parameters
parameters {
  matrix[10,10] b0_matrix; // parameters we are monitoring - in this case it's the intercept matrix
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // Set the priors for each of the non-zero elements of the b0 matrix
  //  for (i in 1:nmovements){
  //    b0_matrix[movements[i,1], movements[i,2]] ~ normal(0,100); // In stan we are now writing it as mean, SD (rather than precision)
  //  }
  
  // Loop through the detection matrices for each individual
  for (i in 1:n_ind){
    // Loop through each of the observations, stopping at the loss column
    for (j in 1:(n_obs[i])){
      
      // I think that stan does let you exponentiate a vector, so this will be more compact
      // I think that the categorical_logit distribution works here
      // Create a vector where it's c(alive transitions, loss)
      y[,j+1,i] ~ categorical_logit(c(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]) / 
      (1 + sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]))),
      1 - sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]) / 
      (1 + sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]))))));
      
    }
    
  }

}

