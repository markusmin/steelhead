// 01_stan_sim_int_only

data {
  array[10, 41, 1200] int y; // array with dimensions 10 states (rows), 48 columns (maximum number of site visits), 1200 matrix slices (number of fish); has to be int for multinomial logit
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
    // Elements of a vectors can only take real values, and for loops can only use integer values.
    // Let's extract elements of the vector as integers to get over this
    // This is not straightforward - solution from here: https://groups.google.com/g/stan-users/c/fX6RHzrzkTU
    // int int_floor(real x) { 
    //   int i; 
    //   i <- 0; 
    //   if (x >= 0) while (x > i + 1) i <- i + 1; 
    //   else while (x < i) i <- i - 1; 
    //   return i; 
    //   } 
    
    // print("fish # ", i);
    for (j in 1:n_obs[i]){
      
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
      
      // Calculate the denominator
      real denom;
      denom = 1 + sum(num_vec[1:9]);
      
      // Populate p_vec
      for (l in 1:9){
        p_vec[l] = num_vec[l]/denom;
      }
      
      // Assign the loss probability
      p_vec[10] = 1 - sum(p_vec[1:9]);
      
      
      // Evaluate the categorical logit
      // From stan reference manual: As of Stan 2.18, the categorical-logit distribution is not 
      // vectorized for parameter arguments, so the loop is required.
      // y[,j+1,i] ~ categorical_logit(p_vec);
      for (z in 1:10){
        // y[z,j+1,i] ~ categorical_logit(p_vec[z]);
        // y[1:10,j+1,i] ~ multinomial_logit_lpmf(10 | p_vec)
        target += multinomial_logit_lpmf(y[1:10,j+1,i] | p_vec);
      }
      
      // y[,j+1,i] ~ categorical_logit(c(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]) / 
      // (1 + sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]))),
      // 1 - sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]) / 
      // (1 + sum(possible_states[states_mat[i,j],1:9] * exp(b0_matrix[states_mat[i,j],]))))));
      
    }
    
  }

}

