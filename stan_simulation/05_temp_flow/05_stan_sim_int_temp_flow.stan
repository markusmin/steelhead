// 05_stan_sim_int_temp_flow

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
  array[1200, max_visits] int dates; // a matrix containing dates (as integers) where rows = number of fish and columns = site visits
  array[9] int temp_index; // get the index of the temperature value of interest for each state
  matrix[1095,8] temp_matrix; // matrix containing the temperature value on each of the 1095 possible dates + 8 different temp data sources
  array[9] int flow_index; // get the index of the flow value of interest for each state
  matrix[1095,8] flow_matrix; // matrix containing the flow value on each of the 1095 possible dates + 8 different temp data sources
  int n_notmovements; // an integer value containing the number of non- possible transitions (same as rows in not_movements data)
  array[10, 10] int possible_states; // the transition matrix with 1s for allowable transitions and 0s for non-allowable
  array[1200, 4] int cat_X_mat; // a matrix with rows = individual fish and columns = various categorical covariates (rear and origin)
}

transformed data{
  // Using the possible_states matrix, we will figure out for each state what the possible transitions are
  
  array[9,9] int possible_transitions;
  array[9,9] int impossible_transitions;
        
  
  for (i in 1:9){
    int counter1;
    counter1 = 1;
  
    int counter2;
    counter2 = 1;
  
    int index;
    index = 1;
    for (j in 1:9){
          if (possible_states[i,j] == 1){
            possible_transitions[i, counter1] = index;
            counter1 = counter1 + 1;
          }
          else{
            impossible_transitions[i, counter2] = index;
            counter2 = counter2 + 1;
          }
          index = index + 1;
      }
    
  }
  // print("possible_transitions: ", possible_transitions);
  // print("impossible_transitions: ", impossible_transitions);
}

// The parameters accepted by the model. Our model has only one matrix of parameters
parameters {
  // intercept parameters
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
  
  // temp parameters
  real btemp_matrix_2_1;
  real btemp_matrix_1_2;
  real btemp_matrix_3_2;
  real btemp_matrix_6_2;
  real btemp_matrix_7_2;
  real btemp_matrix_2_3;
  real btemp_matrix_4_3;
  real btemp_matrix_5_3;
  real btemp_matrix_9_3;
  real btemp_matrix_3_4;
  real btemp_matrix_3_5;
  real btemp_matrix_8_5;
  real btemp_matrix_2_6;
  real btemp_matrix_2_7;
  real btemp_matrix_5_8;
  real btemp_matrix_3_9;
  
  // flow parameters
  real bflow_matrix_2_1;
  real bflow_matrix_1_2;
  real bflow_matrix_3_2;
  real bflow_matrix_6_2;
  real bflow_matrix_7_2;
  real bflow_matrix_2_3;
  real bflow_matrix_4_3;
  real bflow_matrix_5_3;
  real bflow_matrix_9_3;
  real bflow_matrix_3_4;
  real bflow_matrix_3_5;
  real bflow_matrix_8_5;
  real bflow_matrix_2_6;
  real bflow_matrix_2_7;
  real bflow_matrix_5_8;
  real bflow_matrix_3_9;
  
  
}

transformed parameters {
  // Declare a matrix to store b0 params
  matrix[10,10] b0_matrix;
  
  // Declare a matrix to store btemp params
  matrix[10,10] btemp_matrix;
  
  // Declare a matrix to store bflow params
  matrix[10,10] bflow_matrix;
  
  // Set all of the elements of each parameter matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  b0_matrix = rep_matrix(-100000, 10, 10);
  btemp_matrix = rep_matrix(-100000, 10, 10);
  bflow_matrix = rep_matrix(-100000, 10, 10);
  
  // Fill in all of the non-zero elements in each matrix
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
  
  btemp_matrix[2,1] = btemp_matrix_2_1;
  btemp_matrix[1,2] = btemp_matrix_1_2;
  btemp_matrix[3,2] = btemp_matrix_3_2;
  btemp_matrix[6,2] = btemp_matrix_6_2;
  btemp_matrix[7,2] = btemp_matrix_7_2;
  btemp_matrix[2,3] = btemp_matrix_2_3;
  btemp_matrix[4,3] = btemp_matrix_4_3;
  btemp_matrix[5,3] = btemp_matrix_5_3;
  btemp_matrix[9,3] = btemp_matrix_9_3;
  btemp_matrix[3,4] = btemp_matrix_3_4;
  btemp_matrix[3,5] = btemp_matrix_3_5;
  btemp_matrix[8,5] = btemp_matrix_8_5;
  btemp_matrix[2,6] = btemp_matrix_2_6;
  btemp_matrix[2,7] = btemp_matrix_2_7;
  btemp_matrix[5,8] = btemp_matrix_5_8;
  btemp_matrix[3,9] = btemp_matrix_3_9;
    
  bflow_matrix[2,1] = bflow_matrix_2_1;
  bflow_matrix[1,2] = bflow_matrix_1_2;
  bflow_matrix[3,2] = bflow_matrix_3_2;
  bflow_matrix[6,2] = bflow_matrix_6_2;
  bflow_matrix[7,2] = bflow_matrix_7_2;
  bflow_matrix[2,3] = bflow_matrix_2_3;
  bflow_matrix[4,3] = bflow_matrix_4_3;
  bflow_matrix[5,3] = bflow_matrix_5_3;
  bflow_matrix[9,3] = bflow_matrix_9_3;
  bflow_matrix[3,4] = bflow_matrix_3_4;
  bflow_matrix[3,5] = bflow_matrix_3_5;
  bflow_matrix[8,5] = bflow_matrix_8_5;
  bflow_matrix[2,6] = bflow_matrix_2_6;
  bflow_matrix[2,7] = bflow_matrix_2_7;
  bflow_matrix[5,8] = bflow_matrix_5_8;
  bflow_matrix[3,9] = bflow_matrix_3_9;
        

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
  
  btemp_matrix_2_1 ~ normal(0,10);
  btemp_matrix_1_2 ~ normal(0,10);
  btemp_matrix_3_2 ~ normal(0,10);
  btemp_matrix_6_2 ~ normal(0,10);
  btemp_matrix_7_2 ~ normal(0,10);
  btemp_matrix_2_3 ~ normal(0,10);
  btemp_matrix_4_3 ~ normal(0,10);
  btemp_matrix_5_3 ~ normal(0,10);
  btemp_matrix_9_3 ~ normal(0,10);
  btemp_matrix_3_4 ~ normal(0,10);
  btemp_matrix_3_5 ~ normal(0,10);
  btemp_matrix_8_5 ~ normal(0,10);
  btemp_matrix_2_6 ~ normal(0,10);
  btemp_matrix_2_7 ~ normal(0,10);
  btemp_matrix_5_8 ~ normal(0,10);
  btemp_matrix_3_9 ~ normal(0,10);
  
  bflow_matrix_2_1 ~ normal(0,10);
  bflow_matrix_1_2 ~ normal(0,10);
  bflow_matrix_3_2 ~ normal(0,10);
  bflow_matrix_6_2 ~ normal(0,10);
  bflow_matrix_7_2 ~ normal(0,10);
  bflow_matrix_2_3 ~ normal(0,10);
  bflow_matrix_4_3 ~ normal(0,10);
  bflow_matrix_5_3 ~ normal(0,10);
  bflow_matrix_9_3 ~ normal(0,10);
  bflow_matrix_3_4 ~ normal(0,10);
  bflow_matrix_3_5 ~ normal(0,10);
  bflow_matrix_8_5 ~ normal(0,10);
  bflow_matrix_2_6 ~ normal(0,10);
  bflow_matrix_2_7 ~ normal(0,10);
  bflow_matrix_5_8 ~ normal(0,10);
  bflow_matrix_3_9 ~ normal(0,10);
  
  // Loop through the detection matrices for each individual
  for (i in 1:n_ind){
    // Loop through each of the observations, stopping at the loss column
    
    // print("fish # ", i);
    for (j in 1:n_obs[i]){

        // vector for logits
        // vector[10] logits;
        vector[10] logits;
        
        // derived proportions
        // simplex[10] p_vec;
        // So it looks like you're not allowed to use a simplex as a local variable, so we'll use a vector instead
        vector[10] p_vec;
        
        // Get the index of the current state
        int current;
        current = states_mat[i,j];
        
        // Get the indices of the possible transitions
        int n_possible_transitions;
        n_possible_transitions = sum(possible_states[current,1:9]);
        // vector[n_possible_transitions] possible_transitions;
        // 
        // int counter;
        // counter = 1;
        // int index;
        // index = 1;
        // print("current: ", current);
        // print("n_possible_transitions: ", n_possible_transitions);
        for (l in 1:n_possible_transitions){
          logits[possible_transitions[current,l]] = b0_matrix[current, possible_transitions[current,l]] + 
          btemp_matrix[current,possible_transitions[current,l]] * // index to btemp value, multiply by the correct temp value
          temp_matrix[dates[i,j],temp_index[current]] + 
          bflow_matrix[current,possible_transitions[current,l]] * // index to bflow value, multiply by the correct flow value
          flow_matrix[dates[i,j],flow_index[current]];
        }
        
        // Populate each of the first nine (non-loss)
        // for (k in 1:9){ // Change this line so it instead loops through only the possible transitions
        // for (k in possible_transitions){
        // // for (k in 1:n_possible_transitions){
        //   // int kk;
        //   // kk = possible_transitions[k];
        //   logits[k] = b0_matrix[current, k] + cat_X_mat[i,2] * borigin1_matrix[current,k] + cat_X_mat[i,3] * borigin2_matrix[current,k];
        // }
        
        // int k = 1;
        // while (k <= n_possible_transitions){
        //   int kk = possible_transitions[k];
        //   logits[kk] = b0_matrix[current, kk] + cat_X_mat[i,2] * borigin1_matrix[current,kk] + cat_X_mat[i,3] * borigin2_matrix[current,kk];
        //   k = k + 1;
        // }
        
        // Get the indices of the impossible transitions
        int n_impossible_transitions;
        n_impossible_transitions = 9 - n_possible_transitions;
        
        // int m = 1;
        // while (m <= n_impossible_transitions){
        //   int mm = impossible_transitions[m];
        //   logits[mm] = -100000;
        // }
        for (m in 1:n_impossible_transitions){
          logits[impossible_transitions[current, m]] = -100000;
        }
        
        // print("possible transitions: ", possible_transitions);
        // print("impossible transitions: ", impossible_transitions);
        
        
        
        // Make a second loop for all of the impossible transitions
        
        // loss param
        logits[10] = 0;
                  // print("logits: ", logits);
        // proper proportion vector
        p_vec = softmax(logits);
                  // print("p_vec: ", p_vec);
            
        // print("p_vec: ", p_vec);
        // print("y[i,j]: ", y[i,j+1]);
        target += categorical_lpmf(y[i,j+1] | p_vec);
      
    }
    
  }
  

}

