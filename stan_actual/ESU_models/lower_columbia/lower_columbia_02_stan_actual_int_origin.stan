// 02_stan_actual_upper_columbia_int_origin

data {
  // array[10, 41, 1200] int y; // array with dimensions 10 states (rows), 48 columns (maximum number of site visits), 1200 matrix slices (number of fish); has to be int for multinomial logit
  int max_visits;
  int n_ind; // number of individuals (this is nfish)
  array[n_ind,max_visits] int y; // array with dimensions nfish (number of fish), 41 (maximum number of site visits)
  array[n_ind] int n_obs; // n_obs is a vector that contains the number of site visits for each individual - but needs to be declared as array to take integer values
  vector[29] possible_movements; // a vector containing the number of possible transitions out of each state
  array[n_ind, max_visits-1] int states_mat; // a matrix (array to take integer values) with rows = each fish and columns = the site visits
  // array[54, 2] int movements; // a matrix that contains the possible transitions out of each state
  array[43, 2] int movements; // a matrix that contains the possible transitions out of each state
  // array[758,2] int not_movements; // a matrix containing all non-allowed state transitions
  array[769,2] int not_movements; // a matrix containing all non-allowed state transitions
  int nmovements; // an integer value containing the number of possible transitions (same as rows in movements data)
  // array[n_ind, 48] int dates; // a matrix containing dates (as integers) where rows = number of fish and columns = site visits
  int n_notmovements; // an integer value containing the number of non- possible transitions (same as rows in not_movements data)
  matrix[29, 29] possible_states; // the transition matrix with 1s for allowable transitions and 0s for non-allowable
  array[n_ind, 5] int cat_X_mat; // a matrix with rows = individual fish and columns = various categorical covariates (rear and origin)
}

// The parameters accepted by the model. Our model has only one matrix of parameters
parameters {
  // matrix[10,10] b0_matrix; // parameters we are monitoring - in this case it's the intercept matrix
  // Let's instead only select certain parameters, then put them in a matrix later
real b0_matrix_2_1;
real b0_matrix_1_2;
real b0_matrix_3_2;
real b0_matrix_10_2;
real b0_matrix_12_2;
real b0_matrix_2_3;
real b0_matrix_4_3;
real b0_matrix_8_3;
real b0_matrix_3_4;
real b0_matrix_5_4;
real b0_matrix_4_5;
real b0_matrix_6_5;
real b0_matrix_17_5;
real b0_matrix_5_6;
real b0_matrix_7_6;
real b0_matrix_18_6;
real b0_matrix_6_7;
real b0_matrix_19_7;
real b0_matrix_20_7;
real b0_matrix_28_7;
real b0_matrix_3_8;
real b0_matrix_9_8;
real b0_matrix_8_9;
real b0_matrix_23_9;
real b0_matrix_2_10;
real b0_matrix_2_11;
real b0_matrix_2_12;
real b0_matrix_2_13;
real b0_matrix_2_14;
real b0_matrix_3_15;
real b0_matrix_3_16;
real b0_matrix_5_17;
real b0_matrix_6_18;
real b0_matrix_7_19;
real b0_matrix_7_20;
real b0_matrix_8_21;
real b0_matrix_9_22;
real b0_matrix_9_23;
real b0_matrix_9_24;
real b0_matrix_9_25;
real b0_matrix_9_26;
real b0_matrix_2_27;
real b0_matrix_7_28;

real borigin1_matrix_4_3;
real borigin1_matrix_3_4;
real borigin1_matrix_5_4;
real borigin1_matrix_4_5;
real borigin1_matrix_6_5;
real borigin1_matrix_17_5;
real borigin1_matrix_5_6;
real borigin1_matrix_7_6;
real borigin1_matrix_18_6;
real borigin1_matrix_6_7;
real borigin1_matrix_19_7;
real borigin1_matrix_20_7;
real borigin1_matrix_28_7;
real borigin1_matrix_5_17;
real borigin1_matrix_6_18;
real borigin1_matrix_7_19;
real borigin1_matrix_7_20;
real borigin1_matrix_7_28;

real borigin2_matrix_4_3;
real borigin2_matrix_3_4;
real borigin2_matrix_5_4;
real borigin2_matrix_4_5;
real borigin2_matrix_6_5;
real borigin2_matrix_17_5;
real borigin2_matrix_5_6;
real borigin2_matrix_7_6;
real borigin2_matrix_18_6;
real borigin2_matrix_6_7;
real borigin2_matrix_19_7;
real borigin2_matrix_20_7;
real borigin2_matrix_28_7;
real borigin2_matrix_5_17;
real borigin2_matrix_6_18;
real borigin2_matrix_7_19;
real borigin2_matrix_7_20;
real borigin2_matrix_7_28;

real borigin3_matrix_4_3;
real borigin3_matrix_3_4;
real borigin3_matrix_5_4;
real borigin3_matrix_4_5;
real borigin3_matrix_6_5;
real borigin3_matrix_17_5;
real borigin3_matrix_5_6;
real borigin3_matrix_7_6;
real borigin3_matrix_18_6;
real borigin3_matrix_6_7;
real borigin3_matrix_19_7;
real borigin3_matrix_20_7;
real borigin3_matrix_28_7;
real borigin3_matrix_5_17;
real borigin3_matrix_6_18;
real borigin3_matrix_7_19;
real borigin3_matrix_7_20;
real borigin3_matrix_7_28;
  
}

transformed parameters {
  // Declare a matrix to store b0 params
  matrix[29,29] b0_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  b0_matrix = rep_matrix(-100000, 29, 29);
  
  // Declare a matrix to store borigin1 params
  matrix[29,29] borigin1_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Use a zero instead, since we already have the -100000 in the b0 term. So we just want these to have no effect
  // Non-zero elements will be overwritten
  // borigin1_matrix = rep_matrix(-100000, 29, 29);
  borigin1_matrix = rep_matrix(0, 29, 29);
  
  // Declare a matrix to store borigin2 params
  matrix[29,29] borigin2_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // borigin2_matrix = rep_matrix(-100000, 29, 29);
  borigin2_matrix = rep_matrix(0, 29, 29);
  
  // Declare a matrix to store borigin3 params
  matrix[29,29] borigin3_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // borigin3_matrix = rep_matrix(-100000, 29, 29);
  borigin3_matrix = rep_matrix(0, 29, 29);
  
  
  // Populate this matrix with betas
  // Fill in all of the non-zero elements into the b0_matrix
b0_matrix[2,1] = b0_matrix_2_1;
b0_matrix[1,2] = b0_matrix_1_2;
b0_matrix[3,2] = b0_matrix_3_2;
b0_matrix[10,2] = b0_matrix_10_2;
b0_matrix[12,2] = b0_matrix_12_2;
b0_matrix[2,3] = b0_matrix_2_3;
b0_matrix[4,3] = b0_matrix_4_3;
b0_matrix[8,3] = b0_matrix_8_3;
b0_matrix[3,4] = b0_matrix_3_4;
b0_matrix[5,4] = b0_matrix_5_4;
b0_matrix[4,5] = b0_matrix_4_5;
b0_matrix[6,5] = b0_matrix_6_5;
b0_matrix[17,5] = b0_matrix_17_5;
b0_matrix[5,6] = b0_matrix_5_6;
b0_matrix[7,6] = b0_matrix_7_6;
b0_matrix[18,6] = b0_matrix_18_6;
b0_matrix[6,7] = b0_matrix_6_7;
b0_matrix[19,7] = b0_matrix_19_7;
b0_matrix[20,7] = b0_matrix_20_7;
b0_matrix[28,7] = b0_matrix_28_7;
b0_matrix[3,8] = b0_matrix_3_8;
b0_matrix[9,8] = b0_matrix_9_8;
b0_matrix[8,9] = b0_matrix_8_9;
b0_matrix[23,9] = b0_matrix_23_9;
b0_matrix[2,10] = b0_matrix_2_10;
b0_matrix[2,11] = b0_matrix_2_11;
b0_matrix[2,12] = b0_matrix_2_12;
b0_matrix[2,13] = b0_matrix_2_13;
b0_matrix[2,14] = b0_matrix_2_14;
b0_matrix[3,15] = b0_matrix_3_15;
b0_matrix[3,16] = b0_matrix_3_16;
b0_matrix[5,17] = b0_matrix_5_17;
b0_matrix[6,18] = b0_matrix_6_18;
b0_matrix[7,19] = b0_matrix_7_19;
b0_matrix[7,20] = b0_matrix_7_20;
b0_matrix[8,21] = b0_matrix_8_21;
b0_matrix[9,22] = b0_matrix_9_22;
b0_matrix[9,23] = b0_matrix_9_23;
b0_matrix[9,24] = b0_matrix_9_24;
b0_matrix[9,25] = b0_matrix_9_25;
b0_matrix[9,26] = b0_matrix_9_26;
b0_matrix[2,27] = b0_matrix_2_27;
b0_matrix[7,28] = b0_matrix_7_28;

borigin1_matrix[4,3] = borigin1_matrix_4_3;
borigin1_matrix[3,4] = borigin1_matrix_3_4;
borigin1_matrix[5,4] = borigin1_matrix_5_4;
borigin1_matrix[4,5] = borigin1_matrix_4_5;
borigin1_matrix[6,5] = borigin1_matrix_6_5;
borigin1_matrix[17,5] = borigin1_matrix_17_5;
borigin1_matrix[5,6] = borigin1_matrix_5_6;
borigin1_matrix[7,6] = borigin1_matrix_7_6;
borigin1_matrix[18,6] = borigin1_matrix_18_6;
borigin1_matrix[6,7] = borigin1_matrix_6_7;
borigin1_matrix[19,7] = borigin1_matrix_19_7;
borigin1_matrix[20,7] = borigin1_matrix_20_7;
borigin1_matrix[28,7] = borigin1_matrix_28_7;
borigin1_matrix[5,17] = borigin1_matrix_5_17;
borigin1_matrix[6,18] = borigin1_matrix_6_18;
borigin1_matrix[7,19] = borigin1_matrix_7_19;
borigin1_matrix[7,20] = borigin1_matrix_7_20;
borigin1_matrix[7,28] = borigin1_matrix_7_28;

borigin2_matrix[4,3] = borigin2_matrix_4_3;
borigin2_matrix[3,4] = borigin2_matrix_3_4;
borigin2_matrix[5,4] = borigin2_matrix_5_4;
borigin2_matrix[4,5] = borigin2_matrix_4_5;
borigin2_matrix[6,5] = borigin2_matrix_6_5;
borigin2_matrix[17,5] = borigin2_matrix_17_5;
borigin2_matrix[5,6] = borigin2_matrix_5_6;
borigin2_matrix[7,6] = borigin2_matrix_7_6;
borigin2_matrix[18,6] = borigin2_matrix_18_6;
borigin2_matrix[6,7] = borigin2_matrix_6_7;
borigin2_matrix[19,7] = borigin2_matrix_19_7;
borigin2_matrix[20,7] = borigin2_matrix_20_7;
borigin2_matrix[28,7] = borigin2_matrix_28_7;
borigin2_matrix[5,17] = borigin2_matrix_5_17;
borigin2_matrix[6,18] = borigin2_matrix_6_18;
borigin2_matrix[7,19] = borigin2_matrix_7_19;
borigin2_matrix[7,20] = borigin2_matrix_7_20;
borigin2_matrix[7,28] = borigin2_matrix_7_28;

borigin3_matrix[4,3] = borigin3_matrix_4_3;
borigin3_matrix[3,4] = borigin3_matrix_3_4;
borigin3_matrix[5,4] = borigin3_matrix_5_4;
borigin3_matrix[4,5] = borigin3_matrix_4_5;
borigin3_matrix[6,5] = borigin3_matrix_6_5;
borigin3_matrix[17,5] = borigin3_matrix_17_5;
borigin3_matrix[5,6] = borigin3_matrix_5_6;
borigin3_matrix[7,6] = borigin3_matrix_7_6;
borigin3_matrix[18,6] = borigin3_matrix_18_6;
borigin3_matrix[6,7] = borigin3_matrix_6_7;
borigin3_matrix[19,7] = borigin3_matrix_19_7;
borigin3_matrix[20,7] = borigin3_matrix_20_7;
borigin3_matrix[28,7] = borigin3_matrix_28_7;
borigin3_matrix[5,17] = borigin3_matrix_5_17;
borigin3_matrix[6,18] = borigin3_matrix_6_18;
borigin3_matrix[7,19] = borigin3_matrix_7_19;
borigin3_matrix[7,20] = borigin3_matrix_7_20;
borigin3_matrix[7,28] = borigin3_matrix_7_28;

}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  // Set the priors for each of the non-zero elements of the b0 matrix
b0_matrix_2_1 ~ normal(0,10);
b0_matrix_1_2 ~ normal(0,10);
b0_matrix_3_2 ~ normal(0,10);
b0_matrix_10_2 ~ normal(0,10);
b0_matrix_12_2 ~ normal(0,10);
b0_matrix_2_3 ~ normal(0,10);
b0_matrix_4_3 ~ normal(0,10);
b0_matrix_8_3 ~ normal(0,10);
b0_matrix_3_4 ~ normal(0,10);
b0_matrix_5_4 ~ normal(0,10);
b0_matrix_4_5 ~ normal(0,10);
b0_matrix_6_5 ~ normal(0,10);
b0_matrix_17_5 ~ normal(0,10);
b0_matrix_5_6 ~ normal(0,10);
b0_matrix_7_6 ~ normal(0,10);
b0_matrix_18_6 ~ normal(0,10);
b0_matrix_6_7 ~ normal(0,10);
b0_matrix_19_7 ~ normal(0,10);
b0_matrix_20_7 ~ normal(0,10);
b0_matrix_28_7 ~ normal(0,10);
b0_matrix_3_8 ~ normal(0,10);
b0_matrix_9_8 ~ normal(0,10);
b0_matrix_8_9 ~ normal(0,10);
b0_matrix_23_9 ~ normal(0,10);
b0_matrix_2_10 ~ normal(0,10);
b0_matrix_2_11 ~ normal(0,10);
b0_matrix_2_12 ~ normal(0,10);
b0_matrix_2_13 ~ normal(0,10);
b0_matrix_2_14 ~ normal(0,10);
b0_matrix_3_15 ~ normal(0,10);
b0_matrix_3_16 ~ normal(0,10);
b0_matrix_5_17 ~ normal(0,10);
b0_matrix_6_18 ~ normal(0,10);
b0_matrix_7_19 ~ normal(0,10);
b0_matrix_7_20 ~ normal(0,10);
b0_matrix_8_21 ~ normal(0,10);
b0_matrix_9_22 ~ normal(0,10);
b0_matrix_9_23 ~ normal(0,10);
b0_matrix_9_24 ~ normal(0,10);
b0_matrix_9_25 ~ normal(0,10);
b0_matrix_9_26 ~ normal(0,10);
b0_matrix_2_27 ~ normal(0,10);
b0_matrix_7_28 ~ normal(0,10);

borigin1_matrix_4_3 ~ normal(0,10);
borigin1_matrix_3_4 ~ normal(0,10);
borigin1_matrix_5_4 ~ normal(0,10);
borigin1_matrix_4_5 ~ normal(0,10);
borigin1_matrix_6_5 ~ normal(0,10);
borigin1_matrix_17_5 ~ normal(0,10);
borigin1_matrix_5_6 ~ normal(0,10);
borigin1_matrix_7_6 ~ normal(0,10);
borigin1_matrix_18_6 ~ normal(0,10);
borigin1_matrix_6_7 ~ normal(0,10);
borigin1_matrix_19_7 ~ normal(0,10);
borigin1_matrix_20_7 ~ normal(0,10);
borigin1_matrix_28_7 ~ normal(0,10);
borigin1_matrix_5_17 ~ normal(0,10);
borigin1_matrix_6_18 ~ normal(0,10);
borigin1_matrix_7_19 ~ normal(0,10);
borigin1_matrix_7_20 ~ normal(0,10);
borigin1_matrix_7_28 ~ normal(0,10);

borigin2_matrix_4_3 ~ normal(0,10);
borigin2_matrix_3_4 ~ normal(0,10);
borigin2_matrix_5_4 ~ normal(0,10);
borigin2_matrix_4_5 ~ normal(0,10);
borigin2_matrix_6_5 ~ normal(0,10);
borigin2_matrix_17_5 ~ normal(0,10);
borigin2_matrix_5_6 ~ normal(0,10);
borigin2_matrix_7_6 ~ normal(0,10);
borigin2_matrix_18_6 ~ normal(0,10);
borigin2_matrix_6_7 ~ normal(0,10);
borigin2_matrix_19_7 ~ normal(0,10);
borigin2_matrix_20_7 ~ normal(0,10);
borigin2_matrix_28_7 ~ normal(0,10);
borigin2_matrix_5_17 ~ normal(0,10);
borigin2_matrix_6_18 ~ normal(0,10);
borigin2_matrix_7_19 ~ normal(0,10);
borigin2_matrix_7_20 ~ normal(0,10);
borigin2_matrix_7_28 ~ normal(0,10);

borigin3_matrix_4_3 ~ normal(0,10);
borigin3_matrix_3_4 ~ normal(0,10);
borigin3_matrix_5_4 ~ normal(0,10);
borigin3_matrix_4_5 ~ normal(0,10);
borigin3_matrix_6_5 ~ normal(0,10);
borigin3_matrix_17_5 ~ normal(0,10);
borigin3_matrix_5_6 ~ normal(0,10);
borigin3_matrix_7_6 ~ normal(0,10);
borigin3_matrix_18_6 ~ normal(0,10);
borigin3_matrix_6_7 ~ normal(0,10);
borigin3_matrix_19_7 ~ normal(0,10);
borigin3_matrix_20_7 ~ normal(0,10);
borigin3_matrix_28_7 ~ normal(0,10);
borigin3_matrix_5_17 ~ normal(0,10);
borigin3_matrix_6_18 ~ normal(0,10);
borigin3_matrix_7_19 ~ normal(0,10);
borigin3_matrix_7_20 ~ normal(0,10);
borigin3_matrix_7_28 ~ normal(0,10);
  
  // Loop through the detection matrices for each individual
  for (i in 1:n_ind){
    // Loop through each of the observations, stopping at the loss column
    
    // print("fish # ", i);
    for (j in 1:n_obs[i]){

        // vector for logits
        vector[29] logits;
        
        // derived proportions
        // simplex[29] p_vec;
        // So it looks like you're not allowed to use a simplex as a local variable. That's annoying.
        vector[29] p_vec;
        

        
        // Get the index of the current state
        int current;
        current = states_mat[i,j];
        

        // Populate each of the first 28 (non-loss)
        for (k in 1:28){
          logits[k] = b0_matrix[current, k]+ 
          // cat_X_mat[i,2] * borigin1_matrix[current,k] + # this is rear, which we are currently not using
          cat_X_mat[i,3] * borigin1_matrix[current,k] + 
          cat_X_mat[i,4] * borigin2_matrix[current,k] + 
          cat_X_mat[i,5] * borigin3_matrix[current,k];
        }
        
        // loss param
        logits[29] = 0;
        // proper proportion vector
        p_vec = softmax(logits);
            
        // print("p_vec: ", p_vec);
        // print("y[i,j]: ", y[i,j+1]);
        target += categorical_lpmf(y[i,j+1] | p_vec);
      
    }
    
  }
  

}

