// 03_parallel_snake_stan_actual_int_origin_deteff
// This model includes a detection efficiency calculation for the tributaries.

functions{
    real partial_sum_lpmf( // int[] slice_n_fish, # I don't think that we need this, given that we are using start and end to index the detection histories
                        array[,] int slice_y,
                        int start, int end,
                        // I think we need to re-declare things for our function that are also already declared in our data?
                        int n_ind,
                        int max_visits,
                        // array[n_ind, 7] int cat_X_mat,
                        array[,] int cat_X_mat,
                        // array[n_ind, max_visits-1] int states_mat,
                        array[,] int states_mat,
                        // array[n_ind,max_visits] int y, # I don't think that we need this for the same reason as above
                        // array[,] int y,
                        array[] int n_obs,
                        // I think we also need to declare the dimensions for our parameters
                        // matrix[29,29] b0_matrix,
                        // matrix[29,29] borigin1_matrix, 
                        // matrix[29,29] borigin2_matrix,
                        // matrix[29,29] borigin3_matrix,
                        // matrix[29,29] borigin4_matrix,
                        // matrix[29,29] borigin5_matrix
                        // matrix[,] b0_matrix,
                        // matrix[,] borigin1_matrix,
                        // matrix[,] borigin2_matrix,
                        // matrix[,] borigin3_matrix,
                        // matrix[,] borigin4_matrix,
                        // matrix[,] borigin5_matrix
                        
                        // For detection efficiency calculation purposes: We need two of each of these one for transitions
                        // where we can calculate detection efficiency, one for where we can't. Most parameters
                        // will be shared between the two, but transitions into tributaries won't.
                        // DE = detection efficiency, NDE = no detection efficiency
                        array[,] real b0_matrix_DE,
                        array[,] real borigin1_matrix_DE,
                        array[,] real borigin2_matrix_DE,
                        array[,] real borigin3_matrix_DE,
                        array[,] real borigin4_matrix_DE,
                        array[,] real borigin5_matrix_DE,
                        
                        array[,] real b0_matrix_NDE,
                        array[,] real borigin1_matrix_NDE,
                        array[,] real borigin2_matrix_NDE,
                        array[,] real borigin3_matrix_NDE,
                        array[,] real borigin4_matrix_NDE,
                        array[,] real borigin5_matrix_NDE,
                        // below is new data for detection efficiency
                        // declare an array that you can detection probabilities in, 
                        // indexed by rows = run years, columns = tributaries.
                        // We will calculate estimated detection probability
                        // This matrix will contain the parameters (both alpha/intercept from era
                        // and beta, for slope of discharge relationship)
                        // array[,] real det_eff_param_matrix;
                        // array[,] real discharge_matrix;
                        
                        // different approach - include it directly in the model
                        // need to declare an array with dimensions (number of tributaries for which we can estimate detection efficiency) x 
                        // (number of run years) x (number of detection efficiency parameters)
                        // this is in effect an array of design matrices
                        array[,,] real tributary_design_matrices_array,
                        
                        // declare a vector for the parameters for the detection efficiency glm
                        vector detection_efficiency_parameters,
                        
                        // declare a vector to store the run years of fish
                        vector fish_run_years,
                        
                        
                        
                        
                        
                        ) { # I don't think we need this either? Since we're just indexing it again with start and end
                          
  // First, declare the total lp (log probability) variable that we will be returning
  real total_lp = 0;
  // Now, loop through the detection matrices for each individual IN THIS SLICE OF THE DATA
  // for (i in 1:slice_n_fish){
    
    // start - end apparently doesn't work because reduce_sum() resets the indices for each slice (confusing) - according to this post: https://discourse.mc-stan.org/t/parallelization-in-cmdstanr-using-reduce-sum-on-a-multinomial-likelihood/24607/7
  // or maybe not, this post seems to contradict this: https://discourse.mc-stan.org/t/help-with-multi-threading-a-simple-ordinal-probit-model-using-reduce-sum/15353
  for (i in start:end){
    // for (i in 1:end-start+1){
    // Loop through each of the observations, stopping at the loss column
    
    // print("fish # ", i);
    
    // Here, create a vector to store the lp at each observation for a fish
    // vector[n_obs[i]] lp_fish;
    // Let's initialize this instead as a real value starting at zero
    real lp_fish = 0;
    for (j in 1:n_obs[i]){
      // for (j in 1:n_obs[i - start + 1]){

        // vector for logits
        vector[43] logits;
        
        // Get the index of the current state
        int current;
        current = states_mat[i,j];
        // current = states_mat[i - start + 1,j];
        

        // Populate each of the first 42 (non-loss)
        for (k in 1:42){
          logits[k] = b0_matrix[current, k]+ 
          // cat_X_mat[i,2] * borigin1_matrix[current,k] + # this is rear, which we are currently not using
          cat_X_mat[i,3] * borigin1_matrix[current,k] +
          cat_X_mat[i,4] * borigin2_matrix[current,k] +
          cat_X_mat[i,5] * borigin3_matrix[current,k] +
          cat_X_mat[i,6] * borigin4_matrix[current,k] +
          cat_X_mat[i,7] * borigin5_matrix[current,k];
          // cat_X_mat[i - start + 1,3] * borigin1_matrix[current,k] + 
          // cat_X_mat[i - start + 1,4] * borigin2_matrix[current,k] + 
          // cat_X_mat[i - start + 1,5] * borigin3_matrix[current,k] + 
          // cat_X_mat[i - start + 1,6] * borigin4_matrix[current,k] + 
          // cat_X_mat[i - start + 1,7] * borigin5_matrix[current,k];
        }
        
        // loss param
        logits[43] = 0;
        
        // Now, we need an if/else statement: If it's in a tributary, correct for detection efficiency, otherwise leave it alone
        
        // if it is in a state that connects to at tributary for which we can estimate detection efficiency:
        if (current %in% mainstem_trib_states){
          // declare vector to store movement probabilities
          vector[43] p_vec;
          
          // declare vector to store uncorrected for detection efficiency movement probabilities
          vector[43] p_vec_uncorrected;
          
          
        // proper proportion vector, uncorrected for detection efficiency
        p_vec_uncorrected = softmax(logits);
        
        // calculate detection efficiency (pseudocode)
        
        // first: we need to determine how many different detection probabilities we need to estimate, and the loop through them
        // declare an integer for the length of the for loop
        int n_deteffs;
        n_deteffs = n_detection_efficiencies[current];
        
        // now declare a vector with that length to store the detection probability values
        vector[n_deteffs] detection_efficiencies;
        
        // loop through all of the detection efficiencies that you need to calculate
        for (l in 1:n_deteffs){
          // step 1: select the appropriate tributary design matrix
          // declare an integer to store the index of the state
          int trib_index;
          // figure out the state of the index
          trib_index = mainstem_trib_connections[[which(current == mainstem_trib_states)]][l];
        
          // pull the appropriate design matrix out
          // tributary_design_matrices_array[trib_index,,];
        
          // step 2: select the correct run year from the tributary design matrix
          int fish_run_year_index;
          fish_run_year_index = fish_run_years[i];
        
          // step 3: Multiply this row of the tributary design matrix by the parameter vector to get eta, the linear predictor
          real eta;
          eta = tributary_design_matrices_array[trib_index,fish_run_year_index,] * det_eff_params;
          
        
          // step 4: get the estimated detection probability using the inverse logit
          detection_efficiencies[l] = exp(eta)/(1 + exp(eta));
        }
        
        // now, store all of the detection efficiencies in the appropriate places in all of the transition probabilities
        
        vector[42] correction_factors;
        
        
        
        
        
        // now incorporate detection efficiency (pseudocode)
        for (m in 1:(n_states-1)){
          p_vec[m] = p_vec_uncorrected[m] * detection_efficiencies[m];
          
        }
        
          // if it's not, just calculate as normal
        } else{
          
          // declare vector to store movement probabilities
          vector[43] p_vec;
          
          // populate that vector from logits
          p_vec = softmax(logits);
          
          
        }
        
        

        
            
        // print("p_vec: ", p_vec);
        // print("i = ",i);
        // print("j = ",j);
        // print("slice_y[i,j]: ", slice_y[i,j+1]);
        // Store the log probability of that individual transition in the vector that we declared
        // lp_fish[j] = categorical_lpmf(slice_y[i,j+1] | p_vec);
        // Changed the indexing here, based on:https://discourse.mc-stan.org/t/help-with-multi-threading-a-simple-ordinal-probit-model-using-reduce-sum/15353
        lp_fish += categorical_lpmf(slice_y[i - start + 1,j+1] | p_vec);
        // lp_fish += categorical_lpmf(slice_y[i,j+1] | p_vec);
      
    }
    
    // Now, increment the log density
    // total_lp += sum(lp_fish);
    total_lp += lp_fish;
    
    
  }
    return total_lp;
}

}



data {
  // array[10, 41, 1200] int y; // array with dimensions 10 states (rows), 48 columns (maximum number of site visits), 1200 matrix slices (number of fish); has to be int for multinomial logit
  int max_visits;
  int n_ind; // number of individuals (this is nfish)
  array[n_ind,max_visits] int y; // array with dimensions nfish (number of fish), 41 (maximum number of site visits)
  array[n_ind] int n_obs; // n_obs is a vector that contains the number of site visits for each individual - but needs to be declared as array to take integer values
  vector[43] possible_movements; // a vector containing the number of possible transitions out of each state
  array[n_ind, max_visits-1] int states_mat; // a matrix (array to take integer values) with rows = each fish and columns = the site visits
  // array[54, 2] int movements; // a matrix that contains the possible transitions out of each state
  int nmovements; // an integer value containing the number of possible transitions (same as rows in movements data)
  array[nmovements, 2] int movements; // a matrix that contains the possible transitions out of each state
  // array[758,2] int not_movements; // a matrix containing all non-allowed state transitions
  int n_notmovements; // an integer value containing the number of non- possible transitions (same as rows in not_movements data)
  array[n_notmovements,2] int not_movements; // a matrix containing all non-allowed state transitions
  // array[n_ind, 48] int dates; // a matrix containing dates (as integers) where rows = number of fish and columns = site visits
  matrix[43, 43] possible_states; // the transition matrix with 1s for allowable transitions and 0s for non-allowable
  array[n_ind, 8] int cat_X_mat; // a matrix with rows = individual fish and columns = various categorical covariates (rear and origin and run year)
  
  // new data for tributary detection efficiency
  array[7] vector[5] mainstem_trib_connections;
  
  // Declare data for parallelization (reduce sum function)
  int<lower=0> N;
  // array[N,max_visits] int slice_y; # doesn't need to be declared as data!
  // array[N,max_visits] int n_fish;
  int<lower=1> grainsize;
}

// The parameters accepted by the model. Our model has only one matrix of parameters
parameters {
  // matrix[10,10] b0_matrix; // parameters we are monitoring - in this case it's the intercept matrix
  // Let's instead only select certain parameters, then put them in a matrix later
real b0_matrix_2_1;
real b0_matrix_1_2;
real b0_matrix_3_2;
real b0_matrix_10_2;
real b0_matrix_11_2;
real b0_matrix_12_2;
real b0_matrix_13_2;
real b0_matrix_14_2;
real b0_matrix_16_2;
real b0_matrix_17_2;
real b0_matrix_18_2;
real b0_matrix_19_2;
real b0_matrix_41_2;
real b0_matrix_2_3;
real b0_matrix_8_3;
real b0_matrix_20_3;
real b0_matrix_21_3;
real b0_matrix_22_3;
real b0_matrix_23_3;
real b0_matrix_3_4;
real b0_matrix_5_4;
real b0_matrix_6_5;
real b0_matrix_24_5;
real b0_matrix_5_6;
real b0_matrix_7_6;
real b0_matrix_26_6;
real b0_matrix_27_6;
real b0_matrix_6_7;
real b0_matrix_28_7;
real b0_matrix_29_7;
real b0_matrix_30_7;
real b0_matrix_31_7;
real b0_matrix_42_7;
real b0_matrix_3_8;
real b0_matrix_9_8;
real b0_matrix_32_8;
real b0_matrix_33_8;
real b0_matrix_8_9;
real b0_matrix_34_9;
real b0_matrix_35_9;
real b0_matrix_36_9;
real b0_matrix_37_9;
real b0_matrix_38_9;
real b0_matrix_39_9;
real b0_matrix_40_9;
real b0_matrix_2_10_DE;
real b0_matrix_2_10_NDE;
real b0_matrix_2_12_DE;
real b0_matrix_2_12_NDE;
real b0_matrix_2_14_DE;
real b0_matrix_2_14_NDE;
real b0_matrix_2_16_DE;
real b0_matrix_2_16_NDE;
real b0_matrix_2_18_DE;
real b0_matrix_2_18_NDE;
real b0_matrix_3_20_DE;
real b0_matrix_3_20_NDE;
real b0_matrix_3_22_DE;
real b0_matrix_3_22_NDE;
real b0_matrix_5_24_DE;
real b0_matrix_5_24_NDE;
real b0_matrix_6_26_DE;
real b0_matrix_6_26_NDE;
real b0_matrix_7_28_DE;
real b0_matrix_7_28_NDE;
real b0_matrix_7_30_DE;
real b0_matrix_7_30_NDE;
real b0_matrix_8_32_DE;
real b0_matrix_8_32_NDE;
real b0_matrix_9_34_DE;
real b0_matrix_9_34_NDE;
real b0_matrix_9_36;
real b0_matrix_9_37;
real b0_matrix_9_38;
real b0_matrix_2_41;
real b0_matrix_7_42;

real borigin1_matrix_8_3;
real borigin1_matrix_3_8;
real borigin1_matrix_9_8;
real borigin1_matrix_32_8;
real borigin1_matrix_33_8;
real borigin1_matrix_8_9;
real borigin1_matrix_34_9;
real borigin1_matrix_35_9;
real borigin1_matrix_36_9;
real borigin1_matrix_37_9;
real borigin1_matrix_38_9;
real borigin1_matrix_39_9;
real borigin1_matrix_40_9;
real borigin1_matrix_8_32_DE;
real borigin1_matrix_8_32_NDE;
real borigin1_matrix_9_34_DE;
real borigin1_matrix_9_34_NDE;
real borigin1_matrix_9_36;
real borigin1_matrix_9_37;
real borigin1_matrix_9_38;

real borigin2_matrix_8_3;
real borigin2_matrix_3_8;
real borigin2_matrix_9_8;
real borigin2_matrix_32_8;
real borigin2_matrix_33_8;
real borigin2_matrix_8_9;
real borigin2_matrix_34_9;
real borigin2_matrix_35_9;
real borigin2_matrix_36_9;
real borigin2_matrix_37_9;
real borigin2_matrix_38_9;
real borigin2_matrix_39_9;
real borigin2_matrix_40_9;
real borigin2_matrix_8_32_DE;
real borigin2_matrix_8_32_NDE;
real borigin2_matrix_9_34_DE;
real borigin2_matrix_9_34_NDE;
real borigin2_matrix_9_36;
real borigin2_matrix_9_37;
real borigin2_matrix_9_38;

real borigin3_matrix_8_3;
real borigin3_matrix_3_8;
real borigin3_matrix_9_8;
real borigin3_matrix_32_8;
real borigin3_matrix_33_8;
real borigin3_matrix_8_9;
real borigin3_matrix_34_9;
real borigin3_matrix_35_9;
real borigin3_matrix_36_9;
real borigin3_matrix_37_9;
real borigin3_matrix_38_9;
real borigin3_matrix_39_9;
real borigin3_matrix_40_9;
real borigin3_matrix_8_32_DE;
real borigin3_matrix_8_32_NDE;
real borigin3_matrix_9_34_DE;
real borigin3_matrix_9_34_NDE;
real borigin3_matrix_9_36;
real borigin3_matrix_9_37;
real borigin3_matrix_9_38;

real borigin4_matrix_8_3;
real borigin4_matrix_3_8;
real borigin4_matrix_9_8;
real borigin4_matrix_32_8;
real borigin4_matrix_33_8;
real borigin4_matrix_8_9;
real borigin4_matrix_34_9;
real borigin4_matrix_35_9;
real borigin4_matrix_36_9;
real borigin4_matrix_37_9;
real borigin4_matrix_38_9;
real borigin4_matrix_39_9;
real borigin4_matrix_40_9;
real borigin4_matrix_8_32_DE;
real borigin4_matrix_8_32_NDE;
real borigin4_matrix_9_34_DE;
real borigin4_matrix_9_34_NDE;
real borigin4_matrix_9_36;
real borigin4_matrix_9_37;
real borigin4_matrix_9_38;

real borigin5_matrix_8_3;
real borigin5_matrix_3_8;
real borigin5_matrix_9_8;
real borigin5_matrix_32_8;
real borigin5_matrix_33_8;
real borigin5_matrix_8_9;
real borigin5_matrix_34_9;
real borigin5_matrix_35_9;
real borigin5_matrix_36_9;
real borigin5_matrix_37_9;
real borigin5_matrix_38_9;
real borigin5_matrix_39_9;
real borigin5_matrix_40_9;
real borigin5_matrix_8_32_DE;
real borigin5_matrix_8_32_NDE;
real borigin5_matrix_9_34_DE;
real borigin5_matrix_9_34_NDE;
real borigin5_matrix_9_36;
real borigin5_matrix_9_37;
real borigin5_matrix_9_38;
  
}

transformed parameters {
  // Declare a matrix to store b0 params
  // matrix[43,43] b0_matrix;
  array[43,43] real b0_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // b0_matrix = rep_matrix(-100000, 43, 43);
  b0_matrix = rep_array(-100000, 43, 43);
  
  // Declare a matrix to store borigin1 params
  // matrix[43,43] borigin1_matrix;
  array[43,43] real borigin1_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Use a zero instead, since we already have the -100000 in the b0 term. So we just want these to have no effect
  // Non-zero elements will be overwritten
  // borigin1_matrix = rep_matrix(-100000, 43, 43);
  // borigin1_matrix = rep_matrix(0, 43, 43);
  borigin1_matrix = rep_array(0, 43, 43);
  
  // Declare a matrix to store borigin2 params
  // matrix[43,43] borigin2_matrix;
  array[43,43] real borigin2_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // borigin2_matrix = rep_matrix(-100000, 43, 43);
  // borigin2_matrix = rep_matrix(0, 43, 43);
  borigin2_matrix = rep_array(0, 43, 43);
  
  // Declare a matrix to store borigin3 params
  // matrix[43,43] borigin3_matrix;
  array[43,43] real borigin3_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // borigin3_matrix = rep_matrix(-100000, 43, 43);
  // borigin3_matrix = rep_matrix(0, 43, 43);
  borigin3_matrix = rep_array(0, 43, 43);
  
    // Declare a matrix to store borigin4 params
  // matrix[43,43] borigin4_matrix;
  array[43,43] real borigin4_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // borigin3_matrix = rep_matrix(-100000, 43, 43);
  // borigin4_matrix = rep_matrix(0, 43, 43);
  borigin4_matrix = rep_array(0, 43, 43);
  
    // Declare a matrix to store borigin5 params
  // matrix[43,43] borigin5_matrix;
  array[43,43] real borigin5_matrix;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // borigin3_matrix = rep_matrix(-100000, 43, 43);
  // borigin5_matrix = rep_matrix(0, 43, 43);
  borigin5_matrix = rep_array(0, 43, 43);
  
  // Finally, declare the parameter vector for the detection probability calculation
  det_eff_params = vector[34];
  
  
  
  // Populate this matrix with betas
  // Fill in all of the non-zero elements into the b0_matrix
b0_matrix_DE[2,1] = b0_matrix_2_1;
b0_matrix_NDE[2,1] = b0_matrix_2_1;
b0_matrix_DE[1,2] = b0_matrix_1_2;
b0_matrix_NDE[1,2] = b0_matrix_1_2;
b0_matrix_DE[3,2] = b0_matrix_3_2;
b0_matrix_NDE[3,2] = b0_matrix_3_2;
b0_matrix_DE[10,2] = b0_matrix_10_2;
b0_matrix_NDE[10,2] = b0_matrix_10_2;
b0_matrix_DE[11,2] = b0_matrix_11_2;
b0_matrix_NDE[11,2] = b0_matrix_11_2;
b0_matrix_DE[12,2] = b0_matrix_12_2;
b0_matrix_NDE[12,2] = b0_matrix_12_2;
b0_matrix_DE[13,2] = b0_matrix_13_2;
b0_matrix_NDE[13,2] = b0_matrix_13_2;
b0_matrix_DE[14,2] = b0_matrix_14_2;
b0_matrix_NDE[14,2] = b0_matrix_14_2;
b0_matrix_DE[16,2] = b0_matrix_16_2;
b0_matrix_NDE[16,2] = b0_matrix_16_2;
b0_matrix_DE[17,2] = b0_matrix_17_2;
b0_matrix_NDE[17,2] = b0_matrix_17_2;
b0_matrix_DE[18,2] = b0_matrix_18_2;
b0_matrix_NDE[18,2] = b0_matrix_18_2;
b0_matrix_DE[19,2] = b0_matrix_19_2;
b0_matrix_NDE[19,2] = b0_matrix_19_2;
b0_matrix_DE[41,2] = b0_matrix_41_2;
b0_matrix_NDE[41,2] = b0_matrix_41_2;
b0_matrix_DE[2,3] = b0_matrix_2_3;
b0_matrix_NDE[2,3] = b0_matrix_2_3;
b0_matrix_DE[8,3] = b0_matrix_8_3;
b0_matrix_NDE[8,3] = b0_matrix_8_3;
b0_matrix_DE[20,3] = b0_matrix_20_3;
b0_matrix_NDE[20,3] = b0_matrix_20_3;
b0_matrix_DE[21,3] = b0_matrix_21_3;
b0_matrix_NDE[21,3] = b0_matrix_21_3;
b0_matrix_DE[22,3] = b0_matrix_22_3;
b0_matrix_NDE[22,3] = b0_matrix_22_3;
b0_matrix_DE[23,3] = b0_matrix_23_3;
b0_matrix_NDE[23,3] = b0_matrix_23_3;
b0_matrix_DE[3,4] = b0_matrix_3_4;
b0_matrix_NDE[3,4] = b0_matrix_3_4;
b0_matrix_DE[5,4] = b0_matrix_5_4;
b0_matrix_NDE[5,4] = b0_matrix_5_4;
b0_matrix_DE[6,5] = b0_matrix_6_5;
b0_matrix_NDE[6,5] = b0_matrix_6_5;
b0_matrix_DE[24,5] = b0_matrix_24_5;
b0_matrix_NDE[24,5] = b0_matrix_24_5;
b0_matrix_DE[5,6] = b0_matrix_5_6;
b0_matrix_NDE[5,6] = b0_matrix_5_6;
b0_matrix_DE[7,6] = b0_matrix_7_6;
b0_matrix_NDE[7,6] = b0_matrix_7_6;
b0_matrix_DE[26,6] = b0_matrix_26_6;
b0_matrix_NDE[26,6] = b0_matrix_26_6;
b0_matrix_DE[27,6] = b0_matrix_27_6;
b0_matrix_NDE[27,6] = b0_matrix_27_6;
b0_matrix_DE[6,7] = b0_matrix_6_7;
b0_matrix_NDE[6,7] = b0_matrix_6_7;
b0_matrix_DE[28,7] = b0_matrix_28_7;
b0_matrix_NDE[28,7] = b0_matrix_28_7;
b0_matrix_DE[29,7] = b0_matrix_29_7;
b0_matrix_NDE[29,7] = b0_matrix_29_7;
b0_matrix_DE[30,7] = b0_matrix_30_7;
b0_matrix_NDE[30,7] = b0_matrix_30_7;
b0_matrix_DE[31,7] = b0_matrix_31_7;
b0_matrix_NDE[31,7] = b0_matrix_31_7;
b0_matrix_DE[42,7] = b0_matrix_42_7;
b0_matrix_NDE[42,7] = b0_matrix_42_7;
b0_matrix_DE[3,8] = b0_matrix_3_8;
b0_matrix_NDE[3,8] = b0_matrix_3_8;
b0_matrix_DE[9,8] = b0_matrix_9_8;
b0_matrix_NDE[9,8] = b0_matrix_9_8;
b0_matrix_DE[32,8] = b0_matrix_32_8;
b0_matrix_NDE[32,8] = b0_matrix_32_8;
b0_matrix_DE[33,8] = b0_matrix_33_8;
b0_matrix_NDE[33,8] = b0_matrix_33_8;
b0_matrix_DE[8,9] = b0_matrix_8_9;
b0_matrix_NDE[8,9] = b0_matrix_8_9;
b0_matrix_DE[34,9] = b0_matrix_34_9;
b0_matrix_NDE[34,9] = b0_matrix_34_9;
b0_matrix_DE[35,9] = b0_matrix_35_9;
b0_matrix_NDE[35,9] = b0_matrix_35_9;
b0_matrix_DE[36,9] = b0_matrix_36_9;
b0_matrix_NDE[36,9] = b0_matrix_36_9;
b0_matrix_DE[37,9] = b0_matrix_37_9;
b0_matrix_NDE[37,9] = b0_matrix_37_9;
b0_matrix_DE[38,9] = b0_matrix_38_9;
b0_matrix_NDE[38,9] = b0_matrix_38_9;
b0_matrix_DE[39,9] = b0_matrix_39_9;
b0_matrix_NDE[39,9] = b0_matrix_39_9;
b0_matrix_DE[40,9] = b0_matrix_40_9;
b0_matrix_NDE[40,9] = b0_matrix_40_9;
b0_matrix_DE[2,10] = b0_matrix_2_10_DE;
b0_matrix_NDE[2,10] = b0_matrix_2_10_NDE;
b0_matrix_DE[11,10] = b0_matrix_11_10;
b0_matrix_NDE[11,10] = b0_matrix_11_10;
b0_matrix_DE[2,11] = b0_matrix_2_11_DE;
b0_matrix_NDE[2,11] = b0_matrix_2_11_NDE;
b0_matrix_DE[10,11] = b0_matrix_10_11;
b0_matrix_NDE[10,11] = b0_matrix_10_11;
b0_matrix_DE[2,12] = b0_matrix_2_12_DE;
b0_matrix_NDE[2,12] = b0_matrix_2_12_NDE;
b0_matrix_DE[13,12] = b0_matrix_13_12;
b0_matrix_NDE[13,12] = b0_matrix_13_12;
b0_matrix_DE[2,13] = b0_matrix_2_13_DE;
b0_matrix_NDE[2,13] = b0_matrix_2_13_NDE;
b0_matrix_DE[12,13] = b0_matrix_12_13;
b0_matrix_NDE[12,13] = b0_matrix_12_13;
b0_matrix_DE[2,14] = b0_matrix_2_14_DE;
b0_matrix_NDE[2,14] = b0_matrix_2_14_NDE;
b0_matrix_DE[2,15] = b0_matrix_2_15_DE;
b0_matrix_NDE[2,15] = b0_matrix_2_15_NDE;
b0_matrix_DE[14,15] = b0_matrix_14_15;
b0_matrix_NDE[14,15] = b0_matrix_14_15;
b0_matrix_DE[2,16] = b0_matrix_2_16_DE;
b0_matrix_NDE[2,16] = b0_matrix_2_16_NDE;
b0_matrix_DE[17,16] = b0_matrix_17_16;
b0_matrix_NDE[17,16] = b0_matrix_17_16;
b0_matrix_DE[2,17] = b0_matrix_2_17_DE;
b0_matrix_NDE[2,17] = b0_matrix_2_17_NDE;
b0_matrix_DE[16,17] = b0_matrix_16_17;
b0_matrix_NDE[16,17] = b0_matrix_16_17;
b0_matrix_DE[2,18] = b0_matrix_2_18_DE;
b0_matrix_NDE[2,18] = b0_matrix_2_18_NDE;
b0_matrix_DE[19,18] = b0_matrix_19_18;
b0_matrix_NDE[19,18] = b0_matrix_19_18;
b0_matrix_DE[2,19] = b0_matrix_2_19_DE;
b0_matrix_NDE[2,19] = b0_matrix_2_19_NDE;
b0_matrix_DE[18,19] = b0_matrix_18_19;
b0_matrix_NDE[18,19] = b0_matrix_18_19;
b0_matrix_DE[3,20] = b0_matrix_3_20_DE;
b0_matrix_NDE[3,20] = b0_matrix_3_20_NDE;
b0_matrix_DE[21,20] = b0_matrix_21_20;
b0_matrix_NDE[21,20] = b0_matrix_21_20;
b0_matrix_DE[3,21] = b0_matrix_3_21_DE;
b0_matrix_NDE[3,21] = b0_matrix_3_21_NDE;
b0_matrix_DE[20,21] = b0_matrix_20_21;
b0_matrix_NDE[20,21] = b0_matrix_20_21;
b0_matrix_DE[3,22] = b0_matrix_3_22_DE;
b0_matrix_NDE[3,22] = b0_matrix_3_22_NDE;
b0_matrix_DE[23,22] = b0_matrix_23_22;
b0_matrix_NDE[23,22] = b0_matrix_23_22;
b0_matrix_DE[3,23] = b0_matrix_3_23_DE;
b0_matrix_NDE[3,23] = b0_matrix_3_23_NDE;
b0_matrix_DE[22,23] = b0_matrix_22_23;
b0_matrix_NDE[22,23] = b0_matrix_22_23;
b0_matrix_DE[5,24] = b0_matrix_5_24_DE;
b0_matrix_NDE[5,24] = b0_matrix_5_24_NDE;
b0_matrix_DE[5,25] = b0_matrix_5_25_DE;
b0_matrix_NDE[5,25] = b0_matrix_5_25_NDE;
b0_matrix_DE[24,25] = b0_matrix_24_25;
b0_matrix_NDE[24,25] = b0_matrix_24_25;
b0_matrix_DE[6,26] = b0_matrix_6_26_DE;
b0_matrix_NDE[6,26] = b0_matrix_6_26_NDE;
b0_matrix_DE[27,26] = b0_matrix_27_26;
b0_matrix_NDE[27,26] = b0_matrix_27_26;
b0_matrix_DE[6,27] = b0_matrix_6_27_DE;
b0_matrix_NDE[6,27] = b0_matrix_6_27_NDE;
b0_matrix_DE[26,27] = b0_matrix_26_27;
b0_matrix_NDE[26,27] = b0_matrix_26_27;
b0_matrix_DE[7,28] = b0_matrix_7_28_DE;
b0_matrix_NDE[7,28] = b0_matrix_7_28_NDE;
b0_matrix_DE[29,28] = b0_matrix_29_28;
b0_matrix_NDE[29,28] = b0_matrix_29_28;
b0_matrix_DE[7,29] = b0_matrix_7_29_DE;
b0_matrix_NDE[7,29] = b0_matrix_7_29_NDE;
b0_matrix_DE[28,29] = b0_matrix_28_29;
b0_matrix_NDE[28,29] = b0_matrix_28_29;
b0_matrix_DE[7,30] = b0_matrix_7_30_DE;
b0_matrix_NDE[7,30] = b0_matrix_7_30_NDE;
b0_matrix_DE[31,30] = b0_matrix_31_30;
b0_matrix_NDE[31,30] = b0_matrix_31_30;
b0_matrix_DE[7,31] = b0_matrix_7_31_DE;
b0_matrix_NDE[7,31] = b0_matrix_7_31_NDE;
b0_matrix_DE[30,31] = b0_matrix_30_31;
b0_matrix_NDE[30,31] = b0_matrix_30_31;
b0_matrix_DE[8,32] = b0_matrix_8_32_DE;
b0_matrix_NDE[8,32] = b0_matrix_8_32_NDE;
b0_matrix_DE[33,32] = b0_matrix_33_32;
b0_matrix_NDE[33,32] = b0_matrix_33_32;
b0_matrix_DE[8,33] = b0_matrix_8_33_DE;
b0_matrix_NDE[8,33] = b0_matrix_8_33_NDE;
b0_matrix_DE[32,33] = b0_matrix_32_33;
b0_matrix_NDE[32,33] = b0_matrix_32_33;
b0_matrix_DE[9,34] = b0_matrix_9_34_DE;
b0_matrix_NDE[9,34] = b0_matrix_9_34_NDE;
b0_matrix_DE[35,34] = b0_matrix_35_34;
b0_matrix_NDE[35,34] = b0_matrix_35_34;
b0_matrix_DE[9,35] = b0_matrix_9_35_DE;
b0_matrix_NDE[9,35] = b0_matrix_9_35_NDE;
b0_matrix_DE[34,35] = b0_matrix_34_35;
b0_matrix_NDE[34,35] = b0_matrix_34_35;
b0_matrix_DE[9,36] = b0_matrix_9_36;
b0_matrix_NDE[9,36] = b0_matrix_9_36;
b0_matrix_DE[9,37] = b0_matrix_9_37;
b0_matrix_NDE[9,37] = b0_matrix_9_37;
b0_matrix_DE[9,38] = b0_matrix_9_38;
b0_matrix_NDE[9,38] = b0_matrix_9_38;
b0_matrix_DE[40,39] = b0_matrix_40_39;
b0_matrix_NDE[40,39] = b0_matrix_40_39;
b0_matrix_DE[39,40] = b0_matrix_39_40;
b0_matrix_NDE[39,40] = b0_matrix_39_40;
b0_matrix_DE[2,41] = b0_matrix_2_41;
b0_matrix_NDE[2,41] = b0_matrix_2_41;
b0_matrix_DE[7,42] = b0_matrix_7_42;
b0_matrix_NDE[7,42] = b0_matrix_7_42;

borigin1_matrix_DE[8,3] = borigin1_matrix_8_3;
borigin1_matrix_NDE[8,3] = borigin1_matrix_8_3;
borigin1_matrix_DE[3,8] = borigin1_matrix_3_8;
borigin1_matrix_NDE[3,8] = borigin1_matrix_3_8;
borigin1_matrix_DE[9,8] = borigin1_matrix_9_8;
borigin1_matrix_NDE[9,8] = borigin1_matrix_9_8;
borigin1_matrix_DE[32,8] = borigin1_matrix_32_8;
borigin1_matrix_NDE[32,8] = borigin1_matrix_32_8;
borigin1_matrix_DE[33,8] = borigin1_matrix_33_8;
borigin1_matrix_NDE[33,8] = borigin1_matrix_33_8;
borigin1_matrix_DE[8,9] = borigin1_matrix_8_9;
borigin1_matrix_NDE[8,9] = borigin1_matrix_8_9;
borigin1_matrix_DE[34,9] = borigin1_matrix_34_9;
borigin1_matrix_NDE[34,9] = borigin1_matrix_34_9;
borigin1_matrix_DE[35,9] = borigin1_matrix_35_9;
borigin1_matrix_NDE[35,9] = borigin1_matrix_35_9;
borigin1_matrix_DE[36,9] = borigin1_matrix_36_9;
borigin1_matrix_NDE[36,9] = borigin1_matrix_36_9;
borigin1_matrix_DE[37,9] = borigin1_matrix_37_9;
borigin1_matrix_NDE[37,9] = borigin1_matrix_37_9;
borigin1_matrix_DE[38,9] = borigin1_matrix_38_9;
borigin1_matrix_NDE[38,9] = borigin1_matrix_38_9;
borigin1_matrix_DE[39,9] = borigin1_matrix_39_9;
borigin1_matrix_NDE[39,9] = borigin1_matrix_39_9;
borigin1_matrix_DE[40,9] = borigin1_matrix_40_9;
borigin1_matrix_NDE[40,9] = borigin1_matrix_40_9;
borigin1_matrix_DE[8,32] = borigin1_matrix_8_32_DE;
borigin1_matrix_NDE[8,32] = borigin1_matrix_8_32_NDE;
borigin1_matrix_DE[33,32] = borigin1_matrix_33_32;
borigin1_matrix_NDE[33,32] = borigin1_matrix_33_32;
borigin1_matrix_DE[8,33] = borigin1_matrix_8_33_DE;
borigin1_matrix_NDE[8,33] = borigin1_matrix_8_33_NDE;
borigin1_matrix_DE[32,33] = borigin1_matrix_32_33;
borigin1_matrix_NDE[32,33] = borigin1_matrix_32_33;
borigin1_matrix_DE[9,34] = borigin1_matrix_9_34_DE;
borigin1_matrix_NDE[9,34] = borigin1_matrix_9_34_NDE;
borigin1_matrix_DE[35,34] = borigin1_matrix_35_34;
borigin1_matrix_NDE[35,34] = borigin1_matrix_35_34;
borigin1_matrix_DE[9,35] = borigin1_matrix_9_35_DE;
borigin1_matrix_NDE[9,35] = borigin1_matrix_9_35_NDE;
borigin1_matrix_DE[34,35] = borigin1_matrix_34_35;
borigin1_matrix_NDE[34,35] = borigin1_matrix_34_35;
borigin1_matrix_DE[9,36] = borigin1_matrix_9_36;
borigin1_matrix_NDE[9,36] = borigin1_matrix_9_36;
borigin1_matrix_DE[9,37] = borigin1_matrix_9_37;
borigin1_matrix_NDE[9,37] = borigin1_matrix_9_37;
borigin1_matrix_DE[9,38] = borigin1_matrix_9_38;
borigin1_matrix_NDE[9,38] = borigin1_matrix_9_38;
borigin1_matrix_DE[40,39] = borigin1_matrix_40_39;
borigin1_matrix_NDE[40,39] = borigin1_matrix_40_39;
borigin1_matrix_DE[39,40] = borigin1_matrix_39_40;
borigin1_matrix_NDE[39,40] = borigin1_matrix_39_40;

borigin2_matrix_DE[8,3] = borigin2_matrix_8_3;
borigin2_matrix_NDE[8,3] = borigin2_matrix_8_3;
borigin2_matrix_DE[3,8] = borigin2_matrix_3_8;
borigin2_matrix_NDE[3,8] = borigin2_matrix_3_8;
borigin2_matrix_DE[9,8] = borigin2_matrix_9_8;
borigin2_matrix_NDE[9,8] = borigin2_matrix_9_8;
borigin2_matrix_DE[32,8] = borigin2_matrix_32_8;
borigin2_matrix_NDE[32,8] = borigin2_matrix_32_8;
borigin2_matrix_DE[33,8] = borigin2_matrix_33_8;
borigin2_matrix_NDE[33,8] = borigin2_matrix_33_8;
borigin2_matrix_DE[8,9] = borigin2_matrix_8_9;
borigin2_matrix_NDE[8,9] = borigin2_matrix_8_9;
borigin2_matrix_DE[34,9] = borigin2_matrix_34_9;
borigin2_matrix_NDE[34,9] = borigin2_matrix_34_9;
borigin2_matrix_DE[35,9] = borigin2_matrix_35_9;
borigin2_matrix_NDE[35,9] = borigin2_matrix_35_9;
borigin2_matrix_DE[36,9] = borigin2_matrix_36_9;
borigin2_matrix_NDE[36,9] = borigin2_matrix_36_9;
borigin2_matrix_DE[37,9] = borigin2_matrix_37_9;
borigin2_matrix_NDE[37,9] = borigin2_matrix_37_9;
borigin2_matrix_DE[38,9] = borigin2_matrix_38_9;
borigin2_matrix_NDE[38,9] = borigin2_matrix_38_9;
borigin2_matrix_DE[39,9] = borigin2_matrix_39_9;
borigin2_matrix_NDE[39,9] = borigin2_matrix_39_9;
borigin2_matrix_DE[40,9] = borigin2_matrix_40_9;
borigin2_matrix_NDE[40,9] = borigin2_matrix_40_9;
borigin2_matrix_DE[8,32] = borigin2_matrix_8_32_DE;
borigin2_matrix_NDE[8,32] = borigin2_matrix_8_32_NDE;
borigin2_matrix_DE[33,32] = borigin2_matrix_33_32;
borigin2_matrix_NDE[33,32] = borigin2_matrix_33_32;
borigin2_matrix_DE[8,33] = borigin2_matrix_8_33_DE;
borigin2_matrix_NDE[8,33] = borigin2_matrix_8_33_NDE;
borigin2_matrix_DE[32,33] = borigin2_matrix_32_33;
borigin2_matrix_NDE[32,33] = borigin2_matrix_32_33;
borigin2_matrix_DE[9,34] = borigin2_matrix_9_34_DE;
borigin2_matrix_NDE[9,34] = borigin2_matrix_9_34_NDE;
borigin2_matrix_DE[35,34] = borigin2_matrix_35_34;
borigin2_matrix_NDE[35,34] = borigin2_matrix_35_34;
borigin2_matrix_DE[9,35] = borigin2_matrix_9_35_DE;
borigin2_matrix_NDE[9,35] = borigin2_matrix_9_35_NDE;
borigin2_matrix_DE[34,35] = borigin2_matrix_34_35;
borigin2_matrix_NDE[34,35] = borigin2_matrix_34_35;
borigin2_matrix_DE[9,36] = borigin2_matrix_9_36;
borigin2_matrix_NDE[9,36] = borigin2_matrix_9_36;
borigin2_matrix_DE[9,37] = borigin2_matrix_9_37;
borigin2_matrix_NDE[9,37] = borigin2_matrix_9_37;
borigin2_matrix_DE[9,38] = borigin2_matrix_9_38;
borigin2_matrix_NDE[9,38] = borigin2_matrix_9_38;
borigin2_matrix_DE[40,39] = borigin2_matrix_40_39;
borigin2_matrix_NDE[40,39] = borigin2_matrix_40_39;
borigin2_matrix_DE[39,40] = borigin2_matrix_39_40;
borigin2_matrix_NDE[39,40] = borigin2_matrix_39_40;

borigin3_matrix_DE[8,3] = borigin3_matrix_8_3;
borigin3_matrix_NDE[8,3] = borigin3_matrix_8_3;
borigin3_matrix_DE[3,8] = borigin3_matrix_3_8;
borigin3_matrix_NDE[3,8] = borigin3_matrix_3_8;
borigin3_matrix_DE[9,8] = borigin3_matrix_9_8;
borigin3_matrix_NDE[9,8] = borigin3_matrix_9_8;
borigin3_matrix_DE[32,8] = borigin3_matrix_32_8;
borigin3_matrix_NDE[32,8] = borigin3_matrix_32_8;
borigin3_matrix_DE[33,8] = borigin3_matrix_33_8;
borigin3_matrix_NDE[33,8] = borigin3_matrix_33_8;
borigin3_matrix_DE[8,9] = borigin3_matrix_8_9;
borigin3_matrix_NDE[8,9] = borigin3_matrix_8_9;
borigin3_matrix_DE[34,9] = borigin3_matrix_34_9;
borigin3_matrix_NDE[34,9] = borigin3_matrix_34_9;
borigin3_matrix_DE[35,9] = borigin3_matrix_35_9;
borigin3_matrix_NDE[35,9] = borigin3_matrix_35_9;
borigin3_matrix_DE[36,9] = borigin3_matrix_36_9;
borigin3_matrix_NDE[36,9] = borigin3_matrix_36_9;
borigin3_matrix_DE[37,9] = borigin3_matrix_37_9;
borigin3_matrix_NDE[37,9] = borigin3_matrix_37_9;
borigin3_matrix_DE[38,9] = borigin3_matrix_38_9;
borigin3_matrix_NDE[38,9] = borigin3_matrix_38_9;
borigin3_matrix_DE[39,9] = borigin3_matrix_39_9;
borigin3_matrix_NDE[39,9] = borigin3_matrix_39_9;
borigin3_matrix_DE[40,9] = borigin3_matrix_40_9;
borigin3_matrix_NDE[40,9] = borigin3_matrix_40_9;
borigin3_matrix_DE[8,32] = borigin3_matrix_8_32_DE;
borigin3_matrix_NDE[8,32] = borigin3_matrix_8_32_NDE;
borigin3_matrix_DE[33,32] = borigin3_matrix_33_32;
borigin3_matrix_NDE[33,32] = borigin3_matrix_33_32;
borigin3_matrix_DE[8,33] = borigin3_matrix_8_33_DE;
borigin3_matrix_NDE[8,33] = borigin3_matrix_8_33_NDE;
borigin3_matrix_DE[32,33] = borigin3_matrix_32_33;
borigin3_matrix_NDE[32,33] = borigin3_matrix_32_33;
borigin3_matrix_DE[9,34] = borigin3_matrix_9_34_DE;
borigin3_matrix_NDE[9,34] = borigin3_matrix_9_34_NDE;
borigin3_matrix_DE[35,34] = borigin3_matrix_35_34;
borigin3_matrix_NDE[35,34] = borigin3_matrix_35_34;
borigin3_matrix_DE[9,35] = borigin3_matrix_9_35_DE;
borigin3_matrix_NDE[9,35] = borigin3_matrix_9_35_NDE;
borigin3_matrix_DE[34,35] = borigin3_matrix_34_35;
borigin3_matrix_NDE[34,35] = borigin3_matrix_34_35;
borigin3_matrix_DE[9,36] = borigin3_matrix_9_36;
borigin3_matrix_NDE[9,36] = borigin3_matrix_9_36;
borigin3_matrix_DE[9,37] = borigin3_matrix_9_37;
borigin3_matrix_NDE[9,37] = borigin3_matrix_9_37;
borigin3_matrix_DE[9,38] = borigin3_matrix_9_38;
borigin3_matrix_NDE[9,38] = borigin3_matrix_9_38;
borigin3_matrix_DE[40,39] = borigin3_matrix_40_39;
borigin3_matrix_NDE[40,39] = borigin3_matrix_40_39;
borigin3_matrix_DE[39,40] = borigin3_matrix_39_40;
borigin3_matrix_NDE[39,40] = borigin3_matrix_39_40;

borigin4_matrix_DE[8,3] = borigin4_matrix_8_3;
borigin4_matrix_NDE[8,3] = borigin4_matrix_8_3;
borigin4_matrix_DE[3,8] = borigin4_matrix_3_8;
borigin4_matrix_NDE[3,8] = borigin4_matrix_3_8;
borigin4_matrix_DE[9,8] = borigin4_matrix_9_8;
borigin4_matrix_NDE[9,8] = borigin4_matrix_9_8;
borigin4_matrix_DE[32,8] = borigin4_matrix_32_8;
borigin4_matrix_NDE[32,8] = borigin4_matrix_32_8;
borigin4_matrix_DE[33,8] = borigin4_matrix_33_8;
borigin4_matrix_NDE[33,8] = borigin4_matrix_33_8;
borigin4_matrix_DE[8,9] = borigin4_matrix_8_9;
borigin4_matrix_NDE[8,9] = borigin4_matrix_8_9;
borigin4_matrix_DE[34,9] = borigin4_matrix_34_9;
borigin4_matrix_NDE[34,9] = borigin4_matrix_34_9;
borigin4_matrix_DE[35,9] = borigin4_matrix_35_9;
borigin4_matrix_NDE[35,9] = borigin4_matrix_35_9;
borigin4_matrix_DE[36,9] = borigin4_matrix_36_9;
borigin4_matrix_NDE[36,9] = borigin4_matrix_36_9;
borigin4_matrix_DE[37,9] = borigin4_matrix_37_9;
borigin4_matrix_NDE[37,9] = borigin4_matrix_37_9;
borigin4_matrix_DE[38,9] = borigin4_matrix_38_9;
borigin4_matrix_NDE[38,9] = borigin4_matrix_38_9;
borigin4_matrix_DE[39,9] = borigin4_matrix_39_9;
borigin4_matrix_NDE[39,9] = borigin4_matrix_39_9;
borigin4_matrix_DE[40,9] = borigin4_matrix_40_9;
borigin4_matrix_NDE[40,9] = borigin4_matrix_40_9;
borigin4_matrix_DE[8,32] = borigin4_matrix_8_32_DE;
borigin4_matrix_NDE[8,32] = borigin4_matrix_8_32_NDE;
borigin4_matrix_DE[33,32] = borigin4_matrix_33_32;
borigin4_matrix_NDE[33,32] = borigin4_matrix_33_32;
borigin4_matrix_DE[8,33] = borigin4_matrix_8_33_DE;
borigin4_matrix_NDE[8,33] = borigin4_matrix_8_33_NDE;
borigin4_matrix_DE[32,33] = borigin4_matrix_32_33;
borigin4_matrix_NDE[32,33] = borigin4_matrix_32_33;
borigin4_matrix_DE[9,34] = borigin4_matrix_9_34_DE;
borigin4_matrix_NDE[9,34] = borigin4_matrix_9_34_NDE;
borigin4_matrix_DE[35,34] = borigin4_matrix_35_34;
borigin4_matrix_NDE[35,34] = borigin4_matrix_35_34;
borigin4_matrix_DE[9,35] = borigin4_matrix_9_35_DE;
borigin4_matrix_NDE[9,35] = borigin4_matrix_9_35_NDE;
borigin4_matrix_DE[34,35] = borigin4_matrix_34_35;
borigin4_matrix_NDE[34,35] = borigin4_matrix_34_35;
borigin4_matrix_DE[9,36] = borigin4_matrix_9_36;
borigin4_matrix_NDE[9,36] = borigin4_matrix_9_36;
borigin4_matrix_DE[9,37] = borigin4_matrix_9_37;
borigin4_matrix_NDE[9,37] = borigin4_matrix_9_37;
borigin4_matrix_DE[9,38] = borigin4_matrix_9_38;
borigin4_matrix_NDE[9,38] = borigin4_matrix_9_38;
borigin4_matrix_DE[40,39] = borigin4_matrix_40_39;
borigin4_matrix_NDE[40,39] = borigin4_matrix_40_39;
borigin4_matrix_DE[39,40] = borigin4_matrix_39_40;
borigin4_matrix_NDE[39,40] = borigin4_matrix_39_40;

borigin5_matrix_DE[8,3] = borigin5_matrix_8_3;
borigin5_matrix_NDE[8,3] = borigin5_matrix_8_3;
borigin5_matrix_DE[3,8] = borigin5_matrix_3_8;
borigin5_matrix_NDE[3,8] = borigin5_matrix_3_8;
borigin5_matrix_DE[9,8] = borigin5_matrix_9_8;
borigin5_matrix_NDE[9,8] = borigin5_matrix_9_8;
borigin5_matrix_DE[32,8] = borigin5_matrix_32_8;
borigin5_matrix_NDE[32,8] = borigin5_matrix_32_8;
borigin5_matrix_DE[33,8] = borigin5_matrix_33_8;
borigin5_matrix_NDE[33,8] = borigin5_matrix_33_8;
borigin5_matrix_DE[8,9] = borigin5_matrix_8_9;
borigin5_matrix_NDE[8,9] = borigin5_matrix_8_9;
borigin5_matrix_DE[34,9] = borigin5_matrix_34_9;
borigin5_matrix_NDE[34,9] = borigin5_matrix_34_9;
borigin5_matrix_DE[35,9] = borigin5_matrix_35_9;
borigin5_matrix_NDE[35,9] = borigin5_matrix_35_9;
borigin5_matrix_DE[36,9] = borigin5_matrix_36_9;
borigin5_matrix_NDE[36,9] = borigin5_matrix_36_9;
borigin5_matrix_DE[37,9] = borigin5_matrix_37_9;
borigin5_matrix_NDE[37,9] = borigin5_matrix_37_9;
borigin5_matrix_DE[38,9] = borigin5_matrix_38_9;
borigin5_matrix_NDE[38,9] = borigin5_matrix_38_9;
borigin5_matrix_DE[39,9] = borigin5_matrix_39_9;
borigin5_matrix_NDE[39,9] = borigin5_matrix_39_9;
borigin5_matrix_DE[40,9] = borigin5_matrix_40_9;
borigin5_matrix_NDE[40,9] = borigin5_matrix_40_9;
borigin5_matrix_DE[8,32] = borigin5_matrix_8_32_DE;
borigin5_matrix_NDE[8,32] = borigin5_matrix_8_32_NDE;
borigin5_matrix_DE[33,32] = borigin5_matrix_33_32;
borigin5_matrix_NDE[33,32] = borigin5_matrix_33_32;
borigin5_matrix_DE[8,33] = borigin5_matrix_8_33_DE;
borigin5_matrix_NDE[8,33] = borigin5_matrix_8_33_NDE;
borigin5_matrix_DE[32,33] = borigin5_matrix_32_33;
borigin5_matrix_NDE[32,33] = borigin5_matrix_32_33;
borigin5_matrix_DE[9,34] = borigin5_matrix_9_34_DE;
borigin5_matrix_NDE[9,34] = borigin5_matrix_9_34_NDE;
borigin5_matrix_DE[35,34] = borigin5_matrix_35_34;
borigin5_matrix_NDE[35,34] = borigin5_matrix_35_34;
borigin5_matrix_DE[9,35] = borigin5_matrix_9_35_DE;
borigin5_matrix_NDE[9,35] = borigin5_matrix_9_35_NDE;
borigin5_matrix_DE[34,35] = borigin5_matrix_34_35;
borigin5_matrix_NDE[34,35] = borigin5_matrix_34_35;
borigin5_matrix_DE[9,36] = borigin5_matrix_9_36;
borigin5_matrix_NDE[9,36] = borigin5_matrix_9_36;
borigin5_matrix_DE[9,37] = borigin5_matrix_9_37;
borigin5_matrix_NDE[9,37] = borigin5_matrix_9_37;
borigin5_matrix_DE[9,38] = borigin5_matrix_9_38;
borigin5_matrix_NDE[9,38] = borigin5_matrix_9_38;
borigin5_matrix_DE[40,39] = borigin5_matrix_40_39;
borigin5_matrix_NDE[40,39] = borigin5_matrix_40_39;
borigin5_matrix_DE[39,40] = borigin5_matrix_39_40;
borigin5_matrix_NDE[39,40] = borigin5_matrix_39_40;

#### Calculate movement probabilities as derived parameters
matrix[43,43] origin1_probs;
origin1_probs = rep_matrix(0, 43, 43);
matrix[43,43] origin2_probs;
origin2_probs = rep_matrix(0, 43, 43);
matrix[43,43] origin3_probs;
origin3_probs = rep_matrix(0, 43, 43);
matrix[43,43] origin4_probs;
origin4_probs = rep_matrix(0, 43, 43);
matrix[43,43] origin5_probs;
origin5_probs = rep_matrix(0, 43, 43);
matrix[43,43] origin6_probs;
origin6_probs = rep_matrix(0, 43, 43);

  # for each of the non-loss states:
for (i in 1:42){
  for (j in 1:42){
    origin1_probs[i,j] = exp(b0_matrix[i,j]+ borigin1_matrix[i,j])/(1 + sum(exp(to_row_vector(b0_matrix[i,]) + to_row_vector(borigin1_matrix[i,]))));
    origin2_probs[i,j] = exp(b0_matrix[i,j]+ borigin2_matrix[i,j])/(1 + sum(exp(to_row_vector(b0_matrix[i,]) + to_row_vector(borigin2_matrix[i,]))));
    origin3_probs[i,j] = exp(b0_matrix[i,j]+ borigin3_matrix[i,j])/(1 + sum(exp(to_row_vector(b0_matrix[i,]) + to_row_vector(borigin3_matrix[i,]))));
    origin4_probs[i,j] = exp(b0_matrix[i,j]+ borigin4_matrix[i,j])/(1 + sum(exp(to_row_vector(b0_matrix[i,]) + to_row_vector(borigin4_matrix[i,]))));
    origin5_probs[i,j] = exp(b0_matrix[i,j]+ borigin5_matrix[i,j])/(1 + sum(exp(to_row_vector(b0_matrix[i,]) + to_row_vector(borigin5_matrix[i,]))));
    origin6_probs[i,j] = exp(b0_matrix[i,j]- borigin1_matrix[i,j] - borigin2_matrix[i,j] - borigin3_matrix[i,j] - borigin4_matrix[i,j] - borigin5_matrix[i,j])/
    (1 + sum(exp(to_row_vector(b0_matrix[i,]) - to_row_vector(borigin1_matrix[i,]) - to_row_vector(borigin2_matrix[i,]) - to_row_vector(borigin3_matrix[i,]) - 
    to_row_vector(borigin4_matrix[i,]) - to_row_vector(borigin5_matrix[i,]))));
  }
}

# calculate loss states
for (i in 1:42){
  origin1_probs[i,43] = 1 - sum(origin1_probs[i,1:42]);
  origin2_probs[i,43] = 1 - sum(origin2_probs[i,1:42]);
  origin3_probs[i,43] = 1 - sum(origin3_probs[i,1:42]);
  origin4_probs[i,43] = 1 - sum(origin4_probs[i,1:42]);
  origin5_probs[i,43] = 1 - sum(origin5_probs[i,1:42]);
  origin6_probs[i,43] = 1 - sum(origin6_probs[i,1:42]);
  
}

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
b0_matrix_11_2 ~ normal(0,10);
b0_matrix_12_2 ~ normal(0,10);
b0_matrix_13_2 ~ normal(0,10);
b0_matrix_14_2 ~ normal(0,10);
b0_matrix_16_2 ~ normal(0,10);
b0_matrix_17_2 ~ normal(0,10);
b0_matrix_18_2 ~ normal(0,10);
b0_matrix_19_2 ~ normal(0,10);
b0_matrix_41_2 ~ normal(0,10);
b0_matrix_2_3 ~ normal(0,10);
b0_matrix_8_3 ~ normal(0,10);
b0_matrix_20_3 ~ normal(0,10);
b0_matrix_21_3 ~ normal(0,10);
b0_matrix_22_3 ~ normal(0,10);
b0_matrix_23_3 ~ normal(0,10);
b0_matrix_3_4 ~ normal(0,10);
b0_matrix_5_4 ~ normal(0,10);
b0_matrix_6_5 ~ normal(0,10);
b0_matrix_24_5 ~ normal(0,10);
b0_matrix_5_6 ~ normal(0,10);
b0_matrix_7_6 ~ normal(0,10);
b0_matrix_26_6 ~ normal(0,10);
b0_matrix_27_6 ~ normal(0,10);
b0_matrix_6_7 ~ normal(0,10);
b0_matrix_28_7 ~ normal(0,10);
b0_matrix_29_7 ~ normal(0,10);
b0_matrix_30_7 ~ normal(0,10);
b0_matrix_31_7 ~ normal(0,10);
b0_matrix_42_7 ~ normal(0,10);
b0_matrix_3_8 ~ normal(0,10);
b0_matrix_9_8 ~ normal(0,10);
b0_matrix_32_8 ~ normal(0,10);
b0_matrix_33_8 ~ normal(0,10);
b0_matrix_8_9 ~ normal(0,10);
b0_matrix_34_9 ~ normal(0,10);
b0_matrix_35_9 ~ normal(0,10);
b0_matrix_36_9 ~ normal(0,10);
b0_matrix_37_9 ~ normal(0,10);
b0_matrix_38_9 ~ normal(0,10);
b0_matrix_39_9 ~ normal(0,10);
b0_matrix_40_9 ~ normal(0,10);
b0_matrix_2_10_DE ~ normal(0,10);
b0_matrix_2_10_NDE ~ normal(0,10);
b0_matrix_11_10 ~ normal(0,10);
b0_matrix_2_11_DE ~ normal(0,10);
b0_matrix_2_11_NDE ~ normal(0,10);
b0_matrix_10_11 ~ normal(0,10);
b0_matrix_2_12_DE ~ normal(0,10);
b0_matrix_2_12_NDE ~ normal(0,10);
b0_matrix_13_12 ~ normal(0,10);
b0_matrix_2_13_DE ~ normal(0,10);
b0_matrix_2_13_NDE ~ normal(0,10);
b0_matrix_12_13 ~ normal(0,10);
b0_matrix_2_14_DE ~ normal(0,10);
b0_matrix_2_14_NDE ~ normal(0,10);
b0_matrix_2_15_DE ~ normal(0,10);
b0_matrix_2_15_NDE ~ normal(0,10);
b0_matrix_14_15 ~ normal(0,10);
b0_matrix_2_16_DE ~ normal(0,10);
b0_matrix_2_16_NDE ~ normal(0,10);
b0_matrix_17_16 ~ normal(0,10);
b0_matrix_2_17_DE ~ normal(0,10);
b0_matrix_2_17_NDE ~ normal(0,10);
b0_matrix_16_17 ~ normal(0,10);
b0_matrix_2_18_DE ~ normal(0,10);
b0_matrix_2_18_NDE ~ normal(0,10);
b0_matrix_19_18 ~ normal(0,10);
b0_matrix_2_19_DE ~ normal(0,10);
b0_matrix_2_19_NDE ~ normal(0,10);
b0_matrix_18_19 ~ normal(0,10);
b0_matrix_3_20_DE ~ normal(0,10);
b0_matrix_3_20_NDE ~ normal(0,10);
b0_matrix_21_20 ~ normal(0,10);
b0_matrix_3_21_DE ~ normal(0,10);
b0_matrix_3_21_NDE ~ normal(0,10);
b0_matrix_20_21 ~ normal(0,10);
b0_matrix_3_22_DE ~ normal(0,10);
b0_matrix_3_22_NDE ~ normal(0,10);
b0_matrix_23_22 ~ normal(0,10);
b0_matrix_3_23_DE ~ normal(0,10);
b0_matrix_3_23_NDE ~ normal(0,10);
b0_matrix_22_23 ~ normal(0,10);
b0_matrix_5_24_DE ~ normal(0,10);
b0_matrix_5_24_NDE ~ normal(0,10);
b0_matrix_5_25_DE ~ normal(0,10);
b0_matrix_5_25_NDE ~ normal(0,10);
b0_matrix_24_25 ~ normal(0,10);
b0_matrix_6_26_DE ~ normal(0,10);
b0_matrix_6_26_NDE ~ normal(0,10);
b0_matrix_27_26 ~ normal(0,10);
b0_matrix_6_27_DE ~ normal(0,10);
b0_matrix_6_27_NDE ~ normal(0,10);
b0_matrix_26_27 ~ normal(0,10);
b0_matrix_7_28_DE ~ normal(0,10);
b0_matrix_7_28_NDE ~ normal(0,10);
b0_matrix_29_28 ~ normal(0,10);
b0_matrix_7_29_DE ~ normal(0,10);
b0_matrix_7_29_NDE ~ normal(0,10);
b0_matrix_28_29 ~ normal(0,10);
b0_matrix_7_30_DE ~ normal(0,10);
b0_matrix_7_30_NDE ~ normal(0,10);
b0_matrix_31_30 ~ normal(0,10);
b0_matrix_7_31_DE ~ normal(0,10);
b0_matrix_7_31_NDE ~ normal(0,10);
b0_matrix_30_31 ~ normal(0,10);
b0_matrix_8_32_DE ~ normal(0,10);
b0_matrix_8_32_NDE ~ normal(0,10);
b0_matrix_33_32 ~ normal(0,10);
b0_matrix_8_33_DE ~ normal(0,10);
b0_matrix_8_33_NDE ~ normal(0,10);
b0_matrix_32_33 ~ normal(0,10);
b0_matrix_9_34_DE ~ normal(0,10);
b0_matrix_9_34_NDE ~ normal(0,10);
b0_matrix_35_34 ~ normal(0,10);
b0_matrix_9_35_DE ~ normal(0,10);
b0_matrix_9_35_NDE ~ normal(0,10);
b0_matrix_34_35 ~ normal(0,10);
b0_matrix_9_36 ~ normal(0,10);
b0_matrix_9_37 ~ normal(0,10);
b0_matrix_9_38 ~ normal(0,10);
b0_matrix_40_39 ~ normal(0,10);
b0_matrix_39_40 ~ normal(0,10);
b0_matrix_2_41 ~ normal(0,10);
b0_matrix_7_42 ~ normal(0,10);

borigin1_matrix_8_3 ~ normal(0,10);
borigin1_matrix_3_8 ~ normal(0,10);
borigin1_matrix_9_8 ~ normal(0,10);
borigin1_matrix_32_8 ~ normal(0,10);
borigin1_matrix_33_8 ~ normal(0,10);
borigin1_matrix_8_9 ~ normal(0,10);
borigin1_matrix_34_9 ~ normal(0,10);
borigin1_matrix_35_9 ~ normal(0,10);
borigin1_matrix_36_9 ~ normal(0,10);
borigin1_matrix_37_9 ~ normal(0,10);
borigin1_matrix_38_9 ~ normal(0,10);
borigin1_matrix_39_9 ~ normal(0,10);
borigin1_matrix_40_9 ~ normal(0,10);
borigin1_matrix_8_32_DE ~ normal(0,10);
borigin1_matrix_8_32_NDE ~ normal(0,10);
borigin1_matrix_33_32 ~ normal(0,10);
borigin1_matrix_8_33_DE ~ normal(0,10);
borigin1_matrix_8_33_NDE ~ normal(0,10);
borigin1_matrix_32_33 ~ normal(0,10);
borigin1_matrix_9_34_DE ~ normal(0,10);
borigin1_matrix_9_34_NDE ~ normal(0,10);
borigin1_matrix_35_34 ~ normal(0,10);
borigin1_matrix_9_35_DE ~ normal(0,10);
borigin1_matrix_9_35_NDE ~ normal(0,10);
borigin1_matrix_34_35 ~ normal(0,10);
borigin1_matrix_9_36 ~ normal(0,10);
borigin1_matrix_9_37 ~ normal(0,10);
borigin1_matrix_9_38 ~ normal(0,10);
borigin1_matrix_40_39 ~ normal(0,10);
borigin1_matrix_39_40 ~ normal(0,10);

borigin2_matrix_8_3 ~ normal(0,10);
borigin2_matrix_3_8 ~ normal(0,10);
borigin2_matrix_9_8 ~ normal(0,10);
borigin2_matrix_32_8 ~ normal(0,10);
borigin2_matrix_33_8 ~ normal(0,10);
borigin2_matrix_8_9 ~ normal(0,10);
borigin2_matrix_34_9 ~ normal(0,10);
borigin2_matrix_35_9 ~ normal(0,10);
borigin2_matrix_36_9 ~ normal(0,10);
borigin2_matrix_37_9 ~ normal(0,10);
borigin2_matrix_38_9 ~ normal(0,10);
borigin2_matrix_39_9 ~ normal(0,10);
borigin2_matrix_40_9 ~ normal(0,10);
borigin2_matrix_8_32_DE ~ normal(0,10);
borigin2_matrix_8_32_NDE ~ normal(0,10);
borigin2_matrix_33_32 ~ normal(0,10);
borigin2_matrix_8_33_DE ~ normal(0,10);
borigin2_matrix_8_33_NDE ~ normal(0,10);
borigin2_matrix_32_33 ~ normal(0,10);
borigin2_matrix_9_34_DE ~ normal(0,10);
borigin2_matrix_9_34_NDE ~ normal(0,10);
borigin2_matrix_35_34 ~ normal(0,10);
borigin2_matrix_9_35_DE ~ normal(0,10);
borigin2_matrix_9_35_NDE ~ normal(0,10);
borigin2_matrix_34_35 ~ normal(0,10);
borigin2_matrix_9_36 ~ normal(0,10);
borigin2_matrix_9_37 ~ normal(0,10);
borigin2_matrix_9_38 ~ normal(0,10);
borigin2_matrix_40_39 ~ normal(0,10);
borigin2_matrix_39_40 ~ normal(0,10);

borigin3_matrix_8_3 ~ normal(0,10);
borigin3_matrix_3_8 ~ normal(0,10);
borigin3_matrix_9_8 ~ normal(0,10);
borigin3_matrix_32_8 ~ normal(0,10);
borigin3_matrix_33_8 ~ normal(0,10);
borigin3_matrix_8_9 ~ normal(0,10);
borigin3_matrix_34_9 ~ normal(0,10);
borigin3_matrix_35_9 ~ normal(0,10);
borigin3_matrix_36_9 ~ normal(0,10);
borigin3_matrix_37_9 ~ normal(0,10);
borigin3_matrix_38_9 ~ normal(0,10);
borigin3_matrix_39_9 ~ normal(0,10);
borigin3_matrix_40_9 ~ normal(0,10);
borigin3_matrix_8_32_DE ~ normal(0,10);
borigin3_matrix_8_32_NDE ~ normal(0,10);
borigin3_matrix_33_32 ~ normal(0,10);
borigin3_matrix_8_33_DE ~ normal(0,10);
borigin3_matrix_8_33_NDE ~ normal(0,10);
borigin3_matrix_32_33 ~ normal(0,10);
borigin3_matrix_9_34_DE ~ normal(0,10);
borigin3_matrix_9_34_NDE ~ normal(0,10);
borigin3_matrix_35_34 ~ normal(0,10);
borigin3_matrix_9_35_DE ~ normal(0,10);
borigin3_matrix_9_35_NDE ~ normal(0,10);
borigin3_matrix_34_35 ~ normal(0,10);
borigin3_matrix_9_36 ~ normal(0,10);
borigin3_matrix_9_37 ~ normal(0,10);
borigin3_matrix_9_38 ~ normal(0,10);
borigin3_matrix_40_39 ~ normal(0,10);
borigin3_matrix_39_40 ~ normal(0,10);

borigin4_matrix_8_3 ~ normal(0,10);
borigin4_matrix_3_8 ~ normal(0,10);
borigin4_matrix_9_8 ~ normal(0,10);
borigin4_matrix_32_8 ~ normal(0,10);
borigin4_matrix_33_8 ~ normal(0,10);
borigin4_matrix_8_9 ~ normal(0,10);
borigin4_matrix_34_9 ~ normal(0,10);
borigin4_matrix_35_9 ~ normal(0,10);
borigin4_matrix_36_9 ~ normal(0,10);
borigin4_matrix_37_9 ~ normal(0,10);
borigin4_matrix_38_9 ~ normal(0,10);
borigin4_matrix_39_9 ~ normal(0,10);
borigin4_matrix_40_9 ~ normal(0,10);
borigin4_matrix_8_32_DE ~ normal(0,10);
borigin4_matrix_8_32_NDE ~ normal(0,10);
borigin4_matrix_33_32 ~ normal(0,10);
borigin4_matrix_8_33_DE ~ normal(0,10);
borigin4_matrix_8_33_NDE ~ normal(0,10);
borigin4_matrix_32_33 ~ normal(0,10);
borigin4_matrix_9_34_DE ~ normal(0,10);
borigin4_matrix_9_34_NDE ~ normal(0,10);
borigin4_matrix_35_34 ~ normal(0,10);
borigin4_matrix_9_35_DE ~ normal(0,10);
borigin4_matrix_9_35_NDE ~ normal(0,10);
borigin4_matrix_34_35 ~ normal(0,10);
borigin4_matrix_9_36 ~ normal(0,10);
borigin4_matrix_9_37 ~ normal(0,10);
borigin4_matrix_9_38 ~ normal(0,10);
borigin4_matrix_40_39 ~ normal(0,10);
borigin4_matrix_39_40 ~ normal(0,10);

borigin5_matrix_8_3 ~ normal(0,10);
borigin5_matrix_3_8 ~ normal(0,10);
borigin5_matrix_9_8 ~ normal(0,10);
borigin5_matrix_32_8 ~ normal(0,10);
borigin5_matrix_33_8 ~ normal(0,10);
borigin5_matrix_8_9 ~ normal(0,10);
borigin5_matrix_34_9 ~ normal(0,10);
borigin5_matrix_35_9 ~ normal(0,10);
borigin5_matrix_36_9 ~ normal(0,10);
borigin5_matrix_37_9 ~ normal(0,10);
borigin5_matrix_38_9 ~ normal(0,10);
borigin5_matrix_39_9 ~ normal(0,10);
borigin5_matrix_40_9 ~ normal(0,10);
borigin5_matrix_8_32_DE ~ normal(0,10);
borigin5_matrix_8_32_NDE ~ normal(0,10);
borigin5_matrix_33_32 ~ normal(0,10);
borigin5_matrix_8_33_DE ~ normal(0,10);
borigin5_matrix_8_33_NDE ~ normal(0,10);
borigin5_matrix_32_33 ~ normal(0,10);
borigin5_matrix_9_34_DE ~ normal(0,10);
borigin5_matrix_9_34_NDE ~ normal(0,10);
borigin5_matrix_35_34 ~ normal(0,10);
borigin5_matrix_9_35_DE ~ normal(0,10);
borigin5_matrix_9_35_NDE ~ normal(0,10);
borigin5_matrix_34_35 ~ normal(0,10);
borigin5_matrix_9_36 ~ normal(0,10);
borigin5_matrix_9_37 ~ normal(0,10);
borigin5_matrix_9_38 ~ normal(0,10);
borigin5_matrix_40_39 ~ normal(0,10);
borigin5_matrix_39_40 ~ normal(0,10);

// Detection efficiency GLM
// This is calculated 



// PARALLELIZATION EDITS 
// What we will do is modify the incremental log density statement, to bring target into the first loop by individual.
// This should allow us to increment log density across fish, rather than observations within a fish, allowing us to 
// break up the dataset by fish and therefore run different chunks of the dataset in parallel.

  target += reduce_sum(partial_sum_lupmf, y, grainsize, n_ind, max_visits, cat_X_mat, states_mat, n_obs,
  b0_matrix, borigin1_matrix, borigin2_matrix, borigin3_matrix, borigin4_matrix, borigin5_matrix);

}

