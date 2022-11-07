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
                        array[,] real b0_matrix,
                        array[,] real borigin1_matrix,
                        array[,] real borigin2_matrix,
                        array[,] real borigin3_matrix,
                        array[,] real borigin4_matrix,
                        array[,] real borigin5_matrix,
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
real b0_matrix_16_2;
real b0_matrix_17_2;
real b0_matrix_18_2;
real b0_matrix_19_2;
real b0_matrix_41_2;
real b0_matrix_2_3;
real b0_matrix_4_3;
real b0_matrix_8_3;
real b0_matrix_20_3;
real b0_matrix_21_3;
real b0_matrix_22_3;
real b0_matrix_23_3;
real b0_matrix_3_4;
real b0_matrix_5_4;
real b0_matrix_4_5;
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
real b0_matrix_11_10;
real b0_matrix_10_11;
real b0_matrix_13_12;
real b0_matrix_12_13;
real b0_matrix_17_16;
real b0_matrix_16_17;
real b0_matrix_19_18;
real b0_matrix_18_19;
real b0_matrix_21_20;
real b0_matrix_20_21;
real b0_matrix_23_22;
real b0_matrix_22_23;
real b0_matrix_24_25;
real b0_matrix_27_26;
real b0_matrix_26_27;
real b0_matrix_29_28;
real b0_matrix_28_29;
real b0_matrix_31_30;
real b0_matrix_30_31;
real b0_matrix_33_32;
real b0_matrix_32_33;
real b0_matrix_35_34;
real b0_matrix_34_35;
real b0_matrix_9_36;
real b0_matrix_9_37;
real b0_matrix_9_38;
real b0_matrix_40_39;
real b0_matrix_39_40;
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
real borigin1_matrix_33_32;
real borigin1_matrix_32_33;
real borigin1_matrix_35_34;
real borigin1_matrix_34_35;
real borigin1_matrix_9_36;
real borigin1_matrix_9_37;
real borigin1_matrix_9_38;
real borigin1_matrix_40_39;
real borigin1_matrix_39_40;

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
real borigin2_matrix_33_32;
real borigin2_matrix_32_33;
real borigin2_matrix_35_34;
real borigin2_matrix_34_35;
real borigin2_matrix_9_36;
real borigin2_matrix_9_37;
real borigin2_matrix_9_38;
real borigin2_matrix_40_39;
real borigin2_matrix_39_40;

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
real borigin3_matrix_33_32;
real borigin3_matrix_32_33;
real borigin3_matrix_35_34;
real borigin3_matrix_34_35;
real borigin3_matrix_9_36;
real borigin3_matrix_9_37;
real borigin3_matrix_9_38;
real borigin3_matrix_40_39;
real borigin3_matrix_39_40;

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
real borigin4_matrix_33_32;
real borigin4_matrix_32_33;
real borigin4_matrix_35_34;
real borigin4_matrix_34_35;
real borigin4_matrix_9_36;
real borigin4_matrix_9_37;
real borigin4_matrix_9_38;
real borigin4_matrix_40_39;
real borigin4_matrix_39_40;

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
real borigin5_matrix_33_32;
real borigin5_matrix_32_33;
real borigin5_matrix_35_34;
real borigin5_matrix_34_35;
real borigin5_matrix_9_36;
real borigin5_matrix_9_37;
real borigin5_matrix_9_38;
real borigin5_matrix_40_39;
real borigin5_matrix_39_40;
  
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
b0_matrix[2,1] = b0_matrix_2_1;
b0_matrix[1,2] = b0_matrix_1_2;
b0_matrix[3,2] = b0_matrix_3_2;
b0_matrix[10,2] = b0_matrix_10_2;
b0_matrix[11,2] = b0_matrix_11_2;
b0_matrix[12,2] = b0_matrix_12_2;
b0_matrix[13,2] = b0_matrix_13_2;
b0_matrix[16,2] = b0_matrix_16_2;
b0_matrix[17,2] = b0_matrix_17_2;
b0_matrix[18,2] = b0_matrix_18_2;
b0_matrix[19,2] = b0_matrix_19_2;
b0_matrix[41,2] = b0_matrix_41_2;
b0_matrix[2,3] = b0_matrix_2_3;
b0_matrix[4,3] = b0_matrix_4_3;
b0_matrix[8,3] = b0_matrix_8_3;
b0_matrix[20,3] = b0_matrix_20_3;
b0_matrix[21,3] = b0_matrix_21_3;
b0_matrix[22,3] = b0_matrix_22_3;
b0_matrix[23,3] = b0_matrix_23_3;
b0_matrix[3,4] = b0_matrix_3_4;
b0_matrix[5,4] = b0_matrix_5_4;
b0_matrix[4,5] = b0_matrix_4_5;
b0_matrix[6,5] = b0_matrix_6_5;
b0_matrix[24,5] = b0_matrix_24_5;
b0_matrix[5,6] = b0_matrix_5_6;
b0_matrix[7,6] = b0_matrix_7_6;
b0_matrix[26,6] = b0_matrix_26_6;
b0_matrix[27,6] = b0_matrix_27_6;
b0_matrix[6,7] = b0_matrix_6_7;
b0_matrix[28,7] = b0_matrix_28_7;
b0_matrix[29,7] = b0_matrix_29_7;
b0_matrix[30,7] = b0_matrix_30_7;
b0_matrix[31,7] = b0_matrix_31_7;
b0_matrix[42,7] = b0_matrix_42_7;
b0_matrix[3,8] = b0_matrix_3_8;
b0_matrix[9,8] = b0_matrix_9_8;
b0_matrix[32,8] = b0_matrix_32_8;
b0_matrix[33,8] = b0_matrix_33_8;
b0_matrix[8,9] = b0_matrix_8_9;
b0_matrix[34,9] = b0_matrix_34_9;
b0_matrix[35,9] = b0_matrix_35_9;
b0_matrix[36,9] = b0_matrix_36_9;
b0_matrix[37,9] = b0_matrix_37_9;
b0_matrix[38,9] = b0_matrix_38_9;
b0_matrix[39,9] = b0_matrix_39_9;
b0_matrix[40,9] = b0_matrix_40_9;
b0_matrix[11,10] = b0_matrix_11_10;
b0_matrix[10,11] = b0_matrix_10_11;
b0_matrix[13,12] = b0_matrix_13_12;
b0_matrix[12,13] = b0_matrix_12_13;
b0_matrix[17,16] = b0_matrix_17_16;
b0_matrix[16,17] = b0_matrix_16_17;
b0_matrix[19,18] = b0_matrix_19_18;
b0_matrix[18,19] = b0_matrix_18_19;
b0_matrix[21,20] = b0_matrix_21_20;
b0_matrix[20,21] = b0_matrix_20_21;
b0_matrix[23,22] = b0_matrix_23_22;
b0_matrix[22,23] = b0_matrix_22_23;
b0_matrix[24,25] = b0_matrix_24_25;
b0_matrix[27,26] = b0_matrix_27_26;
b0_matrix[26,27] = b0_matrix_26_27;
b0_matrix[29,28] = b0_matrix_29_28;
b0_matrix[28,29] = b0_matrix_28_29;
b0_matrix[31,30] = b0_matrix_31_30;
b0_matrix[30,31] = b0_matrix_30_31;
b0_matrix[33,32] = b0_matrix_33_32;
b0_matrix[32,33] = b0_matrix_32_33;
b0_matrix[35,34] = b0_matrix_35_34;
b0_matrix[34,35] = b0_matrix_34_35;
b0_matrix[9,36] = b0_matrix_9_36;
b0_matrix[9,37] = b0_matrix_9_37;
b0_matrix[9,38] = b0_matrix_9_38;
b0_matrix[40,39] = b0_matrix_40_39;
b0_matrix[39,40] = b0_matrix_39_40;
b0_matrix[2,41] = b0_matrix_2_41;
b0_matrix[7,42] = b0_matrix_7_42;

borigin1_matrix[8,3] = borigin1_matrix_8_3;
borigin1_matrix[3,8] = borigin1_matrix_3_8;
borigin1_matrix[9,8] = borigin1_matrix_9_8;
borigin1_matrix[32,8] = borigin1_matrix_32_8;
borigin1_matrix[33,8] = borigin1_matrix_33_8;
borigin1_matrix[8,9] = borigin1_matrix_8_9;
borigin1_matrix[34,9] = borigin1_matrix_34_9;
borigin1_matrix[35,9] = borigin1_matrix_35_9;
borigin1_matrix[36,9] = borigin1_matrix_36_9;
borigin1_matrix[37,9] = borigin1_matrix_37_9;
borigin1_matrix[38,9] = borigin1_matrix_38_9;
borigin1_matrix[39,9] = borigin1_matrix_39_9;
borigin1_matrix[40,9] = borigin1_matrix_40_9;
borigin1_matrix[33,32] = borigin1_matrix_33_32;
borigin1_matrix[32,33] = borigin1_matrix_32_33;
borigin1_matrix[35,34] = borigin1_matrix_35_34;
borigin1_matrix[34,35] = borigin1_matrix_34_35;
borigin1_matrix[9,36] = borigin1_matrix_9_36;
borigin1_matrix[9,37] = borigin1_matrix_9_37;
borigin1_matrix[9,38] = borigin1_matrix_9_38;
borigin1_matrix[40,39] = borigin1_matrix_40_39;
borigin1_matrix[39,40] = borigin1_matrix_39_40;

borigin2_matrix[8,3] = borigin2_matrix_8_3;
borigin2_matrix[3,8] = borigin2_matrix_3_8;
borigin2_matrix[9,8] = borigin2_matrix_9_8;
borigin2_matrix[32,8] = borigin2_matrix_32_8;
borigin2_matrix[33,8] = borigin2_matrix_33_8;
borigin2_matrix[8,9] = borigin2_matrix_8_9;
borigin2_matrix[34,9] = borigin2_matrix_34_9;
borigin2_matrix[35,9] = borigin2_matrix_35_9;
borigin2_matrix[36,9] = borigin2_matrix_36_9;
borigin2_matrix[37,9] = borigin2_matrix_37_9;
borigin2_matrix[38,9] = borigin2_matrix_38_9;
borigin2_matrix[39,9] = borigin2_matrix_39_9;
borigin2_matrix[40,9] = borigin2_matrix_40_9;
borigin2_matrix[33,32] = borigin2_matrix_33_32;
borigin2_matrix[32,33] = borigin2_matrix_32_33;
borigin2_matrix[35,34] = borigin2_matrix_35_34;
borigin2_matrix[34,35] = borigin2_matrix_34_35;
borigin2_matrix[9,36] = borigin2_matrix_9_36;
borigin2_matrix[9,37] = borigin2_matrix_9_37;
borigin2_matrix[9,38] = borigin2_matrix_9_38;
borigin2_matrix[40,39] = borigin2_matrix_40_39;
borigin2_matrix[39,40] = borigin2_matrix_39_40;

borigin3_matrix[8,3] = borigin3_matrix_8_3;
borigin3_matrix[3,8] = borigin3_matrix_3_8;
borigin3_matrix[9,8] = borigin3_matrix_9_8;
borigin3_matrix[32,8] = borigin3_matrix_32_8;
borigin3_matrix[33,8] = borigin3_matrix_33_8;
borigin3_matrix[8,9] = borigin3_matrix_8_9;
borigin3_matrix[34,9] = borigin3_matrix_34_9;
borigin3_matrix[35,9] = borigin3_matrix_35_9;
borigin3_matrix[36,9] = borigin3_matrix_36_9;
borigin3_matrix[37,9] = borigin3_matrix_37_9;
borigin3_matrix[38,9] = borigin3_matrix_38_9;
borigin3_matrix[39,9] = borigin3_matrix_39_9;
borigin3_matrix[40,9] = borigin3_matrix_40_9;
borigin3_matrix[33,32] = borigin3_matrix_33_32;
borigin3_matrix[32,33] = borigin3_matrix_32_33;
borigin3_matrix[35,34] = borigin3_matrix_35_34;
borigin3_matrix[34,35] = borigin3_matrix_34_35;
borigin3_matrix[9,36] = borigin3_matrix_9_36;
borigin3_matrix[9,37] = borigin3_matrix_9_37;
borigin3_matrix[9,38] = borigin3_matrix_9_38;
borigin3_matrix[40,39] = borigin3_matrix_40_39;
borigin3_matrix[39,40] = borigin3_matrix_39_40;

borigin4_matrix[8,3] = borigin4_matrix_8_3;
borigin4_matrix[3,8] = borigin4_matrix_3_8;
borigin4_matrix[9,8] = borigin4_matrix_9_8;
borigin4_matrix[32,8] = borigin4_matrix_32_8;
borigin4_matrix[33,8] = borigin4_matrix_33_8;
borigin4_matrix[8,9] = borigin4_matrix_8_9;
borigin4_matrix[34,9] = borigin4_matrix_34_9;
borigin4_matrix[35,9] = borigin4_matrix_35_9;
borigin4_matrix[36,9] = borigin4_matrix_36_9;
borigin4_matrix[37,9] = borigin4_matrix_37_9;
borigin4_matrix[38,9] = borigin4_matrix_38_9;
borigin4_matrix[39,9] = borigin4_matrix_39_9;
borigin4_matrix[40,9] = borigin4_matrix_40_9;
borigin4_matrix[33,32] = borigin4_matrix_33_32;
borigin4_matrix[32,33] = borigin4_matrix_32_33;
borigin4_matrix[35,34] = borigin4_matrix_35_34;
borigin4_matrix[34,35] = borigin4_matrix_34_35;
borigin4_matrix[9,36] = borigin4_matrix_9_36;
borigin4_matrix[9,37] = borigin4_matrix_9_37;
borigin4_matrix[9,38] = borigin4_matrix_9_38;
borigin4_matrix[40,39] = borigin4_matrix_40_39;
borigin4_matrix[39,40] = borigin4_matrix_39_40;

borigin5_matrix[8,3] = borigin5_matrix_8_3;
borigin5_matrix[3,8] = borigin5_matrix_3_8;
borigin5_matrix[9,8] = borigin5_matrix_9_8;
borigin5_matrix[32,8] = borigin5_matrix_32_8;
borigin5_matrix[33,8] = borigin5_matrix_33_8;
borigin5_matrix[8,9] = borigin5_matrix_8_9;
borigin5_matrix[34,9] = borigin5_matrix_34_9;
borigin5_matrix[35,9] = borigin5_matrix_35_9;
borigin5_matrix[36,9] = borigin5_matrix_36_9;
borigin5_matrix[37,9] = borigin5_matrix_37_9;
borigin5_matrix[38,9] = borigin5_matrix_38_9;
borigin5_matrix[39,9] = borigin5_matrix_39_9;
borigin5_matrix[40,9] = borigin5_matrix_40_9;
borigin5_matrix[33,32] = borigin5_matrix_33_32;
borigin5_matrix[32,33] = borigin5_matrix_32_33;
borigin5_matrix[35,34] = borigin5_matrix_35_34;
borigin5_matrix[34,35] = borigin5_matrix_34_35;
borigin5_matrix[9,36] = borigin5_matrix_9_36;
borigin5_matrix[9,37] = borigin5_matrix_9_37;
borigin5_matrix[9,38] = borigin5_matrix_9_38;
borigin5_matrix[40,39] = borigin5_matrix_40_39;
borigin5_matrix[39,40] = borigin5_matrix_39_40;

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
b0_matrix_16_2 ~ normal(0,10);
b0_matrix_17_2 ~ normal(0,10);
b0_matrix_18_2 ~ normal(0,10);
b0_matrix_19_2 ~ normal(0,10);
b0_matrix_41_2 ~ normal(0,10);
b0_matrix_2_3 ~ normal(0,10);
b0_matrix_4_3 ~ normal(0,10);
b0_matrix_8_3 ~ normal(0,10);
b0_matrix_20_3 ~ normal(0,10);
b0_matrix_21_3 ~ normal(0,10);
b0_matrix_22_3 ~ normal(0,10);
b0_matrix_23_3 ~ normal(0,10);
b0_matrix_3_4 ~ normal(0,10);
b0_matrix_5_4 ~ normal(0,10);
b0_matrix_4_5 ~ normal(0,10);
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
b0_matrix_11_10 ~ normal(0,10);
b0_matrix_10_11 ~ normal(0,10);
b0_matrix_13_12 ~ normal(0,10);
b0_matrix_12_13 ~ normal(0,10);
b0_matrix_17_16 ~ normal(0,10);
b0_matrix_16_17 ~ normal(0,10);
b0_matrix_19_18 ~ normal(0,10);
b0_matrix_18_19 ~ normal(0,10);
b0_matrix_21_20 ~ normal(0,10);
b0_matrix_20_21 ~ normal(0,10);
b0_matrix_23_22 ~ normal(0,10);
b0_matrix_22_23 ~ normal(0,10);
b0_matrix_24_25 ~ normal(0,10);
b0_matrix_27_26 ~ normal(0,10);
b0_matrix_26_27 ~ normal(0,10);
b0_matrix_29_28 ~ normal(0,10);
b0_matrix_28_29 ~ normal(0,10);
b0_matrix_31_30 ~ normal(0,10);
b0_matrix_30_31 ~ normal(0,10);
b0_matrix_33_32 ~ normal(0,10);
b0_matrix_32_33 ~ normal(0,10);
b0_matrix_35_34 ~ normal(0,10);
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
borigin1_matrix_33_32 ~ normal(0,10);
borigin1_matrix_32_33 ~ normal(0,10);
borigin1_matrix_35_34 ~ normal(0,10);
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
borigin2_matrix_33_32 ~ normal(0,10);
borigin2_matrix_32_33 ~ normal(0,10);
borigin2_matrix_35_34 ~ normal(0,10);
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
borigin3_matrix_33_32 ~ normal(0,10);
borigin3_matrix_32_33 ~ normal(0,10);
borigin3_matrix_35_34 ~ normal(0,10);
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
borigin4_matrix_33_32 ~ normal(0,10);
borigin4_matrix_32_33 ~ normal(0,10);
borigin4_matrix_35_34 ~ normal(0,10);
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
borigin5_matrix_33_32 ~ normal(0,10);
borigin5_matrix_32_33 ~ normal(0,10);
borigin5_matrix_35_34 ~ normal(0,10);
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

