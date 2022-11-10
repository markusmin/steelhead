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
                        
                        // declare a vector to store the run year that each transition occurs in
                        array[] int transition_run_years,
                        
                        // declare the array that says in which run years in which states we calculate DE
                        array[,,] int run_year_DE_array,
                        
                        // declare the matrix that contains the posteriors from the det eff script:
                        matrix det_eff_param_posteriors) { // I don't think we need this either? Since we're just indexing it again with start and end
                          
  // First, declare the total lp (log probability) variable that we will be returning
  real total_lp = 0;
  // Now, loop through the detection matrices for each individual IN THIS SLICE OF THE DATA
  // for (i in 1:slice_n_fish){
    
    // start - end apparently doesn't work because reduce_sum() resets the indices for each slice (confusing) - according to this post: https://discourse.mc-stan.org/t/parallelization-in-cmdstanr-using-reduce-sum-on-a-multinomial-likelihood/24607/7
  // or maybe not, this post seems to contradict this: https://discourse.mc-stan.org/t/help-with-multi-threading-a-simple-ordinal-probit-model-using-reduce-sum/15353
  // i is the number of fish
  for (i in start:end){
    // for (i in 1:end-start+1){
    // Loop through each of the observations, stopping at the loss column
    
    // print("fish # ", i);
    
    // Here, create a vector to store the lp at each observation for a fish
    // vector[n_obs[i]] lp_fish;
    // Let's initialize this instead as a real value starting at zero
    real lp_fish = 0;
    // j is the index of the observation (i.e., each individual state transition)
    for (j in 1:n_obs[i]){
      // for (j in 1:n_obs[i - start + 1]){

        // vector for logits
        vector[43] logits;
        
        // Get the index of the current state
        int current;
        current = states_mat[i,j];
        // current = states_mat[i - start + 1,j];
        
        // here: create a vector with the length of the states, where we will store whether or not they are affected by DE
        vector[42] DE_correction;

        // Populate each of the first 42 (non-loss)
        for (k in 1:42){
          
          // here - add an if else statement. For each transition (combination of current state, and to state),
          // use the run year of the fish to decide whether to pull from the DE or NDE matrix.
          // we need an array that we can index - current (from, rows) x k (to, columns) x run years (slices) - and the values
          // indicate which of the two matrices to use
          
          // what we'll do:
          // create a vector (or array?). Every time we need to select from the DE matrix, note
          // which transition this is. Then, once we're done, we can then index to those and do
          // the correction, where we multiply by a p, and then subtract that from the loss category
          
          // If we are in a run year in a state transition where we are estimating DE, use that matrix;
          // for indexing by run year - some up all n_obs prior to this fish, plus whatever index we're on for the current fish
          // need to make an exception for the first fish, which uses only the index j
          if (i == 1){
                      if (run_year_DE_array[current,k,transition_run_years[j]] == 1){
                      logits[k] = b0_matrix_DE[current, k]+ 
                      // cat_X_mat[i,2] * borigin1_matrix[current,k] + # this is rear, which we are currently not using
                      cat_X_mat[i,3] * borigin1_matrix_DE[current,k] +
                      cat_X_mat[i,4] * borigin2_matrix_DE[current,k] +
                      cat_X_mat[i,5] * borigin3_matrix_DE[current,k] +
                      cat_X_mat[i,6] * borigin4_matrix_DE[current,k] +
                      cat_X_mat[i,7] * borigin5_matrix_DE[current,k];
                      
                      // store in the vector that it's DE
                      DE_correction[k] = 1;
                      
                      // Else, use the NDE matrix
                      } else {
                      logits[k] = b0_matrix_NDE[current, k]+ 
                      // cat_X_mat[i,2] * borigin1_matrix[current,k] + # this is rear, which we are currently not using
                      cat_X_mat[i,3] * borigin1_matrix_NDE[current,k] +
                      cat_X_mat[i,4] * borigin2_matrix_NDE[current,k] +
                      cat_X_mat[i,5] * borigin3_matrix_NDE[current,k] +
                      cat_X_mat[i,6] * borigin4_matrix_NDE[current,k] +
                      cat_X_mat[i,7] * borigin5_matrix_NDE[current,k];
                      
                      // otherwise, store in the vector that it's not DE
                      DE_correction[k] = 0;
                      }
                      
          
          } else {
            // If we're in a state/run year combo that needs DE correction, correct for it
                    if (run_year_DE_array[current,k,transition_run_years[sum(n_obs[1:i-1], j)]] == 1){
                    logits[k] = b0_matrix_DE[current, k]+ 
                    // cat_X_mat[i,2] * borigin1_matrix[current,k] + # this is rear, which we are currently not using
                    cat_X_mat[i,3] * borigin1_matrix_DE[current,k] +
                    cat_X_mat[i,4] * borigin2_matrix_DE[current,k] +
                    cat_X_mat[i,5] * borigin3_matrix_DE[current,k] +
                    cat_X_mat[i,6] * borigin4_matrix_DE[current,k] +
                    cat_X_mat[i,7] * borigin5_matrix_DE[current,k];
                    
                    // store in the vector that it's DE
                    DE_correction[k] = 1;
                    
                    // Else, use the NDE matrix, and don't correct for detection efficiency
                    } else {
                    logits[k] = b0_matrix_NDE[current, k]+ 
                    // cat_X_mat[i,2] * borigin1_matrix[current,k] + # this is rear, which we are currently not using
                    cat_X_mat[i,3] * borigin1_matrix_NDE[current,k] +
                    cat_X_mat[i,4] * borigin2_matrix_NDE[current,k] +
                    cat_X_mat[i,5] * borigin3_matrix_NDE[current,k] +
                    cat_X_mat[i,6] * borigin4_matrix_NDE[current,k] +
                    cat_X_mat[i,7] * borigin5_matrix_NDE[current,k];
                    
                    // otherwise, store in the vector that it's not DE
                    DE_correction[k] = 0;
                    }
                    
          }
          
          
        }
          // loss param
          logits[43] = 0;
                    
          // declare vector to store true movement probabilities
          vector[43] p_vec_actual;
          
          // declare vector to store observed movement probabilities, which are a product of the true probabilities, mulitplied by detection efficiency
          vector[43] p_vec_observed;
          
          // proportion vector, uncorrected for detection efficiency
          p_vec_actual = softmax(logits);
          
          // now, loop through transitions again, and calculate detection efficiency for each where DE_correction == 1
          // declare one vector for linear predictors and one for actual detection efficiency
          vector[42] det_eff_eta;
          vector[42] det_eff;
        for (k in 1:42){
          if (DE_correction[k] == 1) {
            // the exception for if i = 1 for the indexing:
            if (i == 1){
            
            det_eff_eta[k] = tributary_design_matrices_array[transition_run_years[j],,k] * det_eff_param_vector; 
            det_eff[k] = exp(det_eff_eta[k])/(1 + exp(det_eff_eta[k]));
              
            } else {
              
            // to calculate detection efficiency by indexing:
            // The tributary design matrices array will have the same number of slices as states (42, for 43 - loss).
            // Only 14 of these slices (the 14 tributaries that have DE calculations) will have non-zero values, but this will make the indexing simpler.
            // we will index slices using k.
            // Rows will be run years, indexed using transition_run_years[sum(n_obs[1:i-1], j)], same as above
            // the columns will then be the appropriate row of the design matrix, for that state and tributary.
            // That will then be multiplied by the full, 34 length parameter vector. This will be written out in
            // the transformed parameters section, and will have 20 alpha terms and 14 beta terms.
            det_eff_eta[k] = tributary_design_matrices_array[transition_run_years[sum(n_obs[1:i-1], j)],,k] * det_eff_param_vector; 
            det_eff[k] = exp(det_eff_eta[k])/(1 + exp(det_eff_eta[k]));
              
            }
            


            
          } else {
            // If we don't have to calculate a detection efficiency, don't do anything
            
          }
          
          
        }
          
          
          // now loop through transitions and modify p_vec to account for this
        for (k in 1:42){
          if (DE_correction[k] == 1) {
            p_vec_observed[k] = p_vec_actual[k] * det_eff[k];
            
            // each time you modify a term, modify the loss term by the same amount
            p_vec_actual[43] = p_vec_actual[43] + p_vec_actual[k] * (1 - det_eff[k]);
            
          } else {
            p_vec_observed[k] = p_vec[k];
            
          }
          
          
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
    
    
    return total_lp;
}

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
real b0_matrix_12_2;
real b0_matrix_14_2;
real b0_matrix_16_2;
real b0_matrix_18_2;
real b0_matrix_41_2;
real b0_matrix_2_3;
real b0_matrix_8_3;
real b0_matrix_20_3;
real b0_matrix_22_3;
real b0_matrix_3_4;
real b0_matrix_5_4;
real b0_matrix_6_5;
real b0_matrix_24_5;
real b0_matrix_5_6;
real b0_matrix_7_6;
real b0_matrix_26_6;
real b0_matrix_6_7;
real b0_matrix_28_7;
real b0_matrix_30_7;
real b0_matrix_42_7;
real b0_matrix_3_8;
real b0_matrix_9_8;
real b0_matrix_32_8;
real b0_matrix_8_9;
real b0_matrix_34_9;
real b0_matrix_36_9;
real b0_matrix_37_9;
real b0_matrix_38_9;
real b0_matrix_39_9;
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
real borigin1_matrix_8_9;
real borigin1_matrix_34_9;
real borigin1_matrix_36_9;
real borigin1_matrix_37_9;
real borigin1_matrix_38_9;
real borigin1_matrix_39_9;
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
real borigin2_matrix_8_9;
real borigin2_matrix_34_9;
real borigin2_matrix_36_9;
real borigin2_matrix_37_9;
real borigin2_matrix_38_9;
real borigin2_matrix_39_9;
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
real borigin3_matrix_8_9;
real borigin3_matrix_34_9;
real borigin3_matrix_36_9;
real borigin3_matrix_37_9;
real borigin3_matrix_38_9;
real borigin3_matrix_39_9;
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
real borigin4_matrix_8_9;
real borigin4_matrix_34_9;
real borigin4_matrix_36_9;
real borigin4_matrix_37_9;
real borigin4_matrix_38_9;
real borigin4_matrix_39_9;
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
real borigin5_matrix_8_9;
real borigin5_matrix_34_9;
real borigin5_matrix_36_9;
real borigin5_matrix_37_9;
real borigin5_matrix_38_9;
real borigin5_matrix_39_9;
real borigin5_matrix_8_32_DE;
real borigin5_matrix_8_32_NDE;
real borigin5_matrix_9_34_DE;
real borigin5_matrix_9_34_NDE;
real borigin5_matrix_9_36;
real borigin5_matrix_9_37;
real borigin5_matrix_9_38;

// here, write out all of the parameters for detection efficiency
// twenty terms for intercepts for different eras (configurations of antennas) in the different tributaries
real asotin_alpha1;
real asotin_alpha2;
real deschutes_alpha1;
real entiat_alpha1;
real fifteenmile_alpha1;
real hood_alpha1;
real imnaha_alpha1;
real john_day_alpha1;
real methow_alpha1;
real methow_alpha2;
real okanogan_alpha1;
real tucannon_alpha1;
real tucannon_alpha2;
real umatilla_alpha1;
real umatilla_alpha2;
real walla_walla_alpha1;
real walla_walla_alpha2;
real walla_walla_alpha3;
real wenatchee_alpha1;
real yakima_alpha1;

// 14 terms for discharge relationship, one for each tributary
real asotin_beta;
real deschutes_beta;
real entiat_beta;
real fifteenmile_beta;
real hood_beta;
real imnaha_beta;
real john_day_beta;
real methow_beta;
real okanogan_beta;
real tucannon_beta;
real umatilla_beta;
real walla_walla_beta;
real wenatchee_beta;
real yakima_beta;

  
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
  vector[34] det_eff_params;
  
  
  
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
b0_matrix_DE[11,2] = b0_matrix_10_2_DE;
b0_matrix_NDE[11,2] = b0_matrix_10_2_NDE;
b0_matrix_DE[12,2] = b0_matrix_12_2;
b0_matrix_NDE[12,2] = b0_matrix_12_2;
b0_matrix_DE[13,2] = b0_matrix_12_2_DE;
b0_matrix_NDE[13,2] = b0_matrix_12_2_NDE;
b0_matrix_DE[14,2] = b0_matrix_14_2;
b0_matrix_NDE[14,2] = b0_matrix_14_2;
b0_matrix_DE[16,2] = b0_matrix_16_2;
b0_matrix_NDE[16,2] = b0_matrix_16_2;
b0_matrix_DE[17,2] = b0_matrix_16_2_DE;
b0_matrix_NDE[17,2] = b0_matrix_16_2_NDE;
b0_matrix_DE[18,2] = b0_matrix_18_2;
b0_matrix_NDE[18,2] = b0_matrix_18_2;
b0_matrix_DE[19,2] = b0_matrix_18_2_DE;
b0_matrix_NDE[19,2] = b0_matrix_18_2_NDE;
b0_matrix_DE[41,2] = b0_matrix_41_2;
b0_matrix_NDE[41,2] = b0_matrix_41_2;
b0_matrix_DE[2,3] = b0_matrix_2_3;
b0_matrix_NDE[2,3] = b0_matrix_2_3;
b0_matrix_DE[8,3] = b0_matrix_8_3;
b0_matrix_NDE[8,3] = b0_matrix_8_3;
b0_matrix_DE[20,3] = b0_matrix_20_3;
b0_matrix_NDE[20,3] = b0_matrix_20_3;
b0_matrix_DE[21,3] = b0_matrix_20_3_DE;
b0_matrix_NDE[21,3] = b0_matrix_20_3_NDE;
b0_matrix_DE[22,3] = b0_matrix_22_3;
b0_matrix_NDE[22,3] = b0_matrix_22_3;
b0_matrix_DE[23,3] = b0_matrix_22_3_DE;
b0_matrix_NDE[23,3] = b0_matrix_22_3_NDE;
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
b0_matrix_DE[27,6] = b0_matrix_26_6_DE;
b0_matrix_NDE[27,6] = b0_matrix_26_6_NDE;
b0_matrix_DE[6,7] = b0_matrix_6_7;
b0_matrix_NDE[6,7] = b0_matrix_6_7;
b0_matrix_DE[28,7] = b0_matrix_28_7;
b0_matrix_NDE[28,7] = b0_matrix_28_7;
b0_matrix_DE[29,7] = b0_matrix_28_7_DE;
b0_matrix_NDE[29,7] = b0_matrix_28_7_NDE;
b0_matrix_DE[30,7] = b0_matrix_30_7;
b0_matrix_NDE[30,7] = b0_matrix_30_7;
b0_matrix_DE[31,7] = b0_matrix_30_7_DE;
b0_matrix_NDE[31,7] = b0_matrix_30_7_NDE;
b0_matrix_DE[42,7] = b0_matrix_42_7;
b0_matrix_NDE[42,7] = b0_matrix_42_7;
b0_matrix_DE[3,8] = b0_matrix_3_8;
b0_matrix_NDE[3,8] = b0_matrix_3_8;
b0_matrix_DE[9,8] = b0_matrix_9_8;
b0_matrix_NDE[9,8] = b0_matrix_9_8;
b0_matrix_DE[32,8] = b0_matrix_32_8;
b0_matrix_NDE[32,8] = b0_matrix_32_8;
b0_matrix_DE[33,8] = b0_matrix_32_8_DE;
b0_matrix_NDE[33,8] = b0_matrix_32_8_NDE;
b0_matrix_DE[8,9] = b0_matrix_8_9;
b0_matrix_NDE[8,9] = b0_matrix_8_9;
b0_matrix_DE[34,9] = b0_matrix_34_9;
b0_matrix_NDE[34,9] = b0_matrix_34_9;
b0_matrix_DE[35,9] = b0_matrix_34_9_DE;
b0_matrix_NDE[35,9] = b0_matrix_34_9_NDE;
b0_matrix_DE[36,9] = b0_matrix_36_9;
b0_matrix_NDE[36,9] = b0_matrix_36_9;
b0_matrix_DE[37,9] = b0_matrix_37_9;
b0_matrix_NDE[37,9] = b0_matrix_37_9;
b0_matrix_DE[38,9] = b0_matrix_38_9;
b0_matrix_NDE[38,9] = b0_matrix_38_9;
b0_matrix_DE[39,9] = b0_matrix_39_9;
b0_matrix_NDE[39,9] = b0_matrix_39_9;
b0_matrix_DE[40,9] = b0_matrix_39_9_DE;
b0_matrix_NDE[40,9] = b0_matrix_39_9_NDE;
b0_matrix_DE[2,10] = b0_matrix_2_10_DE;
b0_matrix_NDE[2,10] = b0_matrix_2_10_NDE;
b0_matrix_NDE[2,11] = b0_matrix_11_10_NDE;
b0_matrix_DE[2,12] = b0_matrix_2_12_DE;
b0_matrix_NDE[2,12] = b0_matrix_2_12_NDE;
b0_matrix_NDE[2,13] = b0_matrix_13_12_NDE;
b0_matrix_DE[2,14] = b0_matrix_2_14_DE;
b0_matrix_NDE[2,14] = b0_matrix_2_14_NDE;
b0_matrix_NDE[2,15] = b0_matrix_2_14_NDE;
b0_matrix_DE[2,16] = b0_matrix_2_16_DE;
b0_matrix_NDE[2,16] = b0_matrix_2_16_NDE;
b0_matrix_NDE[2,17] = b0_matrix_17_16_NDE;
b0_matrix_DE[2,18] = b0_matrix_2_18_DE;
b0_matrix_NDE[2,18] = b0_matrix_2_18_NDE;
b0_matrix_NDE[2,19] = b0_matrix_19_18_NDE;
b0_matrix_DE[3,20] = b0_matrix_3_20_DE;
b0_matrix_NDE[3,20] = b0_matrix_3_20_NDE;
b0_matrix_NDE[3,21] = b0_matrix_21_20_NDE;
b0_matrix_DE[3,22] = b0_matrix_3_22_DE;
b0_matrix_NDE[3,22] = b0_matrix_3_22_NDE;
b0_matrix_NDE[3,23] = b0_matrix_23_22_NDE;
b0_matrix_DE[5,24] = b0_matrix_5_24_DE;
b0_matrix_NDE[5,24] = b0_matrix_5_24_NDE;
b0_matrix_NDE[5,25] = b0_matrix_5_24_NDE;
b0_matrix_DE[6,26] = b0_matrix_6_26_DE;
b0_matrix_NDE[6,26] = b0_matrix_6_26_NDE;
b0_matrix_NDE[6,27] = b0_matrix_27_26_NDE;
b0_matrix_DE[7,28] = b0_matrix_7_28_DE;
b0_matrix_NDE[7,28] = b0_matrix_7_28_NDE;
b0_matrix_NDE[7,29] = b0_matrix_29_28_NDE;
b0_matrix_DE[7,30] = b0_matrix_7_30_DE;
b0_matrix_NDE[7,30] = b0_matrix_7_30_NDE;
b0_matrix_NDE[7,31] = b0_matrix_31_30_NDE;
b0_matrix_DE[8,32] = b0_matrix_8_32_DE;
b0_matrix_NDE[8,32] = b0_matrix_8_32_NDE;
b0_matrix_NDE[8,33] = b0_matrix_33_32_NDE;
b0_matrix_DE[9,34] = b0_matrix_9_34_DE;
b0_matrix_NDE[9,34] = b0_matrix_9_34_NDE;
b0_matrix_NDE[9,35] = b0_matrix_35_34_NDE;
b0_matrix_DE[9,36] = b0_matrix_9_36;
b0_matrix_NDE[9,36] = b0_matrix_9_36;
b0_matrix_DE[9,37] = b0_matrix_9_37;
b0_matrix_NDE[9,37] = b0_matrix_9_37;
b0_matrix_DE[9,38] = b0_matrix_9_38;
b0_matrix_NDE[9,38] = b0_matrix_9_38;
b0_matrix_DE[2,41] = b0_matrix_2_41;
b0_matrix_NDE[2,41] = b0_matrix_2_41;
b0_matrix_DE[7,42] = b0_matrix_7_42;
b0_matrix_NDE[7,42] = b0_matrix_7_42;

borigin1_matrix_DE[8,3] = borigin1_matrix_8_3_DE;
borigin1_matrix_NDE[8,3] = borigin1_matrix_8_3_NDE;
borigin1_matrix_DE[3,8] = borigin1_matrix_3_8_DE;
borigin1_matrix_NDE[3,8] = borigin1_matrix_3_8_NDE;
borigin1_matrix_DE[9,8] = borigin1_matrix_9_8_DE;
borigin1_matrix_NDE[9,8] = borigin1_matrix_9_8_NDE;
borigin1_matrix_DE[32,8] = borigin1_matrix_32_8_DE;
borigin1_matrix_NDE[32,8] = borigin1_matrix_32_8_NDE;
borigin1_matrix_DE[33,8] = borigin1_matrix_32_8_DE;
borigin1_matrix_NDE[33,8] = borigin1_matrix_32_8_NDE;
borigin1_matrix_DE[8,9] = borigin1_matrix_8_9_DE;
borigin1_matrix_NDE[8,9] = borigin1_matrix_8_9_NDE;
borigin1_matrix_DE[34,9] = borigin1_matrix_34_9_DE;
borigin1_matrix_NDE[34,9] = borigin1_matrix_34_9_NDE;
borigin1_matrix_DE[35,9] = borigin1_matrix_34_9_DE;
borigin1_matrix_NDE[35,9] = borigin1_matrix_34_9_NDE;
borigin1_matrix_DE[36,9] = borigin1_matrix_36_9_DE;
borigin1_matrix_NDE[36,9] = borigin1_matrix_36_9_NDE;
borigin1_matrix_DE[37,9] = borigin1_matrix_37_9_DE;
borigin1_matrix_NDE[37,9] = borigin1_matrix_37_9_NDE;
borigin1_matrix_DE[38,9] = borigin1_matrix_38_9_DE;
borigin1_matrix_NDE[38,9] = borigin1_matrix_38_9_NDE;
borigin1_matrix_DE[39,9] = borigin1_matrix_39_9_DE;
borigin1_matrix_NDE[39,9] = borigin1_matrix_39_9_NDE;
borigin1_matrix_DE[40,9] = borigin1_matrix_39_9_DE;
borigin1_matrix_NDE[40,9] = borigin1_matrix_39_9_NDE;
borigin1_matrix_DE[8,32] = borigin1_matrix_8_32_DE;
borigin1_matrix_NDE[8,32] = borigin1_matrix_8_32_NDE;
borigin1_matrix_NDE[8,33] = borigin1_matrix_33_32_NDE;
borigin1_matrix_DE[9,34] = borigin1_matrix_9_34_DE;
borigin1_matrix_NDE[9,34] = borigin1_matrix_9_34_NDE;
borigin1_matrix_NDE[9,35] = borigin1_matrix_35_34_NDE;
borigin1_matrix_DE[9,36] = borigin1_matrix_9_36_DE;
borigin1_matrix_NDE[9,36] = borigin1_matrix_9_36_NDE;
borigin1_matrix_DE[9,37] = borigin1_matrix_9_37_DE;
borigin1_matrix_NDE[9,37] = borigin1_matrix_9_37_NDE;
borigin1_matrix_DE[9,38] = borigin1_matrix_9_38_DE;
borigin1_matrix_NDE[9,38] = borigin1_matrix_9_38_NDE;

borigin2_matrix_DE[8,3] = borigin2_matrix_8_3_DE;
borigin2_matrix_NDE[8,3] = borigin2_matrix_8_3_NDE;
borigin2_matrix_DE[3,8] = borigin2_matrix_3_8_DE;
borigin2_matrix_NDE[3,8] = borigin2_matrix_3_8_NDE;
borigin2_matrix_DE[9,8] = borigin2_matrix_9_8_DE;
borigin2_matrix_NDE[9,8] = borigin2_matrix_9_8_NDE;
borigin2_matrix_DE[32,8] = borigin2_matrix_32_8_DE;
borigin2_matrix_NDE[32,8] = borigin2_matrix_32_8_NDE;
borigin2_matrix_DE[33,8] = borigin2_matrix_32_8_DE;
borigin2_matrix_NDE[33,8] = borigin2_matrix_32_8_NDE;
borigin2_matrix_DE[8,9] = borigin2_matrix_8_9_DE;
borigin2_matrix_NDE[8,9] = borigin2_matrix_8_9_NDE;
borigin2_matrix_DE[34,9] = borigin2_matrix_34_9_DE;
borigin2_matrix_NDE[34,9] = borigin2_matrix_34_9_NDE;
borigin2_matrix_DE[35,9] = borigin2_matrix_34_9_DE;
borigin2_matrix_NDE[35,9] = borigin2_matrix_34_9_NDE;
borigin2_matrix_DE[36,9] = borigin2_matrix_36_9_DE;
borigin2_matrix_NDE[36,9] = borigin2_matrix_36_9_NDE;
borigin2_matrix_DE[37,9] = borigin2_matrix_37_9_DE;
borigin2_matrix_NDE[37,9] = borigin2_matrix_37_9_NDE;
borigin2_matrix_DE[38,9] = borigin2_matrix_38_9_DE;
borigin2_matrix_NDE[38,9] = borigin2_matrix_38_9_NDE;
borigin2_matrix_DE[39,9] = borigin2_matrix_39_9_DE;
borigin2_matrix_NDE[39,9] = borigin2_matrix_39_9_NDE;
borigin2_matrix_DE[40,9] = borigin2_matrix_39_9_DE;
borigin2_matrix_NDE[40,9] = borigin2_matrix_39_9_NDE;
borigin2_matrix_DE[8,32] = borigin2_matrix_8_32_DE;
borigin2_matrix_NDE[8,32] = borigin2_matrix_8_32_NDE;
borigin2_matrix_NDE[8,33] = borigin2_matrix_33_32_NDE;
borigin2_matrix_DE[9,34] = borigin2_matrix_9_34_DE;
borigin2_matrix_NDE[9,34] = borigin2_matrix_9_34_NDE;
borigin2_matrix_NDE[9,35] = borigin2_matrix_35_34_NDE;
borigin2_matrix_DE[9,36] = borigin2_matrix_9_36_DE;
borigin2_matrix_NDE[9,36] = borigin2_matrix_9_36_NDE;
borigin2_matrix_DE[9,37] = borigin2_matrix_9_37_DE;
borigin2_matrix_NDE[9,37] = borigin2_matrix_9_37_NDE;
borigin2_matrix_DE[9,38] = borigin2_matrix_9_38_DE;
borigin2_matrix_NDE[9,38] = borigin2_matrix_9_38_NDE;

borigin3_matrix_DE[8,3] = borigin3_matrix_8_3_DE;
borigin3_matrix_NDE[8,3] = borigin3_matrix_8_3_NDE;
borigin3_matrix_DE[3,8] = borigin3_matrix_3_8_DE;
borigin3_matrix_NDE[3,8] = borigin3_matrix_3_8_NDE;
borigin3_matrix_DE[9,8] = borigin3_matrix_9_8_DE;
borigin3_matrix_NDE[9,8] = borigin3_matrix_9_8_NDE;
borigin3_matrix_DE[32,8] = borigin3_matrix_32_8_DE;
borigin3_matrix_NDE[32,8] = borigin3_matrix_32_8_NDE;
borigin3_matrix_DE[33,8] = borigin3_matrix_32_8_DE;
borigin3_matrix_NDE[33,8] = borigin3_matrix_32_8_NDE;
borigin3_matrix_DE[8,9] = borigin3_matrix_8_9_DE;
borigin3_matrix_NDE[8,9] = borigin3_matrix_8_9_NDE;
borigin3_matrix_DE[34,9] = borigin3_matrix_34_9_DE;
borigin3_matrix_NDE[34,9] = borigin3_matrix_34_9_NDE;
borigin3_matrix_DE[35,9] = borigin3_matrix_34_9_DE;
borigin3_matrix_NDE[35,9] = borigin3_matrix_34_9_NDE;
borigin3_matrix_DE[36,9] = borigin3_matrix_36_9_DE;
borigin3_matrix_NDE[36,9] = borigin3_matrix_36_9_NDE;
borigin3_matrix_DE[37,9] = borigin3_matrix_37_9_DE;
borigin3_matrix_NDE[37,9] = borigin3_matrix_37_9_NDE;
borigin3_matrix_DE[38,9] = borigin3_matrix_38_9_DE;
borigin3_matrix_NDE[38,9] = borigin3_matrix_38_9_NDE;
borigin3_matrix_DE[39,9] = borigin3_matrix_39_9_DE;
borigin3_matrix_NDE[39,9] = borigin3_matrix_39_9_NDE;
borigin3_matrix_DE[40,9] = borigin3_matrix_39_9_DE;
borigin3_matrix_NDE[40,9] = borigin3_matrix_39_9_NDE;
borigin3_matrix_DE[8,32] = borigin3_matrix_8_32_DE;
borigin3_matrix_NDE[8,32] = borigin3_matrix_8_32_NDE;
borigin3_matrix_NDE[8,33] = borigin3_matrix_33_32_NDE;
borigin3_matrix_DE[9,34] = borigin3_matrix_9_34_DE;
borigin3_matrix_NDE[9,34] = borigin3_matrix_9_34_NDE;
borigin3_matrix_NDE[9,35] = borigin3_matrix_35_34_NDE;
borigin3_matrix_DE[9,36] = borigin3_matrix_9_36_DE;
borigin3_matrix_NDE[9,36] = borigin3_matrix_9_36_NDE;
borigin3_matrix_DE[9,37] = borigin3_matrix_9_37_DE;
borigin3_matrix_NDE[9,37] = borigin3_matrix_9_37_NDE;
borigin3_matrix_DE[9,38] = borigin3_matrix_9_38_DE;
borigin3_matrix_NDE[9,38] = borigin3_matrix_9_38_NDE;

borigin4_matrix_DE[8,3] = borigin4_matrix_8_3_DE;
borigin4_matrix_NDE[8,3] = borigin4_matrix_8_3_NDE;
borigin4_matrix_DE[3,8] = borigin4_matrix_3_8_DE;
borigin4_matrix_NDE[3,8] = borigin4_matrix_3_8_NDE;
borigin4_matrix_DE[9,8] = borigin4_matrix_9_8_DE;
borigin4_matrix_NDE[9,8] = borigin4_matrix_9_8_NDE;
borigin4_matrix_DE[32,8] = borigin4_matrix_32_8_DE;
borigin4_matrix_NDE[32,8] = borigin4_matrix_32_8_NDE;
borigin4_matrix_DE[33,8] = borigin4_matrix_32_8_DE;
borigin4_matrix_NDE[33,8] = borigin4_matrix_32_8_NDE;
borigin4_matrix_DE[8,9] = borigin4_matrix_8_9_DE;
borigin4_matrix_NDE[8,9] = borigin4_matrix_8_9_NDE;
borigin4_matrix_DE[34,9] = borigin4_matrix_34_9_DE;
borigin4_matrix_NDE[34,9] = borigin4_matrix_34_9_NDE;
borigin4_matrix_DE[35,9] = borigin4_matrix_34_9_DE;
borigin4_matrix_NDE[35,9] = borigin4_matrix_34_9_NDE;
borigin4_matrix_DE[36,9] = borigin4_matrix_36_9_DE;
borigin4_matrix_NDE[36,9] = borigin4_matrix_36_9_NDE;
borigin4_matrix_DE[37,9] = borigin4_matrix_37_9_DE;
borigin4_matrix_NDE[37,9] = borigin4_matrix_37_9_NDE;
borigin4_matrix_DE[38,9] = borigin4_matrix_38_9_DE;
borigin4_matrix_NDE[38,9] = borigin4_matrix_38_9_NDE;
borigin4_matrix_DE[39,9] = borigin4_matrix_39_9_DE;
borigin4_matrix_NDE[39,9] = borigin4_matrix_39_9_NDE;
borigin4_matrix_DE[40,9] = borigin4_matrix_39_9_DE;
borigin4_matrix_NDE[40,9] = borigin4_matrix_39_9_NDE;
borigin4_matrix_DE[8,32] = borigin4_matrix_8_32_DE;
borigin4_matrix_NDE[8,32] = borigin4_matrix_8_32_NDE;
borigin4_matrix_NDE[8,33] = borigin4_matrix_33_32_NDE;
borigin4_matrix_DE[9,34] = borigin4_matrix_9_34_DE;
borigin4_matrix_NDE[9,34] = borigin4_matrix_9_34_NDE;
borigin4_matrix_NDE[9,35] = borigin4_matrix_35_34_NDE;
borigin4_matrix_DE[9,36] = borigin4_matrix_9_36_DE;
borigin4_matrix_NDE[9,36] = borigin4_matrix_9_36_NDE;
borigin4_matrix_DE[9,37] = borigin4_matrix_9_37_DE;
borigin4_matrix_NDE[9,37] = borigin4_matrix_9_37_NDE;
borigin4_matrix_DE[9,38] = borigin4_matrix_9_38_DE;
borigin4_matrix_NDE[9,38] = borigin4_matrix_9_38_NDE;

borigin5_matrix_DE[8,3] = borigin5_matrix_8_3_DE;
borigin5_matrix_NDE[8,3] = borigin5_matrix_8_3_NDE;
borigin5_matrix_DE[3,8] = borigin5_matrix_3_8_DE;
borigin5_matrix_NDE[3,8] = borigin5_matrix_3_8_NDE;
borigin5_matrix_DE[9,8] = borigin5_matrix_9_8_DE;
borigin5_matrix_NDE[9,8] = borigin5_matrix_9_8_NDE;
borigin5_matrix_DE[32,8] = borigin5_matrix_32_8_DE;
borigin5_matrix_NDE[32,8] = borigin5_matrix_32_8_NDE;
borigin5_matrix_DE[33,8] = borigin5_matrix_32_8_DE;
borigin5_matrix_NDE[33,8] = borigin5_matrix_32_8_NDE;
borigin5_matrix_DE[8,9] = borigin5_matrix_8_9_DE;
borigin5_matrix_NDE[8,9] = borigin5_matrix_8_9_NDE;
borigin5_matrix_DE[34,9] = borigin5_matrix_34_9_DE;
borigin5_matrix_NDE[34,9] = borigin5_matrix_34_9_NDE;
borigin5_matrix_DE[35,9] = borigin5_matrix_34_9_DE;
borigin5_matrix_NDE[35,9] = borigin5_matrix_34_9_NDE;
borigin5_matrix_DE[36,9] = borigin5_matrix_36_9_DE;
borigin5_matrix_NDE[36,9] = borigin5_matrix_36_9_NDE;
borigin5_matrix_DE[37,9] = borigin5_matrix_37_9_DE;
borigin5_matrix_NDE[37,9] = borigin5_matrix_37_9_NDE;
borigin5_matrix_DE[38,9] = borigin5_matrix_38_9_DE;
borigin5_matrix_NDE[38,9] = borigin5_matrix_38_9_NDE;
borigin5_matrix_DE[39,9] = borigin5_matrix_39_9_DE;
borigin5_matrix_NDE[39,9] = borigin5_matrix_39_9_NDE;
borigin5_matrix_DE[40,9] = borigin5_matrix_39_9_DE;
borigin5_matrix_NDE[40,9] = borigin5_matrix_39_9_NDE;
borigin5_matrix_DE[8,32] = borigin5_matrix_8_32_DE;
borigin5_matrix_NDE[8,32] = borigin5_matrix_8_32_NDE;
borigin5_matrix_NDE[8,33] = borigin5_matrix_33_32_NDE;
borigin5_matrix_DE[9,34] = borigin5_matrix_9_34_DE;
borigin5_matrix_NDE[9,34] = borigin5_matrix_9_34_NDE;
borigin5_matrix_NDE[9,35] = borigin5_matrix_35_34_NDE;
borigin5_matrix_DE[9,36] = borigin5_matrix_9_36_DE;
borigin5_matrix_NDE[9,36] = borigin5_matrix_9_36_NDE;
borigin5_matrix_DE[9,37] = borigin5_matrix_9_37_DE;
borigin5_matrix_NDE[9,37] = borigin5_matrix_9_37_NDE;
borigin5_matrix_DE[9,38] = borigin5_matrix_9_38_DE;
borigin5_matrix_NDE[9,38] = borigin5_matrix_9_38_NDE;

// detection efficiency - create a vector that stores all parameters
vector[43] det_eff_param_vector;

// populate the vector
det_eff_param_vector[1] = asotin_alpha1;
det_eff_param_vector[2] = asotin_alpha2;
det_eff_param_vector[3] = deschutes_alpha1;
det_eff_param_vector[4] = entiat_alpha1;
det_eff_param_vector[5] = fifteenmile_alpha1;
det_eff_param_vector[6] = hood_alpha1;
det_eff_param_vector[7] = imnaha_alpha1;
det_eff_param_vector[8] = john_day_alpha1;
det_eff_param_vector[9] = methow_alpha1;
det_eff_param_vector[10] = methow_alpha2;
det_eff_param_vector[11] = okanogan_alpha1;
det_eff_param_vector[12] = tucannon_alpha1;
det_eff_param_vector[13] = tucannon_alpha2;
det_eff_param_vector[14] = umatilla_alpha1;
det_eff_param_vector[15] = umatilla_alpha2;
det_eff_param_vector[16] = walla_walla_alpha1;
det_eff_param_vector[17] = walla_walla_alpha2;
det_eff_param_vector[18] = walla_walla_alpha3;
det_eff_param_vector[19] = wenatchee_alpha1;
det_eff_param_vector[20] = yakima_alpha1;

// 14 terms for discharge relationship, one for each tributary
det_eff_param_vector[21] = asotin_beta;
det_eff_param_vector[22] = deschutes_beta;
det_eff_param_vector[23] = entiat_beta;
det_eff_param_vector[24] = fifteenmile_beta;
det_eff_param_vector[25] = hood_beta;
det_eff_param_vector[26] = imnaha_beta;
det_eff_param_vector[27] = john_day_beta;
det_eff_param_vector[28] = methow_beta;
det_eff_param_vector[29] = okanogan_beta;
det_eff_param_vector[30] = tucannon_beta;
det_eff_param_vector[31] = umatilla_beta;
det_eff_param_vector[32] = walla_walla_beta;
det_eff_param_vector[33] = wenatchee_beta;
det_eff_param_vector[34] = yakima_beta;



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
b0_matrix_12_2 ~ normal(0,10);
b0_matrix_14_2 ~ normal(0,10);
b0_matrix_16_2 ~ normal(0,10);
b0_matrix_18_2 ~ normal(0,10);
b0_matrix_41_2 ~ normal(0,10);
b0_matrix_2_3 ~ normal(0,10);
b0_matrix_8_3 ~ normal(0,10);
b0_matrix_20_3 ~ normal(0,10);
b0_matrix_22_3 ~ normal(0,10);
b0_matrix_3_4 ~ normal(0,10);
b0_matrix_5_4 ~ normal(0,10);
b0_matrix_6_5 ~ normal(0,10);
b0_matrix_24_5 ~ normal(0,10);
b0_matrix_5_6 ~ normal(0,10);
b0_matrix_7_6 ~ normal(0,10);
b0_matrix_26_6 ~ normal(0,10);
b0_matrix_6_7 ~ normal(0,10);
b0_matrix_28_7 ~ normal(0,10);
b0_matrix_30_7 ~ normal(0,10);
b0_matrix_42_7 ~ normal(0,10);
b0_matrix_3_8 ~ normal(0,10);
b0_matrix_9_8 ~ normal(0,10);
b0_matrix_32_8 ~ normal(0,10);
b0_matrix_8_9 ~ normal(0,10);
b0_matrix_34_9 ~ normal(0,10);
b0_matrix_36_9 ~ normal(0,10);
b0_matrix_37_9 ~ normal(0,10);
b0_matrix_38_9 ~ normal(0,10);
b0_matrix_39_9 ~ normal(0,10);
b0_matrix_2_10_DE ~ normal(0,10);
b0_matrix_2_10_NDE ~ normal(0,10);
b0_matrix_2_12_DE ~ normal(0,10);
b0_matrix_2_12_NDE ~ normal(0,10);
b0_matrix_2_14_DE ~ normal(0,10);
b0_matrix_2_14_NDE ~ normal(0,10);
b0_matrix_2_16_DE ~ normal(0,10);
b0_matrix_2_16_NDE ~ normal(0,10);
b0_matrix_2_18_DE ~ normal(0,10);
b0_matrix_2_18_NDE ~ normal(0,10);
b0_matrix_3_20_DE ~ normal(0,10);
b0_matrix_3_20_NDE ~ normal(0,10);
b0_matrix_3_22_DE ~ normal(0,10);
b0_matrix_3_22_NDE ~ normal(0,10);
b0_matrix_5_24_DE ~ normal(0,10);
b0_matrix_5_24_NDE ~ normal(0,10);
b0_matrix_6_26_DE ~ normal(0,10);
b0_matrix_6_26_NDE ~ normal(0,10);
b0_matrix_7_28_DE ~ normal(0,10);
b0_matrix_7_28_NDE ~ normal(0,10);
b0_matrix_7_30_DE ~ normal(0,10);
b0_matrix_7_30_NDE ~ normal(0,10);
b0_matrix_8_32_DE ~ normal(0,10);
b0_matrix_8_32_NDE ~ normal(0,10);
b0_matrix_9_34_DE ~ normal(0,10);
b0_matrix_9_34_NDE ~ normal(0,10);
b0_matrix_9_36 ~ normal(0,10);
b0_matrix_9_37 ~ normal(0,10);
b0_matrix_9_38 ~ normal(0,10);
b0_matrix_2_41 ~ normal(0,10);
b0_matrix_7_42 ~ normal(0,10);

borigin1_matrix_8_3 ~ normal(0,10);
borigin1_matrix_3_8 ~ normal(0,10);
borigin1_matrix_9_8 ~ normal(0,10);
borigin1_matrix_32_8 ~ normal(0,10);
borigin1_matrix_8_9 ~ normal(0,10);
borigin1_matrix_34_9 ~ normal(0,10);
borigin1_matrix_36_9 ~ normal(0,10);
borigin1_matrix_37_9 ~ normal(0,10);
borigin1_matrix_38_9 ~ normal(0,10);
borigin1_matrix_39_9 ~ normal(0,10);
borigin1_matrix_8_32_DE ~ normal(0,10);
borigin1_matrix_8_32_NDE ~ normal(0,10);
borigin1_matrix_9_34_DE ~ normal(0,10);
borigin1_matrix_9_34_NDE ~ normal(0,10);
borigin1_matrix_9_36 ~ normal(0,10);
borigin1_matrix_9_37 ~ normal(0,10);
borigin1_matrix_9_38 ~ normal(0,10);

borigin2_matrix_8_3 ~ normal(0,10);
borigin2_matrix_3_8 ~ normal(0,10);
borigin2_matrix_9_8 ~ normal(0,10);
borigin2_matrix_32_8 ~ normal(0,10);
borigin2_matrix_8_9 ~ normal(0,10);
borigin2_matrix_34_9 ~ normal(0,10);
borigin2_matrix_36_9 ~ normal(0,10);
borigin2_matrix_37_9 ~ normal(0,10);
borigin2_matrix_38_9 ~ normal(0,10);
borigin2_matrix_39_9 ~ normal(0,10);
borigin2_matrix_8_32_DE ~ normal(0,10);
borigin2_matrix_8_32_NDE ~ normal(0,10);
borigin2_matrix_9_34_DE ~ normal(0,10);
borigin2_matrix_9_34_NDE ~ normal(0,10);
borigin2_matrix_9_36 ~ normal(0,10);
borigin2_matrix_9_37 ~ normal(0,10);
borigin2_matrix_9_38 ~ normal(0,10);

borigin3_matrix_8_3 ~ normal(0,10);
borigin3_matrix_3_8 ~ normal(0,10);
borigin3_matrix_9_8 ~ normal(0,10);
borigin3_matrix_32_8 ~ normal(0,10);
borigin3_matrix_8_9 ~ normal(0,10);
borigin3_matrix_34_9 ~ normal(0,10);
borigin3_matrix_36_9 ~ normal(0,10);
borigin3_matrix_37_9 ~ normal(0,10);
borigin3_matrix_38_9 ~ normal(0,10);
borigin3_matrix_39_9 ~ normal(0,10);
borigin3_matrix_8_32_DE ~ normal(0,10);
borigin3_matrix_8_32_NDE ~ normal(0,10);
borigin3_matrix_9_34_DE ~ normal(0,10);
borigin3_matrix_9_34_NDE ~ normal(0,10);
borigin3_matrix_9_36 ~ normal(0,10);
borigin3_matrix_9_37 ~ normal(0,10);
borigin3_matrix_9_38 ~ normal(0,10);

borigin4_matrix_8_3 ~ normal(0,10);
borigin4_matrix_3_8 ~ normal(0,10);
borigin4_matrix_9_8 ~ normal(0,10);
borigin4_matrix_32_8 ~ normal(0,10);
borigin4_matrix_8_9 ~ normal(0,10);
borigin4_matrix_34_9 ~ normal(0,10);
borigin4_matrix_36_9 ~ normal(0,10);
borigin4_matrix_37_9 ~ normal(0,10);
borigin4_matrix_38_9 ~ normal(0,10);
borigin4_matrix_39_9 ~ normal(0,10);
borigin4_matrix_8_32_DE ~ normal(0,10);
borigin4_matrix_8_32_NDE ~ normal(0,10);
borigin4_matrix_9_34_DE ~ normal(0,10);
borigin4_matrix_9_34_NDE ~ normal(0,10);
borigin4_matrix_9_36 ~ normal(0,10);
borigin4_matrix_9_37 ~ normal(0,10);
borigin4_matrix_9_38 ~ normal(0,10);

borigin5_matrix_8_3 ~ normal(0,10);
borigin5_matrix_3_8 ~ normal(0,10);
borigin5_matrix_9_8 ~ normal(0,10);
borigin5_matrix_32_8 ~ normal(0,10);
borigin5_matrix_8_9 ~ normal(0,10);
borigin5_matrix_34_9 ~ normal(0,10);
borigin5_matrix_36_9 ~ normal(0,10);
borigin5_matrix_37_9 ~ normal(0,10);
borigin5_matrix_38_9 ~ normal(0,10);
borigin5_matrix_39_9 ~ normal(0,10);
borigin5_matrix_8_32_DE ~ normal(0,10);
borigin5_matrix_8_32_NDE ~ normal(0,10);
borigin5_matrix_9_34_DE ~ normal(0,10);
borigin5_matrix_9_34_NDE ~ normal(0,10);
borigin5_matrix_9_36 ~ normal(0,10);
borigin5_matrix_9_37 ~ normal(0,10);
borigin5_matrix_9_38 ~ normal(0,10);

// Prior on detection efficiency parameters - from the other stan script for detection efficiency
// here, write out all of the parameters for detection efficiency
// twenty terms for intercepts for different eras (configurations of antennas) in the different tributaries
// the outputs from that model will be the first column [,1] containing central tendency (mean)
// and the second column [,2] containing the standard deviation
asotin_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
asotin_alpha2 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
deschutes_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
entiat_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
fifteenmile_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
hood_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
imnaha_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
john_day_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
methow_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
methow_alpha2 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
okanogan_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
tucannon_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
tucannon_alpha2 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
umatilla_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
umatilla_alpha2 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
walla_walla_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
walla_walla_alpha2 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
walla_walla_alpha3 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
wenatchee_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
yakima_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);

// 14 terms for discharge relationship, one for each tributary
asotin_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
deschutes_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
entiat_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
fifteenmile_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
hood_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
imnaha_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
john_day_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
methow_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
okanogan_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
tucannon_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
umatilla_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
walla_walla_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
wenatchee_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
yakima_beta ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);



// PARALLELIZATION EDITS 
// What we will do is modify the incremental log density statement, to bring target into the first loop by individual.
// This should allow us to increment log density across fish, rather than observations within a fish, allowing us to 
// break up the dataset by fish and therefore run different chunks of the dataset in parallel.

  target += reduce_sum(partial_sum_lupmf, y, grainsize, n_ind, max_visits, cat_X_mat, states_mat, n_obs,
  b0_matrix, borigin1_matrix, borigin2_matrix, borigin3_matrix, borigin4_matrix, borigin5_matrix);

}

