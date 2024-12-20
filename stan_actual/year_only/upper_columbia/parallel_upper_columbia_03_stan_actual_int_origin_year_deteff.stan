// 03_parallel_snake_stan_actual_int_origin_wild_deteff
// wild, upper Columbia Steelhead
// This model includes a detection efficiency calculation for the tributaries
// this model has as covariate temperature and origin and a random effect of year

functions{
    real partial_sum_lpmf( // int[] slice_n_fish, // I don't think that we need this, given that we are using start and end to index the detection histories
                        array[,] int slice_y,
                        int start, int end,
                        // I think we need to re-declare things for our function that are also already declared in our data?
                        int n_ind,
                        int max_visits,
                        int nmovements,
                        
                        // declare an array with the date of transition for fish
                        // this will have the same dimensions/structure as states_mat
                        // rows = individual fish, columns = transitions
                        // dates will not be formatted as actual dates, but instead as integers
                        // because this is necessary to use this to index covariates
                        array[,] int transition_dates,
                        
                        // transition seasons (to index temperature parameters)
                        array[] int transition_seasons_vector,

                        
                        // categorical covariates (in this case, just origin parameters)
                        // array[n_ind, 7] int cat_X_mat,
                        array[,] int cat_X_mat,
                        
                        
                        
                        // array[n_ind, max_visits-1] int states_mat,
                        array[,] int states_mat,
                        // array[n_ind,max_visits] int y, // I don't think that we need this for the same reason as above
                        // array[,] int y,
                        array[] int n_obs,
                        
                        // For detection efficiency calculation purposes: We need two of each of these one for transitions
                        // where we can calculate detection efficiency, one for where we can't. Most parameters
                        // will be shared between the two, but transitions into tributaries won't.
                        // DE = detection efficiency, NDE = no detection efficiency
                        array[,] real b0_matrix_DE,
                        array[,] real borigin1_matrix_DE,
                        array[,] real borigin2_matrix_DE,
                        
                        array[,] real b0_matrix_NDE,
                        array[,] real borigin1_matrix_NDE,
                        array[,] real borigin2_matrix_NDE,
                        
                        
                        // indices for sigma_year
                        array[,] int sigma_year_indices,
                        // arrays to store sigma_year_matrices (DE and NDE)
                        array[,] real sigma_year_matrix_DE,
                        array[,] real sigma_year_matrix_NDE,
                        
                        
                        // arrays to store byear parameters
                        // first, create a temporary df, with each of the individual byear parameters (drawn from normal, using sigma_year)
                        // create an empty array with 3 dimensions: matrix for parameters, slices for years
                        // matt trick, transform these outside
                        array [,,] real byear_raw_parameters_array_DE,
                        array [,,] real byear_raw_parameters_array_NDE,
                        // these are also arguments, since they're transformed parameters
                        array [,,] real byear_actual_parameters_array_DE,
                        array [,,] real byear_actual_parameters_array_NDE,
                        
                        
                        // below is new data for detection efficiency
                        // declare an array that you can put detection probabilities in, 
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
                        
                        // declare a vector to store the run year that each transition occurs in
                        array[] int transition_run_years, // or does this need to be array[] int?
                        
                        // declare how many run years there are
                        int nyears,
                        
                        // declare the array that says in which run years in which states we calculate DE
                        array[,,] int run_year_DE_array,
                        
                        // declare the vector that contains the parameters for detection efficiency
                        vector det_eff_param_vector
                        ) { // I don't think we need this either? Since we're just indexing it again with start and end
                          
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
          // print("i = ",i);
            // print("j = ",j);

        // vector for logits
        vector[43] logits;
        
        // Get the index of the current state
        int current;
        current = states_mat[i,j];
        // current = states_mat[i - start + 1,j];
        
        // Get the current date
        int date;
        date = transition_dates[i,j];
        
        // here: create a vector with the length of the states, where we will store whether or not they are affected by DE
        vector[42] DE_correction;
        
            // Get the current run year
            
            // indexing is different if i == 1
              int run_year_index;      
              int run_year_actual;
              // index for transition season as well
              int season_actual;
              
              
            if (i == 1) {
              // now, sum the vector to get a single integer value for indexing the run year vector
              run_year_index = j;
              // get the actual run year
              run_year_actual = transition_run_years[run_year_index];
              // index for transition season as well
              season_actual = transition_seasons_vector[run_year_index];
              
            } else {
                // indexing tweaks
                // first: declare a vector that stores all of the indices that we need to sum to get the right run year
                array[i] int run_year_indices_vector;
                // then populate: first elements are all the number of transitions from all previous fish
                run_year_indices_vector[1:(i-1)] = n_obs[1:(i-1)];
                // last element is the transition number of the current fish  
                run_year_indices_vector[i] = j;
                // now, sum the vector to get a single integer value for indexing the run year vector
                run_year_index = sum(run_year_indices_vector);
                // get the actual run year
                run_year_actual = transition_run_years[run_year_index];
                
                // index for transition season as well
              season_actual = transition_seasons_vector[run_year_index];

            }
          

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

            // If we're in a state/run year combo that needs DE correction, correct for it
                    // if (run_year_DE_array[current,k,transition_run_years[sum(n_obs[1:i-1], j)]] == 1){
                      // need to declare another intermediate value: an integer that takes either 0 or 1, for whether it's DE or NDE
                      int DE_index;
                      // now, determine if it's 0 or 1
                      DE_index = run_year_DE_array[current,k,run_year_actual];
                      // if (run_year_DE_array[current,k,transition_run_years[run_year_indices_vector]] == 1){
                  if (DE_index == 1){
                    logits[k] = b0_matrix_DE[current, k]+ 
                    cat_X_mat[i,2] * borigin1_matrix_DE[current,k] +
                    cat_X_mat[i,3] * borigin2_matrix_DE[current,k] +
                    // add year as a random effect
                    byear_actual_parameters_array_DE[current, k, run_year_actual];
                    
                    // store in the vector that it's DE
                    DE_correction[k] = 1;
                    
                    // Else, use the NDE matrix, and don't correct for detection efficiency
                  } else {
                    logits[k] = b0_matrix_NDE[current, k]+
                    cat_X_mat[i,2] * borigin1_matrix_NDE[current,k] +
                    cat_X_mat[i,3] * borigin2_matrix_NDE[current,k] + 
                    // add year as a random effect
                    byear_actual_parameters_array_NDE[current, k, run_year_actual];
                    
                    
                    // otherwise, store in the vector that it's not DE
                    DE_correction[k] = 0;
                    }
                    
                    
          
          
        }
          // loss param
          logits[43] = 0;
                    
          // declare vector to store true movement probabilities
          vector[43] p_vec_actual;
          
          // declare vector to store observed movement probabilities, which are a product of the true probabilities, mulitplied by detection efficiency
          vector[43] p_vec_observed;
          
          // print("i!=1; logits = ",logits);
          
          // print("logits = ", logits);
          
          // proportion vector, uncorrected for detection efficiency
          p_vec_actual = softmax(logits);
          
          // print("p_vec_actual: ", p_vec_actual);
        // print("i = ",i);
        // print("j = ",j);
          
          // now, loop through transitions again, and calculate detection efficiency for each where DE_correction == 1
          // declare one vector for linear predictors and one for actual detection efficiency
          
          vector[42] det_eff_eta;
          vector[42] det_eff;
        
        for (k in 1:42){
          if (DE_correction[k] == 1) {

            // to calculate detection efficiency by indexing:
            // The tributary design matrices array will have the same number of slices as states (42, for 43 - loss).
            // Only 14 of these slices (the 14 tributaries that have DE calculations) will have non-zero values, but this will make the indexing simpler.
            // we will index slices using k.
            // Rows will be run years, indexed using transition_run_years[sum(n_obs[1:i-1], j)], same as above
            // the columns will then be the appropriate row of the design matrix, for that state and tributary.
            // That will then be multiplied by the full, 34 length parameter vector. This will be written out in
            // the transformed parameters section, and will have 20 alpha terms and 14 beta terms.

            // vector[34] trib_design_matrix_result;
            // trib_design_matrix_result = to_row_vector(tributary_design_matrices_array[run_year_actual,,k]) * det_eff_param_vector;
            // det_eff_eta[k] = sum(trib_design_matrix_result);
            det_eff_eta[k] = to_row_vector(tributary_design_matrices_array[run_year_actual,,k]) * to_vector(det_eff_param_vector);
            det_eff[k] = exp(det_eff_eta[k])/(1 + exp(det_eff_eta[k]));


          } else {
            // If we don't have to calculate a detection efficiency, don't do anything

          }


        }
          
          
          // now loop through transitions and modify p_vec to account for this
        
          
          // // Create a vector to store all of the loss corrections for detection efficiency
          vector[42] loss_term_DE_corrections;
          
          // for testing: create a vector of p_vec_observed_test that we can fill
          // vector[42] p_vec_observed_test;
       
        for (k in 1:42){
          if (DE_correction[k] == 1) {
            p_vec_observed[k] = p_vec_actual[k] * det_eff[k];
            // p_vec_observed[k] = p_vec_actual[k]; // again this is just for testing - remove det eff correction for now

            // each time you modify a term, modify the loss term by the same amount
            // p_vec_observed[43] = p_vec_observed[43] + p_vec_actual[k] * (1 - det_eff[k]);
            loss_term_DE_corrections[k] = p_vec_actual[k] * (1 - det_eff[k]);

          } else {
            p_vec_observed[k] = p_vec_actual[k];

            // Just put a zero there, otherwise it'll be NAs
            loss_term_DE_corrections[k] = 0;

          }


        }
        
        // Once you've looped through all states, correct loss term
        p_vec_observed[43] = p_vec_actual[43] + sum(loss_term_DE_corrections);
        
        // p_vec_observed[43] = p_vec_actual[43]; // this line is just for testing - not actually the right loss term, but I just need to make sure that the
        // above line isn't the cause of the errors
        // print("i = ",i);
        // print("j = ",j);
        // print("current state: ", current);
        // print("run year actual = ", run_year_actual);
        // print("DE_correction: ", DE_correction);
        // print("det eff eta: ", det_eff_eta);
        // print("det eff: ", det_eff);
        // print("b0_matrix_DE[current,]: ", b0_matrix_DE[current,]);
        // print("b0_matrix_NDE[current,]: ", b0_matrix_NDE[current,]);
        // print("logits: ", logits);
        // print("p_vec_actual: ", p_vec_actual);
        // print("p_vec_observed: ", p_vec_observed);
        // print("sum p_vec_observed = ", sum(p_vec_observed));
        
        
        // 
        // print("slice_y[i - start + 1,j+1]: ", slice_y[i - start + 1,j+1]);
        // print("next state: ", slice_y[i - start + 1,j+1], "; prob of next state = ", p_vec_observed[slice_y[i - start + 1,j+1]]);
        // Store the log probability of that individual transition in the vector that we declared
        // lp_fish[j] = categorical_lpmf(slice_y[i,j+1] | p_vec);
        // Changed the indexing here, based on:https://discourse.mc-stan.org/t/help-with-multi-threading-a-simple-ordinal-probit-model-using-reduce-sum/15353
        
        // inspect the data vs. the vector of probabilities
        // print("data: ", slice_y[i - start + 1,j+1]);
        
        
        // 2022-11-29 edits for performance:
        // We are going to try to vectorize this and sum it outside of the loop, rather than repeatedly updating a variable.
        
        // print - what the logits are 
        // print("logits = ", logits);
        // print what the byear params are
        // print("byear_actual_parameters_array_NDE = ", byear_actual_parameters_array_NDE);
        // print("byear_actual_parameters_array_DE = ", byear_actual_parameters_array_DE);
        
        
        lp_fish += categorical_lpmf(slice_y[i - start + 1,j+1] | p_vec_observed);
      
    }
    
    // Now, increment the log density
    total_lp += lp_fish;
    
    
} // end looping through fish in this slice

return total_lp;

} // end partial_sum_lpmf function

} // end functions block

data {
  // array[10, 41, 1200] int y; // array with dimensions 10 states (rows), 48 columns (maximum number of site visits), 1200 matrix slices (number of fish); has to be int for multinomial logit
  int max_visits;
  int n_ind; // number of individuals (this is nfish)
  array[n_ind,max_visits] int y; // array with dimensions nfish (number of fish), 41 (maximum number of site visits)
  array[n_ind] int n_obs; // n_obs is a vector that contains the number of site visits for each individual - but needs to be declared as array to take integer values
  vector[43] possible_movements; // a vector containing the number of possible transitions out of each state
  array[n_ind, max_visits] int transition_dates; // a matrix of the dates (as numeric string with 2005-06-01 as 1) that each transition took place
  array[n_ind, max_visits-1] int states_mat; // a matrix (array to take integer values) with rows = each fish and columns = the site visits
  // array[54, 2] int movements; // a matrix that contains the possible transitions out of each state
  int nmovements; // an integer value containing the number of possible transitions (same as rows in movements data)
  array[nmovements, 2] int movements; // a matrix that contains the possible transitions out of each state
  // array[758,2] int not_movements; // a matrix containing all non-allowed state transitions
  int n_notmovements; // an integer value containing the number of non- possible transitions (same as rows in not_movements data)
  array[n_notmovements,2] int not_movements; // a matrix containing all non-allowed state transitions
  // array[n_ind, 48] int dates; // a matrix containing dates (as integers) where rows = number of fish and columns = site visits
  matrix[43, 43] possible_states; // the transition matrix with 1s for allowable transitions and 0s for non-allowable
  array[n_ind, 3] int cat_X_mat; // a matrix with rows = individual fish and columns = various categorical covariates (origin)
  
  matrix[35,2] det_eff_param_posteriors; // declare the matrix that contains the posteriors from the det eff script:
  
  // array that contains design matrices for each tributary
  array[18,35,42] int tributary_design_matrices_array;
  
  // an array that matches the sigma_year values to the correct states
  array[nmovements, 2] int sigma_year_indices;
  
  // vector that contains the run years in which each individual transition occurred
  int ntransitions;
  array[ntransitions] int transition_run_years;
  
  // a vector of the seasons that each transition took place (0 = winter/spring, before June 1; 1 = summer/fall, June 1 of later)
  array[ntransitions] int transition_seasons_vector; 
  
  int nyears;
  
  // array that contains which run years in which transitions need DE correction
  array[43,43,18] int run_year_DE_array;
  
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
real b0_matrix_1_2;
real b0_matrix_2_1;
real b0_matrix_2_3;
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
real b0_matrix_2_41;
real b0_matrix_3_2;
real b0_matrix_3_4;
real b0_matrix_3_8;
real b0_matrix_3_20_DE;
real b0_matrix_3_20_NDE;
real b0_matrix_3_22_DE;
real b0_matrix_3_22_NDE;
real b0_matrix_4_3;
real b0_matrix_4_5;
real b0_matrix_5_4;
real b0_matrix_5_6;
real b0_matrix_5_24_DE;
real b0_matrix_5_24_NDE;
real b0_matrix_6_5;
real b0_matrix_6_7;
real b0_matrix_6_26_DE;
real b0_matrix_6_26_NDE;
real b0_matrix_7_6;
real b0_matrix_7_28_DE;
real b0_matrix_7_28_NDE;
real b0_matrix_7_30_DE;
real b0_matrix_7_30_NDE;
real b0_matrix_7_42;
real b0_matrix_8_3;
real b0_matrix_8_9;
real b0_matrix_8_32_DE;
real b0_matrix_8_32_NDE;
real b0_matrix_9_8;
real b0_matrix_9_34_DE;
real b0_matrix_9_34_NDE;
real b0_matrix_9_36;
real b0_matrix_9_37;
real b0_matrix_9_38;
real b0_matrix_9_39_DE;
real b0_matrix_9_39_NDE;
real b0_matrix_10_2;
real b0_matrix_12_2;
real b0_matrix_14_2;
real b0_matrix_16_2;
real b0_matrix_18_2;
real b0_matrix_20_3;
real b0_matrix_22_3;
real b0_matrix_24_5;
real b0_matrix_26_6;
real b0_matrix_28_7;
real b0_matrix_30_7;
real b0_matrix_32_8;
real b0_matrix_34_9;
real b0_matrix_36_9;
real b0_matrix_37_9;
real b0_matrix_38_9;
real b0_matrix_39_9;
real b0_matrix_41_2;
real b0_matrix_42_7;

// random effect of year parameters
// array [nyears] real byear;
// real<lower=0> sigma_year;
real<lower=0>sigma_year_matrix_1_2;
real<lower=0>sigma_year_matrix_2_1;
real<lower=0>sigma_year_matrix_2_3;
real<lower=0>sigma_year_matrix_2_10_DE;
real<lower=0>sigma_year_matrix_2_10_NDE;
real<lower=0>sigma_year_matrix_2_12_DE;
real<lower=0>sigma_year_matrix_2_12_NDE;
real<lower=0>sigma_year_matrix_2_14_DE;
real<lower=0>sigma_year_matrix_2_14_NDE;
real<lower=0>sigma_year_matrix_2_16_DE;
real<lower=0>sigma_year_matrix_2_16_NDE;
real<lower=0>sigma_year_matrix_2_18_DE;
real<lower=0>sigma_year_matrix_2_18_NDE;
real<lower=0>sigma_year_matrix_2_41;
real<lower=0>sigma_year_matrix_3_2;
real<lower=0>sigma_year_matrix_3_4;
real<lower=0>sigma_year_matrix_3_8;
real<lower=0>sigma_year_matrix_3_20_DE;
real<lower=0>sigma_year_matrix_3_20_NDE;
real<lower=0>sigma_year_matrix_3_22_DE;
real<lower=0>sigma_year_matrix_3_22_NDE;
real<lower=0>sigma_year_matrix_4_3;
real<lower=0>sigma_year_matrix_4_5;
real<lower=0>sigma_year_matrix_5_4;
real<lower=0>sigma_year_matrix_5_6;
real<lower=0>sigma_year_matrix_5_24_DE;
real<lower=0>sigma_year_matrix_5_24_NDE;
real<lower=0>sigma_year_matrix_6_5;
real<lower=0>sigma_year_matrix_6_7;
real<lower=0>sigma_year_matrix_6_26_DE;
real<lower=0>sigma_year_matrix_6_26_NDE;
real<lower=0>sigma_year_matrix_7_6;
real<lower=0>sigma_year_matrix_7_28_DE;
real<lower=0>sigma_year_matrix_7_28_NDE;
real<lower=0>sigma_year_matrix_7_30_DE;
real<lower=0>sigma_year_matrix_7_30_NDE;
real<lower=0>sigma_year_matrix_7_42;
real<lower=0>sigma_year_matrix_8_3;
real<lower=0>sigma_year_matrix_8_9;
real<lower=0>sigma_year_matrix_8_32_DE;
real<lower=0>sigma_year_matrix_8_32_NDE;
real<lower=0>sigma_year_matrix_9_8;
real<lower=0>sigma_year_matrix_9_34_DE;
real<lower=0>sigma_year_matrix_9_34_NDE;
real<lower=0>sigma_year_matrix_9_36;
real<lower=0>sigma_year_matrix_9_37;
real<lower=0>sigma_year_matrix_9_38;
real<lower=0>sigma_year_matrix_9_39_DE;
real<lower=0>sigma_year_matrix_9_39_NDE;
real<lower=0>sigma_year_matrix_10_2;
real<lower=0>sigma_year_matrix_14_2;
real<lower=0>sigma_year_matrix_24_5;
real<lower=0>sigma_year_matrix_26_6;
real<lower=0>sigma_year_matrix_28_7;
real<lower=0>sigma_year_matrix_30_7;
real<lower=0>sigma_year_matrix_36_9;
real<lower=0>sigma_year_matrix_42_7;

// vectors to store byear parameters for each movement
vector[nyears] byear_raw_vector_1_2;
vector[nyears] byear_raw_vector_2_1;
vector[nyears] byear_raw_vector_2_3;
vector[nyears] byear_raw_vector_2_10_DE;
vector[nyears] byear_raw_vector_2_10_NDE;
vector[nyears] byear_raw_vector_2_12_DE;
vector[nyears] byear_raw_vector_2_12_NDE;
vector[nyears] byear_raw_vector_2_14_DE;
vector[nyears] byear_raw_vector_2_14_NDE;
vector[nyears] byear_raw_vector_2_16_DE;
vector[nyears] byear_raw_vector_2_16_NDE;
vector[nyears] byear_raw_vector_2_18_DE;
vector[nyears] byear_raw_vector_2_18_NDE;
vector[nyears] byear_raw_vector_2_41;
vector[nyears] byear_raw_vector_3_2;
vector[nyears] byear_raw_vector_3_4;
vector[nyears] byear_raw_vector_3_8;
vector[nyears] byear_raw_vector_3_20_DE;
vector[nyears] byear_raw_vector_3_20_NDE;
vector[nyears] byear_raw_vector_3_22_DE;
vector[nyears] byear_raw_vector_3_22_NDE;
vector[nyears] byear_raw_vector_4_3;
vector[nyears] byear_raw_vector_4_5;
vector[nyears] byear_raw_vector_5_4;
vector[nyears] byear_raw_vector_5_6;
vector[nyears] byear_raw_vector_5_24_DE;
vector[nyears] byear_raw_vector_5_24_NDE;
vector[nyears] byear_raw_vector_6_5;
vector[nyears] byear_raw_vector_6_7;
vector[nyears] byear_raw_vector_6_26_DE;
vector[nyears] byear_raw_vector_6_26_NDE;
vector[nyears] byear_raw_vector_7_6;
vector[nyears] byear_raw_vector_7_28_DE;
vector[nyears] byear_raw_vector_7_28_NDE;
vector[nyears] byear_raw_vector_7_30_DE;
vector[nyears] byear_raw_vector_7_30_NDE;
vector[nyears] byear_raw_vector_7_42;
vector[nyears] byear_raw_vector_8_3;
vector[nyears] byear_raw_vector_8_9;
vector[nyears] byear_raw_vector_8_32_DE;
vector[nyears] byear_raw_vector_8_32_NDE;
vector[nyears] byear_raw_vector_9_8;
vector[nyears] byear_raw_vector_9_34_DE;
vector[nyears] byear_raw_vector_9_34_NDE;
vector[nyears] byear_raw_vector_9_36;
vector[nyears] byear_raw_vector_9_37;
vector[nyears] byear_raw_vector_9_38;
vector[nyears] byear_raw_vector_9_39_DE;
vector[nyears] byear_raw_vector_9_39_NDE;
vector[nyears] byear_raw_vector_10_2;
vector[nyears] byear_raw_vector_14_2;
vector[nyears] byear_raw_vector_24_5;
vector[nyears] byear_raw_vector_26_6;
vector[nyears] byear_raw_vector_28_7;
vector[nyears] byear_raw_vector_30_7;
vector[nyears] byear_raw_vector_36_9;
vector[nyears] byear_raw_vector_42_7;

// origin parameters
real borigin1_matrix_3_4;
real borigin1_matrix_4_3;
real borigin1_matrix_4_5;
real borigin1_matrix_5_4;
real borigin1_matrix_5_6;
real borigin1_matrix_5_24_DE;
real borigin1_matrix_5_24_NDE;
real borigin1_matrix_6_5;
real borigin1_matrix_6_7;
real borigin1_matrix_6_26_DE;
real borigin1_matrix_6_26_NDE;
real borigin1_matrix_7_6;
real borigin1_matrix_7_28_DE;
real borigin1_matrix_7_28_NDE;
real borigin1_matrix_7_30_DE;
real borigin1_matrix_7_30_NDE;
real borigin1_matrix_7_42;
real borigin1_matrix_24_5;
real borigin1_matrix_26_6;
real borigin1_matrix_28_7;
real borigin1_matrix_30_7;
real borigin1_matrix_42_7;

real borigin2_matrix_3_4;
real borigin2_matrix_4_3;
real borigin2_matrix_4_5;
real borigin2_matrix_5_4;
real borigin2_matrix_5_6;
real borigin2_matrix_5_24_DE;
real borigin2_matrix_5_24_NDE;
real borigin2_matrix_6_5;
real borigin2_matrix_6_7;
real borigin2_matrix_6_26_DE;
real borigin2_matrix_6_26_NDE;
real borigin2_matrix_7_6;
real borigin2_matrix_7_28_DE;
real borigin2_matrix_7_28_NDE;
real borigin2_matrix_7_30_DE;
real borigin2_matrix_7_30_NDE;
real borigin2_matrix_7_42;
real borigin2_matrix_24_5;
real borigin2_matrix_26_6;
real borigin2_matrix_28_7;
real borigin2_matrix_30_7;
real borigin2_matrix_42_7;


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
real walla_walla_alpha4;
real wenatchee_alpha1;
real yakima_alpha1;

// 14 terms for discharge relationship, one for each tributary
real asotin_beta;
real deschutes_beta;
real entiat_beta;
// real fifteenmile_beta;
real hood_beta;
// real imnaha_beta;
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
  // We now need two matrices always - one for DE params, and one for NDE params
  
  
  // Declare a matrix to store b0 params
  // matrix[43,43] b0_matrix;
  array[43,43] real b0_matrix_DE;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // b0_matrix = rep_matrix(-100000, 43, 43);
  b0_matrix_DE = rep_array(-100000, 43, 43);
  
  array[43,43] real b0_matrix_NDE;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // b0_matrix = rep_matrix(-100000, 43, 43);
  b0_matrix_NDE = rep_array(-100000, 43, 43);
  
    // Declare matrices to store sigma_year params 
    // 2023-03-60 - I don't think we need these anymore, given updated syntax
    // 2023-07-05 - I do think we need these, because you need to be able
    // to pass them as arguments to reduce_sum so it knows to change them (and the 
    // parameters within) when finding the minimum NLL
  array[43,43] real sigma_year_matrix_DE;
  sigma_year_matrix_DE = rep_array(0, 43, 43);
  // 
  array[43,43] real sigma_year_matrix_NDE;
  sigma_year_matrix_NDE = rep_array(0, 43, 43);
  
  // Declare matrices to store byear actual parameters (DE and NDE)
  array [43,43,nyears] real byear_actual_parameters_array_DE;
  byear_actual_parameters_array_DE = rep_array(0, 43, 43, nyears);
  array [43,43,nyears] real byear_actual_parameters_array_NDE;
  byear_actual_parameters_array_NDE = rep_array(0, 43, 43, nyears);
  
  // Declare matrices to store byear raw parameters (DE and NDE)
  array [43,43,nyears] real byear_raw_parameters_array_DE;
  byear_raw_parameters_array_DE = rep_array(0, 43, 43, nyears);
  array [43,43,nyears] real byear_raw_parameters_array_NDE;
  byear_raw_parameters_array_NDE = rep_array(0, 43, 43, nyears);
  
  // store the byear_raw_vectors in the array
byear_raw_parameters_array_DE[1,2, ] = to_array_1d(byear_raw_vector_1_2);
byear_raw_parameters_array_NDE[1,2, ] = to_array_1d(byear_raw_vector_1_2);
byear_raw_parameters_array_DE[2,1, ] = to_array_1d(byear_raw_vector_2_1);
byear_raw_parameters_array_NDE[2,1, ] = to_array_1d(byear_raw_vector_2_1);
byear_raw_parameters_array_DE[2,3, ] = to_array_1d(byear_raw_vector_2_3);
byear_raw_parameters_array_NDE[2,3, ] = to_array_1d(byear_raw_vector_2_3);
byear_raw_parameters_array_DE[2,10, ] = to_array_1d(byear_raw_vector_2_10_DE);
byear_raw_parameters_array_NDE[2,10, ] = to_array_1d(byear_raw_vector_2_10_NDE);
byear_raw_parameters_array_DE[2,11, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[2,11, ] = to_array_1d(byear_raw_vector_2_10_NDE);
byear_raw_parameters_array_DE[2,12, ] = to_array_1d(byear_raw_vector_2_12_DE);
byear_raw_parameters_array_NDE[2,12, ] = to_array_1d(byear_raw_vector_2_12_NDE);
byear_raw_parameters_array_DE[2,13, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[2,13, ] = to_array_1d(byear_raw_vector_2_12_NDE);
byear_raw_parameters_array_DE[2,14, ] = to_array_1d(byear_raw_vector_2_14_DE);
byear_raw_parameters_array_NDE[2,14, ] = to_array_1d(byear_raw_vector_2_14_NDE);
byear_raw_parameters_array_DE[2,15, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[2,15, ] = to_array_1d(byear_raw_vector_2_14_NDE);
byear_raw_parameters_array_DE[2,16, ] = to_array_1d(byear_raw_vector_2_16_DE);
byear_raw_parameters_array_NDE[2,16, ] = to_array_1d(byear_raw_vector_2_16_NDE);
byear_raw_parameters_array_DE[2,17, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[2,17, ] = to_array_1d(byear_raw_vector_2_16_NDE);
byear_raw_parameters_array_DE[2,18, ] = to_array_1d(byear_raw_vector_2_18_DE);
byear_raw_parameters_array_NDE[2,18, ] = to_array_1d(byear_raw_vector_2_18_NDE);
byear_raw_parameters_array_DE[2,19, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[2,19, ] = to_array_1d(byear_raw_vector_2_18_NDE);
byear_raw_parameters_array_DE[2,41, ] = to_array_1d(byear_raw_vector_2_41);
byear_raw_parameters_array_NDE[2,41, ] = to_array_1d(byear_raw_vector_2_41);
byear_raw_parameters_array_DE[3,2, ] = to_array_1d(byear_raw_vector_3_2);
byear_raw_parameters_array_NDE[3,2, ] = to_array_1d(byear_raw_vector_3_2);
byear_raw_parameters_array_DE[3,4, ] = to_array_1d(byear_raw_vector_3_4);
byear_raw_parameters_array_NDE[3,4, ] = to_array_1d(byear_raw_vector_3_4);
byear_raw_parameters_array_DE[3,8, ] = to_array_1d(byear_raw_vector_3_8);
byear_raw_parameters_array_NDE[3,8, ] = to_array_1d(byear_raw_vector_3_8);
byear_raw_parameters_array_DE[3,20, ] = to_array_1d(byear_raw_vector_3_20_DE);
byear_raw_parameters_array_NDE[3,20, ] = to_array_1d(byear_raw_vector_3_20_NDE);
byear_raw_parameters_array_DE[3,21, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[3,21, ] = to_array_1d(byear_raw_vector_3_20_NDE);
byear_raw_parameters_array_DE[3,22, ] = to_array_1d(byear_raw_vector_3_22_DE);
byear_raw_parameters_array_NDE[3,22, ] = to_array_1d(byear_raw_vector_3_22_NDE);
byear_raw_parameters_array_DE[3,23, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[3,23, ] = to_array_1d(byear_raw_vector_3_22_NDE);
byear_raw_parameters_array_DE[4,3, ] = to_array_1d(byear_raw_vector_4_3);
byear_raw_parameters_array_NDE[4,3, ] = to_array_1d(byear_raw_vector_4_3);
byear_raw_parameters_array_DE[4,5, ] = to_array_1d(byear_raw_vector_4_5);
byear_raw_parameters_array_NDE[4,5, ] = to_array_1d(byear_raw_vector_4_5);
byear_raw_parameters_array_DE[5,4, ] = to_array_1d(byear_raw_vector_5_4);
byear_raw_parameters_array_NDE[5,4, ] = to_array_1d(byear_raw_vector_5_4);
byear_raw_parameters_array_DE[5,6, ] = to_array_1d(byear_raw_vector_5_6);
byear_raw_parameters_array_NDE[5,6, ] = to_array_1d(byear_raw_vector_5_6);
byear_raw_parameters_array_DE[5,24, ] = to_array_1d(byear_raw_vector_5_24_DE);
byear_raw_parameters_array_NDE[5,24, ] = to_array_1d(byear_raw_vector_5_24_NDE);
byear_raw_parameters_array_DE[5,25, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[5,25, ] = to_array_1d(byear_raw_vector_5_24_NDE);
byear_raw_parameters_array_DE[6,5, ] = to_array_1d(byear_raw_vector_6_5);
byear_raw_parameters_array_NDE[6,5, ] = to_array_1d(byear_raw_vector_6_5);
byear_raw_parameters_array_DE[6,7, ] = to_array_1d(byear_raw_vector_6_7);
byear_raw_parameters_array_NDE[6,7, ] = to_array_1d(byear_raw_vector_6_7);
byear_raw_parameters_array_DE[6,26, ] = to_array_1d(byear_raw_vector_6_26_DE);
byear_raw_parameters_array_NDE[6,26, ] = to_array_1d(byear_raw_vector_6_26_NDE);
byear_raw_parameters_array_DE[6,27, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[6,27, ] = to_array_1d(byear_raw_vector_6_26_NDE);
byear_raw_parameters_array_DE[7,6, ] = to_array_1d(byear_raw_vector_7_6);
byear_raw_parameters_array_NDE[7,6, ] = to_array_1d(byear_raw_vector_7_6);
byear_raw_parameters_array_DE[7,28, ] = to_array_1d(byear_raw_vector_7_28_DE);
byear_raw_parameters_array_NDE[7,28, ] = to_array_1d(byear_raw_vector_7_28_NDE);
byear_raw_parameters_array_DE[7,29, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[7,29, ] = to_array_1d(byear_raw_vector_7_28_NDE);
byear_raw_parameters_array_DE[7,30, ] = to_array_1d(byear_raw_vector_7_30_DE);
byear_raw_parameters_array_NDE[7,30, ] = to_array_1d(byear_raw_vector_7_30_NDE);
byear_raw_parameters_array_DE[7,31, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[7,31, ] = to_array_1d(byear_raw_vector_7_30_NDE);
byear_raw_parameters_array_DE[7,42, ] = to_array_1d(byear_raw_vector_7_42);
byear_raw_parameters_array_NDE[7,42, ] = to_array_1d(byear_raw_vector_7_42);
byear_raw_parameters_array_DE[8,3, ] = to_array_1d(byear_raw_vector_8_3);
byear_raw_parameters_array_NDE[8,3, ] = to_array_1d(byear_raw_vector_8_3);
byear_raw_parameters_array_DE[8,9, ] = to_array_1d(byear_raw_vector_8_9);
byear_raw_parameters_array_NDE[8,9, ] = to_array_1d(byear_raw_vector_8_9);
byear_raw_parameters_array_DE[8,32, ] = to_array_1d(byear_raw_vector_8_32_DE);
byear_raw_parameters_array_NDE[8,32, ] = to_array_1d(byear_raw_vector_8_32_NDE);
byear_raw_parameters_array_DE[8,33, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[8,33, ] = to_array_1d(byear_raw_vector_8_32_NDE);
byear_raw_parameters_array_DE[9,8, ] = to_array_1d(byear_raw_vector_9_8);
byear_raw_parameters_array_NDE[9,8, ] = to_array_1d(byear_raw_vector_9_8);
byear_raw_parameters_array_DE[9,34, ] = to_array_1d(byear_raw_vector_9_34_DE);
byear_raw_parameters_array_NDE[9,34, ] = to_array_1d(byear_raw_vector_9_34_NDE);
byear_raw_parameters_array_DE[9,35, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[9,35, ] = to_array_1d(byear_raw_vector_9_34_NDE);
byear_raw_parameters_array_DE[9,36, ] = to_array_1d(byear_raw_vector_9_36);
byear_raw_parameters_array_NDE[9,36, ] = to_array_1d(byear_raw_vector_9_36);
byear_raw_parameters_array_DE[9,37, ] = to_array_1d(byear_raw_vector_9_37);
byear_raw_parameters_array_NDE[9,37, ] = to_array_1d(byear_raw_vector_9_37);
byear_raw_parameters_array_DE[9,38, ] = to_array_1d(byear_raw_vector_9_38);
byear_raw_parameters_array_NDE[9,38, ] = to_array_1d(byear_raw_vector_9_38);
byear_raw_parameters_array_DE[9,39, ] = to_array_1d(byear_raw_vector_9_39_DE);
byear_raw_parameters_array_NDE[9,39, ] = to_array_1d(byear_raw_vector_9_39_NDE);
byear_raw_parameters_array_DE[9,40, ] = rep_array(-100000, nyears);
byear_raw_parameters_array_NDE[9,40, ] = to_array_1d(byear_raw_vector_9_39_NDE);
byear_raw_parameters_array_DE[10,2, ] = to_array_1d(byear_raw_vector_10_2);
byear_raw_parameters_array_NDE[10,2, ] = to_array_1d(byear_raw_vector_10_2);
byear_raw_parameters_array_DE[14,2, ] = to_array_1d(byear_raw_vector_14_2);
byear_raw_parameters_array_NDE[14,2, ] = to_array_1d(byear_raw_vector_14_2);
byear_raw_parameters_array_DE[24,5, ] = to_array_1d(byear_raw_vector_24_5);
byear_raw_parameters_array_NDE[24,5, ] = to_array_1d(byear_raw_vector_24_5);
byear_raw_parameters_array_DE[25,5, ] = to_array_1d(byear_raw_vector_24_5);
byear_raw_parameters_array_NDE[25,5, ] = to_array_1d(byear_raw_vector_24_5);
byear_raw_parameters_array_DE[26,6, ] = to_array_1d(byear_raw_vector_26_6);
byear_raw_parameters_array_NDE[26,6, ] = to_array_1d(byear_raw_vector_26_6);
byear_raw_parameters_array_DE[27,6, ] = to_array_1d(byear_raw_vector_26_6);
byear_raw_parameters_array_NDE[27,6, ] = to_array_1d(byear_raw_vector_26_6);
byear_raw_parameters_array_DE[28,7, ] = to_array_1d(byear_raw_vector_28_7);
byear_raw_parameters_array_NDE[28,7, ] = to_array_1d(byear_raw_vector_28_7);
byear_raw_parameters_array_DE[29,7, ] = to_array_1d(byear_raw_vector_28_7);
byear_raw_parameters_array_NDE[29,7, ] = to_array_1d(byear_raw_vector_28_7);
byear_raw_parameters_array_DE[30,7, ] = to_array_1d(byear_raw_vector_30_7);
byear_raw_parameters_array_NDE[30,7, ] = to_array_1d(byear_raw_vector_30_7);
byear_raw_parameters_array_DE[31,7, ] = to_array_1d(byear_raw_vector_30_7);
byear_raw_parameters_array_NDE[31,7, ] = to_array_1d(byear_raw_vector_30_7);
byear_raw_parameters_array_DE[36,9, ] = to_array_1d(byear_raw_vector_36_9);
byear_raw_parameters_array_NDE[36,9, ] = to_array_1d(byear_raw_vector_36_9);
byear_raw_parameters_array_DE[42,7, ] = to_array_1d(byear_raw_vector_42_7);
byear_raw_parameters_array_NDE[42,7, ] = to_array_1d(byear_raw_vector_42_7);
  
  
  
  // now, take the byear_raw_vector parameters and multiply them by the sigma (scale) to store them)
    // take the raw matrices and transform them to the actual matrices
    // byear_actual_parameters_array_DE[sigma_year_indices[p,1], sigma_year_indices[p,2],] =  to_array_1d(to_row_vector(byear_raw_parameters_array_DE[sigma_year_indices[p,1], sigma_year_indices[p,2],]) * sigma_year_matrix_DE[sigma_year_indices[p,1], sigma_year_indices[p,2]]);
    // byear_actual_parameters_array_NDE[sigma_year_indices[p,1], sigma_year_indices[p,2],] =  to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[sigma_year_indices[p,1], sigma_year_indices[p,2],]) * sigma_year_matrix_NDE[sigma_year_indices[p,1], sigma_year_indices[p,2]]);
    
byear_actual_parameters_array_DE[1,2, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[1,2, ]) * sigma_year_matrix_DE[1,2]);
byear_actual_parameters_array_NDE[1,2, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[1,2, ]) * sigma_year_matrix_NDE[1,2]);
byear_actual_parameters_array_DE[2,1, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[2,1, ]) * sigma_year_matrix_DE[2,1]);
byear_actual_parameters_array_NDE[2,1, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,1, ]) * sigma_year_matrix_NDE[2,1]);
byear_actual_parameters_array_DE[2,3, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[2,3, ]) * sigma_year_matrix_DE[2,3]);
byear_actual_parameters_array_NDE[2,3, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,3, ]) * sigma_year_matrix_NDE[2,3]);
byear_actual_parameters_array_DE[2,10, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[2,10, ]) * sigma_year_matrix_DE[2,10]);
byear_actual_parameters_array_NDE[2,10, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,10, ]) * sigma_year_matrix_NDE[2,10]);
byear_actual_parameters_array_DE[2,11, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[2,11, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,10, ]) * sigma_year_matrix_NDE[2,10]);
byear_actual_parameters_array_DE[2,12, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[2,12, ]) * sigma_year_matrix_DE[2,12]);
byear_actual_parameters_array_NDE[2,12, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,12, ]) * sigma_year_matrix_NDE[2,12]);
byear_actual_parameters_array_DE[2,13, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[2,13, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,12, ]) * sigma_year_matrix_NDE[2,12]);
byear_actual_parameters_array_DE[2,14, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[2,14, ]) * sigma_year_matrix_DE[2,14]);
byear_actual_parameters_array_NDE[2,14, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,14, ]) * sigma_year_matrix_NDE[2,14]);
byear_actual_parameters_array_DE[2,15, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[2,15, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,14, ]) * sigma_year_matrix_NDE[2,14]);
byear_actual_parameters_array_DE[2,16, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[2,16, ]) * sigma_year_matrix_DE[2,16]);
byear_actual_parameters_array_NDE[2,16, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,16, ]) * sigma_year_matrix_NDE[2,16]);
byear_actual_parameters_array_DE[2,17, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[2,17, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,16, ]) * sigma_year_matrix_NDE[2,16]);
byear_actual_parameters_array_DE[2,18, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[2,18, ]) * sigma_year_matrix_DE[2,18]);
byear_actual_parameters_array_NDE[2,18, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,18, ]) * sigma_year_matrix_NDE[2,18]);
byear_actual_parameters_array_DE[2,19, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[2,19, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,18, ]) * sigma_year_matrix_NDE[2,18]);
byear_actual_parameters_array_DE[2,41, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[2,41, ]) * sigma_year_matrix_DE[2,41]);
byear_actual_parameters_array_NDE[2,41, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[2,41, ]) * sigma_year_matrix_NDE[2,41]);
byear_actual_parameters_array_DE[3,2, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[3,2, ]) * sigma_year_matrix_DE[3,2]);
byear_actual_parameters_array_NDE[3,2, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[3,2, ]) * sigma_year_matrix_NDE[3,2]);
byear_actual_parameters_array_DE[3,4, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[3,4, ]) * sigma_year_matrix_DE[3,4]);
byear_actual_parameters_array_NDE[3,4, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[3,4, ]) * sigma_year_matrix_NDE[3,4]);
byear_actual_parameters_array_DE[3,8, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[3,8, ]) * sigma_year_matrix_DE[3,8]);
byear_actual_parameters_array_NDE[3,8, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[3,8, ]) * sigma_year_matrix_NDE[3,8]);
byear_actual_parameters_array_DE[3,20, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[3,20, ]) * sigma_year_matrix_DE[3,20]);
byear_actual_parameters_array_NDE[3,20, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[3,20, ]) * sigma_year_matrix_NDE[3,20]);
byear_actual_parameters_array_DE[3,21, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[3,21, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[3,20, ]) * sigma_year_matrix_NDE[3,20]);
byear_actual_parameters_array_DE[3,22, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[3,22, ]) * sigma_year_matrix_DE[3,22]);
byear_actual_parameters_array_NDE[3,22, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[3,22, ]) * sigma_year_matrix_NDE[3,22]);
byear_actual_parameters_array_DE[3,23, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[3,23, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[3,22, ]) * sigma_year_matrix_NDE[3,22]);
byear_actual_parameters_array_DE[4,3, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[4,3, ]) * sigma_year_matrix_DE[4,3]);
byear_actual_parameters_array_NDE[4,3, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[4,3, ]) * sigma_year_matrix_NDE[4,3]);
byear_actual_parameters_array_DE[4,5, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[4,5, ]) * sigma_year_matrix_DE[4,5]);
byear_actual_parameters_array_NDE[4,5, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[4,5, ]) * sigma_year_matrix_NDE[4,5]);
byear_actual_parameters_array_DE[5,4, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[5,4, ]) * sigma_year_matrix_DE[5,4]);
byear_actual_parameters_array_NDE[5,4, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[5,4, ]) * sigma_year_matrix_NDE[5,4]);
byear_actual_parameters_array_DE[5,6, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[5,6, ]) * sigma_year_matrix_DE[5,6]);
byear_actual_parameters_array_NDE[5,6, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[5,6, ]) * sigma_year_matrix_NDE[5,6]);
byear_actual_parameters_array_DE[5,24, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[5,24, ]) * sigma_year_matrix_DE[5,24]);
byear_actual_parameters_array_NDE[5,24, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[5,24, ]) * sigma_year_matrix_NDE[5,24]);
byear_actual_parameters_array_DE[5,25, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[5,25, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[5,24, ]) * sigma_year_matrix_NDE[5,24]);
byear_actual_parameters_array_DE[6,5, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[6,5, ]) * sigma_year_matrix_DE[6,5]);
byear_actual_parameters_array_NDE[6,5, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[6,5, ]) * sigma_year_matrix_NDE[6,5]);
byear_actual_parameters_array_DE[6,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[6,7, ]) * sigma_year_matrix_DE[6,7]);
byear_actual_parameters_array_NDE[6,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[6,7, ]) * sigma_year_matrix_NDE[6,7]);
byear_actual_parameters_array_DE[6,26, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[6,26, ]) * sigma_year_matrix_DE[6,26]);
byear_actual_parameters_array_NDE[6,26, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[6,26, ]) * sigma_year_matrix_NDE[6,26]);
byear_actual_parameters_array_DE[6,27, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[6,27, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[6,26, ]) * sigma_year_matrix_NDE[6,26]);
byear_actual_parameters_array_DE[7,6, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[7,6, ]) * sigma_year_matrix_DE[7,6]);
byear_actual_parameters_array_NDE[7,6, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[7,6, ]) * sigma_year_matrix_NDE[7,6]);
byear_actual_parameters_array_DE[7,28, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[7,28, ]) * sigma_year_matrix_DE[7,28]);
byear_actual_parameters_array_NDE[7,28, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[7,28, ]) * sigma_year_matrix_NDE[7,28]);
byear_actual_parameters_array_DE[7,29, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[7,29, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[7,28, ]) * sigma_year_matrix_NDE[7,28]);
byear_actual_parameters_array_DE[7,30, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[7,30, ]) * sigma_year_matrix_DE[7,30]);
byear_actual_parameters_array_NDE[7,30, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[7,30, ]) * sigma_year_matrix_NDE[7,30]);
byear_actual_parameters_array_DE[7,31, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[7,31, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[7,30, ]) * sigma_year_matrix_NDE[7,30]);
byear_actual_parameters_array_DE[7,42, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[7,42, ]) * sigma_year_matrix_DE[7,42]);
byear_actual_parameters_array_NDE[7,42, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[7,42, ]) * sigma_year_matrix_NDE[7,42]);
byear_actual_parameters_array_DE[8,3, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[8,3, ]) * sigma_year_matrix_DE[8,3]);
byear_actual_parameters_array_NDE[8,3, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[8,3, ]) * sigma_year_matrix_NDE[8,3]);
byear_actual_parameters_array_DE[8,9, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[8,9, ]) * sigma_year_matrix_DE[8,9]);
byear_actual_parameters_array_NDE[8,9, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[8,9, ]) * sigma_year_matrix_NDE[8,9]);
byear_actual_parameters_array_DE[8,32, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[8,32, ]) * sigma_year_matrix_DE[8,32]);
byear_actual_parameters_array_NDE[8,32, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[8,32, ]) * sigma_year_matrix_NDE[8,32]);
byear_actual_parameters_array_DE[8,33, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[8,33, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[8,32, ]) * sigma_year_matrix_NDE[8,32]);
byear_actual_parameters_array_DE[9,8, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[9,8, ]) * sigma_year_matrix_DE[9,8]);
byear_actual_parameters_array_NDE[9,8, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[9,8, ]) * sigma_year_matrix_NDE[9,8]);
byear_actual_parameters_array_DE[9,34, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[9,34, ]) * sigma_year_matrix_DE[9,34]);
byear_actual_parameters_array_NDE[9,34, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[9,34, ]) * sigma_year_matrix_NDE[9,34]);
byear_actual_parameters_array_DE[9,35, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[9,35, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[9,34, ]) * sigma_year_matrix_NDE[9,34]);
byear_actual_parameters_array_DE[9,36, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[9,36, ]) * sigma_year_matrix_DE[9,36]);
byear_actual_parameters_array_NDE[9,36, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[9,36, ]) * sigma_year_matrix_NDE[9,36]);
byear_actual_parameters_array_DE[9,37, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[9,37, ]) * sigma_year_matrix_DE[9,37]);
byear_actual_parameters_array_NDE[9,37, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[9,37, ]) * sigma_year_matrix_NDE[9,37]);
byear_actual_parameters_array_DE[9,38, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[9,38, ]) * sigma_year_matrix_DE[9,38]);
byear_actual_parameters_array_NDE[9,38, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[9,38, ]) * sigma_year_matrix_NDE[9,38]);
byear_actual_parameters_array_DE[9,39, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[9,39, ]) * sigma_year_matrix_DE[9,39]);
byear_actual_parameters_array_NDE[9,39, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[9,39, ]) * sigma_year_matrix_NDE[9,39]);
byear_actual_parameters_array_DE[9,40, ] = rep_array(-100000, nyears);
byear_actual_parameters_array_NDE[9,40, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[9,39, ]) * sigma_year_matrix_NDE[9,39]);
byear_actual_parameters_array_DE[10,2, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[10,2, ]) * sigma_year_matrix_DE[10,2]);
byear_actual_parameters_array_NDE[10,2, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[10,2, ]) * sigma_year_matrix_NDE[10,2]);
byear_actual_parameters_array_DE[14,2, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[14,2, ]) * sigma_year_matrix_DE[14,2]);
byear_actual_parameters_array_NDE[14,2, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[14,2, ]) * sigma_year_matrix_NDE[14,2]);
byear_actual_parameters_array_DE[24,5, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[24,5, ]) * sigma_year_matrix_DE[24,5]);
byear_actual_parameters_array_NDE[24,5, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[24,5, ]) * sigma_year_matrix_NDE[24,5]);
byear_actual_parameters_array_DE[25,5, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[24,5, ]) * sigma_year_matrix_DE[24,5]);
byear_actual_parameters_array_NDE[25,5, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[24,5, ]) * sigma_year_matrix_NDE[24,5]);
byear_actual_parameters_array_DE[26,6, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[26,6, ]) * sigma_year_matrix_DE[26,6]);
byear_actual_parameters_array_NDE[26,6, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[26,6, ]) * sigma_year_matrix_NDE[26,6]);
byear_actual_parameters_array_DE[27,6, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[26,6, ]) * sigma_year_matrix_DE[26,6]);
byear_actual_parameters_array_NDE[27,6, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[26,6, ]) * sigma_year_matrix_NDE[26,6]);
byear_actual_parameters_array_DE[28,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[28,7, ]) * sigma_year_matrix_DE[28,7]);
byear_actual_parameters_array_NDE[28,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[28,7, ]) * sigma_year_matrix_NDE[28,7]);
byear_actual_parameters_array_DE[29,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[28,7, ]) * sigma_year_matrix_DE[28,7]);
byear_actual_parameters_array_NDE[29,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[28,7, ]) * sigma_year_matrix_NDE[28,7]);
byear_actual_parameters_array_DE[30,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[30,7, ]) * sigma_year_matrix_DE[30,7]);
byear_actual_parameters_array_NDE[30,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[30,7, ]) * sigma_year_matrix_NDE[30,7]);
byear_actual_parameters_array_DE[31,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[30,7, ]) * sigma_year_matrix_DE[30,7]);
byear_actual_parameters_array_NDE[31,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[30,7, ]) * sigma_year_matrix_NDE[30,7]);
byear_actual_parameters_array_DE[36,9, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[36,9, ]) * sigma_year_matrix_DE[36,9]);
byear_actual_parameters_array_NDE[36,9, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[36,9, ]) * sigma_year_matrix_NDE[36,9]);
byear_actual_parameters_array_DE[42,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_DE[42,7, ]) * sigma_year_matrix_DE[42,7]);
byear_actual_parameters_array_NDE[42,7, ] = to_array_1d(to_row_vector(byear_raw_parameters_array_NDE[42,7, ]) * sigma_year_matrix_NDE[42,7]);

  // Declare a matrix to store borigin1 params
  // matrix[43,43] borigin1_matrix;
  array[43,43] real borigin1_matrix_DE;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Use a zero instead, since we already have the -100000 in the b0 term. So we just want these to have no effect
  // Non-zero elements will be overwritten
  // borigin1_matrix = rep_matrix(-100000, 43, 43);
  // borigin1_matrix = rep_matrix(0, 43, 43);
  borigin1_matrix_DE = rep_array(0, 43, 43);
  
  array[43,43] real borigin1_matrix_NDE;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Use a zero instead, since we already have the -100000 in the b0 term. So we just want these to have no effect
  // Non-zero elements will be overwritten
  // borigin1_matrix = rep_matrix(-100000, 43, 43);
  // borigin1_matrix = rep_matrix(0, 43, 43);
  borigin1_matrix_NDE = rep_array(0, 43, 43);
  
  // Declare a matrix to store borigin2 params
  // matrix[43,43] borigin2_matrix;
  array[43,43] real borigin2_matrix_DE;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // borigin2_matrix = rep_matrix(-100000, 43, 43);
  // borigin2_matrix = rep_matrix(0, 43, 43);
  borigin2_matrix_DE = rep_array(0, 43, 43);
  
  array[43,43] real borigin2_matrix_NDE;
  // Set all of the elements of the b0 matrix to -100000 (effectively a 0 in logit space);
  // Non-zero elements will be overwritten
  // borigin2_matrix = rep_matrix(-100000, 43, 43);
  // borigin2_matrix = rep_matrix(0, 43, 43);
  borigin2_matrix_NDE = rep_array(0, 43, 43);

  
  // Populate this matrix with betas
  // Fill in all of the non-zero elements into the b0_matrix
b0_matrix_DE[1,2] = b0_matrix_1_2;
b0_matrix_NDE[1,2] = b0_matrix_1_2;
b0_matrix_DE[2,1] = b0_matrix_2_1;
b0_matrix_NDE[2,1] = b0_matrix_2_1;
b0_matrix_DE[2,3] = b0_matrix_2_3;
b0_matrix_NDE[2,3] = b0_matrix_2_3;
b0_matrix_DE[2,10] = b0_matrix_2_10_DE;
b0_matrix_NDE[2,10] = b0_matrix_2_10_NDE;
b0_matrix_DE[2,11] = -100000;
b0_matrix_NDE[2,11] = b0_matrix_2_10_NDE;
b0_matrix_DE[2,12] = b0_matrix_2_12_DE;
b0_matrix_NDE[2,12] = b0_matrix_2_12_NDE;
b0_matrix_DE[2,13] = -100000;
b0_matrix_NDE[2,13] = b0_matrix_2_12_NDE;
b0_matrix_DE[2,14] = b0_matrix_2_14_DE;
b0_matrix_NDE[2,14] = b0_matrix_2_14_NDE;
b0_matrix_DE[2,15] = -100000;
b0_matrix_NDE[2,15] = b0_matrix_2_14_NDE;
b0_matrix_DE[2,16] = b0_matrix_2_16_DE;
b0_matrix_NDE[2,16] = b0_matrix_2_16_NDE;
b0_matrix_DE[2,17] = -100000;
b0_matrix_NDE[2,17] = b0_matrix_2_16_NDE;
b0_matrix_DE[2,18] = b0_matrix_2_18_DE;
b0_matrix_NDE[2,18] = b0_matrix_2_18_NDE;
b0_matrix_DE[2,19] = -100000;
b0_matrix_NDE[2,19] = b0_matrix_2_18_NDE;
b0_matrix_DE[2,41] = b0_matrix_2_41;
b0_matrix_NDE[2,41] = b0_matrix_2_41;
b0_matrix_DE[3,2] = b0_matrix_3_2;
b0_matrix_NDE[3,2] = b0_matrix_3_2;
b0_matrix_DE[3,4] = b0_matrix_3_4;
b0_matrix_NDE[3,4] = b0_matrix_3_4;
b0_matrix_DE[3,8] = b0_matrix_3_8;
b0_matrix_NDE[3,8] = b0_matrix_3_8;
b0_matrix_DE[3,20] = b0_matrix_3_20_DE;
b0_matrix_NDE[3,20] = b0_matrix_3_20_NDE;
b0_matrix_DE[3,21] = -100000;
b0_matrix_NDE[3,21] = b0_matrix_3_20_NDE;
b0_matrix_DE[3,22] = b0_matrix_3_22_DE;
b0_matrix_NDE[3,22] = b0_matrix_3_22_NDE;
b0_matrix_DE[3,23] = -100000;
b0_matrix_NDE[3,23] = b0_matrix_3_22_NDE;
b0_matrix_DE[4,3] = b0_matrix_4_3;
b0_matrix_NDE[4,3] = b0_matrix_4_3;
b0_matrix_DE[4,5] = b0_matrix_4_5;
b0_matrix_NDE[4,5] = b0_matrix_4_5;
b0_matrix_DE[5,4] = b0_matrix_5_4;
b0_matrix_NDE[5,4] = b0_matrix_5_4;
b0_matrix_DE[5,6] = b0_matrix_5_6;
b0_matrix_NDE[5,6] = b0_matrix_5_6;
b0_matrix_DE[5,24] = b0_matrix_5_24_DE;
b0_matrix_NDE[5,24] = b0_matrix_5_24_NDE;
b0_matrix_DE[5,25] = -100000;
b0_matrix_NDE[5,25] = b0_matrix_5_24_NDE;
b0_matrix_DE[6,5] = b0_matrix_6_5;
b0_matrix_NDE[6,5] = b0_matrix_6_5;
b0_matrix_DE[6,7] = b0_matrix_6_7;
b0_matrix_NDE[6,7] = b0_matrix_6_7;
b0_matrix_DE[6,26] = b0_matrix_6_26_DE;
b0_matrix_NDE[6,26] = b0_matrix_6_26_NDE;
b0_matrix_DE[6,27] = -100000;
b0_matrix_NDE[6,27] = b0_matrix_6_26_NDE;
b0_matrix_DE[7,6] = b0_matrix_7_6;
b0_matrix_NDE[7,6] = b0_matrix_7_6;
b0_matrix_DE[7,28] = b0_matrix_7_28_DE;
b0_matrix_NDE[7,28] = b0_matrix_7_28_NDE;
b0_matrix_DE[7,29] = -100000;
b0_matrix_NDE[7,29] = b0_matrix_7_28_NDE;
b0_matrix_DE[7,30] = b0_matrix_7_30_DE;
b0_matrix_NDE[7,30] = b0_matrix_7_30_NDE;
b0_matrix_DE[7,31] = -100000;
b0_matrix_NDE[7,31] = b0_matrix_7_30_NDE;
b0_matrix_DE[7,42] = b0_matrix_7_42;
b0_matrix_NDE[7,42] = b0_matrix_7_42;
b0_matrix_DE[8,3] = b0_matrix_8_3;
b0_matrix_NDE[8,3] = b0_matrix_8_3;
b0_matrix_DE[8,9] = b0_matrix_8_9;
b0_matrix_NDE[8,9] = b0_matrix_8_9;
b0_matrix_DE[8,32] = b0_matrix_8_32_DE;
b0_matrix_NDE[8,32] = b0_matrix_8_32_NDE;
b0_matrix_DE[8,33] = -100000;
b0_matrix_NDE[8,33] = b0_matrix_8_32_NDE;
b0_matrix_DE[9,8] = b0_matrix_9_8;
b0_matrix_NDE[9,8] = b0_matrix_9_8;
b0_matrix_DE[9,34] = b0_matrix_9_34_DE;
b0_matrix_NDE[9,34] = b0_matrix_9_34_NDE;
b0_matrix_DE[9,35] = -100000;
b0_matrix_NDE[9,35] = b0_matrix_9_34_NDE;
b0_matrix_DE[9,36] = b0_matrix_9_36;
b0_matrix_NDE[9,36] = b0_matrix_9_36;
b0_matrix_DE[9,37] = b0_matrix_9_37;
b0_matrix_NDE[9,37] = b0_matrix_9_37;
b0_matrix_DE[9,38] = b0_matrix_9_38;
b0_matrix_NDE[9,38] = b0_matrix_9_38;
b0_matrix_DE[9,39] = b0_matrix_9_39_DE;
b0_matrix_NDE[9,39] = b0_matrix_9_39_NDE;
b0_matrix_DE[9,40] = -100000;
b0_matrix_NDE[9,40] = b0_matrix_9_39_NDE;
b0_matrix_DE[10,2] = b0_matrix_10_2;
b0_matrix_NDE[10,2] = b0_matrix_10_2;
b0_matrix_DE[11,2] = b0_matrix_10_2;
b0_matrix_NDE[11,2] = b0_matrix_10_2;
b0_matrix_DE[12,2] = b0_matrix_12_2;
b0_matrix_NDE[12,2] = b0_matrix_12_2;
b0_matrix_DE[14,2] = b0_matrix_14_2;
b0_matrix_NDE[14,2] = b0_matrix_14_2;
b0_matrix_DE[16,2] = b0_matrix_16_2;
b0_matrix_NDE[16,2] = b0_matrix_16_2;
b0_matrix_DE[18,2] = b0_matrix_18_2;
b0_matrix_NDE[18,2] = b0_matrix_18_2;
b0_matrix_DE[19,2] = b0_matrix_18_2;
b0_matrix_NDE[19,2] = b0_matrix_18_2;
b0_matrix_DE[20,3] = b0_matrix_20_3;
b0_matrix_NDE[20,3] = b0_matrix_20_3;
b0_matrix_DE[22,3] = b0_matrix_22_3;
b0_matrix_NDE[22,3] = b0_matrix_22_3;
b0_matrix_DE[24,5] = b0_matrix_24_5;
b0_matrix_NDE[24,5] = b0_matrix_24_5;
b0_matrix_DE[26,6] = b0_matrix_26_6;
b0_matrix_NDE[26,6] = b0_matrix_26_6;
b0_matrix_DE[28,7] = b0_matrix_28_7;
b0_matrix_NDE[28,7] = b0_matrix_28_7;
b0_matrix_DE[29,7] = b0_matrix_28_7;
b0_matrix_NDE[29,7] = b0_matrix_28_7;
b0_matrix_DE[30,7] = b0_matrix_30_7;
b0_matrix_NDE[30,7] = b0_matrix_30_7;
b0_matrix_DE[31,7] = b0_matrix_30_7;
b0_matrix_NDE[31,7] = b0_matrix_30_7;
b0_matrix_DE[32,8] = b0_matrix_32_8;
b0_matrix_NDE[32,8] = b0_matrix_32_8;
b0_matrix_DE[34,9] = b0_matrix_34_9;
b0_matrix_NDE[34,9] = b0_matrix_34_9;
b0_matrix_DE[35,9] = b0_matrix_34_9;
b0_matrix_NDE[35,9] = b0_matrix_34_9;
b0_matrix_DE[36,9] = b0_matrix_36_9;
b0_matrix_NDE[36,9] = b0_matrix_36_9;
b0_matrix_DE[37,9] = b0_matrix_37_9;
b0_matrix_NDE[37,9] = b0_matrix_37_9;
b0_matrix_DE[38,9] = b0_matrix_38_9;
b0_matrix_NDE[38,9] = b0_matrix_38_9;
b0_matrix_DE[39,9] = b0_matrix_39_9;
b0_matrix_NDE[39,9] = b0_matrix_39_9;
b0_matrix_DE[41,2] = b0_matrix_41_2;
b0_matrix_NDE[41,2] = b0_matrix_41_2;
b0_matrix_DE[42,7] = b0_matrix_42_7;
b0_matrix_NDE[42,7] = b0_matrix_42_7;

// I think this is necessary so that sigma_year_matrix_DE[,] can be passed
// as an argument to reduce_sum, and therefore those parameters can change
sigma_year_matrix_DE[1,2] = sigma_year_matrix_1_2;
sigma_year_matrix_NDE[1,2] = sigma_year_matrix_1_2;
sigma_year_matrix_DE[2,1] = sigma_year_matrix_2_1;
sigma_year_matrix_NDE[2,1] = sigma_year_matrix_2_1;
sigma_year_matrix_DE[2,3] = sigma_year_matrix_2_3;
sigma_year_matrix_NDE[2,3] = sigma_year_matrix_2_3;
sigma_year_matrix_DE[2,10] = sigma_year_matrix_2_10_DE;
sigma_year_matrix_NDE[2,10] = sigma_year_matrix_2_10_NDE;
sigma_year_matrix_NDE[2,11] = sigma_year_matrix_2_10_NDE;
sigma_year_matrix_DE[2,12] = sigma_year_matrix_2_12_DE;
sigma_year_matrix_NDE[2,12] = sigma_year_matrix_2_12_NDE;
sigma_year_matrix_NDE[2,13] = sigma_year_matrix_2_12_NDE;
sigma_year_matrix_DE[2,14] = sigma_year_matrix_2_14_DE;
sigma_year_matrix_NDE[2,14] = sigma_year_matrix_2_14_NDE;
sigma_year_matrix_NDE[2,15] = sigma_year_matrix_2_14_NDE;
sigma_year_matrix_DE[2,16] = sigma_year_matrix_2_16_DE;
sigma_year_matrix_NDE[2,16] = sigma_year_matrix_2_16_NDE;
sigma_year_matrix_NDE[2,17] = sigma_year_matrix_2_16_NDE;
sigma_year_matrix_DE[2,18] = sigma_year_matrix_2_18_DE;
sigma_year_matrix_NDE[2,18] = sigma_year_matrix_2_18_NDE;
sigma_year_matrix_NDE[2,19] = sigma_year_matrix_2_18_NDE;
sigma_year_matrix_DE[2,41] = sigma_year_matrix_2_41;
sigma_year_matrix_NDE[2,41] = sigma_year_matrix_2_41;
sigma_year_matrix_DE[3,2] = sigma_year_matrix_3_2;
sigma_year_matrix_NDE[3,2] = sigma_year_matrix_3_2;
sigma_year_matrix_DE[3,4] = sigma_year_matrix_3_4;
sigma_year_matrix_NDE[3,4] = sigma_year_matrix_3_4;
sigma_year_matrix_DE[3,8] = sigma_year_matrix_3_8;
sigma_year_matrix_NDE[3,8] = sigma_year_matrix_3_8;
sigma_year_matrix_DE[3,20] = sigma_year_matrix_3_20_DE;
sigma_year_matrix_NDE[3,20] = sigma_year_matrix_3_20_NDE;
sigma_year_matrix_NDE[3,21] = sigma_year_matrix_3_20_NDE;
sigma_year_matrix_DE[3,22] = sigma_year_matrix_3_22_DE;
sigma_year_matrix_NDE[3,22] = sigma_year_matrix_3_22_NDE;
sigma_year_matrix_NDE[3,23] = sigma_year_matrix_3_22_NDE;
sigma_year_matrix_DE[4,3] = sigma_year_matrix_4_3;
sigma_year_matrix_NDE[4,3] = sigma_year_matrix_4_3;
sigma_year_matrix_DE[4,5] = sigma_year_matrix_4_5;
sigma_year_matrix_NDE[4,5] = sigma_year_matrix_4_5;
sigma_year_matrix_DE[5,4] = sigma_year_matrix_5_4;
sigma_year_matrix_NDE[5,4] = sigma_year_matrix_5_4;
sigma_year_matrix_DE[5,6] = sigma_year_matrix_5_6;
sigma_year_matrix_NDE[5,6] = sigma_year_matrix_5_6;
sigma_year_matrix_DE[5,24] = sigma_year_matrix_5_24_DE;
sigma_year_matrix_NDE[5,24] = sigma_year_matrix_5_24_NDE;
sigma_year_matrix_NDE[5,25] = sigma_year_matrix_5_24_NDE;
sigma_year_matrix_DE[6,5] = sigma_year_matrix_6_5;
sigma_year_matrix_NDE[6,5] = sigma_year_matrix_6_5;
sigma_year_matrix_DE[6,7] = sigma_year_matrix_6_7;
sigma_year_matrix_NDE[6,7] = sigma_year_matrix_6_7;
sigma_year_matrix_DE[6,26] = sigma_year_matrix_6_26_DE;
sigma_year_matrix_NDE[6,26] = sigma_year_matrix_6_26_NDE;
sigma_year_matrix_NDE[6,27] = sigma_year_matrix_6_26_NDE;
sigma_year_matrix_DE[7,6] = sigma_year_matrix_7_6;
sigma_year_matrix_NDE[7,6] = sigma_year_matrix_7_6;
sigma_year_matrix_DE[7,28] = sigma_year_matrix_7_28_DE;
sigma_year_matrix_NDE[7,28] = sigma_year_matrix_7_28_NDE;
sigma_year_matrix_NDE[7,29] = sigma_year_matrix_7_28_NDE;
sigma_year_matrix_DE[7,30] = sigma_year_matrix_7_30_DE;
sigma_year_matrix_NDE[7,30] = sigma_year_matrix_7_30_NDE;
sigma_year_matrix_NDE[7,31] = sigma_year_matrix_7_30_NDE;
sigma_year_matrix_DE[7,42] = sigma_year_matrix_7_42;
sigma_year_matrix_NDE[7,42] = sigma_year_matrix_7_42;
sigma_year_matrix_DE[8,3] = sigma_year_matrix_8_3;
sigma_year_matrix_NDE[8,3] = sigma_year_matrix_8_3;
sigma_year_matrix_DE[8,9] = sigma_year_matrix_8_9;
sigma_year_matrix_NDE[8,9] = sigma_year_matrix_8_9;
sigma_year_matrix_DE[8,32] = sigma_year_matrix_8_32_DE;
sigma_year_matrix_NDE[8,32] = sigma_year_matrix_8_32_NDE;
sigma_year_matrix_NDE[8,33] = sigma_year_matrix_8_32_NDE;
sigma_year_matrix_DE[9,8] = sigma_year_matrix_9_8;
sigma_year_matrix_NDE[9,8] = sigma_year_matrix_9_8;
sigma_year_matrix_DE[9,34] = sigma_year_matrix_9_34_DE;
sigma_year_matrix_NDE[9,34] = sigma_year_matrix_9_34_NDE;
sigma_year_matrix_NDE[9,35] = sigma_year_matrix_9_34_NDE;
sigma_year_matrix_DE[9,36] = sigma_year_matrix_9_36;
sigma_year_matrix_NDE[9,36] = sigma_year_matrix_9_36;
sigma_year_matrix_DE[9,37] = sigma_year_matrix_9_37;
sigma_year_matrix_NDE[9,37] = sigma_year_matrix_9_37;
sigma_year_matrix_DE[9,38] = sigma_year_matrix_9_38;
sigma_year_matrix_NDE[9,38] = sigma_year_matrix_9_38;
sigma_year_matrix_DE[9,39] = sigma_year_matrix_9_39_DE;
sigma_year_matrix_NDE[9,39] = sigma_year_matrix_9_39_NDE;
sigma_year_matrix_NDE[9,40] = sigma_year_matrix_9_39_NDE;
sigma_year_matrix_DE[10,2] = sigma_year_matrix_10_2;
sigma_year_matrix_NDE[10,2] = sigma_year_matrix_10_2;
sigma_year_matrix_DE[14,2] = sigma_year_matrix_14_2;
sigma_year_matrix_NDE[14,2] = sigma_year_matrix_14_2;
sigma_year_matrix_DE[24,5] = sigma_year_matrix_24_5;
sigma_year_matrix_NDE[24,5] = sigma_year_matrix_24_5;
sigma_year_matrix_DE[25,5] = sigma_year_matrix_24_5;
sigma_year_matrix_NDE[25,5] = sigma_year_matrix_24_5;
sigma_year_matrix_DE[26,6] = sigma_year_matrix_26_6;
sigma_year_matrix_NDE[26,6] = sigma_year_matrix_26_6;
sigma_year_matrix_DE[27,6] = sigma_year_matrix_26_6;
sigma_year_matrix_NDE[27,6] = sigma_year_matrix_26_6;
sigma_year_matrix_DE[28,7] = sigma_year_matrix_28_7;
sigma_year_matrix_NDE[28,7] = sigma_year_matrix_28_7;
sigma_year_matrix_DE[29,7] = sigma_year_matrix_28_7;
sigma_year_matrix_NDE[29,7] = sigma_year_matrix_28_7;
sigma_year_matrix_DE[30,7] = sigma_year_matrix_30_7;
sigma_year_matrix_NDE[30,7] = sigma_year_matrix_30_7;
sigma_year_matrix_DE[31,7] = sigma_year_matrix_30_7;
sigma_year_matrix_NDE[31,7] = sigma_year_matrix_30_7;
sigma_year_matrix_DE[36,9] = sigma_year_matrix_36_9;
sigma_year_matrix_NDE[36,9] = sigma_year_matrix_36_9;
sigma_year_matrix_DE[42,7] = sigma_year_matrix_42_7;
sigma_year_matrix_NDE[42,7] = sigma_year_matrix_42_7;

borigin1_matrix_DE[3,4] = borigin1_matrix_3_4;
borigin1_matrix_NDE[3,4] = borigin1_matrix_3_4;
borigin1_matrix_DE[4,3] = borigin1_matrix_4_3;
borigin1_matrix_NDE[4,3] = borigin1_matrix_4_3;
borigin1_matrix_DE[4,5] = borigin1_matrix_4_5;
borigin1_matrix_NDE[4,5] = borigin1_matrix_4_5;
borigin1_matrix_DE[5,4] = borigin1_matrix_5_4;
borigin1_matrix_NDE[5,4] = borigin1_matrix_5_4;
borigin1_matrix_DE[5,6] = borigin1_matrix_5_6;
borigin1_matrix_NDE[5,6] = borigin1_matrix_5_6;
borigin1_matrix_DE[5,24] = borigin1_matrix_5_24_DE;
borigin1_matrix_NDE[5,24] = borigin1_matrix_5_24_NDE;
borigin1_matrix_NDE[5,25] = borigin1_matrix_5_24_NDE;
borigin1_matrix_DE[6,5] = borigin1_matrix_6_5;
borigin1_matrix_NDE[6,5] = borigin1_matrix_6_5;
borigin1_matrix_DE[6,7] = borigin1_matrix_6_7;
borigin1_matrix_NDE[6,7] = borigin1_matrix_6_7;
borigin1_matrix_DE[6,26] = borigin1_matrix_6_26_DE;
borigin1_matrix_NDE[6,26] = borigin1_matrix_6_26_NDE;
borigin1_matrix_NDE[6,27] = borigin1_matrix_6_26_NDE;
borigin1_matrix_DE[7,6] = borigin1_matrix_7_6;
borigin1_matrix_NDE[7,6] = borigin1_matrix_7_6;
borigin1_matrix_DE[7,28] = borigin1_matrix_7_28_DE;
borigin1_matrix_NDE[7,28] = borigin1_matrix_7_28_NDE;
borigin1_matrix_NDE[7,29] = borigin1_matrix_7_28_NDE;
borigin1_matrix_DE[7,30] = borigin1_matrix_7_30_DE;
borigin1_matrix_NDE[7,30] = borigin1_matrix_7_30_NDE;
borigin1_matrix_NDE[7,31] = borigin1_matrix_7_30_NDE;
borigin1_matrix_DE[7,42] = borigin1_matrix_7_42;
borigin1_matrix_NDE[7,42] = borigin1_matrix_7_42;
borigin1_matrix_DE[24,5] = borigin1_matrix_24_5;
borigin1_matrix_NDE[24,5] = borigin1_matrix_24_5;
borigin1_matrix_DE[26,6] = borigin1_matrix_26_6;
borigin1_matrix_NDE[26,6] = borigin1_matrix_26_6;
borigin1_matrix_DE[28,7] = borigin1_matrix_28_7;
borigin1_matrix_NDE[28,7] = borigin1_matrix_28_7;
borigin1_matrix_DE[29,7] = borigin1_matrix_28_7;
borigin1_matrix_NDE[29,7] = borigin1_matrix_28_7;
borigin1_matrix_DE[30,7] = borigin1_matrix_30_7;
borigin1_matrix_NDE[30,7] = borigin1_matrix_30_7;
borigin1_matrix_DE[31,7] = borigin1_matrix_30_7;
borigin1_matrix_NDE[31,7] = borigin1_matrix_30_7;
borigin1_matrix_DE[42,7] = borigin1_matrix_42_7;
borigin1_matrix_NDE[42,7] = borigin1_matrix_42_7;

borigin2_matrix_DE[3,4] = borigin2_matrix_3_4;
borigin2_matrix_NDE[3,4] = borigin2_matrix_3_4;
borigin2_matrix_DE[4,3] = borigin2_matrix_4_3;
borigin2_matrix_NDE[4,3] = borigin2_matrix_4_3;
borigin2_matrix_DE[4,5] = borigin2_matrix_4_5;
borigin2_matrix_NDE[4,5] = borigin2_matrix_4_5;
borigin2_matrix_DE[5,4] = borigin2_matrix_5_4;
borigin2_matrix_NDE[5,4] = borigin2_matrix_5_4;
borigin2_matrix_DE[5,6] = borigin2_matrix_5_6;
borigin2_matrix_NDE[5,6] = borigin2_matrix_5_6;
borigin2_matrix_DE[5,24] = borigin2_matrix_5_24_DE;
borigin2_matrix_NDE[5,24] = borigin2_matrix_5_24_NDE;
borigin2_matrix_NDE[5,25] = borigin2_matrix_5_24_NDE;
borigin2_matrix_DE[6,5] = borigin2_matrix_6_5;
borigin2_matrix_NDE[6,5] = borigin2_matrix_6_5;
borigin2_matrix_DE[6,7] = borigin2_matrix_6_7;
borigin2_matrix_NDE[6,7] = borigin2_matrix_6_7;
borigin2_matrix_DE[6,26] = borigin2_matrix_6_26_DE;
borigin2_matrix_NDE[6,26] = borigin2_matrix_6_26_NDE;
borigin2_matrix_NDE[6,27] = borigin2_matrix_6_26_NDE;
borigin2_matrix_DE[7,6] = borigin2_matrix_7_6;
borigin2_matrix_NDE[7,6] = borigin2_matrix_7_6;
borigin2_matrix_DE[7,28] = borigin2_matrix_7_28_DE;
borigin2_matrix_NDE[7,28] = borigin2_matrix_7_28_NDE;
borigin2_matrix_NDE[7,29] = borigin2_matrix_7_28_NDE;
borigin2_matrix_DE[7,30] = borigin2_matrix_7_30_DE;
borigin2_matrix_NDE[7,30] = borigin2_matrix_7_30_NDE;
borigin2_matrix_NDE[7,31] = borigin2_matrix_7_30_NDE;
borigin2_matrix_DE[7,42] = borigin2_matrix_7_42;
borigin2_matrix_NDE[7,42] = borigin2_matrix_7_42;
borigin2_matrix_DE[24,5] = borigin2_matrix_24_5;
borigin2_matrix_NDE[24,5] = borigin2_matrix_24_5;
borigin2_matrix_DE[26,6] = borigin2_matrix_26_6;
borigin2_matrix_NDE[26,6] = borigin2_matrix_26_6;
borigin2_matrix_DE[28,7] = borigin2_matrix_28_7;
borigin2_matrix_NDE[28,7] = borigin2_matrix_28_7;
borigin2_matrix_DE[29,7] = borigin2_matrix_28_7;
borigin2_matrix_NDE[29,7] = borigin2_matrix_28_7;
borigin2_matrix_DE[30,7] = borigin2_matrix_30_7;
borigin2_matrix_NDE[30,7] = borigin2_matrix_30_7;
borigin2_matrix_DE[31,7] = borigin2_matrix_30_7;
borigin2_matrix_NDE[31,7] = borigin2_matrix_30_7;
borigin2_matrix_DE[42,7] = borigin2_matrix_42_7;
borigin2_matrix_NDE[42,7] = borigin2_matrix_42_7;

// detection efficiency - create a vector that stores all parameters
vector[35] det_eff_param_vector; // this is length 35 because that's how many det eff params we have

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
det_eff_param_vector[19] = walla_walla_alpha4;
det_eff_param_vector[20] = wenatchee_alpha1;
det_eff_param_vector[21] = yakima_alpha1;

// det_eff_param_vector[1] = det_eff_param_posteriors[1,1];
// det_eff_param_vector[2] = det_eff_param_posteriors[2,1];
// det_eff_param_vector[3] = det_eff_param_posteriors[3,1];
// det_eff_param_vector[4] = det_eff_param_posteriors[4,1];
// det_eff_param_vector[5] = det_eff_param_posteriors[5,1];
// det_eff_param_vector[6] = det_eff_param_posteriors[6,1];
// det_eff_param_vector[7] = det_eff_param_posteriors[7,1];
// det_eff_param_vector[8] = det_eff_param_posteriors[8,1];
// det_eff_param_vector[9] = det_eff_param_posteriors[9,1];
// det_eff_param_vector[10] = det_eff_param_posteriors[10,1];
// det_eff_param_vector[11] = det_eff_param_posteriors[11,1];
// det_eff_param_vector[12] = det_eff_param_posteriors[12,1];
// det_eff_param_vector[13] = det_eff_param_posteriors[13,1];
// det_eff_param_vector[14] = det_eff_param_posteriors[14,1];
// det_eff_param_vector[15] = det_eff_param_posteriors[15,1];
// det_eff_param_vector[16] = det_eff_param_posteriors[16,1];
// det_eff_param_vector[17] = det_eff_param_posteriors[17,1];
// det_eff_param_vector[18] = det_eff_param_posteriors[18,1];
// det_eff_param_vector[19] = det_eff_param_posteriors[19,1];
// det_eff_param_vector[20] = det_eff_param_posteriors[20,1];

// 14 terms for discharge relationship, one for each tributary
det_eff_param_vector[22] = asotin_beta;
det_eff_param_vector[23] = deschutes_beta;
det_eff_param_vector[24] = entiat_beta;
// det_eff_param_vector[25] = fifteenmile_beta;
// // fix these to zero
det_eff_param_vector[25] = 0;
det_eff_param_vector[26] = hood_beta;
// det_eff_param_vector[27] = imnaha_beta;
det_eff_param_vector[27] = 0;
det_eff_param_vector[28] = john_day_beta;
det_eff_param_vector[29] = methow_beta;
det_eff_param_vector[30] = okanogan_beta;
det_eff_param_vector[31] = tucannon_beta;
det_eff_param_vector[32] = umatilla_beta;
det_eff_param_vector[33] = walla_walla_beta;
det_eff_param_vector[34] = wenatchee_beta;
det_eff_param_vector[35] = yakima_beta;

// what's the deal with this block? it looks like I turned off the discharge
// turns out that I did this for troubleshooting reasons - need to turn this back 
// on when running the actual model
det_eff_param_vector[22] = 0;
det_eff_param_vector[23] = 0;
det_eff_param_vector[24] = 0;
// det_eff_param_vector[25] = fifteenmile_beta;
// fix these to zero
det_eff_param_vector[25] = 0;
det_eff_param_vector[26] = 0;
// det_eff_param_vector[27] = imnaha_beta;
det_eff_param_vector[27] = 0;
det_eff_param_vector[28] = 0;
det_eff_param_vector[29] = 0;
det_eff_param_vector[30] = 0;
det_eff_param_vector[31] = 0;
det_eff_param_vector[32] = 0;
det_eff_param_vector[33] = 0;
det_eff_param_vector[34] = 0;
det_eff_param_vector[35] = 0;

// Calculate each individual byear parameter per movement as a derived parameter
// this is because in actuality, each byear parameter for a movement is 
// calculated as a sigma



// Calculate movement probabilities as derived parameters
// make two - one for DE corrected, one for DE uncorrected
matrix[43,43] origin1_probs_DE;
origin1_probs_DE = rep_matrix(0, 43, 43);
matrix[43,43] origin2_probs_DE;
origin2_probs_DE = rep_matrix(0, 43, 43);
matrix[43,43] origin3_probs_DE;
origin3_probs_DE = rep_matrix(0, 43, 43);

matrix[43,43] origin1_probs_NDE;
origin1_probs_NDE = rep_matrix(0, 43, 43);
matrix[43,43] origin2_probs_NDE;
origin2_probs_NDE = rep_matrix(0, 43, 43);
matrix[43,43] origin3_probs_NDE;
origin3_probs_NDE = rep_matrix(0, 43, 43);

// for each of the non-loss states:
for (i in 1:42){
  for (j in 1:42){
    // first for DE
    
    // origin1_probs_DE[i,j] = exp(b0_matrix_DE[i,j] + btemp_matrix_DE[i,j]*temp + btempxorigin1_matrix_DE[i,j]*temp + borigin1_matrix_DE[i,j])/
    // (1 + sum(exp(to_row_vector(b0_matrix_DE[i,]) + to_row_vector(btemp_matrix_DE[i,]*temp) + to_row_vector(btempxorigin1_matrix_DE[i,]*temp) + to_row_vector(borigin1_matrix_DE[i,]))));
    // origin2_probs_DE[i,j] = exp(b0_matrix_DE[i,j] + btemp_matrix_DE[i,j]*temp + btempxorigin2_matrix_DE[i,j]*temp + borigin2_matrix_DE[i,j])/
    // (1 + sum(exp(to_row_vector(b0_matrix_DE[i,]) + to_row_vector(btemp_matrix_DE[i,]*temp) + to_row_vector(btempxorigin2_matrix_DE[i,]*temp) + to_row_vector(borigin2_matrix_DE[i,]))));
    // origin3_probs_DE[i,j] = exp(b0_matrix_DE[i,j] + btemp_matrix_DE[i,j]*temp - borigin1_matrix_DE[i,j] - borigin2_matrix_DE[i,j])/
    // (1 + sum(exp(to_row_vector(b0_matrix_DE[i,]) + to_row_vector(btemp_matrix_DE[i,]*temp) - to_row_vector(borigin1_matrix_DE[i,]) - to_row_vector(borigin2_matrix_DE[i,]))));
    // 
    // // then for NDE
    // 
    // origin1_probs_NDE[i,j] = exp(b0_matrix_NDE[i,j] + btemp_matrix_NDE[i,j]*temp + btempxorigin1_matrix_NDE[i,j]*temp + borigin1_matrix_NDE[i,j])/
    // (1 + sum(exp(to_row_vector(b0_matrix_NDE[i,]) + to_row_vector(btemp_matrix_NDE[i,j]*temp) + to_row_vector(btempxorigin1_matrix_NDE[i,j]*temp) + to_row_vector(borigin1_matrix_NDE[i,]))));
    // origin2_probs_NDE[i,j] = exp(b0_matrix_NDE[i,j] + btemp_matrix_NDE[i,j]*temp + btempxorigin2_matrix_NDE[i,j]*temp + borigin2_matrix_NDE[i,j])/
    // (1 + sum(exp(to_row_vector(b0_matrix_NDE[i,]) + to_row_vector(btemp_matrix_NDE[i,j]*temp) + to_row_vector(btempxorigin2_matrix_NDE[i,j]*temp) + to_row_vector(borigin2_matrix_NDE[i,]))));
    // origin3_probs_NDE[i,j] = exp(b0_matrix_NDE[i,j] + btemp_matrix_NDE[i,j]*temp - borigin1_matrix_NDE[i,j] - borigin2_matrix_NDE[i,j])/
    // (1 + sum(exp(to_row_vector(b0_matrix_NDE[i,]) + to_row_vector(btemp_matrix_NDE[i,j]*temp) - to_row_vector(borigin1_matrix_NDE[i,]) - to_row_vector(borigin2_matrix_NDE[i,]))));
    
    
    // you can't get a derived probability with temperature in it, because what temperature would you use?
    // here we are getting the movement probabilities at the "average" temperature, BUT keep
    // in mind that that won't be super informative, because the average temperature in
    // the full dataset won't be the average temperature in any one state.
    
    
    // first for DE
    
    origin1_probs_DE[i,j] = exp(b0_matrix_DE[i,j]+ borigin1_matrix_DE[i,j])/(1 + sum(exp(to_row_vector(b0_matrix_DE[i,]) + to_row_vector(borigin1_matrix_DE[i,]))));
    origin2_probs_DE[i,j] = exp(b0_matrix_DE[i,j]+ borigin2_matrix_DE[i,j])/(1 + sum(exp(to_row_vector(b0_matrix_DE[i,]) + to_row_vector(borigin2_matrix_DE[i,]))));
    origin3_probs_DE[i,j] = exp(b0_matrix_DE[i,j]- borigin1_matrix_DE[i,j] - borigin2_matrix_DE[i,j])/
    (1 + sum(exp(to_row_vector(b0_matrix_DE[i,]) - to_row_vector(borigin1_matrix_DE[i,]) - to_row_vector(borigin2_matrix_DE[i,]))));
    
    // then for NDE
    
    origin1_probs_NDE[i,j] = exp(b0_matrix_NDE[i,j]+ borigin1_matrix_NDE[i,j])/(1 + sum(exp(to_row_vector(b0_matrix_NDE[i,]) + to_row_vector(borigin1_matrix_NDE[i,]))));
    origin2_probs_NDE[i,j] = exp(b0_matrix_NDE[i,j]+ borigin2_matrix_NDE[i,j])/(1 + sum(exp(to_row_vector(b0_matrix_NDE[i,]) + to_row_vector(borigin2_matrix_NDE[i,]))));
    origin3_probs_NDE[i,j] = exp(b0_matrix_NDE[i,j]- borigin1_matrix_NDE[i,j] - borigin2_matrix_NDE[i,j])/
    (1 + sum(exp(to_row_vector(b0_matrix_NDE[i,]) - to_row_vector(borigin1_matrix_NDE[i,]) - to_row_vector(borigin2_matrix_NDE[i,]))));
  }
}

// calculate loss states
for (i in 1:42){
  origin1_probs_DE[i,43] = 1 - sum(origin1_probs_DE[i,1:42]);
  origin2_probs_DE[i,43] = 1 - sum(origin2_probs_DE[i,1:42]);
  origin3_probs_DE[i,43] = 1 - sum(origin3_probs_DE[i,1:42]);
  
  origin1_probs_NDE[i,43] = 1 - sum(origin1_probs_NDE[i,1:42]);
  origin2_probs_NDE[i,43] = 1 - sum(origin2_probs_NDE[i,1:42]);
  origin3_probs_NDE[i,43] = 1 - sum(origin3_probs_NDE[i,1:42]);
  
}

}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  // Set the priors for each of the non-zero elements of the b0 matrix
b0_matrix_1_2 ~ normal(0,10);
b0_matrix_2_1 ~ normal(0,10);
b0_matrix_2_3 ~ normal(0,10);
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
b0_matrix_2_41 ~ normal(0,10);
b0_matrix_3_2 ~ normal(0,10);
b0_matrix_3_4 ~ normal(0,10);
b0_matrix_3_8 ~ normal(0,10);
b0_matrix_3_20_DE ~ normal(0,10);
b0_matrix_3_20_NDE ~ normal(0,10);
b0_matrix_3_22_DE ~ normal(0,10);
b0_matrix_3_22_NDE ~ normal(0,10);
b0_matrix_4_3 ~ normal(0,10);
b0_matrix_4_5 ~ normal(0,10);
b0_matrix_5_4 ~ normal(0,10);
b0_matrix_5_6 ~ normal(0,10);
b0_matrix_5_24_DE ~ normal(0,10);
b0_matrix_5_24_NDE ~ normal(0,10);
b0_matrix_6_5 ~ normal(0,10);
b0_matrix_6_7 ~ normal(0,10);
b0_matrix_6_26_DE ~ normal(0,10);
b0_matrix_6_26_NDE ~ normal(0,10);
b0_matrix_7_6 ~ normal(0,10);
b0_matrix_7_28_DE ~ normal(0,10);
b0_matrix_7_28_NDE ~ normal(0,10);
b0_matrix_7_30_DE ~ normal(0,10);
b0_matrix_7_30_NDE ~ normal(0,10);
b0_matrix_7_42 ~ normal(0,10);
b0_matrix_8_3 ~ normal(0,10);
b0_matrix_8_9 ~ normal(0,10);
b0_matrix_8_32_DE ~ normal(0,10);
b0_matrix_8_32_NDE ~ normal(0,10);
b0_matrix_9_8 ~ normal(0,10);
b0_matrix_9_34_DE ~ normal(0,10);
b0_matrix_9_34_NDE ~ normal(0,10);
b0_matrix_9_36 ~ normal(0,10);
b0_matrix_9_37 ~ normal(0,10);
b0_matrix_9_38 ~ normal(0,10);
b0_matrix_9_39_DE ~ normal(0,10);
b0_matrix_9_39_NDE ~ normal(0,10);
b0_matrix_10_2 ~ normal(0,10);
b0_matrix_12_2 ~ normal(0,10);
b0_matrix_14_2 ~ normal(0,10);
b0_matrix_16_2 ~ normal(0,10);
b0_matrix_18_2 ~ normal(0,10);
b0_matrix_20_3 ~ normal(0,10);
b0_matrix_22_3 ~ normal(0,10);
b0_matrix_24_5 ~ normal(0,10);
b0_matrix_26_6 ~ normal(0,10);
b0_matrix_28_7 ~ normal(0,10);
b0_matrix_30_7 ~ normal(0,10);
b0_matrix_32_8 ~ normal(0,10);
b0_matrix_34_9 ~ normal(0,10);
b0_matrix_36_9 ~ normal(0,10);
b0_matrix_37_9 ~ normal(0,10);
b0_matrix_38_9 ~ normal(0,10);
b0_matrix_39_9 ~ normal(0,10);
b0_matrix_41_2 ~ normal(0,10);
b0_matrix_42_7 ~ normal(0,10);

sigma_year_matrix_1_2 ~ cauchy(0,1);
sigma_year_matrix_2_1 ~ cauchy(0,1);
sigma_year_matrix_2_3 ~ cauchy(0,1);
sigma_year_matrix_2_10_DE ~ cauchy(0,1);
sigma_year_matrix_2_10_NDE ~ cauchy(0,1);
sigma_year_matrix_2_12_DE ~ cauchy(0,1);
sigma_year_matrix_2_12_NDE ~ cauchy(0,1);
sigma_year_matrix_2_14_DE ~ cauchy(0,1);
sigma_year_matrix_2_14_NDE ~ cauchy(0,1);
sigma_year_matrix_2_16_DE ~ cauchy(0,1);
sigma_year_matrix_2_16_NDE ~ cauchy(0,1);
sigma_year_matrix_2_18_DE ~ cauchy(0,1);
sigma_year_matrix_2_18_NDE ~ cauchy(0,1);
sigma_year_matrix_2_41 ~ cauchy(0,1);
sigma_year_matrix_3_2 ~ cauchy(0,1);
sigma_year_matrix_3_4 ~ cauchy(0,1);
sigma_year_matrix_3_8 ~ cauchy(0,1);
sigma_year_matrix_3_20_DE ~ cauchy(0,1);
sigma_year_matrix_3_20_NDE ~ cauchy(0,1);
sigma_year_matrix_3_22_DE ~ cauchy(0,1);
sigma_year_matrix_3_22_NDE ~ cauchy(0,1);
sigma_year_matrix_4_3 ~ cauchy(0,1);
sigma_year_matrix_4_5 ~ cauchy(0,1);
sigma_year_matrix_5_4 ~ cauchy(0,1);
sigma_year_matrix_5_6 ~ cauchy(0,1);
sigma_year_matrix_5_24_DE ~ cauchy(0,1);
sigma_year_matrix_5_24_NDE ~ cauchy(0,1);
sigma_year_matrix_6_5 ~ cauchy(0,1);
sigma_year_matrix_6_7 ~ cauchy(0,1);
sigma_year_matrix_6_26_DE ~ cauchy(0,1);
sigma_year_matrix_6_26_NDE ~ cauchy(0,1);
sigma_year_matrix_7_6 ~ cauchy(0,1);
sigma_year_matrix_7_28_DE ~ cauchy(0,1);
sigma_year_matrix_7_28_NDE ~ cauchy(0,1);
sigma_year_matrix_7_30_DE ~ cauchy(0,1);
sigma_year_matrix_7_30_NDE ~ cauchy(0,1);
sigma_year_matrix_7_42 ~ cauchy(0,1);
sigma_year_matrix_8_3 ~ cauchy(0,1);
sigma_year_matrix_8_9 ~ cauchy(0,1);
sigma_year_matrix_8_32_DE ~ cauchy(0,1);
sigma_year_matrix_8_32_NDE ~ cauchy(0,1);
sigma_year_matrix_9_8 ~ cauchy(0,1);
sigma_year_matrix_9_34_DE ~ cauchy(0,1);
sigma_year_matrix_9_34_NDE ~ cauchy(0,1);
sigma_year_matrix_9_36 ~ cauchy(0,1);
sigma_year_matrix_9_37 ~ cauchy(0,1);
sigma_year_matrix_9_38 ~ cauchy(0,1);
sigma_year_matrix_9_39_DE ~ cauchy(0,1);
sigma_year_matrix_9_39_NDE ~ cauchy(0,1);
sigma_year_matrix_10_2 ~ cauchy(0,1);
sigma_year_matrix_14_2 ~ cauchy(0,1);
sigma_year_matrix_24_5 ~ cauchy(0,1);
sigma_year_matrix_26_6 ~ cauchy(0,1);
sigma_year_matrix_28_7 ~ cauchy(0,1);
sigma_year_matrix_30_7 ~ cauchy(0,1);
sigma_year_matrix_36_9 ~ cauchy(0,1);
sigma_year_matrix_42_7 ~ cauchy(0,1);

// prior on the byear_raw parameters arrays
// byear_raw_parameters_array_DE ~ normal(0,1);
// byear_raw_parameters_array_NDE ~ normal(0,1);

// these are now all vectors that get a prior, which is normal(0,1) (and are tranformed using the sigma parameters)
byear_raw_vector_1_2 ~ normal(0,1);
byear_raw_vector_2_1 ~ normal(0,1);
byear_raw_vector_2_3 ~ normal(0,1);
byear_raw_vector_2_10_DE ~ normal(0,1);
byear_raw_vector_2_10_NDE ~ normal(0,1);
byear_raw_vector_2_12_DE ~ normal(0,1);
byear_raw_vector_2_12_NDE ~ normal(0,1);
byear_raw_vector_2_14_DE ~ normal(0,1);
byear_raw_vector_2_14_NDE ~ normal(0,1);
byear_raw_vector_2_16_DE ~ normal(0,1);
byear_raw_vector_2_16_NDE ~ normal(0,1);
byear_raw_vector_2_18_DE ~ normal(0,1);
byear_raw_vector_2_18_NDE ~ normal(0,1);
byear_raw_vector_2_41 ~ normal(0,1);
byear_raw_vector_3_2 ~ normal(0,1);
byear_raw_vector_3_4 ~ normal(0,1);
byear_raw_vector_3_8 ~ normal(0,1);
byear_raw_vector_3_20_DE ~ normal(0,1);
byear_raw_vector_3_20_NDE ~ normal(0,1);
byear_raw_vector_3_22_DE ~ normal(0,1);
byear_raw_vector_3_22_NDE ~ normal(0,1);
byear_raw_vector_4_3 ~ normal(0,1);
byear_raw_vector_4_5 ~ normal(0,1);
byear_raw_vector_5_4 ~ normal(0,1);
byear_raw_vector_5_6 ~ normal(0,1);
byear_raw_vector_5_24_DE ~ normal(0,1);
byear_raw_vector_5_24_NDE ~ normal(0,1);
byear_raw_vector_6_5 ~ normal(0,1);
byear_raw_vector_6_7 ~ normal(0,1);
byear_raw_vector_6_26_DE ~ normal(0,1);
byear_raw_vector_6_26_NDE ~ normal(0,1);
byear_raw_vector_7_6 ~ normal(0,1);
byear_raw_vector_7_28_DE ~ normal(0,1);
byear_raw_vector_7_28_NDE ~ normal(0,1);
byear_raw_vector_7_30_DE ~ normal(0,1);
byear_raw_vector_7_30_NDE ~ normal(0,1);
byear_raw_vector_7_42 ~ normal(0,1);
byear_raw_vector_8_3 ~ normal(0,1);
byear_raw_vector_8_9 ~ normal(0,1);
byear_raw_vector_8_32_DE ~ normal(0,1);
byear_raw_vector_8_32_NDE ~ normal(0,1);
byear_raw_vector_9_8 ~ normal(0,1);
byear_raw_vector_9_34_DE ~ normal(0,1);
byear_raw_vector_9_34_NDE ~ normal(0,1);
byear_raw_vector_9_36 ~ normal(0,1);
byear_raw_vector_9_37 ~ normal(0,1);
byear_raw_vector_9_38 ~ normal(0,1);
byear_raw_vector_9_39_DE ~ normal(0,1);
byear_raw_vector_9_39_NDE ~ normal(0,1);
byear_raw_vector_10_2 ~ normal(0,1);
byear_raw_vector_14_2 ~ normal(0,1);
byear_raw_vector_24_5 ~ normal(0,1);
byear_raw_vector_26_6 ~ normal(0,1);
byear_raw_vector_28_7 ~ normal(0,1);
byear_raw_vector_30_7 ~ normal(0,1);
byear_raw_vector_36_9 ~ normal(0,1);
byear_raw_vector_42_7 ~ normal(0,1);

borigin1_matrix_3_4 ~ normal(0,10);
borigin1_matrix_4_3 ~ normal(0,10);
borigin1_matrix_4_5 ~ normal(0,10);
borigin1_matrix_5_4 ~ normal(0,10);
borigin1_matrix_5_6 ~ normal(0,10);
borigin1_matrix_5_24_DE ~ normal(0,10);
borigin1_matrix_5_24_NDE ~ normal(0,10);
borigin1_matrix_6_5 ~ normal(0,10);
borigin1_matrix_6_7 ~ normal(0,10);
borigin1_matrix_6_26_DE ~ normal(0,10);
borigin1_matrix_6_26_NDE ~ normal(0,10);
borigin1_matrix_7_6 ~ normal(0,10);
borigin1_matrix_7_28_DE ~ normal(0,10);
borigin1_matrix_7_28_NDE ~ normal(0,10);
borigin1_matrix_7_30_DE ~ normal(0,10);
borigin1_matrix_7_30_NDE ~ normal(0,10);
borigin1_matrix_7_42 ~ normal(0,10);
borigin1_matrix_24_5 ~ normal(0,10);
borigin1_matrix_26_6 ~ normal(0,10);
borigin1_matrix_28_7 ~ normal(0,10);
borigin1_matrix_30_7 ~ normal(0,10);
borigin1_matrix_42_7 ~ normal(0,10);

borigin2_matrix_3_4 ~ normal(0,10);
borigin2_matrix_4_3 ~ normal(0,10);
borigin2_matrix_4_5 ~ normal(0,10);
borigin2_matrix_5_4 ~ normal(0,10);
borigin2_matrix_5_6 ~ normal(0,10);
borigin2_matrix_5_24_DE ~ normal(0,10);
borigin2_matrix_5_24_NDE ~ normal(0,10);
borigin2_matrix_6_5 ~ normal(0,10);
borigin2_matrix_6_7 ~ normal(0,10);
borigin2_matrix_6_26_DE ~ normal(0,10);
borigin2_matrix_6_26_NDE ~ normal(0,10);
borigin2_matrix_7_6 ~ normal(0,10);
borigin2_matrix_7_28_DE ~ normal(0,10);
borigin2_matrix_7_28_NDE ~ normal(0,10);
borigin2_matrix_7_30_DE ~ normal(0,10);
borigin2_matrix_7_30_NDE ~ normal(0,10);
borigin2_matrix_7_42 ~ normal(0,10);
borigin2_matrix_24_5 ~ normal(0,10);
borigin2_matrix_26_6 ~ normal(0,10);
borigin2_matrix_28_7 ~ normal(0,10);
borigin2_matrix_30_7 ~ normal(0,10);
borigin2_matrix_42_7 ~ normal(0,10);



// Prior on detection efficiency parameters - from the other stan script for detection efficiency
// here, write out all of the parameters for detection efficiency
// twenty terms for intercepts for different eras (configurations of antennas) in the different tributaries
// the outputs from that model will be the first column [,1] containing central tendency (mean)
// and the second column [,2] containing the standard deviation
asotin_alpha1 ~ normal(det_eff_param_posteriors[1,1], det_eff_param_posteriors[1,2]);
asotin_alpha2 ~ normal(det_eff_param_posteriors[2,1], det_eff_param_posteriors[2,2]);
deschutes_alpha1 ~ normal(det_eff_param_posteriors[3,1], det_eff_param_posteriors[3,2]);
entiat_alpha1 ~ normal(det_eff_param_posteriors[4,1], det_eff_param_posteriors[4,2]);
fifteenmile_alpha1 ~ normal(det_eff_param_posteriors[5,1], det_eff_param_posteriors[5,2]);
hood_alpha1 ~ normal(det_eff_param_posteriors[6,1], det_eff_param_posteriors[6,2]);
imnaha_alpha1 ~ normal(det_eff_param_posteriors[7,1], det_eff_param_posteriors[7,2]);
john_day_alpha1 ~ normal(det_eff_param_posteriors[8,1], det_eff_param_posteriors[8,2]);
methow_alpha1 ~ normal(det_eff_param_posteriors[9,1], det_eff_param_posteriors[9,2]);
methow_alpha2 ~ normal(det_eff_param_posteriors[10,1], det_eff_param_posteriors[10,2]);
okanogan_alpha1 ~ normal(det_eff_param_posteriors[11,1], det_eff_param_posteriors[11,2]);
tucannon_alpha1 ~ normal(det_eff_param_posteriors[12,1], det_eff_param_posteriors[12,2]);
tucannon_alpha2 ~ normal(det_eff_param_posteriors[13,1], det_eff_param_posteriors[13,2]);
umatilla_alpha1 ~ normal(det_eff_param_posteriors[14,1], det_eff_param_posteriors[14,2]);
umatilla_alpha2 ~ normal(det_eff_param_posteriors[15,1], det_eff_param_posteriors[15,2]);
walla_walla_alpha1 ~ normal(det_eff_param_posteriors[16,1], det_eff_param_posteriors[16,2]);
walla_walla_alpha2 ~ normal(det_eff_param_posteriors[17,1], det_eff_param_posteriors[17,2]);
walla_walla_alpha3 ~ normal(det_eff_param_posteriors[17,1], det_eff_param_posteriors[17,2]);
walla_walla_alpha4 ~ normal(det_eff_param_posteriors[19,1], det_eff_param_posteriors[19,2]);
wenatchee_alpha1 ~ normal(det_eff_param_posteriors[20,1], det_eff_param_posteriors[20,2]);
yakima_alpha1 ~ normal(det_eff_param_posteriors[21,1], det_eff_param_posteriors[21,2]);

// 14 terms for discharge relationship, one for each tributary
asotin_beta ~ normal(det_eff_param_posteriors[22,1], det_eff_param_posteriors[22,2]);
deschutes_beta ~ normal(det_eff_param_posteriors[23,1], det_eff_param_posteriors[23,2]);
entiat_beta ~ normal(det_eff_param_posteriors[24,1], det_eff_param_posteriors[24,2]);
// // fifteenmile_beta ~ normal(det_eff_param_posteriors[25,1], det_eff_param_posteriors[25,2]);
// // note that fifteenmile and imnaha don't have any discharge data, so they just get intercepts and the beta term is fixed to zero
// // fifteenmile_beta ~ normal(0,0.000001);
hood_beta ~ normal(det_eff_param_posteriors[26,1], det_eff_param_posteriors[26,2]);
// // imnaha_beta ~ normal(det_eff_param_posteriors[27,1], det_eff_param_posteriors[27,2]);
// // imnaha_beta ~ normal(0,0.000001);
john_day_beta ~ normal(det_eff_param_posteriors[28,1], det_eff_param_posteriors[28,2]);
methow_beta ~ normal(det_eff_param_posteriors[29,1], det_eff_param_posteriors[29,2]);
okanogan_beta ~ normal(det_eff_param_posteriors[30,1], det_eff_param_posteriors[30,2]);
tucannon_beta ~ normal(det_eff_param_posteriors[31,1], det_eff_param_posteriors[31,2]);
umatilla_beta ~ normal(det_eff_param_posteriors[32,1], det_eff_param_posteriors[32,2]);
walla_walla_beta ~ normal(det_eff_param_posteriors[33,1], det_eff_param_posteriors[33,2]);
wenatchee_beta ~ normal(det_eff_param_posteriors[34,1], det_eff_param_posteriors[34,2]);
yakima_beta ~ normal(det_eff_param_posteriors[35,1], det_eff_param_posteriors[35,2]);



// PARALLELIZATION EDITS 
// What we will do is modify the incremental log density statement, to bring target into the first loop by individual.
// This should allow us to increment log density across fish, rather than observations within a fish, allowing us to 
// break up the dataset by fish and therefore run different chunks of the dataset in parallel.

  target += reduce_sum(partial_sum_lupmf, y, grainsize, n_ind, max_visits, nmovements, transition_dates, // arguments 1-7
  transition_seasons_vector, cat_X_mat, states_mat, n_obs, // arguments 8-13
  b0_matrix_DE, borigin1_matrix_DE, borigin2_matrix_DE,  // arguments 14-16
  b0_matrix_NDE, borigin1_matrix_NDE, borigin2_matrix_NDE,  // arguments 17-19
  sigma_year_indices, sigma_year_matrix_DE, sigma_year_matrix_NDE,  // arguments 36-38
  byear_raw_parameters_array_DE, byear_raw_parameters_array_NDE, // arguments 39-40
  byear_actual_parameters_array_DE, byear_actual_parameters_array_NDE, // arguments 41-42
  tributary_design_matrices_array, transition_run_years, nyears, run_year_DE_array, det_eff_param_vector); // arguments 43-47

}

