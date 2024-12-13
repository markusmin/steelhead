# 09-model-temperature-plots

# This script takes the output from the stan model runs in 05-stan-runs and
# plots the effects of temperature

# First, need to load in all of the model runs and all of the packages.
source("analysis/analysis/00-load-model-runs.R")

#### Load helpful stuff copied over from 13...R script
extract_DE_parameters <- function(fit, fit_summary){
  # extract b0 as an array
  parameters <- fit_summary$variable
  parameters[grepl("_alpha|_beta" , parameters)] -> param_subset
  
  # add in fifteenmile_beta and imnaha_beta
  DE_params <- c(param_subset[1:24], "fifteenmile_beta", param_subset[25], "imnaha_beta", param_subset[26:length(param_subset)])
  
  
  # arrange all parameter values together in a df
  DE_param_matrix <- matrix(data = 0, nrow = length(DE_params),
                            length(as.matrix(fit[,,1])))
  
  rownames(DE_param_matrix) <- DE_params
  
  
  
  # extract the draws for each and store in a matrix
  for(i in 1:nrow(DE_param_matrix)){
    if (DE_params[i] %in% c("fifteenmile_beta", "imnaha_beta")){
      DE_param_matrix[i, ] <- 0
    } else {
      DE_param_matrix[i, ] <- as.matrix(fit[,,DE_params[i]])
    }
    
  }
  
  return(DE_param_matrix)
}

calculate_weighted_DE <- function(DE_param_matrix, tributary, tributary_design_matrices_array,
                                  covariate_experiences){
  
  ## calculate estimated detection efficiency by year
  # get the index for the tributary state (needs to be mouth since that's where we're correcting for DE)
  tributary_state <- intersect(grep(tributary, model_states, ignore.case = TRUE), grep("Mouth", model_states, ignore.case = TRUE))
  
  tributary_from_state_df <- data.frame(from = c(rep(2,5),3,3,5,6,7,7,8,9,9),
                                        to = grep("Mouth", model_states, ignore.case = FALSE))
  
  tributary_mainstem_state <- subset(tributary_from_state_df, to == tributary_state)$from
  
  # use the state indexing to get the right design matrix
  trib_design_matrix <- tributary_design_matrices_array[,,tributary_state]
  
  # create a matrix to store DE per year, per iter
  niter <- 4000
  DE_matrix <- matrix(nrow = nrow(trib_design_matrix), ncol = niter)
  
  # for each run year, get a confidence interval around detection efficiency by using the different draws
  for (i in 1:nrow(trib_design_matrix)){
    # if there is no intercept term, that means there is no DE correction for that year - so skip and leave as NA
    if(sum(trib_design_matrix[i,1:21]) == 0){
      
      
    } else {
      for (j in 1:niter){
        eta <- sum(trib_design_matrix[i,] * DE_param_matrix[,j])
        DE_matrix[i,j] <- exp(eta)/(1 + exp(eta))
      }
      
    }
    
    
  }
  
  colnames(DE_matrix) <- paste0("iter", 1:niter)
  
  DE_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("year") %>% 
    mutate(year = as.numeric(year)) %>% 
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
    group_by(year) %>% 
    summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) %>% 
    mutate(rear = "wild") -> DE_matrix_long
  
  DE_matrix_long$year_actual <- 2004:2021
  
  ## calculate number of transitions into that tributary by year
  covariate_experiences %>% 
    mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> covariate_experiences
  
  # get a table of counts by run year, into that tributary
  as.data.frame(table(subset(covariate_experiences, state == tributary_mainstem_state & next_state == tributary_state)$year)) %>% 
    dplyr::rename(index = Var1, count = Freq) %>% 
    mutate(index = as.numeric(as.character(index))) -> trib_entries_by_year
  
  # add in the actual year (not index year). 2004 = year 1
  year_indices <- data.frame(index = 1:19, year_actual = 2004:2022)
  trib_entries_by_year %>% 
    left_join(year_indices, by = "index") -> trib_entries_by_year
  
  DE_matrix_long %>% 
    left_join(dplyr::select(trib_entries_by_year, count, year_actual), by = "year_actual") -> DE_by_year
  
  # create a weighted average
  # make sure to drop any years where DE isn't estimated in the model
  DE_by_year %>% 
    filter(!(is.na(`0.5`))) %>% 
    group_by(rear) %>% 
    mutate(weight = count/sum(count, na.rm = T)) %>% 
    summarise(weighted_DE = sum(`0.5`*weight, na.rm = T)) -> weighted_DE
  
  return(list(DE_by_year = DE_by_year, weighted_average = weighted_DE$weighted_DE))
}

### extract all params ###
UCW_DE_param_matrix <- extract_DE_parameters(fit = UCW_fit, fit_summary = UCW_fit_summary)
UCH_DE_param_matrix <- extract_DE_parameters(fit = UCH_fit, fit_summary = UCH_fit_summary)
MCW_DE_param_matrix <- extract_DE_parameters(fit = MCW_fit, fit_summary = MCW_fit_summary)
MCH_DE_param_matrix <- extract_DE_parameters(fit = MCH_fit, fit_summary = MCH_fit_summary)
SRW_DE_param_matrix <- extract_DE_parameters(fit = SRW_fit, fit_summary = SRW_fit_summary)
SRH_DE_param_matrix <- extract_DE_parameters(fit = SRH_fit, fit_summary = SRH_fit_summary)

#### Temperature plot ####

# Here, we will probably want to look at only certain movements (?)
# For example, overshoot, holding in Deschutes, straying
# But we could plot all of them

# We have DPS-wide temperature effects and origin-specific ones, 
# and we have a winter/spring and summer/fall param (temp0 or temp1)

# to capture posterior correlations, we'll need to select that same iteration
# for all parameters that affect movement


# create a list that maps origin numbers (params) to what they actually are
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(grepl("other tributaries", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

# Use the parameter map to index the right effects
origin_param_map <- data.frame(
  natal_origin = natal_origins,
  hatchery = c(NA, NA, NA, NA, 1, NA, 2, # MC
               1,NA,2,3, # UC
               5,NA,1,4,2,3), # SR,
  wild = c(1,3,NA,2,4,6,5, # MC
           1,2,NA,3, # UC
           6,1,2,5,3,4)) # SR

# Add the DPS
DPS_map <- data.frame(
  natal_origin = c("Middle Columbia", "Upper Columbia", "Snake River"),
  hatchery = rep(999, 3),
  wild = rep(999, 3)
)
origin_param_map %>% 
  bind_rows(DPS_map) -> origin_param_map



# extract_covariate_experiences <- function(envir, rear, origin_select){
#   origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
#   for(i in 1:nrow(envir$data$cat_X_mat)){
#     origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
#   }
#   
#   pop_states_dates <- data.frame(state = as.vector(envir$data$y),
#                                  date = as.vector(envir$data$transition_dates),
#                                  year = ceiling(as.vector(envir$data$transition_dates)/365.25)+1,
#                                  origin = rep(origin_vector, ncol(envir$data$y)))
#   # add mainstem dam for each state
#   dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
#                           state = seq(2,9))
#   pop_states_dates %>% 
#     left_join(., dam_index, by = "state") -> pop_states_dates
#   
#   
#   # reformat covariates so that they can be joined
#   as.data.frame(envir$data$spill_window_data) %>% 
#     rownames_to_column("date") %>% 
#     mutate(date = as.numeric(date)) %>% 
#     pivot_longer(cols = -c(date), names_to = "dam", values_to = "spill_window") -> spill_window_long
#   
#   as.data.frame(envir$data$temperature_data) %>% 
#     rownames_to_column("date") %>% 
#     mutate(date = as.numeric(date)) %>% 
#     pivot_longer(cols = -c(date), names_to = "dam", values_to = "temperature") -> temp_long
#   
#   as.data.frame(envir$data$winter_spill_days_data) %>% 
#     rownames_to_column("year") %>% 
#     mutate(year = as.numeric(year)) %>% 
#     pivot_longer(cols = -c(year), names_to = "dam", values_to = "winter_spill") -> spill_days_long
#   
#   
#   # add temperature, spill window, winter spill days
#   pop_states_dates %>% 
#     left_join(., spill_window_long, by = c("date", "dam")) %>% 
#     left_join(., temp_long, by = c("date", "dam")) %>% 
#     left_join(., spill_days_long, by = c("year", "dam")) -> pop_states_dates_covariates
#   
#   # Now, keep only the origin selected
#   if(rear == "wild"){
#     origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
#   } else {
#     origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
#   }
#   
#   subset(pop_states_dates_covariates, origin == origin_numeric) -> origin_states_dates_covariates
#   
#   return(origin_states_dates_covariates)
#   
# }

extract_covariate_experiences <- function(envir, rear, origin_select = NULL){
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  
  
  # for spill days - include the winter post-overshoot vector, which contains
  # info on whether they could have experienced winter spill conditions or not
  # add a new fish_ID column, which is not the same as tag code but will allow us to differentiate between fish
  pop_states_dates <- data.frame(fish_ID = rep(1:length(origin_vector), each = ncol(envir$data$y)),
    state = as.vector(t(envir$data$y)),
                                 date = as.vector(t(envir$data$transition_dates)),
                                 year = ceiling(as.vector(t(envir$data$transition_dates))/365.25)+1,
                                 origin = rep(origin_vector, each = ncol(envir$data$y)))
  
  
  # add mainstem dam for each state
  # WEL and LGR because we have slow v fast there
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                          state = seq(2,9))
  
  pop_states_dates %>% 
    left_join(., dam_index, by = "state") -> pop_states_dates
  
  
  # reformat covariates so that they can be joined
  as.data.frame(envir$data$spill_window_data) %>% 
    rownames_to_column("date") %>% 
    dplyr::rename(LGR = LGR_quick, WEL = WEL_quick) %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(cols = -c(date), names_to = "dam", values_to = "spill_window") -> spill_window_long
  
  as.data.frame(envir$data$temperature_data) %>% 
    rownames_to_column("date") %>% 
    dplyr::rename(LGR = LGR_quick, WEL = WEL_quick) %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(cols = -c(date), names_to = "dam", values_to = "temperature") -> temp_long
  
  as.data.frame(envir$data$winter_spill_days_data) %>% 
    rownames_to_column("year") %>% 
    mutate(year = as.numeric(year)) %>% 
    pivot_longer(cols = -c(year), names_to = "dam", values_to = "winter_spill") -> spill_days_long
  
  
  # add temperature, spill window, winter spill days
  pop_states_dates %>% 
    left_join(., spill_window_long, by = c("date", "dam")) %>% 
    left_join(., temp_long, by = c("date", "dam")) %>% 
    left_join(., spill_days_long, by = c("year", "dam")) -> pop_states_dates_covariates
  
  # drop observations in the loss state and with index 0
  pop_states_dates_covariates %>% 
    filter(!(state %in% c(0,43))) -> pop_states_dates_covariates
  
  # now add winter post-overshoot vector
  pop_states_dates_covariates$winter_post_overshoot_vector <- as.vector(envir$data$winter_post_overshoot_vector)
  
  
  # Now, keep only the origin selected, unless you want the whole DPS
  if (!(is.null(origin_select))){
    if(rear == "wild"){
      origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
    } else {
      origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
    }
    
    subset(pop_states_dates_covariates, origin == origin_numeric) -> origin_states_dates_covariates
  } else {
    # if we leave origin select as null, return the full DPS
    pop_states_dates_covariates -> origin_states_dates_covariates
  }


  return(origin_states_dates_covariates)
  
}



# Tell it the movements for which you want to estimate temperature effects
# movements are formatted as matrix, with column for from and column for to
estimate_temp_effect_UCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                                     date = as.vector(UCW_envir$data$transition_dates))
      spillwindow_data <- UCW_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(UCW_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- UCW_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    

    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter, i] <- exp(b0_array_UCW[from,to,iter] +
                                                       # btemp0_array_UCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_UCW[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_UCW[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_UCW[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_UCW[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_UCW[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_UCW[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_UCW[from,to,iter]*temp_predict[j]*origin2 + 
                                                       # btemp0xorigin3_array_UCW[from,to,iter]*origin3 +
                                                       btemp1xorigin3_array_UCW[from,to,iter]*temp_predict[j]*origin3 +
                                                       borigin1_array_UCW[from,to,iter]*origin1 +
                                                       borigin2_array_UCW[from,to,iter]*origin2 +
                                                       borigin3_array_UCW[from,to,iter]*origin3)/
          sum(exp(b0_array_UCW[from,possible_movements,iter] +
                    # btemp0_array_UCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_UCW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_UCW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_UCW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_UCW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    borigin1_array_UCW[from,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[from,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[from,possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_UCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                                     date = as.vector(UCH_envir$data$transition_dates))
      spillwindow_data <- UCH_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(UCH_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- UCH_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter, i] <- exp(b0_array_UCH[from,to,iter] +
                                                       # btemp0_array_UCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_UCH[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_UCH[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_UCH[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_UCH[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_UCH[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_UCH[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_UCH[from,to,iter]*temp_predict[j]*origin2 + 
                                                       # btemp0xorigin3_array_UCH[from,to,iter]*origin3 +
                                                       btemp1xorigin3_array_UCH[from,to,iter]*temp_predict[j]*origin3 +
                                                       borigin1_array_UCH[from,to,iter]*origin1 +
                                                       borigin2_array_UCH[from,to,iter]*origin2 +
                                                       borigin3_array_UCH[from,to,iter]*origin3)/
          sum(exp(b0_array_UCH[from,possible_movements,iter] +
                    # btemp0_array_UCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_UCH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_UCH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_UCH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_UCH[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    borigin1_array_UCH[from,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[from,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[from,possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_MCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  origin4 = wild_origin_params[4]
  origin5 = wild_origin_params[5]
  origin6 = wild_origin_params[6]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                                     date = as.vector(MCW_envir$data$transition_dates))
      spillwindow_data <- MCW_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(MCW_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- MCW_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter,i] <- exp(b0_array_MCW[from,to,iter] +
                                                       # btemp0_array_MCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_MCW[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_MCW[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_MCW[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_MCW[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_MCW[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_MCW[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_MCW[from,to,iter]*temp_predict[j]*origin2 + 
                                                       # btemp0xorigin3_array_MCW[from,to,iter]*origin3 +
                                                       btemp1xorigin3_array_MCW[from,to,iter]*temp_predict[j]*origin3 +
                                                       # btemp0xorigin4_array_MCW[from,to,iter]*origin4 +
                                                       btemp1xorigin4_array_MCW[from,to,iter]*temp_predict[j]*origin4 +
                                                       # btemp0xorigin5_array_MCW[from,to,iter]*origin5 +
                                                       btemp1xorigin5_array_MCW[from,to,iter]*temp_predict[j]*origin5 +
                                                       # btemp0xorigin6_array_MCW[from,to,iter]*origin6 +
                                                       btemp1xorigin6_array_MCW[from,to,iter]*temp_predict[j]*origin6 +
                                                       borigin1_array_MCW[from,to,iter]*origin1 +
                                                       borigin2_array_MCW[from,to,iter]*origin2 +
                                                       borigin3_array_MCW[from,to,iter]*origin3 +
                                                       borigin4_array_MCW[from,to,iter]*origin4 +
                                                       borigin5_array_MCW[from,to,iter]*origin5 +
                                                       borigin6_array_MCW[from,to,iter]*origin6)/
          sum(exp(b0_array_MCW[from,possible_movements,iter] +
                    # btemp0_array_MCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_MCW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_MCW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_MCW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_MCW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_MCW[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_MCW[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    # btemp0xorigin6_array_MCW[from,possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin6 +
                    borigin1_array_MCW[from,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[from,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[from,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[from,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[from,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[from,possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_MCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                                     date = as.vector(MCH_envir$data$transition_dates))
      spillwindow_data <- MCH_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(MCH_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- MCH_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j,iter,i] <- exp(b0_array_MCH[from,to,iter] +
                                                       # btemp0_array_MCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_MCH[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_MCH[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_MCH[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_MCH[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_MCH[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_MCH[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_MCH[from,to,iter]*temp_predict[j]*origin2 +
                                                       borigin1_array_MCH[from,to,iter]*origin1 +
                                                       borigin2_array_MCH[from,to,iter]*origin2)/
          sum(exp(b0_array_MCH[from,possible_movements,iter] +
                    # btemp0_array_MCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_MCH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_MCH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_MCH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    borigin1_array_MCH[from,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[from,possible_movements,iter]*origin2))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_SRW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  origin4 = wild_origin_params[4]
  origin5 = wild_origin_params[5]
  origin6 = wild_origin_params[6]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
                                     date = as.vector(SRW_envir$data$transition_dates))
      spillwindow_data <- SRW_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(SRW_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- SRW_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
  
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j,iter,i] <- exp(b0_array_SRW[from,to,iter] +
              # btemp0_array_SRW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
              btemp1_array_SRW[from,to,iter]*temp_predict[j] + 
              bspillwindow_array_SRW[from,to,iter]*med_spillwindow +
              bwinterspill_array_SRW[from,to,iter]*med_winterspill +
              # btemp0xorigin1_array_SRW[from,to,iter]*origin1 +
              btemp1xorigin1_array_SRW[from,to,iter]*temp_predict[j]*origin1 +
              # btemp0xorigin2_array_SRW[from,to,iter]*origin2 + 
              btemp1xorigin2_array_SRW[from,to,iter]*temp_predict[j]*origin2 + 
              # btemp0xorigin3_array_SRW[from,to,iter]*origin3 +
              btemp1xorigin3_array_SRW[from,to,iter]*temp_predict[j]*origin3 +
              # btemp0xorigin4_array_SRW[from,to,iter]*origin4 +
              btemp1xorigin4_array_SRW[from,to,iter]*temp_predict[j]*origin4 +
              # btemp0xorigin5_array_SRW[from,to,iter]*origin5 +
              btemp1xorigin5_array_SRW[from,to,iter]*temp_predict[j]*origin5 +
              # btemp0xorigin6_array_SRW[from,to,iter]*origin6 +
              btemp1xorigin6_array_SRW[from,to,iter]*temp_predict[j]*origin6 +
              borigin1_array_SRW[from,to,iter]*origin1 +
              borigin2_array_SRW[from,to,iter]*origin2 +
              borigin3_array_SRW[from,to,iter]*origin3 +
              borigin4_array_SRW[from,to,iter]*origin4 +
              borigin5_array_SRW[from,to,iter]*origin5 +
              borigin6_array_SRW[from,to,iter]*origin6)/
          sum(exp(b0_array_SRW[from,possible_movements,iter] +
                    # btemp0_array_SRW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_SRW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_SRW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_SRW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_SRW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_SRW[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_SRW[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    # btemp0xorigin6_array_SRW[from,possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin6 +
                    borigin1_array_SRW[from,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[from,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[from,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[from,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[from,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[from,possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_SRH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,5)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  origin4 = hatchery_origin_params[4]
  origin5 = hatchery_origin_params[5]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from %in% c(1:9)){
      # get median winter spill for this state
      SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
                                     date = as.vector(SRH_envir$data$transition_dates))
      spillwindow_data <- SRH_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(SRH_states_dates, state == from)$date,from])
      
      # get median spill window for this state
      winterspill_data <- SRH_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from])      
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter,i] <- exp(b0_array_SRH[from,to,iter] +
                                                       # btemp0_array_SRH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_SRH[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_SRH[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_SRH[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_SRH[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_SRH[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_SRH[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_SRH[from,to,iter]*temp_predict[j]*origin2 + 
                                                       # btemp0xorigin3_array_SRH[from,to,iter]*origin3 +
                                                       btemp1xorigin3_array_SRH[from,to,iter]*temp_predict[j]*origin3 +
                                                       # btemp0xorigin4_array_SRH[from,to,iter]*origin4 +
                                                       btemp1xorigin4_array_SRH[from,to,iter]*temp_predict[j]*origin4 +
                                                       # btemp0xorigin5_array_SRH[from,to,iter]*origin5 +
                                                       btemp1xorigin5_array_SRH[from,to,iter]*temp_predict[j]*origin5 +
                                                       borigin1_array_SRH[from,to,iter]*origin1 +
                                                       borigin2_array_SRH[from,to,iter]*origin2 +
                                                       borigin3_array_SRH[from,to,iter]*origin3 +
                                                       borigin4_array_SRH[from,to,iter]*origin4 +
                                                       borigin5_array_SRH[from,to,iter]*origin5)/
          sum(exp(b0_array_SRH[from,possible_movements,iter] +
                    # btemp0_array_SRH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_SRH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_SRH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_SRH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_SRH[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_SRH[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_SRH[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    borigin1_array_SRH[from,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[from,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[from,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[from,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[from,possible_movements,iter]*origin5))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

# plot temp effect function

# create df to index to right dam temp for plotting
dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                        state = seq(2,9))


plot_temp_effect <- function(move_prob_array, from, to, plot_title = NULL){
  temp_move_prob <- as.data.frame(move_prob_array[from, to,,])
  
  colnames(temp_move_prob) <- paste0("iter", 1:niter) 
  temp_predict <- seq(-2,2,length = 100)
  temp_move_prob$temp <- temp_predict
  
  # Add a column with the actual temperatures
  temp_move_prob %>% 
    mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
             window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> temp_move_prob
  
  temp_move_prob %>% 
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
    group_by(temp_actual) %>% 
    summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> temp_move_prob_quantiles
  
  temp_move_prob_plot <- ggplot(temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
    xlab(expression(~"Temperature" ~ ("°C"))) +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(temp_move_prob_plot)
}

# another plot option to compare between hatchery and wild
rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
plot_compare_rear_temp_effect <- function(origin_select,
                                          wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                          wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                          movements_evaluated,
                                          from, to, plot_title = NULL){
  
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  niter <- 4000 # this is the number of draws we have
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])

    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    temp_predict <- seq(-2,2,length = 100)
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      wild_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    } else {
      wild_temp_move_prob %>% 
        mutate(temp_actual = 1) -> wild_temp_move_prob
    }
    
    
    
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    wild_temp_move_prob_quantiles  -> rear_temp_move_prob_quantiles

  } 
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    } else {
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = 1) -> hatchery_temp_move_prob
    }
    
    
    
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    # combine wild and hatchery
    hatchery_temp_move_prob_quantiles -> rear_temp_move_prob_quantiles
    
  } 
  # else run both
  else {
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    temp_predict <- seq(-2,2,length = 100)
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      wild_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    } else {
      wild_temp_move_prob %>% 
        mutate(temp_actual = 1) -> wild_temp_move_prob
    }
    
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> wild_covariate_experiences
    
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict

    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    } else {
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = 1) -> hatchery_temp_move_prob
    }
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> hatchery_covariate_experiences
    
    # combine wild and hatchery
    wild_temp_move_prob_quantiles %>% 
      bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) -> covariate_experiences
    
  }
  
  # convert temp to temp_actual in covariate experiences
  if (from %in% c(1:9)){
    covariate_experiences %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temperature) -> covariate_experiences
  } else {
    covariate_experiences %>% 
      mutate(temp_actual = 1) -> covariate_experiences
  }
  
  
  # Remove any instances for rug plot of fish that would have received temp0 covariate
  covariate_experiences %>% 
    dplyr::rename(date_numeric = date) %>% 
    # keep only jan/feb/mar 
    mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
    mutate(month = month(date)) %>% 
    filter(!(month %in% c(1,2,3,4,5))) -> covariate_experiences
  
  
  rear_temp_move_prob_plot <- ggplot(rear_temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                        color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_rug(data = subset(covariate_experiences, rear == "wild"), aes(x = temp_actual, color = rear), inherit.aes = FALSE,
             sides = "t", length = unit(0.3, "cm")) +
    geom_rug(data = subset(covariate_experiences, rear == "hatchery"), aes(x = temp_actual, color = rear), inherit.aes = FALSE,
             sides = "b", length = unit(0.3, "cm")) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab(expression(~"Temperature" ~ ("°C"))) +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(rear_temp_move_prob_plot)
}


#### Plot rear type comparison overshoot for all ####

### Upper Columbia ###
# Wenatchee
WEN_movements <- data.frame(from = c(5), to = c(6))
WEN_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Wenatchee River", movements = WEN_movements)
WEN_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Wenatchee River", movements = WEN_movements)
# Wenatchee - get covariate experiences
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")


WEN_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Wenatchee River",
                                                            wild_move_prob_array = WEN_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = WEN_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = WEN_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
                                                            movements_evaluated = WEN_movements,
                                                            from = WEN_movements$from, to = WEN_movements$to, plot_title = "Wenatchee - overshoot RRE")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WEN_compare_overshoot_temp.png"), WEN_compare_overshoot_temp, height = 8, width = 8)

# Entiat
ENT_movements <- data.frame(from = c(6), to = c(7))
ENT_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Entiat River", movements = ENT_movements)
ENT_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Entiat River", movements = ENT_movements)

# Entiat - get covariate experiences
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
ENT_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Entiat River")


ENT_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Entiat River",
                                                            wild_move_prob_array = ENT_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = ENT_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = ENT_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = ENT_hatchery_covariate_experiences,
                                                            movements_evaluated = ENT_movements,
                                                            from = ENT_movements$from, to = ENT_movements$to, plot_title = "Entiat - overshoot WEL")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "ENT_compare_overshoot_temp.png"), ENT_compare_overshoot_temp, height = 8, width = 8)






### Middle Columbia ###

# Deschutes River
DES_movements <- data.frame(from = c(2), to = c(3))
DES_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Deschutes River", movements = DES_movements)
DES_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Deschutes River", movements = DES_movements)

# Deschutes - get covariate experiences
DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")
DES_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Deschutes River")


DES_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Deschutes River",
                                                            wild_move_prob_array = DES_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = DES_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = DES_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = DES_hatchery_covariate_experiences,
                                                            movements_evaluated = DES_movements,
                                                            from = DES_movements$from, to = DES_movements$to, plot_title = "Deschutes - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "DES_compare_overshoot_temp.png"), DES_compare_overshoot_temp, height = 8, width = 8)


# John Day River
JDR_movements <- data.frame(from = c(2), to = c(3))
JDR_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "John Day River", movements = JDR_movements)
JDR_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "John Day River", movements = JDR_movements)

# John Day - get covariate experiences
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")
JDR_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "John Day River")


JDR_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "John Day River",
                                                            wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = JDR_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = JDR_hatchery_covariate_experiences,
                                                            movements_evaluated = JDR_movements,
                                                            from = JDR_movements$from, to = JDR_movements$to, plot_title = "John Day - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "JDR_compare_overshoot_temp.png"), JDR_compare_overshoot_temp, height = 8, width = 8)


# Fifteenmile Creek
FIF_movements <- data.frame(from = c(2), to = c(3))
FIF_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Fifteenmile Creek", movements = FIF_movements)
FIF_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = FIF_movements)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
FIF_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Fifteenmile Creek")


FIF_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Fifteenmile Creek",
                                                            wild_move_prob_array = FIF_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = FIF_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = FIF_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = FIF_hatchery_covariate_experiences,
                                                            movements_evaluated = FIF_movements,
                                                            from = FIF_movements$from, to = FIF_movements$to, plot_title = "Fifteenmile - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "FIF_compare_overshoot_temp.png"), FIF_compare_overshoot_temp, height = 8, width = 8)


# Fifteenmile Creek - movement into the Deschutes
into_deschutes <- data.frame(from = c(2), to = c(10))
FIF_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Fifteenmile Creek", movements = into_deschutes)
FIF_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = into_deschutes)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
FIF_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Fifteenmile Creek")


FIF_compare_move_deschutes_temp <- plot_compare_rear_temp_effect(origin_select = "Fifteenmile Creek",
                                                            wild_move_prob_array = FIF_wild_temp_move_deschutes_prob_array,
                                                            hatchery_move_prob_array = FIF_hatchery_temp_move_deschutes_prob_array,
                                                            wild_covariate_experiences = FIF_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = FIF_hatchery_covariate_experiences,
                                                            movements_evaluated = into_deschutes,
                                                            from = into_deschutes$from, to = into_deschutes$to, plot_title = "Fifteenmile - move into Deschutes")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "FIF_compare_move_deschutes_temp.png"), FIF_compare_move_deschutes_temp, height = 8, width = 8)


# Fifteenmile Creek - movement out of the Deschutes
out_deschutes <- data.frame(from = c(10), to = c(2))
FIF_hatchery_temp_move_out_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Fifteenmile Creek", movements = out_deschutes)
FIF_wild_temp_move_out_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = out_deschutes)

# create simple plot of probabilities of moving out
# pick the first row; all rows are the same because there is no impact of temperature
hist(FIF_wild_temp_move_out_deschutes_prob_array[1,,1])
# this is in line with the data. The problem seems to be an overestimation of movement into the Deschutes.
# Data says that 11.4% of fish move into the Deschutes from BON to MCN; based on 
# relationship with temperature, I'd estimate that about 40% are moving into Deschutes.
# but why?



# Umatilla River
UMA_movements <- data.frame(from = c(2), to = c(3))
UMA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = UMA_movements)
UMA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = UMA_movements)

# Umatilla - get covariate experiences
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")


UMA_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Umatilla River",
                                                            wild_move_prob_array = UMA_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = UMA_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                            movements_evaluated = UMA_movements,
                                                            from = UMA_movements$from, to = UMA_movements$to, plot_title = "Umatilla - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "UMA_compare_overshoot_temp.png"), UMA_compare_overshoot_temp, height = 8, width = 8)


# Umatilla River - movement into the Deschutes
into_deschutes <- data.frame(from = c(2), to = c(10))
UMA_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = into_deschutes)
UMA_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = into_deschutes)


UMA_compare_move_deschutes_temp <- plot_compare_rear_temp_effect(origin_select = "Umatilla River",
                                                                 wild_move_prob_array = UMA_wild_temp_move_deschutes_prob_array,
                                                                 hatchery_move_prob_array = UMA_hatchery_temp_move_deschutes_prob_array,
                                                                 wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                 hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                 movements_evaluated = into_deschutes,
                                                                 from = into_deschutes$from, to = into_deschutes$to, plot_title = "Umatilla - move into Deschutes")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "UMA_compare_move_deschutes_temp.png"), UMA_compare_move_deschutes_temp, height = 8, width = 8)


# Umatilla River - movement out of the Deschutes
out_deschutes <- data.frame(from = c(10), to = c(2))
UMA_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = out_deschutes)
UMA_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = out_deschutes)

hist(UMA_wild_temp_move_deschutes_prob_array[1,,1])
hist(UMA_hatchery_temp_move_deschutes_prob_array[1,,1])

# Yakima River
YAK_movements <- data.frame(from = c(3), to = c(4))
YAK_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Yakima River", movements = YAK_movements)
YAK_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = YAK_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
YAK_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Yakima River")


YAK_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Yakima River",
                                                            wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = YAK_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = YAK_hatchery_covariate_experiences,
                                                            movements_evaluated = YAK_movements,
                                                            from = YAK_movements$from, to = YAK_movements$to, plot_title = "Yakima - overshoot PRA")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "YAK_compare_overshoot_PRA_temp.png"), YAK_compare_overshoot_temp, height = 8, width = 8)

YAK_movements <- data.frame(from = c(3), to = c(8))
YAK_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Yakima River", movements = YAK_movements)
YAK_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = YAK_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
YAK_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Yakima River")


YAK_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Yakima River",
                                                            wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = YAK_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = YAK_hatchery_covariate_experiences,
                                                            movements_evaluated = YAK_movements,
                                                            from = YAK_movements$from, to = YAK_movements$to, plot_title = "Yakima - overshoot ICH")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "YAK_compare_overshoot_ICH_temp.png"), YAK_compare_overshoot_temp, height = 8, width = 8)


# Walla Walla River
WAWA_movements <- data.frame(from = c(3), to = c(4))
WAWA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")


WAWA_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Walla Walla River",
                                                            wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                            movements_evaluated = WAWA_movements,
                                                            from = WAWA_movements$from, to = WAWA_movements$to, plot_title = "Walla Walla - overshoot PRA")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WAWA_compare_overshoot_PRA_temp.png"), WAWA_compare_overshoot_temp, height = 8, width = 8)


WAWA_movements <- data.frame(from = c(3), to = c(8))
WAWA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")


WAWA_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Walla Walla River",
                                                            wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                            movements_evaluated = WAWA_movements,
                                                             from = WAWA_movements$from, to = WAWA_movements$to, plot_title = "Walla Walla - overshoot ICH")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WAWA_compare_overshoot_ICH_temp.png"), WAWA_compare_overshoot_temp, height = 8, width = 8)








### Snake River ###

## Tucannon River
TUC_movements <- data.frame(from = c(8), to = c(9))
TUC_hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = "Tucannon River", movements = TUC_movements)
TUC_wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = "Tucannon River", movements = TUC_movements)

# Tucannon - get covariate experiences
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")


TUC_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Tucannon River",
                                                            wild_move_prob_array = TUC_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = TUC_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = TUC_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
                                                            movements_evaluated = TUC_movements,
                                                            from = TUC_movements$from, to = TUC_movements$to, plot_title = "Tucannon - overshoot LGR")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "TUC_compare_overshoot_temp.png"), TUC_compare_overshoot_temp, height = 8, width = 8)



# Snake only has one overshoot - for Tucannon


#### New version of temperature relationship plot - now with data! ####

# for testing
into_deschutes <- data.frame(from = c(2), to = c(10))
UMA_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = into_deschutes)
UMA_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = into_deschutes)
# Umatilla - get covariate experiences
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")

origin_select <- "Umatilla River"
wild_move_prob_array = UMA_wild_temp_move_deschutes_prob_array
hatchery_move_prob_array = UMA_hatchery_temp_move_deschutes_prob_array
wild_covariate_experiences = UMA_wild_covariate_experiences
hatchery_covariate_experiences = UMA_hatchery_covariate_experiences
movements_evaluated = into_deschutes
from <- 2
to <- 10
plot_title = NULL

plot_compare_rear_temp_effect_fit_to_data <- function(origin_select,
                                          wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                          wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                          movements_evaluated,
                                          from, to, plot_title = NULL){
  
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  niter <- 4000 # this is the number of draws we have
  
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    # modify covariate experiences to show what the next state is
    wild_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> wild_covariate_experiences
    
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    temp_predict <- seq(-2,2,length = 100)
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      wild_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    } else {
      wild_temp_move_prob %>% 
        mutate(temp_actual = 1) -> wild_temp_move_prob
    }
    
    
    
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    wild_temp_move_prob_quantiles  -> rear_temp_move_prob_quantiles
    
  } else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    # If wild is NA, run hatchery only
    
    
    # modify covariate experiences to show what the next state is
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    } else {
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = 1) -> hatchery_temp_move_prob
    }
    
    
    
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    # combine wild and hatchery
    hatchery_temp_move_prob_quantiles -> rear_temp_move_prob_quantiles
    
  } else { # else run both
    
    # modify covariate experiences to show what the next state is
    wild_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> wild_covariate_experiences
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    temp_predict <- seq(-2,2,length = 100)
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      wild_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    } else {
      wild_temp_move_prob %>% 
        mutate(temp_actual = 1) -> wild_temp_move_prob
    }
    
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> wild_covariate_experiences
    
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    } else {
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = 1) -> hatchery_temp_move_prob
    }
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> hatchery_covariate_experiences
    
    # combine wild and hatchery
    wild_temp_move_prob_quantiles %>% 
      bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) -> covariate_experiences
    
  }
  
  # convert temp to temp_actual in covariate experiences
  if (from %in% c(1:9)){
    covariate_experiences %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temperature) -> covariate_experiences
  } else {
    covariate_experiences %>% 
      mutate(temp_actual = 1) -> covariate_experiences
  }
  
  
  # Remove any instances for rug plot of fish that would have received temp0 covariate
  covariate_experiences %>% 
    dplyr::rename(date_numeric = date) %>% 
    # keep only jan/feb/mar 
    mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
    mutate(month = month(date)) %>% 
    filter(!(month %in% c(1,2,3,4,5))) -> covariate_experiences
  
  # create a new DF to show bubbles of points
  covariate_experiences %>% 
    mutate(temp_round = round(temp_actual, 0)) %>% 
    group_by(temp_round, rear) %>% 
    summarise(total = n()) -> temp_rear_total_counts
  
  covariate_experiences %>% 
    mutate(temp_round = round(temp_actual, 0)) %>% 
    group_by(temp_round, rear, next_state) %>% 
    summarise(count = n()) %>% 
    ungroup() %>% 
    complete(temp_round, nesting(rear, next_state), fill = list(count = 0)) -> temp_rear_move_counts
  
  temp_rear_move_counts %>% 
    left_join(temp_rear_total_counts, by = c("temp_round", "rear")) %>% 
    mutate(prop = count/total) -> temp_move_props
  
  rear_temp_move_prob_plot <- ggplot(rear_temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                        color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_point(data = subset(temp_move_props, next_state == to & rear == "wild"), 
               aes(x = temp_round, y = prop, color = rear, size = total),
               inherit.aes = FALSE) +
    geom_point(data = subset(temp_move_props, next_state == to & rear == "hatchery"), 
               aes(x = temp_round, y = prop, color = rear, size = total),
               inherit.aes = FALSE) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0),
                       sec.axis = sec_axis(~ ., name = "Proportion of movements")) +
    scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab(expression(~"Temperature" ~ ("°C"))) +
    ylab("Movement probability") +
    ggtitle(plot_title) +
    # scale_size(range=c(2,8),breaks=c(10, 50, 100, 250, 500, Inf),labels=c("1-9","10-49","50-99","100-249","250-499","500+"),
    scale_size(range=c(0.1,8),breaks=c(1, 10, 50, 100, 250, 500),
               limits = c(0, 1000), guide="legend")
  
  return(rear_temp_move_prob_plot)
}

#### Plot rear type comparison overshoot for all, with data plotted as well ####

### Upper Columbia ###
# Wenatchee
WEN_movements <- data.frame(from = c(5), to = c(6))
WEN_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Wenatchee River", movements = WEN_movements)
WEN_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Wenatchee River", movements = WEN_movements)
# Wenatchee - get covariate experiences
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")


WEN_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Wenatchee River",
                                                                        wild_move_prob_array = WEN_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = WEN_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = WEN_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
                                                                        movements_evaluated = WEN_movements,
                                                                        from = WEN_movements$from, to = WEN_movements$to, plot_title = "Wenatchee - overshoot RRE")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WEN_compare_overshoot_temp_fit_to_data.png"), WEN_compare_overshoot_temp, height = 8, width = 8)

# Entiat
ENT_movements <- data.frame(from = c(6), to = c(7))
ENT_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Entiat River", movements = ENT_movements)
ENT_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Entiat River", movements = ENT_movements)

# Entiat - get covariate experiences
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
ENT_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Entiat River")


ENT_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Entiat River",
                                                                        wild_move_prob_array = ENT_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = ENT_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = ENT_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = ENT_hatchery_covariate_experiences,
                                                                        movements_evaluated = ENT_movements,
                                                                        from = ENT_movements$from, to = ENT_movements$to, plot_title = "Entiat - overshoot WEL")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "ENT_compare_overshoot_temp_fit_to_data.png"), ENT_compare_overshoot_temp, height = 8, width = 8)






### Middle Columbia ###

# Deschutes River
DES_movements <- data.frame(from = c(2), to = c(3))
DES_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Deschutes River", movements = DES_movements)
DES_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Deschutes River", movements = DES_movements)

# Deschutes - get covariate experiences
DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")
DES_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Deschutes River")


DES_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Deschutes River",
                                                                        wild_move_prob_array = DES_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = DES_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = DES_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = DES_hatchery_covariate_experiences,
                                                                        movements_evaluated = DES_movements,
                                                                        from = DES_movements$from, to = DES_movements$to, plot_title = "Deschutes - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "DES_compare_overshoot_temp_fit_to_data.png"), DES_compare_overshoot_temp, height = 8, width = 8)


# John Day River
JDR_movements <- data.frame(from = c(2), to = c(3))
JDR_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "John Day River", movements = JDR_movements)
JDR_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "John Day River", movements = JDR_movements)

# John Day - get covariate experiences
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")
JDR_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "John Day River")


JDR_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "John Day River",
                                                                        wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = JDR_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = JDR_hatchery_covariate_experiences,
                                                                        movements_evaluated = JDR_movements,
                                                                        from = JDR_movements$from, to = JDR_movements$to, plot_title = "John Day - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "JDR_compare_overshoot_temp_fit_to_data.png"), JDR_compare_overshoot_temp, height = 8, width = 8)


# Fifteenmile Creek
FIF_movements <- data.frame(from = c(2), to = c(3))
FIF_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Fifteenmile Creek", movements = FIF_movements)
FIF_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = FIF_movements)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
FIF_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Fifteenmile Creek")


FIF_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Fifteenmile Creek",
                                                                        wild_move_prob_array = FIF_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = FIF_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = FIF_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = FIF_hatchery_covariate_experiences,
                                                                        movements_evaluated = FIF_movements,
                                                                        from = FIF_movements$from, to = FIF_movements$to, plot_title = "Fifteenmile - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "FIF_compare_overshoot_temp_fit_to_data.png"), FIF_compare_overshoot_temp, height = 8, width = 8)


# Fifteenmile Creek - movement into the Deschutes
into_deschutes <- data.frame(from = c(2), to = c(10))
FIF_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Fifteenmile Creek", movements = into_deschutes)
FIF_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = into_deschutes)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
FIF_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Fifteenmile Creek")


FIF_compare_move_deschutes_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Fifteenmile Creek",
                                                                             wild_move_prob_array = FIF_wild_temp_move_deschutes_prob_array,
                                                                             hatchery_move_prob_array = FIF_hatchery_temp_move_deschutes_prob_array,
                                                                             wild_covariate_experiences = FIF_wild_covariate_experiences,
                                                                             hatchery_covariate_experiences = FIF_hatchery_covariate_experiences,
                                                                             movements_evaluated = into_deschutes,
                                                                             from = into_deschutes$from, to = into_deschutes$to, plot_title = "Fifteenmile - move into Deschutes")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "FIF_compare_move_deschutes_temp_fit_to_data.png"), FIF_compare_move_deschutes_temp, height = 8, width = 8)


# Fifteenmile Creek - movement out of the Deschutes
out_deschutes <- data.frame(from = c(10), to = c(2))
FIF_hatchery_temp_move_out_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Fifteenmile Creek", movements = out_deschutes)
FIF_wild_temp_move_out_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = out_deschutes)

# create simple plot of probabilities of moving out
# pick the first row; all rows are the same because there is no impact of temperature
hist(FIF_wild_temp_move_out_deschutes_prob_array[1,,1])
# this is in line with the data. The problem seems to be an overestimation of movement into the Deschutes.
# Data says that 11.4% of fish move into the Deschutes from BON to MCN; based on 
# relationship with temperature, I'd estimate that about 40% are moving into Deschutes.
# but why?



# Umatilla River
UMA_movements <- data.frame(from = c(2), to = c(3))
UMA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = UMA_movements)
UMA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = UMA_movements)

# Umatilla - get covariate experiences
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")


UMA_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Umatilla River",
                                                                        wild_move_prob_array = UMA_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = UMA_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                        movements_evaluated = UMA_movements,
                                                                        from = UMA_movements$from, to = UMA_movements$to, plot_title = "Umatilla - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "UMA_compare_overshoot_temp_fit_to_data.png"), UMA_compare_overshoot_temp, height = 8, width = 8)


# Umatilla River - movement into the Deschutes
into_deschutes <- data.frame(from = c(2), to = c(10))
UMA_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = into_deschutes)
UMA_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = into_deschutes)


UMA_compare_move_deschutes_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Umatilla River",
                                                                             wild_move_prob_array = UMA_wild_temp_move_deschutes_prob_array,
                                                                             hatchery_move_prob_array = UMA_hatchery_temp_move_deschutes_prob_array,
                                                                             wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                             hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                             movements_evaluated = into_deschutes,
                                                                             from = into_deschutes$from, to = into_deschutes$to, plot_title = "Umatilla - move into Deschutes")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "UMA_compare_move_deschutes_temp_fit_to_data.png"), UMA_compare_move_deschutes_temp, height = 8, width = 8)


# Umatilla River - movement out of the Deschutes
out_deschutes <- data.frame(from = c(10), to = c(2))
UMA_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = out_deschutes)
UMA_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = out_deschutes)

hist(UMA_wild_temp_move_deschutes_prob_array[1,,1])
hist(UMA_hatchery_temp_move_deschutes_prob_array[1,,1])

# Yakima River
YAK_movements <- data.frame(from = c(3), to = c(4))
YAK_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Yakima River", movements = YAK_movements)
YAK_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = YAK_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
YAK_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Yakima River")


YAK_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Yakima River",
                                                                        wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = YAK_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = YAK_hatchery_covariate_experiences,
                                                                        movements_evaluated = YAK_movements,
                                                                        from = YAK_movements$from, to = YAK_movements$to, plot_title = "Yakima - overshoot PRA")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "YAK_compare_overshoot_PRA_temp_fit_to_data.png"), YAK_compare_overshoot_temp, height = 8, width = 8)

YAK_movements <- data.frame(from = c(3), to = c(8))
YAK_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Yakima River", movements = YAK_movements)
YAK_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = YAK_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
YAK_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Yakima River")


YAK_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Yakima River",
                                                                        wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = YAK_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = YAK_hatchery_covariate_experiences,
                                                                        movements_evaluated = YAK_movements,
                                                                        from = YAK_movements$from, to = YAK_movements$to, plot_title = "Yakima - overshoot ICH")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "YAK_compare_overshoot_ICH_temp_fit_to_data.png"), YAK_compare_overshoot_temp, height = 8, width = 8)


# Walla Walla River
WAWA_movements <- data.frame(from = c(3), to = c(4))
WAWA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")


WAWA_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Walla Walla River",
                                                                         wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                                         hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                                         wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                                         hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                         movements_evaluated = WAWA_movements,
                                                                         from = WAWA_movements$from, to = WAWA_movements$to, plot_title = "Walla Walla - overshoot PRA")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WAWA_compare_overshoot_PRA_temp_fit_to_data.png"), WAWA_compare_overshoot_temp, height = 8, width = 8)


WAWA_movements <- data.frame(from = c(3), to = c(8))
WAWA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")


WAWA_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Walla Walla River",
                                                                         wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                                         hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                                         wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                                         hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                         movements_evaluated = WAWA_movements,
                                                                         from = WAWA_movements$from, to = WAWA_movements$to, plot_title = "Walla Walla - overshoot ICH")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WAWA_compare_overshoot_ICH_temp_fit_to_data.png"), WAWA_compare_overshoot_temp, height = 8, width = 8)








### Snake River ###

## Tucannon River
TUC_movements <- data.frame(from = c(8), to = c(9))
TUC_hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = "Tucannon River", movements = TUC_movements)
TUC_wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = "Tucannon River", movements = TUC_movements)

# Tucannon - get covariate experiences
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")


TUC_compare_overshoot_temp <- plot_compare_rear_temp_effect_fit_to_data(origin_select = "Tucannon River",
                                                                        wild_move_prob_array = TUC_wild_temp_move_prob_array,
                                                                        hatchery_move_prob_array = TUC_hatchery_temp_move_prob_array,
                                                                        wild_covariate_experiences = TUC_wild_covariate_experiences,
                                                                        hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
                                                                        movements_evaluated = TUC_movements,
                                                                        from = TUC_movements$from, to = TUC_movements$to, plot_title = "Tucannon - overshoot LGR")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "TUC_compare_overshoot_temp_fit_to_data.png"), TUC_compare_overshoot_temp, height = 8, width = 8)



#### New version of temperature relationship plot - now with data and a detection efficiency correction for data! ####
# Before you run this, you will have to source the detection efficiency script to get those functions
# to plot DE
# source(here::here("analysis", "analysis", "13-model-detection-efficiency-plots.R"))


# for testing
into_deschutes <- data.frame(from = c(2), to = c(10))
UMA_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = into_deschutes)
UMA_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = into_deschutes)
# Umatilla - get covariate experiences
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")
wild_covariate_experiences <- UMA_wild_covariate_experiences
hatchery_covariate_experiences <- UMA_hatchery_covariate_experiences

origin_select <- "Umatilla River"
wild_move_prob_array = UMA_wild_temp_move_deschutes_prob_array
hatchery_move_prob_array = UMA_hatchery_temp_move_deschutes_prob_array
wild_covariate_experiences = UMA_wild_covariate_experiences
hatchery_covariate_experiences = UMA_hatchery_covariate_experiences
movements_evaluated = into_deschutes
uma_hatchery_des_average_DE <- calculate_weighted_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array,
                                                     covariate_experiences = UMA_hatchery_covariate_experiences)
uma_wild_des_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                                 covariate_experiences = UMA_wild_covariate_experiences)
wild_DE_correction = uma_wild_des_average_DE
hatchery_DE_correction = uma_hatchery_des_average_DE
from <- 2
to <- 10
plot_title = NULL


plot_compare_rear_temp_effect_fit_to_data_DE <- function(origin_select,
                                                         wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                                         wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                                         wild_DE_correction = NULL, hatchery_DE_correction = NULL,
                                                         movements_evaluated,
                                                         from, to, plot_title = NULL){
  
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  niter <- 4000 # this is the number of draws we have
  
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    # modify covariate experiences to show what the next state is
    wild_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> wild_covariate_experiences
    
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    temp_predict <- seq(-2,2,length = 100)
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      wild_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    } else {
      wild_temp_move_prob %>% 
        mutate(temp_actual = 1) -> wild_temp_move_prob
    }
    
    
    
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    wild_temp_move_prob_quantiles  -> rear_temp_move_prob_quantiles
    
  } else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    # If wild is NA, run hatchery only
    
    
    # modify covariate experiences to show what the next state is
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    } else {
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = 1) -> hatchery_temp_move_prob
    }
    
    
    
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    # combine wild and hatchery
    hatchery_temp_move_prob_quantiles -> rear_temp_move_prob_quantiles
    
  } else { # else run both
    
    # modify covariate experiences to show what the next state is
    wild_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> wild_covariate_experiences
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    temp_predict <- seq(-2,2,length = 100)
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      wild_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    } else {
      wild_temp_move_prob %>% 
        mutate(temp_actual = 1) -> wild_temp_move_prob
    }
    
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> wild_covariate_experiences
    
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(1:9)){
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    } else {
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = 1) -> hatchery_temp_move_prob
    }
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> hatchery_covariate_experiences
    
    # combine wild and hatchery
    wild_temp_move_prob_quantiles %>% 
      bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) -> covariate_experiences
    
  }
  
  # convert temp to temp_actual in covariate experiences
  if (from %in% c(2:9)){
    covariate_experiences %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temperature) -> covariate_experiences
  } else {
    covariate_experiences %>% 
      mutate(temp_actual = 1) -> covariate_experiences
  }
  
  
  # Remove any instances for rug plot of fish that would have received temp0 covariate
  covariate_experiences %>% 
    dplyr::rename(date_numeric = date) %>% 
    # keep only jan/feb/mar 
    mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
    mutate(month = month(date)) %>% 
    filter(!(month %in% c(1,2,3,4,5))) -> covariate_experiences
  
  # create a new DF to show bubbles of points
  covariate_experiences %>% 
    mutate(temp_round = round(temp_actual, 0)) %>% 
    group_by(temp_round, rear) %>% 
    summarise(total = n()) -> temp_rear_total_counts
  
  covariate_experiences %>% 
    mutate(temp_round = round(temp_actual, 0)) %>% 
    group_by(temp_round, rear, next_state) %>% 
    summarise(count = n()) %>% 
    ungroup() %>% 
    complete(temp_round, nesting(rear, next_state), fill = list(count = 0)) -> temp_rear_move_counts
  
  temp_rear_move_counts %>% 
    left_join(temp_rear_total_counts, by = c("temp_round", "rear")) %>% 
    mutate(prop = count/total) -> temp_move_props
  
  # correct for DE
  DE_correction_df <- data.frame(rear = c("wild", "hatchery"),
                                 DE = c(wild_DE_correction$weighted_average, hatchery_DE_correction$weighted_average))
  
  rear_temp_move_prob_quantiles %>% 
    left_join(DE_correction_df, by = "rear") %>% 
    mutate(`0.5*DE` = `0.5` * DE,
           `0.025*DE` = `0.025` * DE,
           `0.975*DE` = `0.975` * DE) -> rear_temp_move_prob_quantiles
  
  # make a long form of this
  rear_temp_move_prob_quantiles %>% 
    dplyr::select(temp_actual, `0.025`, rear, DE, `0.025*DE`) %>% 
    pivot_longer(cols = c(`0.025`, `0.025*DE`), values_to = "0.025", names_to = "Model Fit") %>% 
    mutate(`Model Fit` = ifelse(`Model Fit` == "0.025", "estimated", "observed")) -> rear_temp_move_prob_quantiles_0.025
  
  rear_temp_move_prob_quantiles %>% 
    dplyr::select(temp_actual, `0.5`, rear, DE, `0.5*DE`) %>% 
    pivot_longer(cols = c(`0.5`, `0.5*DE`), values_to = "0.5", names_to = "Model Fit") %>% 
    mutate(`Model Fit` = ifelse(`Model Fit` == "0.5", "estimated", "observed")) -> rear_temp_move_prob_quantiles_0.5
  
  rear_temp_move_prob_quantiles %>% 
    dplyr::select(temp_actual, `0.975`, rear, DE, `0.975*DE`) %>% 
    pivot_longer(cols = c(`0.975`, `0.975*DE`), values_to = "0.975", names_to = "Model Fit") %>% 
    mutate(`Model Fit` = ifelse(`Model Fit` == "0.975", "estimated", "observed")) -> rear_temp_move_prob_quantiles_0.975
  
  rear_temp_move_prob_quantiles_0.025 %>% 
    left_join(., rear_temp_move_prob_quantiles_0.5, by = c("temp_actual", "rear", "Model Fit", "DE")) %>% 
    left_join(., rear_temp_move_prob_quantiles_0.975, by = c("temp_actual", "rear", "Model Fit", "DE")) -> rear_temp_move_prob_quantiles_long
  
  rear_temp_move_prob_plot <- ggplot(rear_temp_move_prob_quantiles_long, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                             color = rear, fill = rear, lty = `Model Fit`, alpha = `Model Fit`)) +
    # the "actual" value
    geom_line(alpha = 1) +
    geom_ribbon(color = NA) +
    geom_point(data = subset(temp_move_props, next_state == to & rear == "wild"), 
               aes(x = temp_round, y = prop, color = rear, size = total),
               inherit.aes = FALSE) +
    geom_point(data = subset(temp_move_props, next_state == to & rear == "hatchery"), 
               aes(x = temp_round, y = prop, color = rear, size = total),
               inherit.aes = FALSE) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0),
                       sec.axis = sec_axis(~ ., name = "Proportion of movements")) +
    scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    scale_alpha_manual(values = c(0.15, 0.35)) +
    xlab(expression(~"Temperature" ~ ("°C"))) +
    ylab("Movement probability") +
    ggtitle(plot_title) +
    # scale_size(range=c(2,8),breaks=c(10, 50, 100, 250, 500, Inf),labels=c("1-9","10-49","50-99","100-249","250-499","500+"),
    scale_size(range=c(0.1,8),breaks=c(1, 10, 50, 100, 250, 500),
               limits = c(0, 1000), guide="legend")
  
  return(rear_temp_move_prob_plot)
}

# Umatilla River - movement into the Deschutes
into_deschutes <- data.frame(from = c(2), to = c(10))
UMA_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = into_deschutes)
UMA_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = into_deschutes)
uma_hatchery_des_average_DE <- calculate_weighted_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array,
                                                     covariate_experiences = UMA_hatchery_covariate_experiences)
uma_wild_des_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                                 covariate_experiences = UMA_wild_covariate_experiences)

UMA_compare_move_deschutes_temp_DE <- plot_compare_rear_temp_effect_fit_to_data_DE(origin_select = "Umatilla River",
                                                                             wild_move_prob_array = UMA_wild_temp_move_deschutes_prob_array,
                                                                             hatchery_move_prob_array = UMA_hatchery_temp_move_deschutes_prob_array,
                                                                             wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                             hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                             wild_DE_correction = uma_wild_des_average_DE,
                                                                             hatchery_DE_correction = uma_hatchery_des_average_DE,
                                                                             movements_evaluated = into_deschutes,
                                                                             from = into_deschutes$from, to = into_deschutes$to, plot_title = "Umatilla - move into Deschutes")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "UMA_compare_move_deschutes_temp_fit_to_data_DE.png"), UMA_compare_move_deschutes_temp_DE, height = 8, width = 8)

# Fifteenmile Creek - movement into the Deschutes
into_deschutes <- data.frame(from = c(2), to = c(10))
FIF_hatchery_temp_move_deschutes_prob_array <- estimate_temp_effect_MCH(origin_select = "Fifteenmile Creek", movements = into_deschutes)
FIF_wild_temp_move_deschutes_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = into_deschutes)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
FIF_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Fifteenmile Creek")
# Get DE
FIF_wild_des_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                                 covariate_experiences = FIF_wild_covariate_experiences)


FIF_compare_move_deschutes_temp_DE <- plot_compare_rear_temp_effect_fit_to_data_DE(origin_select = "Fifteenmile Creek",
                                                                                wild_move_prob_array = FIF_wild_temp_move_deschutes_prob_array,
                                                                                hatchery_move_prob_array = FIF_hatchery_temp_move_deschutes_prob_array,
                                                                                wild_covariate_experiences = FIF_wild_covariate_experiences,
                                                                                hatchery_covariate_experiences = FIF_hatchery_covariate_experiences,
                                                                                movements_evaluated = into_deschutes,
                                                                                wild_DE_correction = FIF_wild_des_average_DE,
                                                                                from = into_deschutes$from, to = into_deschutes$to, plot_title = "Fifteenmile - move into Deschutes")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "FIF_compare_move_deschutes_temp_fit_to_data_DE.png"), FIF_compare_move_deschutes_temp_DE, height = 8, width = 8)


#### Plot fit to data, get goodness of fit ####


plot_GOF_compare_rear_temp_effect_fit_to_data_DE <- function(origin_select,
                                                         wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                                         wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                                         wild_DE_correction = NULL, hatchery_DE_correction = NULL,
                                                         movements_evaluated,
                                                         from, to, plot_title = NULL){
  # create df to index to right dam temp joining; here we need to change names of 
  # WEL and LGR because we have slow v fast there
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL_quick", "ICH", "LGR_quick"),
                          state = seq(2,9))
  
  
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  niter <- 4000 # this is the number of draws we have
  
  # create vector of temps to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  # If origin_select is null - then we run the whole DPS comparison. And every DPS
  # has both hatchery and wild
  
  
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    # modify covariate experiences to show what the next state is
    wild_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> wild_covariate_experiences
    
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(2:9)){
      wild_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    } else {
      wild_temp_move_prob %>% 
        mutate(temp_actual = 1) -> wild_temp_move_prob
    }
    
    
    
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    wild_temp_move_prob_quantiles  -> rear_temp_move_prob_quantiles
    
  } else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    # If wild is NA, run hatchery only
    
    
    # modify covariate experiences to show what the next state is
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(2:9)){
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    } else {
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = 1) -> hatchery_temp_move_prob
    }
    
    
    
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    # combine wild and hatchery
    hatchery_temp_move_prob_quantiles -> rear_temp_move_prob_quantiles
    
  } else { # else run both
    
    # modify covariate experiences to show what the next state is
    wild_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> wild_covariate_experiences
    hatchery_covariate_experiences %>% 
      mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> hatchery_covariate_experiences
    
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(2:9)){
      wild_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    } else {
      wild_temp_move_prob %>% 
        mutate(temp_actual = 1) -> wild_temp_move_prob
    }
    
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> wild_covariate_experiences
    
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    if (from %in% c(2:9)){
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
                 window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    } else {
      hatchery_temp_move_prob %>% 
        mutate(temp_actual = 1) -> hatchery_temp_move_prob
    }
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> hatchery_covariate_experiences
    
    # combine wild and hatchery
    wild_temp_move_prob_quantiles %>% 
      bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) -> covariate_experiences
    
  }
  
  # convert temp to temp_actual in covariate experiences
  if (from %in% c(2:9)){
    covariate_experiences %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temperature) -> covariate_experiences
  } else {
    covariate_experiences %>% 
      mutate(temp_actual = 1) -> covariate_experiences
  }
  
  
  # Remove any instances for rug plot of fish that would have received temp0 covariate
  # run this only if the fish actually visited that state
  if (nrow(covariate_experiences) > 0){
    covariate_experiences %>% 
      dplyr::rename(date_numeric = date) %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) %>% 
      filter(!(month %in% c(1,2,3,4,5))) -> covariate_experiences
    
    # create a new DF to show bubbles of points
    covariate_experiences %>% 
      mutate(temp_round = round(temp_actual, 0)) %>% 
      group_by(temp_round, rear) %>% 
      summarise(total = n()) -> temp_rear_total_counts
    
    covariate_experiences %>% 
      mutate(temp_round = round(temp_actual, 0)) %>% 
      group_by(temp_round, rear, next_state) %>% 
      summarise(count = n()) %>% 
      ungroup() %>% 
      complete(temp_round, nesting(rear, next_state), fill = list(count = 0)) -> temp_rear_move_counts
    
    temp_rear_move_counts %>% 
      left_join(temp_rear_total_counts, by = c("temp_round", "rear")) %>% 
      mutate(prop = count/total) -> temp_move_props
    
  } else {
    temp_move_props <- data.frame(temp_round = as.numeric(NA), rear = as.character(NA), 
                                  next_state = as.character(NA), count  = as.numeric(NA), 
                                  total = as.numeric(NA), prop = as.numeric(NA))
  }
  

  
  # correct for DE - if it's a movement into a tributary
  if (is.null(wild_DE_correction) & is.null(hatchery_DE_correction)){
    # if they're both NULL, then there is no DE correction - so DE is 1
    DE_correction_df <- data.frame(rear = c("wild", "hatchery"),
                                   DE = c(1,1))
  } else {
    DE_correction_df <- data.frame(rear = c("wild", "hatchery"),
                                   DE = c(wild_DE_correction$weighted_average, hatchery_DE_correction$weighted_average))
  }
  
  
  
  rear_temp_move_prob_quantiles %>% 
    left_join(DE_correction_df, by = "rear") %>% 
    mutate(`0.5*DE` = `0.5` * DE,
           `0.025*DE` = `0.025` * DE,
           `0.975*DE` = `0.975` * DE) -> rear_temp_move_prob_quantiles
  
  # make a long form of this
  rear_temp_move_prob_quantiles %>% 
    dplyr::select(temp_actual, `0.025`, rear, DE, `0.025*DE`) %>% 
    pivot_longer(cols = c(`0.025`, `0.025*DE`), values_to = "0.025", names_to = "Model Fit") %>% 
    mutate(`Model Fit` = ifelse(`Model Fit` == "0.025", "estimated", "observed")) -> rear_temp_move_prob_quantiles_0.025
  
  rear_temp_move_prob_quantiles %>% 
    dplyr::select(temp_actual, `0.5`, rear, DE, `0.5*DE`) %>% 
    pivot_longer(cols = c(`0.5`, `0.5*DE`), values_to = "0.5", names_to = "Model Fit") %>% 
    mutate(`Model Fit` = ifelse(`Model Fit` == "0.5", "estimated", "observed")) -> rear_temp_move_prob_quantiles_0.5
  
  rear_temp_move_prob_quantiles %>% 
    dplyr::select(temp_actual, `0.975`, rear, DE, `0.975*DE`) %>% 
    pivot_longer(cols = c(`0.975`, `0.975*DE`), values_to = "0.975", names_to = "Model Fit") %>% 
    mutate(`Model Fit` = ifelse(`Model Fit` == "0.975", "estimated", "observed")) -> rear_temp_move_prob_quantiles_0.975
  
  rear_temp_move_prob_quantiles_0.025 %>% 
    left_join(., rear_temp_move_prob_quantiles_0.5, by = c("temp_actual", "rear", "Model Fit", "DE")) %>% 
    left_join(., rear_temp_move_prob_quantiles_0.975, by = c("temp_actual", "rear", "Model Fit", "DE")) -> rear_temp_move_prob_quantiles_long
  
  # Calculate the credible interval coverage
  # also calculate the average bias (from model estimated median)
  rear_temp_move_prob_quantiles_long %>% 
    subset(`Model Fit` == "observed") %>% 
    mutate(temp_round = round(temp_actual, 0)) %>% 
    # only keep the closest prediction per rounded temperature
    mutate(temp_diff = abs(temp_round-temp_actual)) %>% 
    group_by(temp_round) %>% 
    top_n(-1, temp_diff) %>% 
    left_join(subset(temp_move_props, next_state == to), by = c("temp_round", "rear")) %>% 
    # round the model estimate prop to the nearest 0.001 - otherwise we are getting points that are
    # falling outside of CI because points are zero and model estimates a probability that is basically
    # zero but just above
    mutate(CI_coverage = ifelse(prop >= round(`0.025`, 3) & prop <= round(`0.975`,3), 1, 0)) %>% 
    mutate(CI_coverage_weight = CI_coverage*total) %>% 
    # calculate bias
    mutate(bias = `0.5` - prop) %>% 
    mutate(bias_weight = bias*total) %>% 
    mutate(expected_count = `0.5` * total) -> rear_temp_move_prob_quantiles_CIcov
  
  rear_temp_move_prob_quantiles_CIcov %>% 
    group_by(rear) %>% 
    summarise(points_within_CI = sum(CI_coverage_weight, na.rm = T),
              total_bias = sum(bias_weight, na.rm = T),
              total = sum(total, na.rm = T),
              actual_movements = sum(count, na.rm = T),
              total_expected_count = sum(expected_count, na.rm = T)) %>% 
    mutate(CI_coverage = points_within_CI/total) %>% 
    mutate(overall_bias = total_bias/total) %>% 
    mutate(model_v_data = total_expected_count/actual_movements) %>% 
    # add information on origin and movement
    mutate(origin = origin_select, from = from, to = to) %>% 
    relocate(origin) %>% 
    relocate(from, .after = rear) %>% 
    relocate(to, .after = from) -> goodness_of_fit
  
  
  rear_temp_move_prob_plot <- ggplot(rear_temp_move_prob_quantiles_long, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                             color = rear, fill = rear, lty = `Model Fit`)) +
    # the "actual" value
    geom_line(alpha = 1) +
    geom_ribbon(data = subset(rear_temp_move_prob_quantiles_long, `Model Fit` == "observed"), aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                color = rear, fill = rear, lty = `Model Fit`), color = NA, inherit.aes = FALSE, alpha = 0.2) +
    geom_point(data = subset(temp_move_props, next_state == to & rear == "wild"),
               aes(x = temp_round, y = prop, color = rear, size = total), alpha = 0.5,
               inherit.aes = FALSE) +
    geom_point(data = subset(temp_move_props, next_state == to & rear == "hatchery"),
               aes(x = temp_round, y = prop, color = rear, size = total), alpha = 0.5,
               inherit.aes = FALSE) +
    scale_y_continuous(lim = c(-0.01,1.01), expand = c(0,0),
                       sec.axis = sec_axis(~ ., name = "Proportion of movements")) +
    scale_linetype_manual(values = c("observed" = 1, "estimated" = 2)) +
    scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    # scale_alpha_manual(values = c(0.15, 0.35)) +
    xlab(expression(~"Temperature" ~ ("°C"))) +
    ylab("Movement probability") +
    ggtitle(plot_title) +
    # scale_size(range=c(2,8),breaks=c(10, 50, 100, 250, 500, Inf),labels=c("1-9","10-49","50-99","100-249","250-499","500+"),
    scale_size(range=c(0.1,8),breaks=c(1, 10, 50, 100, 250, 500), guide="legend")
  
  return(list(plot = rear_temp_move_prob_plot,
              goodness_of_fit = goodness_of_fit))
}

## Loop through every movement that has a temperature effect and plot it and log goodness of fit
# Wild and hatchery have the same set up for DPS-wide vs. origin-specific effects


#### Set up tribtuary DE information ####

### DE information
# create a run_year_numeric field
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_numeric = seq(4, 22, 1)
run_year_index = 1:19
run_year_df <- data.frame(run_year,run_year_numeric, run_year_index)

# create a detection efficiency capability DF, for each TRIBUTARY - they're not natal
# origin specific (although they would be expected to primarily impact fish
# of that natal origin)
deschutes_river_trib_det_eff_capability <- data.frame(natal_origin = "Deschutes River",
                                                      run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                      DE = c(rep(0,9), rep(1,6), rep(0,3)))

john_day_river_trib_det_eff_capability <- data.frame(natal_origin = "John Day River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,8), rep(1,10)))

hood_river_trib_det_eff_capability <- data.frame(natal_origin = "Hood River",
                                                 run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                 DE = c(rep(0,9), rep(1,9)))

fifteenmile_creek_trib_det_eff_capability <- data.frame(natal_origin = "Fifteenmile Creek",
                                                        run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                        DE = c(rep(0,8), rep(1,10)))

umatilla_river_trib_det_eff_capability <- data.frame(natal_origin = "Umatilla River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,3), rep(1,15)))

yakima_river_trib_det_eff_capability <- data.frame(natal_origin = "Yakima River",     
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   # DE = c(rep(1,18)))
                                                   # For now - remove 04/05 run year
                                                   DE = c(0, rep(1,17)))

walla_walla_river_trib_det_eff_capability <- data.frame(natal_origin = "Walla Walla River",
                                                        run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                        DE = c(rep(0,1), rep(1,17)))

wenatchee_river_trib_det_eff_capability <- data.frame(natal_origin = "Wenatchee River",
                                                      run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                      DE = c(rep(0,7), rep(1,11)))

entiat_river_trib_det_eff_capability <- data.frame(natal_origin = "Entiat River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,4), rep(1,14)))

okanogan_river_trib_det_eff_capability <- data.frame(natal_origin = "Okanogan River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,9), rep(1,9)))

methow_river_trib_det_eff_capability <- data.frame(natal_origin = "Methow River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,5), rep(1,13)))

tucannon_river_trib_det_eff_capability <- data.frame(natal_origin = "Tucannon River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,7), rep(1,11)))

asotin_creek_trib_det_eff_capability <- data.frame(natal_origin = "Asotin Creek",      
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,7), rep(1,11)))

imnaha_river_trib_det_eff_capability <- data.frame(natal_origin = "Imnaha River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,7), rep(1,11)))


deschutes_river_trib_det_eff_capability %>% 
  bind_rows(., john_day_river_trib_det_eff_capability) %>% 
  bind_rows(., hood_river_trib_det_eff_capability) %>% 
  bind_rows(., fifteenmile_creek_trib_det_eff_capability) %>% 
  bind_rows(., umatilla_river_trib_det_eff_capability) %>% 
  bind_rows(., yakima_river_trib_det_eff_capability) %>% 
  bind_rows(., walla_walla_river_trib_det_eff_capability) %>% 
  bind_rows(., wenatchee_river_trib_det_eff_capability) %>% 
  bind_rows(., entiat_river_trib_det_eff_capability) %>% 
  bind_rows(., okanogan_river_trib_det_eff_capability) %>% 
  bind_rows(., methow_river_trib_det_eff_capability) %>% 
  bind_rows(., tucannon_river_trib_det_eff_capability) %>% 
  bind_rows(., asotin_creek_trib_det_eff_capability) %>% 
  bind_rows(., imnaha_river_trib_det_eff_capability) -> trib_det_eff_capability

trib_det_eff_capability %>% 
  left_join(run_year_df, by = "run_year") -> trib_det_eff_capability

# dplyr::rename(trib_det_eff_capability, tributary = natal_origin) -> trib_det_eff_capability


#### Middle Columbia ####
MC_DPS_temp_movements <- as.data.frame(which((btemp1_array_MCW[,,1] != 0), arr.ind = T))
MC_origin_temp_movements <- as.data.frame(which((btemp1xorigin1_array_MCW[,,1] != 0), arr.ind = T))
MC_DPS_temp_movements %>% 
  bind_rows(., MC_origin_temp_movements) -> MC_temp_movements

# drop any that are out of the mouth state
MC_temp_movements <- subset(MC_temp_movements, row != 1)

MC_origin_rear <- data.frame(origin = c("Deschutes River", "John Day River", "Fifteenmile Creek", 
                 "Umatilla River", "Yakima River", "Walla Walla River"),
                 hatchery = c(0, 0, 0, 1, 0, 1),
                 wild = rep(1, 6))

# Initialize the goodness of fit table
MC_temp_GOF_table <- data.frame()

i <- 4
j <- which(MC_temp_movements$col == 10)

# for each origin
for (i in 1:nrow(MC_origin_rear)){
  print(paste0("Origin: ", MC_origin_rear$origin[i]))
  # get covariate experiences
  if(MC_origin_rear[i,"hatchery"] == 1){
    hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = MC_origin_rear$origin[i])
  } else {
    hatchery_covariate_experiences <- NULL
  }
  
  if(MC_origin_rear[i,"wild"] == 1){
    wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = MC_origin_rear$origin[i])
  } else {
    wild_covariate_experiences <- NULL
  }

  # Loop through all movements
  for (j in 1:nrow(MC_temp_movements)){
  movement <- data.frame(from = MC_temp_movements$row[j],
                         to = MC_temp_movements$col[j])
  
  print(paste0("Movement: ", movement$from, " to ", movement$to))
  
  # if it's a movement into a tributary, correct for detection efficiency AND drop non-DE
  # years from the data; otherwise leave as null
  hatchery_trib_average_DE <- NULL
  wild_trib_average_DE <- NULL
  if (MC_temp_movements$col[j] %in% grep("Mouth", model_states)) {
    tributary_name <- tolower(word(model_states[MC_temp_movements$col[j]],1))
    if(MC_origin_rear[i,"hatchery"] == 1){
      # this transition also must be observed
      if(nrow(subset(hatchery_covariate_experiences, state == movement$to)) > 0){
        hatchery_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array,
                                                          covariate_experiences = hatchery_covariate_experiences)
      }
    }
    
    if(MC_origin_rear[i,"wild"] == 1){
      if (nrow(subset(wild_covariate_experiences, state == movement$to)) > 0){
        # this transition also must be observed
        wild_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                                      covariate_experiences = wild_covariate_experiences)
      }
    }
    
    # drop the non-DE years from the data for a fair comparison
    if(MC_origin_rear[i,"hatchery"] == 1){
      # this transition also must be observed
      if(nrow(subset(hatchery_covariate_experiences, state == movement$to)) > 0){
      hatchery_covariate_experiences %>% 
        mutate(natal_origin = MC_origin_rear$origin[i]) %>% 
        filter(year %in% subset(hatchery_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> hatchery_covariate_experiences_movement_specific
      }  else {
        # this needs to happen so that we don't keep overwriting
        hatchery_covariate_experiences -> hatchery_covariate_experiences_movement_specific
      }
    } 
    
    if(MC_origin_rear[i,"wild"] == 1){
      if (nrow(subset(wild_covariate_experiences, state == movement$to)) > 0){
        # this transition also must be observed
      wild_covariate_experiences %>% 
        mutate(natal_origin = MC_origin_rear$origin[i]) %>% 
        filter(year %in% subset(wild_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> wild_covariate_experiences_movement_specific
      }  else {
        # this needs to happen so that we don't keep overwriting
        wild_covariate_experiences -> wild_covariate_experiences_movement_specific
      }
    } 

    
    
  } else {
    # this needs to happen so that we don't keep overwriting
      hatchery_covariate_experiences -> hatchery_covariate_experiences_movement_specific
      wild_covariate_experiences -> wild_covariate_experiences_movement_specific
  }
  
  # get movement probabilities
    if(MC_origin_rear[i,"hatchery"] == 1){
      hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = MC_origin_rear$origin[i], movements = movement)
    } else {
      hatchery_temp_move_prob_array <- NULL
    }
    
    if(MC_origin_rear[i,"wild"] == 1){
      wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = MC_origin_rear$origin[i], movements = movement)
    } else {
      wild_temp_move_prob_array <- NULL
    }
    
  output <- plot_GOF_compare_rear_temp_effect_fit_to_data_DE(origin_select = MC_origin_rear[i,"origin"],
                                                                                     wild_move_prob_array = wild_temp_move_prob_array,
                                                                                     hatchery_move_prob_array = hatchery_temp_move_prob_array,
                                                                                     wild_covariate_experiences = wild_covariate_experiences_movement_specific,
                                                                                     hatchery_covariate_experiences = hatchery_covariate_experiences_movement_specific,
                                                                                     wild_DE_correction = wild_trib_average_DE,
                                                                                     hatchery_DE_correction = hatchery_trib_average_DE,
                                                                                     movements_evaluated = movement,
                                                                                     from = movement$from, to = movement$to, plot_title = paste0(MC_origin_rear[i,"origin"],": ",
                                                                                                                                                   model_states[movement$from], " to ",
                                                                                                                                                 model_states[movement$to]))
  # save the plot
  ggsave(here::here("stan_actual", "output", "fit_to_data", "temperature", paste0(word(MC_origin_rear[i,"origin"],1),"_",
                                                                                  movement$from, "_to_",
                                                                                  movement$to, ".png")), output$plot, height = 8, width = 8)
  
  # rbind the table
  MC_temp_GOF_table %>% 
    bind_rows(., output$goodness_of_fit) -> MC_temp_GOF_table
  }
}

# Now, run the DPS-wide comparisons
MC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = NULL)
MC_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = NULL)
# Loop through all movements
for (j in 1:nrow(MC_DPS_temp_movements)){
  movement <- data.frame(from = MC_DPS_temp_movements$row[j],
                         to = MC_DPS_temp_movements$col[j])
  
  print(paste0("Movement: ", movement$from, " to ", movement$to))
  
  # if it's a movement into a tributary, correct for detection efficiency AND drop non-DE
  # years from the data; otherwise leave as null
  hatchery_trib_average_DE <- NULL
  wild_trib_average_DE <- NULL
  if (MC_DPS_temp_movements$col[j] %in% grep("Mouth", model_states)) {
    tributary_name <- tolower(word(model_states[MC_DPS_temp_movements$col[j]],1))
      # this transition also must be observed
      if(nrow(subset(MC_hatchery_covariate_experiences, state == movement$to)) > 0){
        hatchery_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array,
                                                          covariate_experiences = MC_hatchery_covariate_experiences)
      }

      if (nrow(subset(MC_wild_covariate_experiences, state == movement$to)) > 0){
        # this transition also must be observed
        wild_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                                      covariate_experiences = MC_wild_covariate_experiences)
      }
    
    # drop the non-DE years from the data for a fair comparison
      # this transition also must be observed
      if(nrow(subset(MC_hatchery_covariate_experiences, state == movement$to)) > 0){
        MC_hatchery_covariate_experiences %>% 
          mutate(natal_origin = "Middle Columbia") %>% 
          filter(year %in% subset(hatchery_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> MC_hatchery_covariate_experiences_movement_specific
      }  else {
        # this needs to happen so that we don't keep overwriting
        MC_hatchery_covariate_experiences -> MC_hatchery_covariate_experiences_movement_specific
      }

    
      # this transition also must be observed
      if(nrow(subset(MC_wild_covariate_experiences, state == movement$to)) > 0){
        MC_wild_covariate_experiences %>% 
          mutate(natal_origin = "Middle Columbia") %>% 
          filter(year %in% subset(wild_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> MC_wild_covariate_experiences_movement_specific
      }  else {
        # this needs to happen so that we don't keep overwriting
        MC_wild_covariate_experiences -> MC_wild_covariate_experiences_movement_specific
      }
    
    
    
  } else {
    # this needs to happen so that we don't keep overwriting
    MC_hatchery_covariate_experiences -> MC_hatchery_covariate_experiences_movement_specific
    MC_wild_covariate_experiences -> MC_wild_covariate_experiences_movement_specific
  }
  
  # get movement probabilities
  # Note that because these are DPS-wide, you can pick any origin (since it doesn't matter)
    hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = MC_origin_rear$origin[1], movements = movement)
    wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = MC_origin_rear$origin[1], movements = movement)
  
  output <- plot_GOF_compare_rear_temp_effect_fit_to_data_DE(origin_select = "Middle Columbia",
                                                             wild_move_prob_array = wild_temp_move_prob_array,
                                                             hatchery_move_prob_array = hatchery_temp_move_prob_array,
                                                             wild_covariate_experiences = MC_wild_covariate_experiences_movement_specific,
                                                             hatchery_covariate_experiences = MC_hatchery_covariate_experiences_movement_specific,
                                                             wild_DE_correction = wild_trib_average_DE,
                                                             hatchery_DE_correction = hatchery_trib_average_DE,
                                                             movements_evaluated = movement,
                                                             from = movement$from, to = movement$to, plot_title = paste0("Middle Columbia fish: ",
                                                                                                                         model_states[movement$from], " to ",
                                                                                                                         model_states[movement$to]))
  # save the plot
  ggsave(here::here("stan_actual", "output", "fit_to_data", "temperature", paste0("MC_",
                                                                                  movement$from, "_to_",
                                                                                  movement$to, ".png")), output$plot, height = 8, width = 8)
  
  # rbind the table
  MC_temp_GOF_table %>% 
    bind_rows(., output$goodness_of_fit) -> MC_temp_GOF_table
}

# save this table
write.csv(MC_temp_GOF_table, file = here::here("stan_actual", "output", "fit_to_data", "temperature", "MC_temp_GOF_table.csv")) 




#### Upper Columbia ####
UC_DPS_temp_movements <- as.data.frame(which((btemp1_array_UCW[,,1] != 0), arr.ind = T))
UC_origin_temp_movements <- as.data.frame(which((btemp1xorigin1_array_UCW[,,1] != 0), arr.ind = T))
UC_DPS_temp_movements %>% 
  bind_rows(., UC_origin_temp_movements) -> UC_temp_movements

# drop any that are out of the mouth state
UC_temp_movements <- subset(UC_temp_movements, row != 1)

UC_origin_rear <- data.frame(origin = c("Wenatchee River", "Entiat River",
                                        "Okanogan River", "Methow River"),
                             hatchery = c(1, 0, 1, 1),
                             wild = c(1, 1, 0, 1))

# Initialize the goodness of fit table
UC_temp_GOF_table <- data.frame()

# for each origin
for (i in 1:nrow(UC_origin_rear)){
  print(paste0("Origin: ", UC_origin_rear$origin[i]))
  # get covariate experiences
  if(UC_origin_rear[i,"hatchery"] == 1){
    hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = UC_origin_rear$origin[i])
  } else {
    hatchery_covariate_experiences <- NULL
  }
  
  if(UC_origin_rear[i,"wild"] == 1){
    wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = UC_origin_rear$origin[i])
  } else {
    wild_covariate_experiences <- NULL
  }
  
  # Loop through all movements
  for (j in 1:nrow(UC_temp_movements)){
    movement <- data.frame(from = UC_temp_movements$row[j],
                           to = UC_temp_movements$col[j])
    
    print(paste0("Movement: ", movement$from, " to ", movement$to))
    
    # if it's a movement into a tributary, correct for detection efficiency AND drop non-DE
    # years from the data; otherwise leave as null
    hatchery_trib_average_DE <- NULL
    wild_trib_average_DE <- NULL
    if (UC_temp_movements$col[j] %in% grep("Mouth", model_states)) {
      tributary_name <- tolower(word(model_states[UC_temp_movements$col[j]],1))
      if(UC_origin_rear[i,"hatchery"] == 1){
        # this transition also must be observed
        if(nrow(subset(hatchery_covariate_experiences, state == movement$to)) > 0){
          hatchery_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array,
                                                            covariate_experiences = hatchery_covariate_experiences)
        }
      }
      
      if(UC_origin_rear[i,"wild"] == 1){
        if (nrow(subset(wild_covariate_experiences, state == movement$to)) > 0){
          # this transition also must be observed
          wild_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array,
                                                        covariate_experiences = wild_covariate_experiences)
        }
      }
      
      # drop the non-DE years from the data for a fair comparison
      if(UC_origin_rear[i,"hatchery"] == 1){
        # this transition also must be observed
        if(nrow(subset(hatchery_covariate_experiences, state == movement$to)) > 0){
        hatchery_covariate_experiences %>% 
          mutate(natal_origin = UC_origin_rear$origin[i]) %>% 
          filter(year %in% subset(hatchery_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> hatchery_covariate_experiences_movement_specific
        }  else {
          # this needs to happen so that we don't keep overwriting
          hatchery_covariate_experiences -> hatchery_covariate_experiences_movement_specific
        }
      } 
      
      if(UC_origin_rear[i,"wild"] == 1){
        # this transition also must be observed
        if(nrow(subset(wild_covariate_experiences, state == movement$to)) > 0){
        wild_covariate_experiences %>% 
          mutate(natal_origin = UC_origin_rear$origin[i]) %>% 
          filter(year %in% subset(wild_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> wild_covariate_experiences_movement_specific
        }  else {
          # this needs to happen so that we don't keep overwriting
          wild_covariate_experiences -> wild_covariate_experiences_movement_specific
        }
      } 
      
      
      
    } else {
      # this needs to happen so that we don't keep overwriting
      hatchery_covariate_experiences -> hatchery_covariate_experiences_movement_specific
      wild_covariate_experiences -> wild_covariate_experiences_movement_specific
    }
    
    # get movement probabilities
    if(UC_origin_rear[i,"hatchery"] == 1){
      hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = UC_origin_rear$origin[i], movements = movement)
    } else {
      hatchery_temp_move_prob_array <- NULL
    }
    
    if(UC_origin_rear[i,"wild"] == 1){
      wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = UC_origin_rear$origin[i], movements = movement)
    } else {
      wild_temp_move_prob_array <- NULL
    }
    
    output <- plot_GOF_compare_rear_temp_effect_fit_to_data_DE(origin_select = UC_origin_rear[i,"origin"],
                                                               wild_move_prob_array = wild_temp_move_prob_array,
                                                               hatchery_move_prob_array = hatchery_temp_move_prob_array,
                                                               wild_covariate_experiences = wild_covariate_experiences_movement_specific,
                                                               hatchery_covariate_experiences = hatchery_covariate_experiences_movement_specific,
                                                               wild_DE_correction = wild_trib_average_DE,
                                                               hatchery_DE_correction = hatchery_trib_average_DE,
                                                               movements_evaluated = movement,
                                                               from = movement$from, to = movement$to, plot_title = paste0(UC_origin_rear[i,"origin"],": ",
                                                                                                                           model_states[movement$from], " to ",
                                                                                                                           model_states[movement$to]))
    # save the plot
    ggsave(here::here("stan_actual", "output", "fit_to_data", "temperature", paste0(word(UC_origin_rear[i,"origin"],1),"_",
                                                                                    movement$from, "_to_",
                                                                                    movement$to, ".png")), output$plot, height = 8, width = 8)
    
    # rbind the table
    UC_temp_GOF_table %>% 
      bind_rows(., output$goodness_of_fit) -> UC_temp_GOF_table
  }
}

# Now, run the DPS-wide comparisons
UC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = NULL)
UC_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = NULL)
# Loop through all movements
for (j in 1:nrow(UC_DPS_temp_movements)){
  movement <- data.frame(from = UC_DPS_temp_movements$row[j],
                         to = UC_DPS_temp_movements$col[j])
  
  print(paste0("Movement: ", movement$from, " to ", movement$to))
  
  # if it's a movement into a tributary, correct for detection efficiency AND drop non-DE
  # years from the data; otherwise leave as null
  hatchery_trib_average_DE <- NULL
  wild_trib_average_DE <- NULL
  if (UC_DPS_temp_movements$col[j] %in% grep("Mouth", model_states)) {
    tributary_name <- tolower(word(model_states[UC_DPS_temp_movements$col[j]],1))
      # this transition also must be observed
      if(nrow(subset(UC_hatchery_covariate_experiences, state == movement$to)) > 0){
        hatchery_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array,
                                                          covariate_experiences = UC_hatchery_covariate_experiences)
      }
    
      if (nrow(subset(UC_wild_covariate_experiences, state == movement$to)) > 0){
        # this transition also must be observed
        wild_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array,
                                                      covariate_experiences = UC_wild_covariate_experiences)
      }

    
    # drop the non-DE years from the data for a fair comparison
      # this transition also must be observed
      if(nrow(subset(UC_hatchery_covariate_experiences, state == movement$to)) > 0){
        UC_hatchery_covariate_experiences %>% 
          mutate(natal_origin = "Upper Columbia") %>% 
          filter(year %in% subset(hatchery_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> UC_hatchery_covariate_experiences_movement_specific
      }  else {
        # this needs to happen so that we don't keep overwriting
        UC_hatchery_covariate_experiences -> UC_hatchery_covariate_experiences_movement_specific
      }
    
      # this transition also must be observed
      if(nrow(subset(UC_wild_covariate_experiences, state == movement$to)) > 0){
        UC_wild_covariate_experiences %>% 
          mutate(natal_origin = "Upper Columbia") %>% 
          filter(year %in% subset(wild_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> UC_wild_covariate_experiences_movement_specific
      }  else {
        # this needs to happen so that we don't keep overwriting
        UC_wild_covariate_experiences -> UC_wild_covariate_experiences_movement_specific
      }
    
    
    
  } else {
    # this needs to happen so that we don't keep overwriting
    UC_hatchery_covariate_experiences -> UC_hatchery_covariate_experiences_movement_specific
    UC_wild_covariate_experiences -> UC_wild_covariate_experiences_movement_specific
  }
  
  # get movement probabilities
  # Note that because these are DPS-wide, you can pick any origin (since it doesn't matter)
    hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = UC_origin_rear$origin[1], movements = movement)
    wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = UC_origin_rear$origin[1], movements = movement)
  
  output <- plot_GOF_compare_rear_temp_effect_fit_to_data_DE(origin_select = "Upper Columbia",
                                                             wild_move_prob_array = wild_temp_move_prob_array,
                                                             hatchery_move_prob_array = hatchery_temp_move_prob_array,
                                                             wild_covariate_experiences = UC_wild_covariate_experiences_movement_specific,
                                                             hatchery_covariate_experiences = UC_hatchery_covariate_experiences_movement_specific,
                                                             wild_DE_correction = wild_trib_average_DE,
                                                             hatchery_DE_correction = hatchery_trib_average_DE,
                                                             movements_evaluated = movement,
                                                             from = movement$from, to = movement$to, plot_title = paste0("Upper Columbia fish: ",
                                                                                                                         model_states[movement$from], " to ",
                                                                                                                         model_states[movement$to]))
  # save the plot
  ggsave(here::here("stan_actual", "output", "fit_to_data", "temperature", paste0("UC_",
                                                                                  movement$from, "_to_",
                                                                                  movement$to, ".png")), output$plot, height = 8, width = 8)
  
  # rbind the table
  UC_temp_GOF_table %>% 
    bind_rows(., output$goodness_of_fit) -> UC_temp_GOF_table
}

# save this table
write.csv(UC_temp_GOF_table, file = here::here("stan_actual", "output", "fit_to_data", "temperature", "UC_temp_GOF_table.csv")) 



#### Snake River ####
SR_DPS_temp_movements <- as.data.frame(which((btemp1_array_SRW[,,1] != 0), arr.ind = T))
SR_origin_temp_movements <- as.data.frame(which((btemp1xorigin1_array_SRW[,,1] != 0), arr.ind = T))
SR_DPS_temp_movements %>% 
  bind_rows(., SR_origin_temp_movements) -> SR_temp_movements

# drop any that are out of the mouth state
SR_DPS_temp_movements <- subset(SR_DPS_temp_movements, row != 1)
SR_temp_movements <- subset(SR_temp_movements, row != 1)

SR_origin_rear <- data.frame(origin = c("Tucannon River", "Asotin Creek", 
                                        "Clearwater River", "Salmon River", 
                                        "Grande Ronde River", "Imnaha River"),
                             hatchery = c(1, 0, 1, 1, 1, 1),
                             wild = rep(1, 6))

# Initialize the goodness of fit table
SR_temp_GOF_table <- data.frame()

# for each origin
for (i in 1:nrow(SR_origin_rear)){
  print(paste0("Origin: ", SR_origin_rear$origin[i]))
  # get covariate experiences
  if(SR_origin_rear[i,"hatchery"] == 1){
    hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = SR_origin_rear$origin[i])
  } else {
    hatchery_covariate_experiences <- NULL
  }
  
  if(SR_origin_rear[i,"wild"] == 1){
    wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = SR_origin_rear$origin[i])
  } else {
    wild_covariate_experiences <- NULL
  }
  
  # Loop through all movements
  for (j in 1:nrow(SR_temp_movements)){
    movement <- data.frame(from = SR_temp_movements$row[j],
                           to = SR_temp_movements$col[j])
    
    print(paste0("Movement: ", movement$from, " to ", movement$to))
    
    # if it's a movement into a tributary, correct for detection efficiency AND drop non-DE
    # years from the data; otherwise leave as null
    hatchery_trib_average_DE <- NULL
    wild_trib_average_DE <- NULL
    if (SR_temp_movements$col[j] %in% grep("Mouth", model_states)) {
      tributary_name <- tolower(word(model_states[SR_temp_movements$col[j]],1))
      if(SR_origin_rear[i,"hatchery"] == 1){
        # this transition also must be observed
        if(nrow(subset(hatchery_covariate_experiences, state == movement$to)) > 0){
          hatchery_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = SRH_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = SRH_envir$data$tributary_design_matrices_array,
                                                            covariate_experiences = hatchery_covariate_experiences)
        }
      }
      
      if(SR_origin_rear[i,"wild"] == 1){
        if (nrow(subset(wild_covariate_experiences, state == movement$to)) > 0){
          # this transition also must be observed
          wild_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array,
                                                        covariate_experiences = wild_covariate_experiences)
        }
      }
      
      # drop the non-DE years from the data for a fair comparison
      if(SR_origin_rear[i,"hatchery"] == 1){
        # this transition also must be observed
        if(nrow(subset(hatchery_covariate_experiences, state == movement$to)) > 0){
        hatchery_covariate_experiences %>% 
          mutate(natal_origin = SR_origin_rear$origin[i]) %>% 
          filter(year %in% subset(hatchery_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> hatchery_covariate_experiences_movement_specific
        }  else {
          # this needs to happen so that we don't keep overwriting
          hatchery_covariate_experiences -> hatchery_covariate_experiences_movement_specific
        }
      } 
      
      if(SR_origin_rear[i,"wild"] == 1){
        # this transition also must be observed
        if(nrow(subset(wild_covariate_experiences, state == movement$to)) > 0){
        wild_covariate_experiences %>% 
          mutate(natal_origin = SR_origin_rear$origin[i]) %>% 
          filter(year %in% subset(wild_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> wild_covariate_experiences_movement_specific
        }  else {
          # this needs to happen so that we don't keep overwriting
          wild_covariate_experiences -> wild_covariate_experiences_movement_specific
        }
      } 
      
      
      
    } else {
      # this needs to happen so that we don't keep overwriting
      hatchery_covariate_experiences -> hatchery_covariate_experiences_movement_specific
      wild_covariate_experiences -> wild_covariate_experiences_movement_specific
    }
    
    # get movement probabilities
    if(SR_origin_rear[i,"hatchery"] == 1){
      hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = SR_origin_rear$origin[i], movements = movement)
    } else {
      hatchery_temp_move_prob_array <- NULL
    }
    
    if(SR_origin_rear[i,"wild"] == 1){
      wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = SR_origin_rear$origin[i], movements = movement)
    } else {
      wild_temp_move_prob_array <- NULL
    }
    
    output <- plot_GOF_compare_rear_temp_effect_fit_to_data_DE(origin_select = SR_origin_rear[i,"origin"],
                                                               wild_move_prob_array = wild_temp_move_prob_array,
                                                               hatchery_move_prob_array = hatchery_temp_move_prob_array,
                                                               wild_covariate_experiences = wild_covariate_experiences_movement_specific,
                                                               hatchery_covariate_experiences = hatchery_covariate_experiences_movement_specific,
                                                               wild_DE_correction = wild_trib_average_DE,
                                                               hatchery_DE_correction = hatchery_trib_average_DE,
                                                               movements_evaluated = movement,
                                                               from = movement$from, to = movement$to, plot_title = paste0(SR_origin_rear[i,"origin"],": ",
                                                                                                                           model_states[movement$from], " to ",
                                                                                                                           model_states[movement$to]))
    # save the plot
    ggsave(here::here("stan_actual", "output", "fit_to_data", "temperature", paste0(word(SR_origin_rear[i,"origin"],1),"_",
                                                                                    movement$from, "_to_",
                                                                                    movement$to, ".png")), output$plot, height = 8, width = 8)
    
    # rbind the table
    SR_temp_GOF_table %>% 
      bind_rows(., output$goodness_of_fit) -> SR_temp_GOF_table
  }
}




# Now, run the DPS-wide comparisons
SR_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = NULL)
SR_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = NULL)
  # Loop through all movements
  for (j in 1:nrow(SR_DPS_temp_movements)){
    movement <- data.frame(from = SR_DPS_temp_movements$row[j],
                           to = SR_DPS_temp_movements$col[j])
    
    print(paste0("Movement: ", movement$from, " to ", movement$to))
    
    # if it's a movement into a tributary, correct for detection efficiency AND drop non-DE
    # years from the data; otherwise leave as null
    hatchery_trib_average_DE <- NULL
    wild_trib_average_DE <- NULL
    if (SR_DPS_temp_movements$col[j] %in% grep("Mouth", model_states)) {
      tributary_name <- tolower(word(model_states[SR_DPS_temp_movements$col[j]],1))
        # this transition also must be observed
        if(nrow(subset(SR_hatchery_covariate_experiences, state == movement$to)) > 0){
          hatchery_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = SRH_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = SRH_envir$data$tributary_design_matrices_array,
                                                            covariate_experiences = SR_hatchery_covariate_experiences)
        }
      
        if (nrow(subset(SR_wild_covariate_experiences, state == movement$to)) > 0){
          # this transition also must be observed
          wild_trib_average_DE <- calculate_weighted_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = tributary_name, tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array,
                                                        covariate_experiences = SR_wild_covariate_experiences)
        }
      
      # drop the non-DE years from the data for a fair comparison
        # this transition also must be observed
        if(nrow(subset(SR_hatchery_covariate_experiences, state == movement$to)) > 0){
          SR_hatchery_covariate_experiences %>% 
            mutate(natal_origin = "Snake River") %>% 
            filter(year %in% subset(hatchery_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> SR_hatchery_covariate_experiences_movement_specific
        }  else {
          # this needs to happen so that we don't keep overwriting
          SR_hatchery_covariate_experiences -> SR_hatchery_covariate_experiences_movement_specific
        }
      
        # this transition also must be observed
        if(nrow(subset(SR_wild_covariate_experiences, state == movement$to)) > 0){
          SR_wild_covariate_experiences %>% 
            mutate(natal_origin = "Snake River") %>% 
            filter(year %in% subset(wild_trib_average_DE$DE_by_year, !(is.na(`0.5`)))$year) -> SR_wild_covariate_experiences_movement_specific
        }  else {
          # this needs to happen so that we don't keep overwriting
          SR_wild_covariate_experiences -> SR_wild_covariate_experiences_movement_specific
        }
      
      
      
    } else {
      # this needs to happen so that we don't keep overwriting
      SR_hatchery_covariate_experiences -> SR_hatchery_covariate_experiences_movement_specific
      SR_wild_covariate_experiences -> SR_wild_covariate_experiences_movement_specific
    }
    
    # get movement probabilities
    # Note that because these are DPS-wide, you can pick any origin (since it doesn't matter)
      hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = SR_origin_rear$origin[1], movements = movement)
      wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = SR_origin_rear$origin[1], movements = movement)
    
    output <- plot_GOF_compare_rear_temp_effect_fit_to_data_DE(origin_select = "Snake River",
                                                               wild_move_prob_array = wild_temp_move_prob_array,
                                                               hatchery_move_prob_array = hatchery_temp_move_prob_array,
                                                               wild_covariate_experiences = SR_wild_covariate_experiences_movement_specific,
                                                               hatchery_covariate_experiences = SR_hatchery_covariate_experiences_movement_specific,
                                                               wild_DE_correction = wild_trib_average_DE,
                                                               hatchery_DE_correction = hatchery_trib_average_DE,
                                                               movements_evaluated = movement,
                                                               from = movement$from, to = movement$to, plot_title = paste0("Snake River fish: ",
                                                                                                                           model_states[movement$from], " to ",
                                                                                                                           model_states[movement$to]))
    # save the plot
    ggsave(here::here("stan_actual", "output", "fit_to_data", "temperature", paste0("SR_",
                                                                                    movement$from, "_to_",
                                                                                    movement$to, ".png")), output$plot, height = 8, width = 8)
    
    # rbind the table
    SR_temp_GOF_table %>% 
      bind_rows(., output$goodness_of_fit) -> SR_temp_GOF_table
  }

# save this table
write.csv(SR_temp_GOF_table, file = here::here("stan_actual", "output", "fit_to_data", "temperature", "SR_temp_GOF_table.csv")) 

#### Alternative plots for 2-3, Deschutes focus ####

# hm, I don't think this is worth the effort, because really you'd have to also 
# correct for DE into the Deschutes, not just select specific years

# Needs to be run for all origins (MC) and DPS-wide (UC and SR)
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_numeric = seq(4, 22, 1)
run_year_index = 1:19
run_year_df <- data.frame(run_year,run_year_numeric, run_year_index)

UMA_H_transitions_temps %>% 
  mutate(run_year_index = ceiling(entry_date/365.25)+1) %>% 
  left_join(., run_year_df, by = "run_year_index") %>% 
  filter(run_year %in% c("13/14", "14/15", "15/16", "16/17", "17/18", "18/19")) -> UMA_H_transitions_temps_Deschutes_DE_years

### Snake River

movement <- data.frame(from = 2,
                       to = 3)

print(paste0("Movement: ", movement$from, " to ", movement$to))

# Drop non-DE years for Deschutes
SR_hatchery_covariate_experiences %>% 
  mutate(run_year_index = ceiling(date/365.25)+1) %>% 
  left_join(., run_year_df, by = "run_year_index") %>% 
  filter(run_year %in% c("13/14", "14/15", "15/16", "16/17", "17/18", "18/19")) -> SR_hatchery_covariate_experiences_Deschutes_DE

SR_wild_covariate_experiences %>% 
  mutate(run_year_index = ceiling(date/365.25)+1) %>% 
  left_join(., run_year_df, by = "run_year_index") %>% 
  filter(run_year %in% c("13/14", "14/15", "15/16", "16/17", "17/18", "18/19")) -> SR_wild_covariate_experiences_Deschutes_DE


# get movement probabilities
# Note that because these are DPS-wide, you can pick any origin (since it doesn't matter)
hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = SR_origin_rear$origin[1], movements = movement)
wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = SR_origin_rear$origin[1], movements = movement)

output_SR_2_3_Deschutes_DE <- plot_GOF_compare_rear_temp_effect_fit_to_data_DE(origin_select = "Snake River",
                                                           wild_move_prob_array = wild_temp_move_prob_array,
                                                           hatchery_move_prob_array = hatchery_temp_move_prob_array,
                                                           wild_covariate_experiences = SR_wild_covariate_experiences_Deschutes_DE,
                                                           hatchery_covariate_experiences = SR_hatchery_covariate_experiences_Deschutes_DE,
                                                           wild_DE_correction = wild_trib_average_DE,
                                                           hatchery_DE_correction = hatchery_trib_average_DE,
                                                           movements_evaluated = movement,
                                                           from = movement$from, to = movement$to, plot_title = paste0("Snake River fish: ",
                                                                                                                       model_states[movement$from], " to ",
                                                                                                                       model_states[movement$to]))
# save the plot
ggsave(here::here("stan_actual", "output", "fit_to_data", "temperature", "SR_2_3_Deschutes_DE_years.png"), 
                  output_SR_2_3_Deschutes_DE$plot, height = 8, width = 8)








#### Temperature vs. year ####

### Inspect year vs. temperature effects - Umatilla Hatchery for example

# select the right origin
subset(origin_param_map, natal_origin == "Umatilla River")
# origin1

# get the temperature effect
UMA_2_10_temp1 <- MCH_fit_summary$variable[grep("btemp1xorigin1_matrix_2_10_DE", MCH_fit_summary$variable)]
subset(MCH_fit_summary, variable %in% UMA_2_10_temp1)

# get the year effects - raw and matt trick scaling
UMA_2_10_raw_year_variables <- MCH_fit_summary$variable[grep("byearxorigin1_raw_vector_2_10", MCH_fit_summary$variable)]
subset(MCH_fit_summary, variable %in% UMA_2_10_raw_year_variables) -> UMA_2_10_raw_year_variables_summary
UMA_2_10_scale_year_variables <- MCH_fit_summary$variable[grep("sigma_yearxorigin1_vector_2_10", MCH_fit_summary$variable)]
subset(MCH_fit_summary, variable %in% UMA_2_10_scale_year_variables)

# Ok, so we don't see a huge year effect here. The effects for individual years are centered around 0.

# So what else could it be?
UMA_2_10_params <- MCH_fit_summary$variable[grep("2_10", MCH_fit_summary$variable)]
UMA_2_10_params <- UMA_2_10_params[!(grepl("origin2", UMA_2_10_params))]


#### Multiple movements on same plot ####

# For some origins, there are multiple interesting and deleterious movements - 
# for example, Middle Columbia fish could go to loss, Deschutes, or overshoot

# plotting function
plot_compare_rear_temp_effect_multiple_movements <- function(origin_select,
                                          wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                          wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                          movements_evaluated,
                                          plot_title = NULL, plot_legend = FALSE){
  
  # create df to index to right dam temp joining; here we need to change names of 
  # WEL and LGR because we have slow v fast there
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL_quick", "ICH", "LGR_quick"),
                          state = seq(2,9))
  
  niter <- 4000 # this is the number of draws we have
  
  # Initialize a df to store every movement, then loop through each movement and bind
  # rear_temp_move_prob_quantiles <- data.frame(temp_actual = NA, `0.025` = NA, `0.5` = NA, `0.975` = NA,
  #                                             rear = NA, from = NA, to = NA)
  rear_temp_move_prob_quantiles <- data.frame()
  for (i in 1:nrow(movements_evaluated)){
    # First, determine if this origin has both a hatchery and a wild population
    
    # If hatchery is NA, run wild only
    if (is.null(hatchery_covariate_experiences)){
      wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,i])
      
      colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
      temp_predict <- seq(-2,2,length = 100)
      wild_temp_move_prob$temp <- temp_predict
      
      # Add a column with the actual temperatures
      if (movements_evaluated$from[i] %in% c(1:9)){
        wild_temp_move_prob %>% 
          mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
                   window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temp) -> wild_temp_move_prob
      } else {
        wild_temp_move_prob %>% 
          mutate(temp_actual = 1) -> wild_temp_move_prob
      }
      
      
      
      
      wild_temp_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(temp_actual) %>% 
        dplyr::summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "wild") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> wild_temp_move_prob_quantiles
      
      # get data organized for rug plot
      wild_covariate_experiences %>% 
        mutate(rear = "wild") -> wild_covariate_experiences
      
      wild_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> covariate_experiences
      
      rear_temp_move_prob_quantiles %>% 
        bind_rows(., wild_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
      
    } 
    # If wild is NA, run hatchery only
    else if (is.null(wild_covariate_experiences)){
      hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,i])
      
      colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
      temp_predict <- seq(-2,2,length = 100)
      hatchery_temp_move_prob$temp <- temp_predict
      
      # Add a column with the actual temperatures
      if (movements_evaluated$from[i] %in% c(1:9)){
        hatchery_temp_move_prob %>% 
          mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
                   window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temp) -> hatchery_temp_move_prob
      } else {
        hatchery_temp_move_prob %>% 
          mutate(temp_actual = 1) -> hatchery_temp_move_prob
      }
      
      
      
      
      hatchery_temp_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(temp_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "hatchery") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> hatchery_temp_move_prob_quantiles
      
      # get data organized for rug plot
      hatchery_covariate_experiences %>% 
        mutate(rear = "hatchery") -> hatchery_covariate_experiences
      
      hatchery_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> covariate_experiences
      
      # combine wild and hatchery
      rear_temp_move_prob_quantiles %>% 
        bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
      
    } 
    # else run both
    else {
      wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,i])
      
      colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
      temp_predict <- seq(-2,2,length = 100)
      wild_temp_move_prob$temp <- temp_predict
      
      # Add a column with the actual temperatures
      if (movements_evaluated$from[i] %in% c(1:9)){
        wild_temp_move_prob %>% 
          mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
                   window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temp) -> wild_temp_move_prob
      } else {
        wild_temp_move_prob %>% 
          mutate(temp_actual = 1) -> wild_temp_move_prob
      }
      
      
      wild_temp_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(temp_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "wild") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> wild_temp_move_prob_quantiles
      
      # get data organized for rug plot
      wild_covariate_experiences %>% 
        mutate(rear = "wild") -> wild_covariate_experiences
      
      wild_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> wild_covariate_experiences
      
      hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,i])
      
      colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
      hatchery_temp_move_prob$temp <- temp_predict
      
      # Add a column with the actual temperatures
      if (movements_evaluated$from[i] %in% c(1:9)){
        hatchery_temp_move_prob %>% 
          mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
                   window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temp) -> hatchery_temp_move_prob
      } else {
        hatchery_temp_move_prob %>% 
          mutate(temp_actual = 1) -> hatchery_temp_move_prob
      }
      
      hatchery_temp_move_prob %>% 
        pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
        group_by(temp_actual) %>% 
        summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
        pivot_wider(names_from = q, values_from = prob) %>% 
        mutate(rear = "hatchery") %>% 
        mutate(from = movements_evaluated$from[i],
               to = movements_evaluated$to[i]) -> hatchery_temp_move_prob_quantiles
      
      # get data organized for rug plot
      hatchery_covariate_experiences %>% 
        mutate(rear = "hatchery") -> hatchery_covariate_experiences
      
      hatchery_covariate_experiences %>% 
        # keep only the state that is our from state
        filter(state == movements_evaluated$from[i]) -> hatchery_covariate_experiences
      
      # combine wild and hatchery
      rear_temp_move_prob_quantiles %>% 
        bind_rows(., wild_temp_move_prob_quantiles) %>% 
        bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
      
      wild_covariate_experiences %>% 
        bind_rows(., hatchery_covariate_experiences) -> covariate_experiences
      
    }
    
  }
  
  

  
  # convert temp to temp_actual in covariate experiences
  if (movements_evaluated$from[i] %in% c(1:9)){
    covariate_experiences %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == movements_evaluated$from[i])$dam, "_sd")]*temperature) -> covariate_experiences
  } else {
    covariate_experiences %>% 
      mutate(temp_actual = 1) -> covariate_experiences
  }
  
  
  # Remove any instances for rug plot of fish that would have received temp0 covariate
  covariate_experiences %>% 
    dplyr::rename(date_numeric = date) %>%
    # keep only jan/feb/mar 
    mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
    mutate(month = month(date)) %>% 
    filter(!(month %in% c(1,2,3,4,5))) -> covariate_experiences
  
  # rear_temp_move_prob_quantiles$to <- as.character(rear_temp_move_prob_quantiles$to)
  
  
  movement_colors <- c("Overshoot" = "#ff7f00", "Overshoot - ICH" = "#ff7f00",
                       "Overshoot - PRA" = "#ff7f00",
                       "Deschutes River" = "#a6cee3",
                       "Home" = "#1f78b4", "Loss" = "#e31a1c")
  
  rear_temp_move_prob_quantiles %>% 
    left_join(., movements_evaluated, by = "to") -> rear_temp_move_prob_quantiles
  
  if (plot_legend == TRUE){
    
    combined_plot <- ggplot(rear_temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = Movement, fill = Movement)) +
      geom_line(linewidth = 2.5) +
      geom_ribbon(alpha = 0.2, color = NA) +
      # geom_rug(data = covariate_experiences, aes(x = temp_actual), inherit.aes = FALSE,
      #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
      scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
      # scale_x_continuous(lim = c(0,ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
      # common x-axis scale across all populations
      scale_x_continuous(lim = c(0, 22.5), expand = c(0,0)) + 
      scale_color_manual(values = movement_colors) +
      scale_fill_manual(values =  movement_colors) +
      xlab(expression(~"Temperature" ~ ("°C"))) +
      ylab("Movement probability") +
      coord_cartesian(clip = "off") +
      theme(panel.grid.major = element_line(color = "gray90"),
            panel.background = element_rect(fill = "white", color = NA),
            panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
            legend.key.height = unit(1.25, "cm"),
            legend.key.width = unit(1.25, "cm"),
            legend.title = element_text(size = 25),
            legend.text = element_text(size = 15),
            # these plot margins are to leave space for the population name on the big figure
            plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
    # for testing
    temp_legend <- ggpubr::get_legend(combined_plot)
    # temp_legend$widths[[1]] <- unit(0, "cm")
    # temp_legend$widths[[5]] <- unit(0, "cm")
    # temp_legend$heights[[1]] <- unit(0, "cm")
    # temp_legend$heights[[5]] <- unit(0, "cm")
    temp_plot_legend_gg <- as_ggplot(temp_legend)
    
    # for testing
    # ggsave(here::here("stan_actual", "output", "paper_figures", "01_legend_test.png"), temp_plot_legend_gg, height = 4, width = 4)
    
  } else {
    # suppress common legend - for combined plot
    rear_temp_move_prob_plot <- ggplot(rear_temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = Movement, fill = Movement)) +
      geom_line() +
      geom_ribbon(alpha = 0.2, color = NA) +
      # geom_rug(data = covariate_experiences, aes(x = temp_actual), inherit.aes = FALSE,
      #          sides = "t", length = unit(0.3, "cm"), outside = TRUE) +
      scale_y_continuous(lim = c(0,1)) +
      # common x-axis scale across all populations
      coord_cartesian(xlim = c(4, 22.5), expand = FALSE) +
      # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
      scale_color_manual(values = movement_colors) +
      scale_fill_manual(values =  movement_colors) +
      xlab(expression(~"Temperature" ~ ("°C"))) +
      ylab("Movement probability") +
      theme(legend.position = "none",
            panel.grid.major = element_line(color = "gray90"),
            panel.background = element_rect(fill = "white", color = "black"),
            panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
            # turn off the axis titles on each individual plot and just show one for whole plot
            axis.title = element_blank(),
            # these plot margins are to leave space for the population name on the big figure
            plot.margin = unit(c(0, 0.2, 0.2, 0.2),"cm"))
    
    density_plot <- ggplot(data = covariate_experiences, aes(temp_actual))+
      # geom_density(alpha = 0.1) +
      # geom_rug(sides = "t", length = unit(0.2, "cm"), outside = FALSE) +
      # geom_histogram(alpha = 0.5, bins = 60) +
      geom_histogram(aes(y=..count../sum(..count..)), alpha = 0.5, bins = 60) +
      ylab("Density") +
      # scale_x_continuous(lim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual))), expand = c(0,0)) +
      # scale_x_continuous(lim = c(0, 23), expand = c(0,0)) +
      # scale_y_continuous(lim = c(0,0.75), expand = c(0,0),
      #                    breaks = c(0, 0.25, 0.50)) +
      # these labels aren't real, they just need to be the same size as the labels from the temp effect plot
      scale_y_continuous(n.breaks = 2, labels = c("0.00", "1.00")) +
      # coord_cartesian(xlim = c(floor(min(covariate_experiences$temp_actual)),ceiling(max(covariate_experiences$temp_actual)))) +
      coord_cartesian(xlim = c(4, 22.5), expand = FALSE) +
      theme(axis.text.x = element_blank(),
            axis.text.y = element_text(color = "white"),
            axis.ticks.x = element_blank(),
            # axis.ticks.length.x=unit(.1, "cm"),
            axis.ticks.length.x=unit(0, "cm"),
            axis.ticks.y = element_line(color = "white"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(color = "white"),
            panel.background = element_blank(),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            # these plot margins are to leave space for the population name on the big figure
            # we actually don't need this anymore I think, because of the histogram
            plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))
    # testing theme
    # theme(
    #       # these plot margins are to leave space for the population name on the big figure
    #       plot.margin = unit(c(0.2, 0.2, 0, 0.2),"cm"))
    
    combined_plot <- ggarrange(density_plot, rear_temp_move_prob_plot, nrow = 2, ncol = 1,
              heights = c(2,6))
    
    # for testing
    # ggsave(here::here("stan_actual", "output", "paper_figures", "01_test.png"), combined_plot, height = 6, width = 6)
    
    
  }

  
  
  return(combined_plot)
}


#### Prepare data to create figures for individual origins: ####
# John Day River (W), Umatilla River (H, W), Yakima River (W), Walla Walla River (H, W),
# Wenatchee River (H, W), Entiat River (W), Tucannon River (H, W), Imnaha River (H, W)

### John Day River ###
JDR_movements <- data.frame(from = c(2, 2, 2, 2), to = c(3, 10, 12, 43),
                            Movement = c("Overshoot", "Deschutes River",
                                         "Home", "Loss"))
JDR_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "John Day River", movements = JDR_movements)
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")

### Umatilla River ###
UMA_movements <- data.frame(from = c(2,2,2,2), to = c(3, 10, 18, 43),
                            Movement = c("Overshoot", "Deschutes River",
                                         "Home", "Loss"))
UMA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = UMA_movements)
UMA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = UMA_movements)
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")

### Yakima River ###
YAK_movements <- data.frame(from = c(3,3,3,3), to = c(8, 4, 20, 43),
                            Movement = c("Overshoot - ICH", "Overshoot - PRA",
                                         "Home", "Loss"))

YAK_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = YAK_movements)
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")

### Walla Walla River ###
WAWA_movements <- data.frame(from = c(3,3,3,3), to = c(8, 4, 22, 43),
                             Movement = c("Overshoot - ICH", "Overshoot - PRA",
                                          "Home", "Loss"))
WAWA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")

### Wenatchee River ###
WEN_movements <- data.frame(from = c(5,5,5), to = c(6, 24, 43),
                            Movement = c("Overshoot",
                                         "Home", "Loss"))
WEN_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Wenatchee River", movements = WEN_movements)
WEN_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Wenatchee River", movements = WEN_movements)
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")

### Entiat River ###
ENT_movements <- data.frame(from = c(6,6,6), to = c(7, 26, 43),
                            Movement = c("Overshoot",
                                         "Home", "Loss"))
ENT_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Entiat River", movements = ENT_movements)
ENT_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Entiat River", movements = ENT_movements)
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
ENT_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Entiat River")

### Tucannon River ###
TUC_movements <- data.frame(from = c(8,8,8), to = c(9, 32, 43),
                            Movement = c("Overshoot", 
                                         "Home", "Loss"))
TUC_wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = "Tucannon River", movements = TUC_movements)
TUC_hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = "Tucannon River", movements = TUC_movements)
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")

### Imnaha River ###
IMN_movements <- data.frame(from = c(9,9), to = c(39, 43),
                            Movement = c("Home", "Loss"))
IMN_wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = "Imnaha River", movements = IMN_movements)
IMN_hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = "Imnaha River", movements = IMN_movements)
IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Imnaha River")
IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Imnaha River")




#### Create figures for individual origins: ####
# John Day River (W), Umatilla River (H, W), Yakima River (W), Walla Walla River (H, W),
# Wenatchee River (H, W), Entiat River (W), Tucannon River (H, W), Imnaha River (H, W)
### John Day River ###
JDR_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "John Day River",
                                                            wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                            
                                                            wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                            
                                                            movements_evaluated = JDR_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "JDR_wild_compare_movement_temp.png"), JDR_wild_compare_movement_temp, height = 8, width = 8)

### Umatilla River ###

# UMA wild plot
UMA_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Umatilla River",
                                                                                   wild_move_prob_array = UMA_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                                                   movements_evaluated = UMA_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "UMA_wild_compare_movement_temp.png"), UMA_wild_compare_movement_temp, height = 8, width = 8)

# UMA hatchery plot
UMA_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Umatilla River",
                                                                                       hatchery_move_prob_array = UMA_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                                       movements_evaluated = UMA_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "UMA_hatchery_compare_movement_temp.png"), UMA_hatchery_compare_movement_temp, height = 8, width = 8)

### Yakima River ###

YAK_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Yakima River",
                                                                                   wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                                                   
                                                                                   wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                                                   
                                                                                   movements_evaluated = YAK_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "YAK_wild_compare_movement_temp.png"), YAK_wild_compare_movement_temp, height = 8, width = 8)



### Walla Walla River ###

# WAWA wild plot
WAWA_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Walla Walla River",
                                                                               wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                                               wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                                               movements_evaluated = WAWA_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WAWA_wild_compare_movement_temp.png"), WAWA_wild_compare_movement_temp, height = 8, width = 8)

# WAWA hatchery plot
WAWA_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Walla Walla River",
                                                                               hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                                               hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                               movements_evaluated = WAWA_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WAWA_hatchery_compare_movement_temp.png"), WAWA_hatchery_compare_movement_temp, height = 8, width = 8)

### Wenatchee River ###

# WEN wild plot
WEN_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Wenatchee River",
                                                                                   wild_move_prob_array = WEN_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = WEN_wild_covariate_experiences,
                                                                                   movements_evaluated = WEN_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WEN_wild_compare_movement_temp.png"), WEN_wild_compare_movement_temp, height = 8, width = 8)

# WEN hatchery plot
WEN_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Wenatchee River",
                                                                                       hatchery_move_prob_array = WEN_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
                                                                                       movements_evaluated = WEN_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WEN_hatchery_compare_movement_temp.png"), WEN_hatchery_compare_movement_temp, height = 8, width = 8)

### Entiat River ###

# ENT wild plot
ENT_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Entiat River",
                                                                                   wild_move_prob_array = ENT_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = ENT_wild_covariate_experiences,
                                                                                   movements_evaluated = ENT_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "ENT_wild_compare_movement_temp.png"), ENT_wild_compare_movement_temp, height = 8, width = 8)


### Tucannon River ###

# TUC wild plot
TUC_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Tucannon River",
                                                                                   wild_move_prob_array = TUC_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = TUC_wild_covariate_experiences,
                                                                                   movements_evaluated = TUC_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "TUC_wild_compare_movement_temp.png"), TUC_wild_compare_movement_temp, height = 8, width = 8)

# TUC hatchery plot
TUC_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Tucannon River",
                                                                                       hatchery_move_prob_array = TUC_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
                                                                                       movements_evaluated = TUC_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "TUC_hatchery_compare_movement_temp.png"), TUC_hatchery_compare_movement_temp, height = 8, width = 8)

### Imnaha River ###

# IMN wild plot
IMN_wild_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Imnaha River",
                                                                                   wild_move_prob_array = IMN_wild_temp_move_prob_array,
                                                                                   wild_covariate_experiences = IMN_wild_covariate_experiences,
                                                                                   movements_evaluated = IMN_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "IMN_wild_compare_movement_temp.png"), IMN_wild_compare_movement_temp, height = 8, width = 8)

# IMN hatchery plot
IMN_hatchery_compare_movement_temp <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "Imnaha River",
                                                                                       hatchery_move_prob_array = IMN_hatchery_temp_move_prob_array,
                                                                                       hatchery_covariate_experiences = IMN_hatchery_covariate_experiences,
                                                                                       movements_evaluated = IMN_movements)

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "IMN_hatchery_compare_movement_temp.png"), IMN_hatchery_compare_movement_temp, height = 8, width = 8)


# Create the legend figure by itself
temp_plot_for_legend <- plot_compare_rear_temp_effect_multiple_movements(origin_select = "John Day River",
                                                                   wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                                   
                                                                   wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                                   
                                                                   movements_evaluated = JDR_movements,
                                                                   plot_legend= TRUE)
temp_legend <- ggpubr::get_legend(temp_plot_for_legend)
temp_plot_legend_gg <- as_ggplot(temp_legend) + theme(panel.background = element_rect(fill = "white", color = "white"))


#### Generate the figure for the paper, using the figures above ####

# combined_movement_temp_plot <- ggarrange(JDR_wild_compare_movement_temp, UMA_wild_compare_movement_temp, UMA_hatchery_compare_movement_temp,
#           YAK_wild_compare_movement_temp, WAWA_wild_compare_movement_temp, WAWA_hatchery_compare_movement_temp,
#           WEN_wild_compare_movement_temp, WEN_hatchery_compare_movement_temp, ENT_wild_compare_movement_temp,
#           TUC_wild_compare_movement_temp, TUC_hatchery_compare_movement_temp, IMN_wild_compare_movement_temp, 
#           IMN_hatchery_compare_movement_temp, plot_legend_gg, nrow = 4, ncol = 4,
#           labels = c("(A)", "(B)", "(C)", "(D)", "(E)", "(F)", "(G)", "(H)", "(I)",
#                      "(J)", "(K)", "(L)", "(M)"),
#           label.x = 0.00, label.y = 0.95, font.label = list(size = 10, face = "plain"),
#           hjust = 0, vjust = 0)
# 
# ggsave(here::here("stan_actual", "output", "paper_figures", "combined_movement_temp_plot_v1.png"), combined_movement_temp_plot, height = 12, width = 16)

# without imnaha
combined_movement_temp_plot <- ggarrange(JDR_wild_compare_movement_temp, UMA_wild_compare_movement_temp, UMA_hatchery_compare_movement_temp,
                                         YAK_wild_compare_movement_temp, WAWA_wild_compare_movement_temp, WAWA_hatchery_compare_movement_temp,
                                         WEN_wild_compare_movement_temp, WEN_hatchery_compare_movement_temp, ENT_wild_compare_movement_temp,
                                         TUC_wild_compare_movement_temp, TUC_hatchery_compare_movement_temp,
                                         temp_plot_legend_gg, nrow = 3, ncol = 4,
                                         labels = c("(A) JDR, Natural", "(B) UMA, Natural", 
                                                    "(C) UMA, Hatchery", "(D) YAK, Natural", 
                                                    "(E) WAWA, Natural", "(F) WAWA, Hatchery", 
                                                    "(G) WEN, Natural", "(H) WEN, Hatchery", 
                                                    "(I) ENT, Natural",
                                                    "(J) TUC, Natural", "(K) TUC, Hatchery"),
                                         label.x = 0.05, label.y = 0.925, font.label = list(size = 14, face = "plain"),
                                         hjust = 0, vjust = 0)

# combined_movement_temp_plot <- annotate_figure(combined_movement_temp_plot,
#                                                bottom = textGrob(expression(~"Temperature" ~ ("°C")), gp = gpar(cex = 1.3)),
#                                                  left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 1.3))) + bgcolor("white")

# Let's try this again, this time using a cowplot solution since ggpubr is
# struggling with the background color
combined_movement_temp_plot <- cowplot::ggdraw(annotate_figure(combined_movement_temp_plot,
                                               bottom = textGrob(expression(~"Temperature" ~ ("°C")), gp = gpar(cex = 1.3)),
                                               left = textGrob("Movement probability", rot = 90, gp = gpar(cex = 1.3)))) +
  theme(plot.background = element_rect(fill="white", color = NA))


ggsave(here::here("stan_actual", "output", "paper_figures", "combined_movement_temp_plot_v2.png"), combined_movement_temp_plot, height = 12, width = 16)

