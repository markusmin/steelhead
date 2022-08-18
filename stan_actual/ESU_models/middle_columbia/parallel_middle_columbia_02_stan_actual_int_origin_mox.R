# middle columbia, parallel, origin + intercept

# this model will fit to the following populations/origins:
# "Deschutes_River", "Fifteenmile_Creek", "John_Day_River", "Umatilla_River", "Yakima_River", "Walla_Walla_River"


# This script will fit an intercept + origin model using stan to the actual dataset

# FOR TESTING: setwd
# setwd("/Users/markusmin/Documents/CBR/steelhead/stan_actual/ESU_models/middle_columbia/")

# library("rstan")
library(cmdstanr)
library(posterior)
library(tidyverse)
library(lubridate)
# library(bayesplot)

##### Get all of the data in order #####

# Create a transition matrix of 1s and 0s for movements from (rows) to (columns)

# Create a starting matrix of sites, with the first occasion populated (BON to MCN)
model_states = c(
  # Mainstem states (9)
  "mainstem, mouth to BON",
  "mainstem, BON to MCN",
  "mainstem, MCN to ICH or PRA",
  "mainstem, PRA to RIS",
  "mainstem, RIS to RRE",
  "mainstem, RRE to WEL",
  "mainstem, upstream of WEL",
  "mainstem, ICH to LGR",
  "mainstem, upstream of LGR",
  
  # Tributary sites (19)
  "Deschutes River",
  "John Day River",
  "Hood River",
  "Fifteenmile Creek",
  "Umatilla River",
  "Yakima River",
  "Walla Walla River",
  "Wenatchee River",
  "Entiat River",
  "Okanogan River",
  "Methow River",
  "Tucannon River",
  "Asotin Creek",
  "Clearwater River",
  "Salmon River",
  "Grande Ronde River",
  "Imnaha River",
  "BON to MCN other tributaries",
  "Upstream WEL other tributaries",
  
  # Loss
  "loss"
)

nstates <- length(model_states)

# Create a rows = from, columns = to matrix for movement probabilities


transition_matrix <- matrix(0, nrow = nstates, ncol = nstates)
rownames(transition_matrix) <- model_states
colnames(transition_matrix) <- model_states

# Populate every possible option with a 1
# 1: mainstem, mouth to BON
transition_matrix["mainstem, mouth to BON", "loss"] <- 1
transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- 1


# 2: mainstem, BON to MCN
transition_matrix["mainstem, BON to MCN", "loss"] <- 1
transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- 1
transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- 1
transition_matrix["mainstem, BON to MCN", "John Day River"] <- 1
transition_matrix["mainstem, BON to MCN", "Hood River"] <- 1
transition_matrix["mainstem, BON to MCN", "Fifteenmile Creek"] <- 1
transition_matrix["mainstem, BON to MCN", "Umatilla River"] <- 1
transition_matrix["mainstem, BON to MCN", "BON to MCN other tributaries"] <- 1


# 3: mainstem, MCN to ICH or PRA
transition_matrix["mainstem, MCN to ICH or PRA", "loss"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, BON to MCN"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Walla Walla River"] <- 1


# 4: mainstem, PRA to RIS
transition_matrix["mainstem, PRA to RIS", "loss"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, RIS to RRE"] <- 1


# 5: mainstem, RIS to RRE
transition_matrix["mainstem, RIS to RRE", "loss"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, RRE to WEL"] <- 1
transition_matrix["mainstem, RIS to RRE", "Wenatchee River"] <- 1


# 6: mainstem, RRE to WEL
transition_matrix["mainstem, RRE to WEL", "loss"] <- 1
transition_matrix["mainstem, RRE to WEL", "mainstem, RIS to RRE"] <- 1
transition_matrix["mainstem, RRE to WEL", "mainstem, upstream of WEL"] <- 1
transition_matrix["mainstem, RRE to WEL", "Entiat River"] <- 1


# 7: mainstem, upstream of WEL
transition_matrix["mainstem, upstream of WEL", "loss"] <- 1
transition_matrix["mainstem, upstream of WEL", "mainstem, RRE to WEL"] <- 1
transition_matrix["mainstem, upstream of WEL", "Okanogan River"] <- 1
transition_matrix["mainstem, upstream of WEL", "Methow River"] <- 1
transition_matrix["mainstem, upstream of WEL", "Upstream WEL other tributaries"] <- 1


# 8: mainstem, ICH to LGR
transition_matrix["mainstem, ICH to LGR", "loss"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, upstream of LGR"] <- 1
transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- 1


# 9: mainstem, upstream of LGR
transition_matrix["mainstem, upstream of LGR", "loss"] <- 1
transition_matrix["mainstem, upstream of LGR", "mainstem, ICH to LGR"] <- 1
transition_matrix["mainstem, upstream of LGR", "Asotin Creek"] <- 1
transition_matrix["mainstem, upstream of LGR", "Clearwater River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Salmon River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Grande Ronde River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Imnaha River"] <- 1


# 10: Deschutes River
transition_matrix["Deschutes River", "loss"] <- 1
transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- 1


# 11: John Day River
transition_matrix["John Day River", "loss"] <- 1
transition_matrix["John Day River", "mainstem, BON to MCN"] <- 1


# 12: Hood River
transition_matrix["Hood River", "loss"] <- 1
transition_matrix["Hood River", "mainstem, BON to MCN"] <- 1


# 13: Fifteenmile Creek
transition_matrix["Fifteenmile Creek", "loss"] <- 1
transition_matrix["Fifteenmile Creek", "mainstem, BON to MCN"] <- 1


# 14: Umatilla River
transition_matrix["Umatilla River", "loss"] <- 1
transition_matrix["Umatilla River", "mainstem, BON to MCN"] <- 1


# 15: Yakima River
transition_matrix["Yakima River", "loss"] <- 1
transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- 1


# 16: Walla Walla River
transition_matrix["Walla Walla River", "loss"] <- 1
transition_matrix["Walla Walla River", "mainstem, MCN to ICH or PRA"] <- 1


# 17: Wenatchee River
transition_matrix["Wenatchee River", "loss"] <- 1
transition_matrix["Wenatchee River", "mainstem, RIS to RRE"] <- 1


# 18: Entiat River
transition_matrix["Entiat River", "loss"] <- 1
transition_matrix["Entiat River", "mainstem, RRE to WEL"] <- 1


# 19: Okanogan River
transition_matrix["Okanogan River", "loss"] <- 1
transition_matrix["Okanogan River", "mainstem, upstream of WEL"] <- 1


# 20: Methow River
transition_matrix["Methow River", "loss"] <- 1
transition_matrix["Methow River", "mainstem, upstream of WEL"] <- 1


# 21: Tucannon River
transition_matrix["Tucannon River", "loss"] <- 1
transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- 1


# 22: Asotin Creek
transition_matrix["Asotin Creek", "loss"] <- 1
transition_matrix["Asotin Creek", "mainstem, upstream of LGR"] <- 1


# 23: Clearwater River
transition_matrix["Clearwater River", "loss"] <- 1
transition_matrix["Clearwater River", "mainstem, upstream of LGR"] <- 1


# 24: Salmon River
transition_matrix["Salmon River", "loss"] <- 1
transition_matrix["Salmon River", "mainstem, upstream of LGR"] <- 1


# 25: Grande Ronde River
transition_matrix["Grande Ronde River", "loss"] <- 1
transition_matrix["Grande Ronde River", "mainstem, upstream of LGR"] <- 1


# 26: Imnaha River
transition_matrix["Imnaha River", "loss"] <- 1
transition_matrix["Imnaha River", "mainstem, upstream of LGR"] <- 1

# 27: BON to MCN other tributaries
transition_matrix["BON to MCN other tributaries", "loss"] <- 1
transition_matrix["BON to MCN other tributaries", "mainstem, BON to MCN"] <- 1

# 28: Upstream WEL other tributaries
transition_matrix["Upstream WEL other tributaries", "loss"] <- 1
transition_matrix["Upstream WEL other tributaries", "mainstem, upstream of WEL"] <- 1


##### for this ESU - cut out state transitions for which we have no data #####
# Remove the states that no fish ever enter - see last lines of this script

# 28: Upstream WEL other tributaries
transition_matrix["Upstream WEL other tributaries", "loss"] <- 0
transition_matrix["Upstream WEL other tributaries", "mainstem, upstream of WEL"] <- 0



##### Data continued #####
# Use the transition matrix to calculate the possible movements
possible_movements <- rowSums(transition_matrix)

# Get the indices that are 1s (except loss, since that's 1 minus the others)
movements <- which(transition_matrix[,1:(nstates-1)] == 1, arr.ind = TRUE)
nmovements <- dim(movements)[1]
# Now get all of the movements which are fixed to zero
not_movements <- which(transition_matrix[,1:(nstates-1)] == 0, arr.ind = TRUE)
n_notmovements <- dim(not_movements)[1]

# Paste movements to name parameters
movements %>% 
  as.data.frame() %>% 
  mutate(b0_matrix_name = paste0("b0_matrix_", row, "_", col)) -> b0_matrix_names

movements %>% 
  as.data.frame() -> movements_df

middle_columbia_movements <- subset(movements_df, row %in% c(2,3,10,11,12,13,14,15,16,27) |
                                     col %in% c(2,3,10,11,12,13,14,15,16,27))

non_middle_columbia_movements  <- subset(movements_df, !(row %in% c(2,3,10,11,12,13,14,15,16,27)) &
                                          !(col %in% c(2,3,10,11,12,13,14,15,16,27)))

for (i in 1:(nrow(b0_matrix_names))){
  # Paste the numerator
  cat("real ", b0_matrix_names$b0_matrix_name[i],";", "\n", sep = "")
}

# There are six origins in this ESU, so we will have five origin parameters.
# We will only allow an origin effect within the ESU (so after PRA). Before and after they will only have an intercept term
for (i in 1:5){
  for (j in 1:nrow(middle_columbia_movements)){
    cat("real ", "borigin", i, "_matrix_", middle_columbia_movements$row[j], "_", middle_columbia_movements$col[j], ";", "\n", sep = "")
  }
  cat ("\n")
}


# write it out for the transformed parameters section
for (i in 1:(nrow(b0_matrix_names))){
  # Paste the numerator
  cat("b0_matrix[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i],";", "\n", sep = "")
}

# There are four origins in this model, so we will have three origin parameters.
for (i in 1:5){
  for (j in 1:nrow(middle_columbia_movements)){
    cat("borigin", i, "_matrix[", middle_columbia_movements$row[j], ",", middle_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", middle_columbia_movements$row[j], "_", middle_columbia_movements$col[j], ";", "\n", sep = "")
  }
  cat ("\n")
}

# Note - we need to change the other elements of the origin matrix to 0, not -100000, otherwise they'll totally change the probability
for (i in 1:5){
  for (j in 1:nrow(non_middle_columbia_movements)){
    cat("borigin", i, "_matrix[", non_middle_columbia_movements$row[j], ",", non_middle_columbia_movements$col[j], "]", " = ", 0, ";", "\n", sep = "")
  }
  cat ("\n")
}


# write out the priors
for (i in 1:(nrow(b0_matrix_names))){
  # Paste the numerator
  cat(b0_matrix_names$b0_matrix_name[i], " ~ normal(0,10);", "\n", sep = "")
}

for (i in 1:5){
  for (j in 1:nrow(middle_columbia_movements)){
    cat("borigin", i, "_matrix_", middle_columbia_movements$row[j], "_", middle_columbia_movements$col[j], " ~ normal(0,10);", "\n", sep = "")
  }
  cat ("\n")
}

##### LOAD and REFORMAT DATA #####

# Load states complete
# states_complete <- read.csv(here::here("from_hyak_transfer", "2022-07-21-complete_det_hist", "states_complete.csv"), row.names = 1)
states_complete <- read.csv("middle_columbia_adults_states_complete.csv", row.names = 1)

unique_tag_code_2s <- unique(states_complete$tag_code_2)

# first, get the maximum number of visits by any fish
states_complete %>% 
  # group_by(tag_code) %>% 
  # We need to use tag_code_2 because this splits apart repeat spawners
  group_by(tag_code_2) %>% 
  count() %>% 
  as.data.frame() -> site_visits_by_tag_code

# + 1 here because we don't yet have a loss state
max_visits <- max(site_visits_by_tag_code$n) + 1

unique_tag_codes <- unique(states_complete$tag_code_2)
nfish <- length(unique_tag_codes)

# second, create an array with the dimensions number of fish, maximum number of site visits, number of states
state_data <- array(data = 0, dim = c(nstates, max_visits, nfish))

# third, populate the state data

# Add a column for the observation
states_complete %>% 
  # group_by(tag_code) %>% 
  group_by(tag_code_2) %>% 
  mutate(order = row_number()) -> states_complete

# add a column to index by fish
states_complete %>% 
  # group_by(tag_code) %>% 
  group_by(tag_code_2) %>% 
  mutate(tag_code_number = cur_group_id()) -> states_complete



# Loop through to convert df to array
print(Sys.time())
for (i in 1:nrow(states_complete)){
# for (i in 1:100){
  # index by 1) which state, 2) which visit, 3) which fish
  state_data[which(model_states == states_complete[i, "state", drop = TRUE]),states_complete[i, "order", drop = TRUE], states_complete[i, "tag_code_number", drop = TRUE]] <- 1 
  
}
print(Sys.time())

# Now, add the loss state
# First, get the number of non-loss states
n_not_loss_states <- vector(length = nfish)
# Now, get the indices of the fish that already have a loss state
states_complete %>% 
  subset(state == "loss") -> trapped_fish

trapped_fish_tags <- trapped_fish$tag_code_2
# Get the indices of these
which(unique_tag_codes %in% trapped_fish_tags) -> trapped_fish_indices

for (i in 1:nfish){
  # 29 is the loss state
  n_not_loss_states[i] <- sum(state_data[,,i])
  state_data[29, n_not_loss_states[i]+1,i] <- 1
}

# remove the extra loss state for the trapped fish
for (i in 1:length(trapped_fish_indices)){
  j <- trapped_fish_indices[i]
  state_data[29, n_not_loss_states[j]+1,j] <-0
}

saveRDS(state_data, "middle_columbia_state_data.csv")

##### Convert the dates into a numeric index #####
# for now, comment this whole block out

# 
# # First, populate the dates for implicit site visits
# # Create a column to indicate if the date_time was interpolated
# # Note the arrival date
# states_complete %>% 
#   mutate(date_source = ifelse(!is.na(date_time), "known_arrival", "interpolated")) %>% 
#   mutate(date = NA) -> states_complete
# 
# # I think we have to loop this
# 
# # First, figure out indices of missing date_time
# missing_date_time <- is.na(states_complete$date_time)
# date_times <- states_complete$date_time
# 
# # For efficiency, only loop through the missing dates
# missing_date_time_indices <- which(missing_date_time == TRUE)
# 
# # Extract date for all known date times
# states_complete %>% 
#   mutate(date =  format(as.Date(date(states_complete[i,"date_time", drop = TRUE]), origin = "1970-01-01"))) -> states_complete
# 
# # Loop to get arrival dates in state
# # We currently don't need this, because we interpolated the time earlier
# # takes about 90 minutes to run
# # for (j in 1:length(missing_date_time_indices)){
# # # for (i in 1:100){
# #   i <- missing_date_time_indices[j]
# # 
# #     # Figure out how far back the last known time was
# #     # Truncate the missing date times to only ones prior to current state
# #     missing_date_time_subset <- missing_date_time[1:(i-1)]
# #     prev_time_index <- max(which(missing_date_time_subset %in% FALSE))
# #     
# #     # Truncate the missing date times to only ones after current state
# #     missing_date_time_subset <- missing_date_time[(i+1):nrow(states_complete)]
# #     next_time_index <- min(which(missing_date_time_subset %in% FALSE)) + i
# #     
# #     # Now, interpolate the missing time
# #     # First, figure out how long between the two known times
# #     prev_time <- ymd_hms(date_times[prev_time_index])
# #     next_time <- ymd_hms(date_times[next_time_index])
# #     time_diff <- next_time - prev_time
# #     
# #     # Get the missing time - add the time difference divided by the number 
# #     # of missing steps plus 1, multiply by which number step it is
# #     missing_time <- prev_time + (time_diff/(next_time_index - prev_time_index) * (i - prev_time_index))
# #     
# #     # Extract just the date
# #     missing_date <- date(missing_time)
# #     
# #     # populate the missing date_time
# #     states_complete[i, "date"] <- format(as.Date(missing_date, origin = "1970-01-01"))
# # }
#   
# 
# # create an empty matrix
# transition_date_matrix <- matrix(data = NA, nrow = nfish, ncol = max_visits)
# # populate it with the dates
# for (i in 1:nrow(states_complete)){
#   transition_date_matrix[states_complete[i,"tag_code_number", drop = TRUE], states_complete[i, "order", drop = TRUE]] <- states_complete[i, "date_time", drop = TRUE]
# }
# 
# 
# # We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
# date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(date(x)) - ymd("2005-05-31"))
# transition_date_matrix %>% 
#   as_tibble() %>% 
#   mutate_all(date_numeric) -> transition_date_numeric


##### Reformat origin and rear information #####
tag_code_metadata <- read.csv("tag_code_metadata.csv")
# keep only the fish that are in the dataset
tag_code_metadata <- subset(tag_code_metadata, tag_code %in% states_complete$tag_code)

# Convert origins into numbers
# All relative to Yakima, since it's the last one alphabetically
origin_numeric <- data.frame(natal_origin = c("Asotin_Creek", 
                                        "Clearwater_River",
                                        "Deschutes_River", 
                                        "Entiat_River", 
                                        "Fifteenmile_Creek", 
                                        "Grande_Ronde_River", 
                                        "Hood_River",
                                        "Imnaha_River",
                                        "John_Day_River", 
                                        "Methow_River", 
                                        "Okanogan_River", 
                                        "Salmon_River", 
                                        "Tucannon_River", 
                                        "Umatilla_River",
                                        "Walla_Walla_River",
                                        "Wenatchee_River", 
                                        "Yakima_River"),
                             natal_origin_numeric = seq(1,17,1))

natal_origin_table <- read.csv("natal_origin_table.csv")
tag_code_metadata %>% 
  left_join(., natal_origin_table, by = "release_site_name") %>% 
  left_join(., origin_numeric, by = "natal_origin") %>% 
  mutate(rear_type_numeric = ifelse(rear_type_code %in% c("H", "U"), 2, 1))-> tag_code_metadata

# At this point, we need to recreate tag_code_metadata but with the tag_code_2 field
states_complete %>% 
  distinct(tag_code_2, .keep_all = TRUE) %>% 
  dplyr::select(tag_code, tag_code_2) -> tag_codes_2

# reformat this into origin_rear info
tag_codes_2 %>%
  left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin_numeric, rear_type_numeric), by = "tag_code") %>% 
  dplyr::rename(natal_origin = natal_origin_numeric, rear_type = rear_type_numeric) %>% 
  dplyr::select(-tag_code) -> origin_rear_actual

write.csv(origin_rear_actual,"middle_columbia_origin_rear_actual.csv")
  

fish_sim_cat_data_actual <- origin_rear_actual
  
  
# Store quantities for loop
# Store the total number of individuals
n.ind <- dim(state_data)[3]
  
  # Store the number of observations per individual
  # -1 because the last state is loss, which isn't actually an observation
  n.obs <- vector(length = n.ind)
  for (i in 1:n.ind){
    n.obs[i] <- sum(state_data[,,i]) - 1
  }
  
  
  
  # Get the state that each fish was in at each n.obs
  
  # First initialize an empty list
  states_list <- list()
  
  for (i in 1:n.ind){
    vec <- vector(length = (n.obs[i]-1))
    
    # Store in list
    states_list[[i]] <- vec
  }
  
  # Store with the states
  for (i in 1:n.ind){
    for (j in 1:(n.obs[i])){
      # states_list[[i]][j] <- rownames(as.data.frame(which(state_data[[i]][,j] == 1)))
      states_list[[i]][j] <- which(state_data[,j,i] == 1) # Get the index of the site instead of the name
    }
  }
  
  # Turn into matrix for stan
  states_mat <- matrix(nrow = n.ind, ncol = max(n.obs))
  for (i in 1:n.ind){
    states_mat[i,1:(n.obs[i])] <- states_list[[i]]
  }
  
  
  # Create the design matrix for categorical variables
  cat_X_mat_actual <- matrix(0, nrow = n.ind, ncol = 7)
  # Start it so that they're all 0s
  # The first column everyone gets a 1 (this is beta 0/grand mean mu)
  cat_X_mat_actual[,1] <- 1
  
  # This is for origin + rear
  for (i in 1:n.ind){
    # Rear type
    if (fish_sim_cat_data_actual$rear_type[i] == 1){
      cat_X_mat_actual[i,2] <- 1
    }
    else {
      cat_X_mat_actual[i,2] <- -1
    }
    # Natal origin
    if (fish_sim_cat_data_actual$natal_origin[i] == 3){ # Deschutes River
      cat_X_mat_actual[i,3] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == 5){ # Fifteenmile Creek
      cat_X_mat_actual[i,4] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == 9){ # John Day River
      cat_X_mat_actual[i,5] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == 14){ # Umatilla River
      cat_X_mat_actual[i,6] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == 17){ # Yakima River
      cat_X_mat_actual[i,7] <- 1
    }
    else { # for Walla Walla River
      cat_X_mat_actual[i,3:7] <- -1
    }
  }
  
  # One tweak: We have to replace all NAs in our input data for stan to accept it
  states_mat[is.na(states_mat)] <- -999
  
  
  # We are going to transform our detection histories to instead be vectors containing
  # the number of each site that was visited, instead of a matrix of sites with 0s for
  # not visited and 1s for visited
  
  # First, create an empty matrix to store the newly formatted detection histories
  state_data_2 <- matrix(0, nrow = dim(state_data)[3], ncol = dim(state_data)[2])
  
  # Now fill it in
  for (i in 1:dim(state_data)[3]){
    det_hist <- state_data[,,i]
    
    # Count the number of site visits
    nsite_visits <- sum(det_hist)
    
    for (j in 1:nsite_visits){
      state_data_2[i,j] <- which(det_hist[,j] == 1, arr.ind = TRUE)
    }
    
  }
    
    
  
  ##### Run stan model #####
  
  # step 0: data in a list #
  data <- list(y = state_data_2, n_ind = n.ind, n_obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, max_visits = dim(state_data_2)[2],
               movements = movements, not_movements = not_movements,
               nmovements = nmovements, # dates = dates,
               n_notmovements = n_notmovements, possible_states = transition_matrix, cat_X_mat = cat_X_mat_actual,
               grainsize = 1, N = dim(state_data_2)[1])
  
  
  print(Sys.time())
  
  # Fit stan model
  
  # fit <- stan(file = '01_stan_sim_int_only.stan', data = data)
  
  # Fit stan model using cmdstan
  # Step 1: load the model
  mod <- cmdstan_model("parallel_middle_columbia_02_stan_actual_int_origin.stan", compile = FALSE)
  
  # Step 2: Compile the model, set up to run in parallel
  mod$compile(cpp_options = list(stan_threads = TRUE))
  
  # Step 3: Run MCMC (HMC)
  fit <- mod$sample(
    data = data, 
    # seed = 123, # this seed gets stuck around 22-24, goes really fast and then at that iteration it slows way down
    # seed = 456,
    # chains = 3, 
    chains = 1,
    parallel_chains = 1,
    # parallel_chains = 3,
    refresh = 10, # print update every iter
    iter_sampling = 1000,
    iter_warmup = 1000,
    # iter_warmup = 10,
    # iter_sampling = 10,
    threads_per_chain = 28,
    init = 1
  )
  
saveRDS(fit, "parallel_middle_columbia_stan_actual_int_origin_stan_fit.rds")

# Troubleshoot our data
# Check to see if every transition in our model is represented
model_states %>% 
  as.data.frame() %>% 
  dplyr::rename(state = ".") %>% 
  mutate(index = row_number()) -> model_states_df

states_complete %>% 
  left_join(., model_states_df, by = "state") %>% 
  mutate(transition = ifelse(tag_code_2 == lag(tag_code_2), paste0(lag(state), " - ", state), NA)) %>% 
  mutate(transition_numeric = ifelse(tag_code_2 == lag(tag_code_2), paste0(lag(index), " - ", index), NA)) -> state_transitions

table(state_transitions$transition)
table(state_transitions$transition_numeric)
setdiff(model_states, unique(state_transitions$state))
