# 01_stan_actual_int_only

# This script will fit an intercept only model using stan to the actual dataset

# FOR TESTING: setwd
setwd("/Users/markusmin/Documents/CBR/steelhead/stan_actual/01_int_only/")

library("rstan")
rstan_options(auto_write = TRUE)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(bayesplot)

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
transition_matrix["mainstem, RIS to RRE", "loss"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, RIS to RRE"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, upstream of WEL"] <- 1
transition_matrix["mainstem, RIS to RRE", "Entiat River"] <- 1


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
transition_matrix["Upstream WEL other tributaries", "mainstem, BON to MCN"] <- 1


# Use the transition matrix to calculate the possible movements
possible_movements <- rowSums(transition_matrix)

# Get the indices that are 1s (except loss, since that's 1 minus the others)
movements <- which(transition_matrix[,1:(nstates-1)] == 1, arr.ind = TRUE)
nmovements <- dim(movements)[1]
# Now get all of the movements which are fixed to zero
not_movements <- which(transition_matrix[,1:(nstates-1)] == 0, arr.ind = TRUE)
n_notmovements <- dim(not_movements)[1]


##### LOAD and REFORMAT DATA #####

# Load states complete
states_complete <- read.csv(here::here("from_hyak_transfer", "2022-07-21-complete_det_hist", "states_complete.csv"), row.names = 1)

# first, get the maximum number of visits by any fish
states_complete %>% 
  group_by(tag_code) %>% 
  count() %>% 
  as.data.frame() -> site_visits_by_tag_code

max_visits <- max(site_visits_by_tag_code$n)

unique_tag_codes <- unique(states_complete$tag_code)
nfish <- length(unique_tag_codes)

# second, create an array with the dimensions number of fish, maximum number of site visits, number of states
state_data <- array(data = 0, dim = c(nstates, max_visits, nfish))

# third, populate the state data

# Add a column for the observation
states_complete %>% 
  group_by(tag_code) %>% 
  mutate(order = row_number()) -> states_complete

# add a column to index by fish
states_complete %>% 
  group_by(tag_code) %>% 
  mutate(tag_code_number = cur_group_id()) -> states_complete



# for (i in 1:nrow(states_complete)){
# Note: this loop will take about 3-4 hours to complete. It's primarily due to the which statement, since it has to search through and find a match
print(Sys.time())
for (i in 1:100){
  # index by 1) which state, 2) which visit, 3) which fish
  state_data[which(model_states == states_complete[i, "state", drop = TRUE]),states_complete[i, "order", drop = TRUE], states_complete[i, "tag_code_number", drop = TRUE]] <- 1 
  
}
print(Sys.time())

# Convert the dates into a numeric index

# First, populate the dates for implicit site visits
# Create a column to indicate if the date_time was interpolated
# Note the arrival date
states_complete %>% 
  mutate(date_source = ifelse(!is.na(date_time), "known_arrival", "interpolated")) %>% 
  mutate(date = NA) -> states_complete

# I think we have to loop this

# First, figure out indices of missing date_time
missing_date_time <- is.na(states_complete$date_time)
date_times <- states_complete$date_time

# For efficiency, only loop through the missing dates
missing_date_time_indices <- which(missing_date_time == TRUE)

# Extract date for all known date times
states_complete %>% 
  mutate(date =  format(as.Date(date(states_complete[i,"date_time", drop = TRUE]), origin = "1970-01-01"))) -> states_complete

# Loop to get arrival dates in state
# takes about 90 minutes to run
for (j in 1:length(missing_date_time_indices)){
# for (i in 1:100){
  i <- missing_date_time_indices[j]

    # Figure out how far back the last known time was
    # Truncate the missing date times to only ones prior to current state
    missing_date_time_subset <- missing_date_time[1:(i-1)]
    prev_time_index <- max(which(missing_date_time_subset %in% FALSE))
    
    # Truncate the missing date times to only ones after current state
    missing_date_time_subset <- missing_date_time[(i+1):nrow(states_complete)]
    next_time_index <- min(which(missing_date_time_subset %in% FALSE)) + i
    
    # Now, interpolate the missing time
    # First, figure out how long between the two known times
    prev_time <- ymd_hms(date_times[prev_time_index])
    next_time <- ymd_hms(date_times[next_time_index])
    time_diff <- next_time - prev_time
    
    # Get the missing time - add the time difference divided by the number 
    # of missing steps plus 1, multiply by which number step it is
    missing_time <- prev_time + (time_diff/(next_time_index - prev_time_index) * (i - prev_time_index))
    
    # Extract just the date
    missing_date <- date(missing_time)
    
    # populate the missing date_time
    states_complete[i, "date"] <- format(as.Date(missing_date, origin = "1970-01-01"))
}
  

# create an empty matrix
transition_date_matrix <- matrix(data = NA, nrow = nfish, ncol = max_visits)
# populate it with the dates
for (i in 1:nrow(states_complete)){
  transition_date_matrix[states_complete[i,"tag_code_number", drop = TRUE], states_complete[i, "order", drop = TRUE]] <- states_complete[i, "date_time", drop = TRUE]
}

# We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2005-05-31"))
transition_date_matrix %>% 
  as_tibble() %>% 
  mutate_all(date_numeric) -> transition_date_numeric



sim_1200_hist_list <- readRDS("sim_1200_hist_list.rds")
sim_1200_dates_list <- readRDS("sim_1200_dates_list.rds")
fish_sim_cat_data_1200 <- as.data.frame(as.matrix(read.csv("origin_rear_1200.csv")))


# Create a list to store JAGS objects
stan_1200_list <- list()
# Loop it
# for (z in 2:length(sim_1200_hist_list)){
for (z in 1:1){
  dates <- sim_1200_dates_list[[z]]
  sim_data <- sim_1200_hist_list[[z]]
  
  
  # Store quantities for loop
  # Store the total number of individuals
  n.ind <- dim(sim_data)[3]
  
  # Store the number of observations per individual
  # -1 because the last state is loss, which isn't actually an observation
  n.obs <- vector(length = n.ind)
  for (i in 1:n.ind){
    n.obs[i] <- sum(sim_data[,,i]) - 1
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
      # states_list[[i]][j] <- rownames(as.data.frame(which(sim_data[[i]][,j] == 1)))
      states_list[[i]][j] <- which(sim_data[,j,i] == 1) # Get the index of the site instead of the name
    }
  }
  
  # Turn into matrix for JAGS
  states_mat <- matrix(nrow = n.ind, ncol = max(n.obs))
  for (i in 1:n.ind){
    states_mat[i,1:(n.obs[i])] <- states_list[[i]]
  }
  
  
  # Create the design matrix for categorical variables
  cat_X_mat_1200 <- matrix(0, nrow = n.ind, ncol = 3)
  # The first column everyone gets a 1 (this is beta 0/grand mean mu)
  cat_X_mat_1200[,1] <- 1
  
  for (i in 1:n.ind){
    # Rear type
    # if (fish_sim_cat_data_1200$rear_type[i] == 1){
    #   cat_X_mat_1200[i,2] <- 1
    # }
    # else {
    #   cat_X_mat_1200[i,2] <- -1
    # }
    
    
    # Natal origin
    if (fish_sim_cat_data_1200$natal_origin[i] == 1){
      cat_X_mat_1200[i,2] <- 1
      cat_X_mat_1200[i,3] <- 0
    }
    else if (fish_sim_cat_data_1200$natal_origin[i] == 2){
      cat_X_mat_1200[i,2] <- 0
      cat_X_mat_1200[i,3] <- 1
    }
    else {
      cat_X_mat_1200[i,2] <- -1
      cat_X_mat_1200[i,3] <- -1
    }
  }
  
  # One tweak: We have to replace all NAs in our input data for stan to accept it
  states_mat[is.na(states_mat)] <- -999
  
  
  # We are going to transform our detection histories to instead be vectors containing
  # the number of each site that was visited, instead of a matrix of sites with 0s for
  # not visited and 1s for visited
  
  # First, create an empty matrix to store the newly formatted detection histories
  sim_data_2 <- matrix(0, nrow = dim(sim_data)[3], ncol = dim(sim_data)[2])
  
  # Now fill it in
  for (i in 1:dim(sim_data)[3]){
    det_hist <- sim_data[,,i]
    
    # Count the number of site visits
    nsite_visits <- sum(det_hist)
    
    for (j in 1:nsite_visits){
      sim_data_2[i,j] <- which(det_hist[,j] == 1, arr.ind = TRUE)
    }
    
  }
    
    
  
  ##### Run stan model #####
  
  # step 0: data in a list #
  data <- list(y = sim_data_2, n_ind = n.ind, n_obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, max_visits = dim(sim_data_2)[2],
               movements = movements, not_movements = not_movements,
               nmovements = nmovements, # dates = dates,
               n_notmovements = n_notmovements, possible_states = transition_matrix, cat_X_mat = cat_X_mat_1200)
  
  
  print(paste0("Run: ", z))
  print(Sys.time())
  
  # Fit stan model
  
  # fit <- stan(file = '01_stan_sim_int_only.stan', data = data)
  
  # Fit stan model using cmdstan
  # Step 1: load the model
  # mod <- cmdstan_model("01_stan_sim_int_only.stan", compile = FALSE)
  mod <- cmdstan_model("01_stan_sim_int_only_v2.stan", compile = FALSE)
  
  # Step 2: Compile the model, set up to run in parallel
  mod$compile(cpp_options = list(stan_threads = TRUE))
  
  # Step 3: Run MCMC (HMC)
  fit <- mod$sample(
    data = data, 
    seed = 123, 
    chains = 1, 
    parallel_chains = 1,
    refresh = 100, # print update every 10 iters
    iter_sampling = 1000,
    iter_warmup = 1000,
    threads_per_chain = 7,
    init = 1,
  )
 
stan_1200_list[[z]] <- fit$summary(variables = c("b0_matrix_2_1", 
                                               "b0_matrix_1_2",
                                               "b0_matrix_3_2",
                                               "b0_matrix_6_2",
                                               "b0_matrix_7_2",
                                               "b0_matrix_2_3",
                                               "b0_matrix_4_3",
                                               "b0_matrix_5_3",
                                               "b0_matrix_9_3",
                                               "b0_matrix_3_4",
                                               "b0_matrix_3_5",
                                               "b0_matrix_8_5",
                                               "b0_matrix_2_6",
                                               "b0_matrix_2_7",
                                               "b0_matrix_5_8",
                                               "b0_matrix_3_9")) 
}



# Inspect chains
draws <- fit$draws()
as_draws_matrix(draws) %>% 
  as.data.frame() -> draws_df

draws_df %>% 
  rownames_to_column("iter") %>% 
  mutate(iter = as.numeric(iter)) %>% 
  pivot_longer(cols = colnames(draws_df[,2:ncol(draws_df)])) -> draws_df_long

ggplot(subset(draws_df_long, name == "b0_matrix_2_1"), aes(x = iter, y = value)) +
  geom_line()

# Looks pretty good to me!

# Plot one example:
mcmc_hist(fit$draws("b0_matrix_2_1"))


##### Plot the same plots as before #####

simulation_plots_nocov <- function(stan_list){
  stan_runs_comp <- data.frame("variable" = NA, "mean" = NA, "q5" = NA, "q95" = NA)
  for (i in 1:length(stan_list)){
    as.data.frame(stan_list[[i]]) %>% 
      dplyr::select(variable, mean, q5, q95) %>% 
      mutate(run = paste0(i)) -> summary
    
    # summary %>% 
    #   mutate(run = factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) -> summary
    
    stan_runs_comp %>% 
      bind_rows(., summary) -> stan_runs_comp
  }
  
  
  # Create the plot
  stan_runs_comp <- subset(stan_runs_comp, !(is.na(variable))) 
  
  
  plot <- ggplot(stan_runs_comp, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q5, ymax = q95)) +
    geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~variable) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(0.5, 1, 1.5), limits = c(0,2)) +
    xlab("Run")
  
  return(plot)
}

# stan_1200_list <- readRDS(here::here("simulation", "stan_nocov_1200_list.rds"))
sim1200_plots <- simulation_plots_nocov(stan_list = stan_1200_list)
ggsave(here::here("stan_simulation", "01_int_only", "stan_origin_only_summary_plot.png"), height = 6, width = 10, sim1200_plots)

