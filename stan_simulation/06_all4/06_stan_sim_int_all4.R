# 06_stan_sim_int_temp

# This script will fit a model with four covariates to the simulated data
# The data is the same as was used for fitting the JAGS model.

# FOR TESTING: setwd
setwd("/Users/markusmin/Documents/CBR/steelhead/stan_simulation/06_all4")

library("rstan")
rstan_options(auto_write = TRUE)
library(cmdstanr)
library(posterior)
library(tidyverse)
library(bayesplot)

##### Get all of the data in order #####

temp_sim <- as.matrix(read.csv("temp_600.csv", row.names = 1))
flow_sim <- as.matrix(read.csv("flow_600.csv", row.names = 1))


# Get number of possible movements from each site
possible_movements <- c("mainstem, mouth to BON" = 2,
                        "mainstem, BON to MCN" = 5,
                        "mainstem, MCN to ICH or PRA" = 5,
                        "mainstem, PRA to RIS" = 2,
                        "mainstem, ICH to LGR" = 3,
                        "Deschutes River" = 2,
                        "John Day River" = 2,
                        "Tucannon River" = 2,
                        "Yakima River" = 2,
                        "loss" = 0)

possible_movements <- c(2,5,5,2,3,2,2,2,2,0)

# Indexing covariate data
cov_index_mat <- matrix(data = c(1, NA, NA, NA,
                                 1, 2, 6, 5,
                                 2, 3, 4, 7,
                                 3, NA, NA, NA,
                                 4, 8, NA, NA,
                                 6, NA, NA, NA,
                                 5, NA, NA, NA,
                                 7, NA, NA, NA,
                                 8, NA, NA, NA), 
                        ncol = 4, nrow = 9, byrow = TRUE)

# For temperature, there is only one temperature per state
# The issue arises again with the branching MCN to ICH or PRA state - which temperature do we choose?
# For now, let's use MCN, since it's a combination of both
# states = same order as possible_movements
temp_index <- c(1, 2, 2, 3, 4, 6, 5, 8, 7)

# Index flow - one flow for state for now, for simplicity
flow_index <- c(1, 2, 2, 3, 4, 6, 5, 8, 7)
# Note that in the simulation, the relationships are currently different, so we'll have to fix that

# Create a transition matrix of 1s and 0s for movements from (rows) to (columns)


# Create a starting matrix of sites, with the first occasion populated (BON to MCN)
sim_states = c("mainstem, mouth to BON", 
               "mainstem, BON to MCN", 
               "mainstem, MCN to ICH or PRA", 
               "mainstem, PRA to RIS",
               "mainstem, ICH to LGR",
               "Deschutes River", 
               "John Day River", 
               "Tucannon River", 
               "Yakima River",
               "loss")

nstates <- length(sim_states)

# Create a rows = from, columns = to matrix for movement probabilities

transition_matrix <- matrix(0, nrow = nstates, ncol = nstates)
rownames(transition_matrix) <- sim_states
colnames(transition_matrix) <- sim_states

# Populate every possible option with a 1
# mouth to BON
transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- 1
transition_matrix["mainstem, mouth to BON", "loss"] <- 1

# BON to MCN
transition_matrix["mainstem, BON to MCN", "loss"] <- 1
transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- 1
transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- 1
transition_matrix["mainstem, BON to MCN", "John Day River"] <- 1


# MCN to ICH or PRA
transition_matrix["mainstem, MCN to ICH or PRA", "loss"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, BON to MCN"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River"] <- 1

# PRA to RIS
transition_matrix["mainstem, PRA to RIS", "loss"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- 1

# ICH to LGR
transition_matrix["mainstem, ICH to LGR", "loss"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- 1

# Deschutes River
transition_matrix["Deschutes River", "loss"] <- 1
transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- 1

# John Day River
transition_matrix["John Day River", "loss"] <- 1
transition_matrix["John Day River", "mainstem, BON to MCN"] <- 1

# Yakima River
transition_matrix["Yakima River", "loss"] <- 1
transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- 1

# Tucannon River
transition_matrix["Tucannon River", "loss"] <- 1
transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- 1

# Get the indices that are 1s (except loss, since that's 1 minus the others)
movements <- which(transition_matrix[,1:9] == 1, arr.ind = TRUE)
nmovements <- dim(movements)[1]
# Now get all of the movements which are fixed to zero
not_movements <- which(transition_matrix[,1:9] == 0, arr.ind = TRUE)
n_notmovements <- dim(not_movements)[1]



sim_1200_hist_list <- readRDS("sim_1200_all4_hist_list.rds")
sim_1200_dates_list <- readRDS("sim_1200_all4_dates_list.rds")
fish_sim_cat_data_1200 <- as.data.frame(as.matrix(read.csv("origin_rear_1200.csv")))


# Create a list to store JAGS objects
stan_1200_all4_list <- list()
# Loop it
for (z in 1:length(sim_1200_hist_list)){
# for (z in 1:1){
  dates <- as.matrix(sim_1200_dates_list[[z]])
  # Can't have NA values in data in stan, so we'll replace them all with -999
  dates[is.na(dates)] <- -999
  
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
  cat_X_mat_1200 <- matrix(0, nrow = n.ind, ncol = 4)
  # The first column everyone gets a 1 (this is beta 0/grand mean mu)
  cat_X_mat_1200[,1] <- 1
  
  for (i in 1:n.ind){
    # Rear type
    if (fish_sim_cat_data_1200$rear_type[i] == 1){
      cat_X_mat_1200[i,2] <- 1
    }
    else {
      cat_X_mat_1200[i,2] <- -1
    }
    
    
    # Natal origin
    if (fish_sim_cat_data_1200$natal_origin[i] == 1){
      cat_X_mat_1200[i,3] <- 1
      cat_X_mat_1200[i,4] <- 0
    }
    else if (fish_sim_cat_data_1200$natal_origin[i] == 2){
      cat_X_mat_1200[i,3] <- 0
      cat_X_mat_1200[i,4] <- 1
    }
    else {
      cat_X_mat_1200[i,3] <- -1
      cat_X_mat_1200[i,4] <- -1
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

  
  ##### Data #####
  data <- list(y = sim_data_2, n_ind = n.ind, n_obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, max_visits = dim(sim_data_2)[2],
               movements = movements, not_movements = not_movements,
               nmovements = nmovements, dates = dates, temp_index = temp_index, temp_matrix = temp_sim,
               flow_index = flow_index, flow_matrix = flow_sim,
               n_notmovements = n_notmovements, possible_states = transition_matrix, cat_X_mat = cat_X_mat_1200)
  
  
  print(paste0("Run: ", z))
  print(Sys.time())
  
  # Fit stan model
  
  # fit <- stan(file = '01_stan_sim_int_only.stan', data = data)
  
  # Fit stan model using cmdstan
  # Step 1: load the model
  mod <- cmdstan_model("06_stan_sim_int_all4.stan", compile = FALSE)
  
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
  
  stan_1200_all4_list[[z]] <- fit$summary(variables = c("b0_matrix_2_1", 
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
                                                   "b0_matrix_3_9",
                                                   "btemp_matrix_2_1", 
                                                   "brear_matrix_2_1", 
                                                   "brear_matrix_1_2",
                                                   "brear_matrix_3_2",
                                                   "brear_matrix_6_2",
                                                   "brear_matrix_7_2",
                                                   "brear_matrix_2_3",
                                                   "brear_matrix_4_3",
                                                   "brear_matrix_5_3",
                                                   "brear_matrix_9_3",
                                                   "brear_matrix_3_4",
                                                   "brear_matrix_3_5",
                                                   "brear_matrix_8_5",
                                                   "brear_matrix_2_6",
                                                   "brear_matrix_2_7",
                                                   "brear_matrix_5_8",
                                                   "brear_matrix_3_9",
                                                   "borigin1_matrix_2_1", 
                                                   "borigin1_matrix_1_2",
                                                   "borigin1_matrix_3_2",
                                                   "borigin1_matrix_6_2",
                                                   "borigin1_matrix_7_2",
                                                   "borigin1_matrix_2_3",
                                                   "borigin1_matrix_4_3",
                                                   "borigin1_matrix_5_3",
                                                   "borigin1_matrix_9_3",
                                                   "borigin1_matrix_3_4",
                                                   "borigin1_matrix_3_5",
                                                   "borigin1_matrix_8_5",
                                                   "borigin1_matrix_2_6",
                                                   "borigin1_matrix_2_7",
                                                   "borigin1_matrix_5_8",
                                                   "borigin1_matrix_3_9",
                                                   "borigin2_matrix_2_1", 
                                                   "borigin2_matrix_1_2",
                                                   "borigin2_matrix_3_2",
                                                   "borigin2_matrix_6_2",
                                                   "borigin2_matrix_7_2",
                                                   "borigin2_matrix_2_3",
                                                   "borigin2_matrix_4_3",
                                                   "borigin2_matrix_5_3",
                                                   "borigin2_matrix_9_3",
                                                   "borigin2_matrix_3_4",
                                                   "borigin2_matrix_3_5",
                                                   "borigin2_matrix_8_5",
                                                   "borigin2_matrix_2_6",
                                                   "borigin2_matrix_2_7",
                                                   "borigin2_matrix_5_8",
                                                   "borigin2_matrix_3_9",
                                                   
                                                   "btemp_matrix_1_2",
                                                   "btemp_matrix_3_2",
                                                   
                                                   
                                                   "btemp_matrix_2_3",
                                                   "btemp_matrix_4_3",
                                                   "btemp_matrix_5_3",
                                                   
                                                   "btemp_matrix_3_4",
                                                   "btemp_matrix_3_5",
                                                   
                                                   "btemp_matrix_2_6",
                                                   "btemp_matrix_2_7",
                                                   "btemp_matrix_5_8",
                                                   "btemp_matrix_3_9",
                                                   
                                                   "bflow_matrix_2_1", 
                                                   "bflow_matrix_1_2",
                                                   "bflow_matrix_3_2",
                                                   
                                                   
                                                   "bflow_matrix_2_3",
                                                   "bflow_matrix_4_3",
                                                   "bflow_matrix_5_3",
                                                   
                                                   "bflow_matrix_3_4",
                                                   "bflow_matrix_3_5",
                                                   
                                                   "bflow_matrix_2_6",
                                                   "bflow_matrix_2_7",
                                                   "bflow_matrix_5_8",
                                                   "bflow_matrix_3_9")) 
}

fit$summary() 


##### Plot the same plots as before #####

simulation_plots_all4 <- function(stan_list){
  stan_runs_comp_b0 <- data.frame("variable" = NA, "mean" = NA, "q5" = NA, "q95" = NA)
  stan_runs_comp_borigin1 <- data.frame("variable" = NA, "mean" = NA, "q5" = NA, "q95" = NA)
  stan_runs_comp_borigin2 <- data.frame("variable" = NA, "mean" = NA, "q5" = NA, "q95" = NA)
  stan_runs_comp_brear <- data.frame("variable" = NA, "mean" = NA, "q5" = NA, "q95" = NA)
  stan_runs_comp_btemp <- data.frame("variable" = NA, "mean" = NA, "q5" = NA, "q95" = NA)
  stan_runs_comp_bflow <- data.frame("variable" = NA, "mean" = NA, "q5" = NA, "q95" = NA)
  
  for (i in 1:length(stan_list)){
    as.data.frame(stan_list[[i]]) %>% 
      dplyr::select(variable, mean, q5, q95) %>% 
      mutate(run = paste0(i)) -> summary
    
    # Change run to a factor for plotting
    summary %>% 
      mutate(run = factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) -> summary
    
    # Split the summary by parameter/variable
    summary %>% 
      filter(grepl("b0", variable)) -> b0_summary
    summary %>% 
      filter(grepl("borigin1", variable))  -> borigin1_summary
    summary %>% 
      filter(grepl("borigin2", variable))  -> borigin2_summary
    summary %>% 
      filter(grepl("brear", variable))  -> brear_summary
    summary %>% 
      filter(grepl("btemp", variable))  -> btemp_summary
    summary %>% 
      filter(grepl("bflow", variable))  -> bflow_summary
    
    
    stan_runs_comp_b0 %>% 
      bind_rows(., b0_summary) -> stan_runs_comp_b0
    
    stan_runs_comp_borigin1 %>% 
      bind_rows(., borigin1_summary) -> stan_runs_comp_borigin1
    
    stan_runs_comp_borigin2 %>% 
      bind_rows(., borigin2_summary) -> stan_runs_comp_borigin2
    
    stan_runs_comp_brear %>% 
      bind_rows(., brear_summary) -> stan_runs_comp_brear
    
    stan_runs_comp_btemp %>% 
      bind_rows(., btemp_summary) -> stan_runs_comp_btemp
    
    stan_runs_comp_bflow %>% 
      bind_rows(., bflow_summary) -> stan_runs_comp_bflow
  }
  
  
  # Create the plot
  stan_runs_comp_b0 <- subset(stan_runs_comp_b0, !(is.na(variable))) 
  stan_runs_comp_borigin1 <- subset(stan_runs_comp_borigin1, !(is.na(variable))) 
  stan_runs_comp_borigin2 <- subset(stan_runs_comp_borigin2, !(is.na(variable))) 
  stan_runs_comp_brear <- subset(stan_runs_comp_brear, !(is.na(variable))) 
  stan_runs_comp_btemp <- subset(stan_runs_comp_btemp, !(is.na(variable))) 
  stan_runs_comp_bflow <- subset(stan_runs_comp_bflow, !(is.na(variable))) 
  
  # Load actual values
  
  actual_values_borigin1 <- data.frame(variable = c("borigin1_matrix_1_2",
                                                    "borigin1_matrix_2_1",
                                                    "borigin1_matrix_2_3",
                                                    "borigin1_matrix_2_6",
                                                    "borigin1_matrix_2_7",
                                                    "borigin1_matrix_3_2",
                                                    "borigin1_matrix_3_4",
                                                    "borigin1_matrix_3_5",
                                                    "borigin1_matrix_3_9",
                                                    "borigin1_matrix_4_3",
                                                    "borigin1_matrix_5_3",
                                                    "borigin1_matrix_5_8",
                                                    "borigin1_matrix_6_2",
                                                    "borigin1_matrix_7_2",
                                                    "borigin1_matrix_8_5",
                                                    "borigin1_matrix_9_3"
  ),
  yint = c(0,0,-0.5,0.25,0.5,0.5,0,-0.5,-0.25,0,0.25,-0.25,0,-0.5,0.25,0.25))
  
  actual_values_borigin2 <- data.frame(variable = c("borigin2_matrix_1_2",
                                                    "borigin2_matrix_2_1",
                                                    "borigin2_matrix_2_3",
                                                    "borigin2_matrix_2_6",
                                                    "borigin2_matrix_2_7",
                                                    "borigin2_matrix_3_2",
                                                    "borigin2_matrix_3_4",
                                                    "borigin2_matrix_3_5",
                                                    "borigin2_matrix_3_9",
                                                    "borigin2_matrix_4_3",
                                                    "borigin2_matrix_5_3",
                                                    "borigin2_matrix_5_8",
                                                    "borigin2_matrix_6_2",
                                                    "borigin2_matrix_7_2",
                                                    "borigin2_matrix_8_5",
                                                    "borigin2_matrix_9_3"
  ),
  yint = c(0,0,0.25,-0.125,-0.25,-0.25,0,0,0.5,0,0.25,-0.25,0,0.25,0.25,-0.5))
  
  actual_values_brear <- data.frame(variable = c("brear_matrix_1_2",
                                                 "brear_matrix_2_1",
                                                 "brear_matrix_2_3",
                                                 "brear_matrix_2_6",
                                                 "brear_matrix_2_7",
                                                 "brear_matrix_3_2",
                                                 "brear_matrix_3_4",
                                                 "brear_matrix_3_5",
                                                 "brear_matrix_3_9",
                                                 "brear_matrix_4_3",
                                                 "brear_matrix_5_3",
                                                 "brear_matrix_5_8",
                                                 "brear_matrix_6_2",
                                                 "brear_matrix_7_2",
                                                 "brear_matrix_8_5",
                                                 "brear_matrix_9_3"
                                                 
  ),
  yint = c(0,0,0,0.5,0,0,0,0,-0.2,0,0,0,-0.5,0,0,0))
  
  # Less values here because no trib temps
  actual_values_btemp <- data.frame(variable = c("btemp_matrix_1_2",
                                            "btemp_matrix_2_1",
                                            "btemp_matrix_2_3",
                                            "btemp_matrix_2_6",
                                            "btemp_matrix_2_7",
                                            "btemp_matrix_3_2",
                                            "btemp_matrix_3_4",
                                            "btemp_matrix_3_5",
                                            "btemp_matrix_3_9",
                                            "btemp_matrix_4_3",
                                            "btemp_matrix_5_3",
                                            "btemp_matrix_5_8"
  ),
  yint = c(0,0,0.5,0,0,-0.5,0.3,1,0,0,0.2,0))
  
  # Less values here because no trib flows
  actual_values_bflow <- data.frame(variable = c("bflow_matrix_1_2",
                                                 "bflow_matrix_2_1",
                                                 "bflow_matrix_2_3",
                                                 "bflow_matrix_2_6",
                                                 "bflow_matrix_2_7",
                                                 "bflow_matrix_3_2",
                                                 "bflow_matrix_3_4",
                                                 "bflow_matrix_3_5",
                                                 "bflow_matrix_3_9",
                                                 "bflow_matrix_4_3",
                                                 "bflow_matrix_5_3",
                                                 "bflow_matrix_5_8"
  ),
  yint = c(0,0.5,0,0,0,0.05,0,-0.3,0,0,0,0))
  

  
  
  b0_plot <- ggplot(stan_runs_comp_b0, aes(y = mean, x = run)) +
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
  
  borigin1_plot <- ggplot(stan_runs_comp_borigin1, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q5, ymax = q95)) +
    geom_hline(data = actual_values_borigin1, aes(yintercept = yint), lty = 2) +
    facet_wrap(~variable) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0 ,1), limits = c(-2,2)) +
    xlab("Run")
  
  borigin2_plot <- ggplot(stan_runs_comp_borigin2, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q5, ymax = q95)) +
    geom_hline(data = actual_values_borigin2, aes(yintercept = yint), lty = 2) +
    facet_wrap(~variable) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0 ,1), limits = c(-2,2)) +
    xlab("Run")
  
  brear_plot <- ggplot(stan_runs_comp_brear, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q5, ymax = q95)) +
    geom_hline(data = actual_values_brear, aes(yintercept = yint), lty = 2) +
    facet_wrap(~variable) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0 ,1), limits = c(-2,2)) +
    xlab("Run")

  btemp_plot <- ggplot(stan_runs_comp_btemp, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q5, ymax = q95)) +
    geom_hline(data = actual_values_btemp, aes(yintercept = yint), lty = 2) +
    facet_wrap(~variable) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0 ,1), limits = c(-2,2)) +
    xlab("Run")
  
  bflow_plot <- ggplot(stan_runs_comp_bflow, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q5, ymax = q95)) +
    geom_hline(data = actual_values_bflow, aes(yintercept = yint), lty = 2) +
    facet_wrap(~variable) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0 ,1), limits = c(-2,2)) +
    xlab("Run")
  
  return(list(b0_plot, brear_plot, borigin1_plot, borigin2_plot, btemp_plot, bflow_plot))
}

# stan_1200_list <- readRDS(here::here("simulation", "stan_nocov_1200_list.rds"))
sim1200_all4_plots <- simulation_plots_all4(stan_list = stan_1200_all4_list)
ggsave(here::here("stan_simulation", "06_all4", "stan_all4_summary_plot_b0.png"), height = 6, width = 10, sim1200_all4_plots[[1]])
ggsave(here::here("stan_simulation", "06_all4", "stan_all4_summary_plot_brear.png"), height = 6, width = 10, sim1200_all4_plots[[2]])
ggsave(here::here("stan_simulation", "06_all4", "stan_all4_summary_plot_borigin1.png"), height = 6, width = 10, sim1200_all4_plots[[3]])
ggsave(here::here("stan_simulation", "06_all4", "stan_all4_summary_plot_borigin2.png"), height = 6, width = 10, sim1200_all4_plots[[4]])
ggsave(here::here("stan_simulation", "06_all4", "stan_all4_summary_plot_btemp.png"), height = 6, width = 10, sim1200_all4_plots[[5]])
ggsave(here::here("stan_simulation", "06_all4", "stan_all4_summary_plot_bflow.png"), height = 6, width = 10, sim1200_all4_plots[[6]])
