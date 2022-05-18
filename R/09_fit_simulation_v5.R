# 09_fit_simulation

# NOTE: In this script, we will be fitting the different simulated datasets.
# The JAGS code is exported from the v4 script, we will just be running it here

# Load libraries
library(tidyverse)
library(lubridate)
library(janitor)
library(rjags)
library(jagsUI)
library(R2jags)
library(coda)

##### Universal data (these are used to parameterize the model) #####
temp_sim <- as.matrix(read.csv(here::here("simulation", "temp_600.csv"), row.names = 1))
flow_sim <- as.matrix(read.csv(here::here("simulation", "flow_600.csv"), row.names = 1))
# flow_sim <- flow_sim_zscore_datenumeric


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


##### 600 fish, no covariates #####
sim_600_hist_list <- readRDS(here::here("simulation", "sim_600_hist_list.rds"))
sim_600_dates_list <- readRDS(here::here("simulation", "sim_600_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_600.csv")))

# Create a list to store JAGS objects
JAGS_600_list <- list()
# Loop it
for (z in 1:length(sim_600_hist_list)){
  dates <- sim_600_dates_list[[z]]
  sim_data <- sim_600_hist_list[[z]]
  
  
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, origin = fish_sim_cat_data[,2], rear = fish_sim_cat_data[,3], 
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix
    ))
  }
  
  
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_nocov.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_600_list[[z]] <- out.jags
}

saveRDS(JAGS_600_list, here::here("simulation", "JAGS_nocov_600_list.rds"))

##### 1200 fish, no covariates #####
sim_1200_hist_list <- readRDS(here::here("simulation", "sim_1200_hist_list.rds"))
sim_1200_dates_list <- readRDS(here::here("simulation", "sim_1200_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_1200.csv")))

# Create a list to store JAGS objects
JAGS_1200_list <- list()
# Loop it
for (z in 1:length(sim_1200_hist_list)){
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, origin = fish_sim_cat_data[,2], rear = fish_sim_cat_data[,3], 
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix
    ))
  }
  
  
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_nocov.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_1200_list[[z]] <- out.jags
}

saveRDS(JAGS_1200_list, here::here("simulation", "JAGS_nocov_1200_list.rds"))

##### 3000 fish, no covariates #####
sim_3000_hist_list <- readRDS(here::here("simulation", "sim_3000_hist_list.rds"))
sim_3000_dates_list <- readRDS(here::here("simulation", "sim_3000_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_3000.csv")))

# Create a list to store JAGS objects
JAGS_3000_list <- list()
# Loop it
for (z in 1:length(sim_3000_hist_list)){
  dates <- sim_3000_dates_list[[z]]
  sim_data <- sim_3000_hist_list[[z]]
  
  
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, origin = fish_sim_cat_data[,2], rear = fish_sim_cat_data[,3], 
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix
    ))
  }
  
  
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_nocov.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_3000_list[[z]] <- out.jags
}

saveRDS(JAGS_3000_list, here::here("simulation", "JAGS_nocov_3000_list.rds"))





##### 600 fish, covariates #####
sim_600_hist_list <- readRDS(here::here("simulation", "sim_600_cov_hist_list.rds"))
sim_600_dates_list <- readRDS(here::here("simulation", "sim_600_cov_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_600.csv")))

# Create a list to store JAGS objects
JAGS_600_list <- list()
# Loop it
for (z in 1:length(sim_600_hist_list)){
  dates <- sim_600_dates_list[[z]]
  sim_data <- sim_600_hist_list[[z]]
  
  
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "bflow_matrix",
    "btemp_matrix",
    "brear_matrix",
    "borigin_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, origin = fish_sim_cat_data[,2], rear = fish_sim_cat_data[,3], 
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    bflow_matrix <- matrix(NA, nrow = 10, ncol = 9)
    btemp_matrix <- matrix(NA, nrow = 10, ncol = 9)
    brear_matrix <- matrix(NA, nrow = 10, ncol = 9)
    borigin_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      bflow_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      btemp_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      brear_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      borigin_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      bflow_matrix = bflow_matrix,
      btemp_matrix = btemp_matrix,
      brear_matrix = brear_matrix,
      borigin_matrix = borigin_matrix
    ))
  }
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_600_list[[z]] <- out.jags
}

saveRDS(JAGS_600_list, here::here("simulation", "JAGS_cov_600_list.rds"))

##### 1200 fish, covariates #####
sim_1200_hist_list <- readRDS(here::here("simulation", "sim_1200_cov_hist_list.rds"))
sim_1200_dates_list <- readRDS(here::here("simulation", "sim_1200_cov_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_1200.csv")))

# Create a list to store JAGS objects
JAGS_1200_list <- list()
# Loop it
# for (z in 1:length(sim_1200_hist_list)){
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "bflow_matrix",
    "btemp_matrix",
    "brear_array",
    "borigin_array"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, origin = fish_sim_cat_data[,2], rear = fish_sim_cat_data[,3], 
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    bflow_matrix <- matrix(NA, nrow = 10, ncol = 9)
    btemp_matrix <- matrix(NA, nrow = 10, ncol = 9)
    brear_array <- array(NA, dim = c(10, 9, 2))
    borigin_array <- array(NA, dim = c(10, 9, 3))
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      bflow_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      btemp_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      brear_array[movements[j,1], movements[j,2],1] <- runif(1,-1,1)
      brear_array[movements[j,1], movements[j,2],2] <- runif(1,-1,1)
      borigin_array[movements[j,1], movements[j,2],1] <- runif(1,-1,1)
      borigin_array[movements[j,1], movements[j,2],2] <- runif(1,-1,1)
      borigin_array[movements[j,1], movements[j,2],3] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      bflow_matrix = bflow_matrix,
      btemp_matrix = btemp_matrix,
      brear_array = brear_array,
      borigin_array = borigin_array
    ))
  }
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 25000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_1200_list[[z]] <- out.jags 
}
saveRDS(out.jags, here::here("simulation", "JAGS_cov_1200_onerun.rds"))
# saveRDS(JAGS_1200_list, here::here("simulation", "JAGS_cov_1200_list.rds"))

##### 3000 fish, covariates #####
sim_3000_hist_list <- readRDS(here::here("simulation", "sim_3000_cov_hist_list.rds"))
sim_3000_dates_list <- readRDS(here::here("simulation", "sim_3000_cov_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_3000.csv")))

# Create a list to store JAGS objects
JAGS_3000_list <- list()
# Loop it
for (z in 1:length(sim_3000_hist_list)){
  dates <- sim_3000_dates_list[[z]]
  sim_data <- sim_3000_hist_list[[z]]
  
  
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "bflow_matrix",
    "btemp_matrix",
    "brear_matrix",
    "borigin_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, origin = fish_sim_cat_data[,2], rear = fish_sim_cat_data[,3], 
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    bflow_matrix <- matrix(NA, nrow = 10, ncol = 9)
    btemp_matrix <- matrix(NA, nrow = 10, ncol = 9)
    brear_matrix <- matrix(NA, nrow = 10, ncol = 9)
    borigin_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      bflow_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      btemp_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      brear_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      borigin_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      bflow_matrix = bflow_matrix,
      btemp_matrix = btemp_matrix,
      brear_matrix = brear_matrix,
      borigin_matrix = borigin_matrix
    ))
  }
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_3000_list[[z]] <- out.jags
}

saveRDS(JAGS_3000_list, here::here("simulation", "JAGS_cov_3000_list.rds"))







##### 600 fish, continuous covariates #####
sim_600_hist_list <- readRDS(here::here("simulation", "sim_600_cov_continuous_hist_list.rds"))
sim_600_dates_list <- readRDS(here::here("simulation", "sim_600_cov_continuous_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_600.csv")))

# Create a list to store JAGS objects
JAGS_600_list <- list()
# Loop it
# for (z in 1:length(sim_600_hist_list)){
for (z in 1:1){
  dates <- sim_600_dates_list[[z]]
  sim_data <- sim_600_hist_list[[z]]
  
  
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "bflow_matrix",
    "btemp_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data, n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat,
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    bflow_matrix <- matrix(NA, nrow = 10, ncol = 9)
    btemp_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      bflow_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      btemp_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      bflow_matrix = bflow_matrix,
      btemp_matrix = btemp_matrix
    ))
  }
  
  print(paste0("Run: ", z))
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov_continuous.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_600_list[[z]] <- out.jags
}

saveRDS(out.jags, here::here("simulation", "JAGS_cov_continuous_600_onerun.rds"))
saveRDS(JAGS_600_list, here::here("simulation", "JAGS_cov_continuous_600_list.rds"))

##### 1200 fish, continuous covariates #####
sim_1200_hist_list <- readRDS(here::here("simulation", "sim_1200_cov_continuous_hist_list.rds"))
sim_1200_dates_list <- readRDS(here::here("simulation", "sim_1200_cov_continuous_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_1200.csv")))

# Create a list to store JAGS objects
JAGS_1200_list <- list()
# Loop it
# for (z in 1:length(sim_1200_hist_list)){
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "bflow_matrix",
    "btemp_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat,
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    bflow_matrix <- matrix(NA, nrow = 10, ncol = 9)
    btemp_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      bflow_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      btemp_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      bflow_matrix = bflow_matrix,
      btemp_matrix = btemp_matrix
    ))
  }
  
  print(paste0("Run: ", z))
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov_continuous.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_1200_list[[z]] <- out.jags
}

saveRDS(out.jags, here::here("simulation", "JAGS_cov_continuous_1200_onerun.rds"))
saveRDS(JAGS_1200_list, here::here("simulation", "JAGS_cov_continuous_1200_list.rds"))


##### 600 fish, temperature only #####
sim_600_hist_list <- readRDS(here::here("simulation", "sim_600_cov_temp_hist_list.rds"))
sim_600_dates_list <- readRDS(here::here("simulation", "sim_600_cov_temp_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_600.csv")))

# Create a list to store JAGS objects
JAGS_600_list <- list()
# Loop it
for (z in 1:length(sim_600_hist_list)){
# for (z in 1:1){
  dates <- sim_600_dates_list[[z]]
  sim_data <- sim_600_hist_list[[z]]
  
  
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
  
  
  
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "btemp_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat,
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    btemp_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      btemp_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      btemp_matrix = btemp_matrix
    ))
  }
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov_temp.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_600_list[[z]] <- out.jags
}

saveRDS(out.jags, here::here("simulation", "JAGS_cov_temp_600_onerun.rds"))
saveRDS(JAGS_600_list, here::here("simulation", "JAGS_cov_temp_600_list.rds"))






##### 600 fish, CATEGORICAL covariates #####
sim_600_hist_list <- readRDS(here::here("simulation", "sim_600_cov_categorical_hist_list.rds"))
sim_600_dates_list <- readRDS(here::here("simulation", "sim_600_cov_categorical_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_600.csv")))

# Create a list to store JAGS objects
JAGS_600_list <- list()
# Loop it
# for (z in 1:length(sim_600_hist_list)){
for (z in 1:3){
  dates <- sim_600_dates_list[[z]]
  sim_data <- sim_600_hist_list[[z]]
  
  
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
  cat_X_mat_600 <- matrix(0, nrow = n.ind, ncol = 4)
  # The first column everyone gets a 1 (this is beta 0/grand mean mu)
  cat_X_mat_600[,1] <- 1
  
  for (i in 1:n.ind){
    # Rear type
    if (fish_sim_cat_data_600$rear_type[i] == 1){
      cat_X_mat_600[i,2] <- 1
    }
    else {
      cat_X_mat_600[i,2] <- -1
    }
    
    
    # Natal origin
    if (fish_sim_cat_data_600$natal_origin[i] == 1){
      cat_X_mat_600[i,3] <- 1
      cat_X_mat_600[i,4] <- 0
    }
    else if (fish_sim_cat_data_600$natal_origin[i] == 2){
      cat_X_mat_600[i,3] <- 0
      cat_X_mat_600[i,4] <- 1
    }
    else {
      cat_X_mat_600[i,3] <- -1
      cat_X_mat_600[i,4] <- -1
    }
  }
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "brear_matrix",
    "borigin1_matrix",
    "borigin2_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat,
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix, cat_X_mat = cat_X_mat_600)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    brear_matrix <- matrix(NA, nrow = 10, ncol = 9)
    borigin1_matrix <- matrix(NA, nrow = 10, ncol = 9)
    borigin2_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      brear_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      borigin1_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      borigin2_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      brear_matrix = brear_matrix,
      borigin1_matrix = borigin1_matrix,
      borigin2_matrix = borigin2_matrix
    ))
  }
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov_categorical.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 10000, n.burnin = 5000,
      n.thin = 10,
      jags.seed = 123
    )
  
  JAGS_600_list[[z]] <- out.jags
}

saveRDS(out.jags, here::here("simulation", "JAGS_cov_categorical_600_onerun.rds"))
saveRDS(JAGS_600_list, here::here("simulation", "JAGS_cov_categorical_600_list.rds"))














##### 600 fish, ORIGIN only #####
sim_600_hist_list <- readRDS(here::here("simulation", "sim_600_cov_origin_hist_list.rds"))
sim_600_dates_list <- readRDS(here::here("simulation", "sim_600_cov_origin_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_600.csv")))

# Create a list to store JAGS objects
JAGS_600_list <- list()
# Loop it
for (z in 1:length(sim_600_hist_list)){
# for (z in 1:1){
  dates <- sim_600_dates_list[[z]]
  sim_data <- sim_600_hist_list[[z]]
  
  
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
  cat_X_mat_600 <- matrix(0, nrow = n.ind, ncol = 3)
  # The first column everyone gets a 1 (this is beta 0/grand mean mu)
  cat_X_mat_600[,1] <- 1
  
  for (i in 1:n.ind){
    # Rear type
    # if (fish_sim_cat_data_600$rear_type[i] == 1){
    #   cat_X_mat_600[i,2] <- 1
    # }
    # else {
    #   cat_X_mat_600[i,2] <- -1
    # }
    
    
    # Natal origin
    if (fish_sim_cat_data_600$natal_origin[i] == 1){
      cat_X_mat_600[i,2] <- 1
      cat_X_mat_600[i,3] <- 0
    }
    else if (fish_sim_cat_data_600$natal_origin[i] == 2){
      cat_X_mat_600[i,2] <- 0
      cat_X_mat_600[i,3] <- 1
    }
    else {
      cat_X_mat_600[i,2] <- -1
      cat_X_mat_600[i,3] <- -1
    }
  }
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "borigin1_matrix",
    "borigin2_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat,
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix, cat_X_mat = cat_X_mat_600)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    borigin1_matrix <- matrix(NA, nrow = 10, ncol = 9)
    borigin2_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      borigin1_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      borigin2_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      borigin1_matrix = borigin1_matrix,
      borigin2_matrix = borigin2_matrix
    ))
  }
  
  print(paste0("Run: ", z))
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov_origin.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 20000, n.burnin = 10000,
      n.thin = 20,
      jags.seed = 123
    )
  
  JAGS_600_list[[z]] <- out.jags
}

saveRDS(out.jags, here::here("simulation", "JAGS_cov_origin_600_onerun.rds"))
saveRDS(JAGS_600_list, here::here("simulation", "JAGS_cov_origin_600_list.rds"))

##### 1200 fish, ORIGIN only #####
sim_1200_hist_list <- readRDS(here::here("simulation", "sim_1200_cov_origin_hist_list.rds"))
sim_1200_dates_list <- readRDS(here::here("simulation", "sim_1200_cov_origin_dates_list.rds"))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_1200.csv")))

# Create a list to store JAGS objects
JAGS_1200_list <- list()
# Loop it
for (z in 1:length(sim_1200_hist_list)){
  # for (z in 1:1){
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
      cat_X_mat_1200[i,2] <- 0
      cat_X_mat_1200[i,3] <- 1
    }
    else if (fish_sim_cat_data_1200$natal_origin[i] == 2){
      cat_X_mat_1200[i,2] <- 1
      cat_X_mat_1200[i,3] <- 0
    }
    else {
      cat_X_mat_1200[i,2] <- -1
      cat_X_mat_1200[i,3] <- -1
    }
  }
  
  
  
  
  ###### Parameters monitored #####
  parameters <- c(
    "b0_matrix",
    "borigin1_matrix",
    "borigin2_matrix"
  )
  
  
  ##### Data #####
  data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat,
               movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
               nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
               n_notmovements = n_notmovements, possible_states = transition_matrix, cat_X_mat = cat_X_mat_1200)
  
  ###### Initial values #####
  # New version, using only b0
  inits <- function(){
    b0_matrix <- matrix(NA, nrow = 10, ncol = 9)
    borigin1_matrix <- matrix(NA, nrow = 10, ncol = 9)
    borigin2_matrix <- matrix(NA, nrow = 10, ncol = 9)
    
    for (j in 1:dim(movements)[1]){
      b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      borigin1_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
      borigin2_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
    }
    
    return(list(
      b0_matrix = b0_matrix,
      borigin1_matrix = borigin1_matrix,
      borigin2_matrix = borigin2_matrix
    ))
  }
  
  print(paste0("Run: ", z))
  
  ##### Run JAGS #####
  out.jags <- 
    jags.parallel(
      data = data,
      inits = inits,
      model.file=here::here("simulation", "sim_model_cov_origin.txt"),
      parameters.to.save = parameters,
      n.chains = 3, n.iter = 20000, n.burnin = 10000,
      n.thin = 20,
      jags.seed = 123
    )
  
  JAGS_1200_list[[z]] <- out.jags
}

saveRDS(out.jags, here::here("simulation", "JAGS_cov_origin_1200_onerun.rds"))
saveRDS(JAGS_1200_list, here::here("simulation", "JAGS_cov_origin_1200_list.rds"))
