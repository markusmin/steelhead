# 09_fit_simulation

# Here we will be fitting our model to our simulated data

# Load libraries
library(tidyverse)
library(lubridate)
library(janitor)
library(rjags)
library(jagsUI)

# Load data to fit model to
sim_data <- readRDS(here::here("simulation", "sim_600_array.rds"))
sim_data_list <- readRDS(here::here("simulation", "sim_600_list.rds"))

dates <- as.matrix(readRDS(here::here("simulation", "dates_600_matrix.rds")))
# Dates are numbers, with 2017-01-01 = 1

# Read in covariates
temp_sim <- as.matrix(read.csv(here::here("simulation", "temp_600.csv"), row.names = 1))
flow_sim <- as.matrix(read.csv(here::here("simulation", "flow_600.csv"), row.names = 1))
fish_sim_cat_data <- as.matrix(read.csv(here::here("simulation", "origin_rear_600.csv")))


# Store quantities for loop
# Store the total number of individuals
n.ind <- length(sim_data_list)

# Store the number of observations per individual
# -1 because the last state is loss, which isn't actually an observation
n.obs <- unlist(lapply(sim_data_list, ncol)) - 1

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





##### JAGS model #####
cat("
model {

#state-space likelihood
for(i in 1:n.ind){ # Loop through the detection matrices for each individual
  for(j in 1:(n.obs[i]-1)){ # Loop through each of the observations, stopping at the loss column (-1)
  
  # Vectorized state transitions
  # Index flow and temperature by date (x) and index of site (y)
  # exp(b0_matrix[states_mat[i,j],] + bflow_matrix[states_mat[i,j],] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  # btemp_matrix[states_mat[i,j],] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] + 
  # brear_matrix[states_mat[i,j],] + borigin_matrix[states_mat[i,j],])/ (1 + sum(exp(b0_matrix[states_mat[i,j],] + bflow_matrix[states_mat[i,j],] + btemp_matrix[states_mat[i,j],] + 
  # brear_matrix[states_mat[i,j],] + borigin_matrix[states_mat[i,j],])))
  # Note here that we also have a 1 - sum() term for the loss probability
  
  # y[,j+1,i] ~ dmulti(c( exp(b0_matrix[states_mat[i,j],] + bflow_matrix[states_mat[i,j],] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  # btemp_matrix[states_mat[i,j],] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] + 
  # brear_matrix[states_mat[i,j],] + borigin_matrix[states_mat[i,j],])/ (1 + sum(exp(b0_matrix[states_mat[i,j],] + bflow_matrix[states_mat[i,j],] + btemp_matrix[states_mat[i,j],] + 
  # brear_matrix[states_mat[i,j],] + borigin_matrix[states_mat[i,j],]))), 1 - sum(  exp(b0_matrix[states_mat[i,j],] + bflow_matrix[states_mat[i,j],] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  # btemp_matrix[states_mat[i,j],] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] + 
  # brear_matrix[states_mat[i,j],] + borigin_matrix[states_mat[i,j],])/ (1 + sum(exp(b0_matrix[states_mat[i,j],] + bflow_matrix[states_mat[i,j],] + btemp_matrix[states_mat[i,j],] + 
  # brear_matrix[states_mat[i,j],] + borigin_matrix[states_mat[i,j],]))))) ,1)
  
  # JAGS won't let you exponentiate a vector (BOO!!) so we'll have to write out each possible transition separately
    ~ dmulti(c(
  # State 1
  exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 2
  exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 3
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 4
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 5
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 6
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 7
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 8
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 9
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
  # Loss
  (1 - sum(  # State 1
  exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 2
  exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 3
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 4
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 5
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 6
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 7
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 8
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])),
  
    # State 9
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9])/(1 +   exp(b0_matrix[states_mat[i,j], 1] + bflow_matrix[states_mat[i,j], 1] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 1] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 1] + borigin_matrix[states_mat[i,j], 1]) + exp(b0_matrix[states_mat[i,j], 2] + bflow_matrix[states_mat[i,j], 2] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 2] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 2] + borigin_matrix[states_mat[i,j], 2])+ 
  exp(b0_matrix[states_mat[i,j], 3] + bflow_matrix[states_mat[i,j], 3] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 3] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 3] + borigin_matrix[states_mat[i,j], 3])+
  exp(b0_matrix[states_mat[i,j], 4] + bflow_matrix[states_mat[i,j], 4] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 4] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 4] + borigin_matrix[states_mat[i,j], 4])+
  exp(b0_matrix[states_mat[i,j], 5] + bflow_matrix[states_mat[i,j], 5] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 5] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 5] + borigin_matrix[states_mat[i,j], 5])+
  exp(b0_matrix[states_mat[i,j], 6] + bflow_matrix[states_mat[i,j], 6] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 6] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 6] + borigin_matrix[states_mat[i,j], 6])+
  exp(b0_matrix[states_mat[i,j], 7] + bflow_matrix[states_mat[i,j], 7] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 7] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 7] + borigin_matrix[states_mat[i,j], 7])+
  exp(b0_matrix[states_mat[i,j], 8] + bflow_matrix[states_mat[i,j], 8] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 8] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 8] + borigin_matrix[states_mat[i,j], 8])+
  exp(b0_matrix[states_mat[i,j], 9] + bflow_matrix[states_mat[i,j], 9] * flow_sim[dates[i,j], flow_index[states_mat[i,j]]] + 
  btemp_matrix[states_mat[i,j], 9] * temp_sim[dates[i,j], temp_index[states_mat[i,j]]] +
  brear_matrix[states_mat[i,j], 9] + borigin_matrix[states_mat[i,j], 9]))))
  ))
  }

}

    ##### PRIORS #####
    
    ### Set the priors by matrix
        for (i in 1:nmovements){
    b0_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)
        }
    
        for (i in 1:nmovements){
    bflow_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)
        }
    
        for (i in 1:nmovements){
    btemp_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)
        }
    
        for (i in 1:nmovements){
    brear_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)
        }
    
        for (i in 1:nmovements){
    borigin_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)
        }
    
    
    
    ### Set every other element to zero
            for (i in 1:n_notmovements){
    b0_matrix[not_movements[i,1], not_movements[i,2]] <- 0
        }
    
        for (i in 1:n_notmovements){
    bflow_matrix[not_movements[i,1], not_movements[i,2]] <- 0
        }
    
        for (i in 1:n_notmovements){
    btemp_matrix[not_movements[i,1], not_movements[i,2]] <- 0
        }
    
        for (i in 1:n_notmovements){
    brear_matrix[not_movements[i,1], not_movements[i,2]] <- 0
        }
    
        for (i in 1:n_notmovements){
    borigin_matrix[not_movements[i,1], not_movements[i,2]] <- 0
        }

}", fill=TRUE, file=here::here("simulation", "sim_model.txt"))


###### Parameters monitored #####
# Here we have a matrix for each term (b0, bflow, btemp, etc.)
# Where each matrix has the term for each from (row) to (column) movement
parameters <- c(
  "b0_matrix",
  "bflow_matrix",
  "btemp_matrix",
  "brear_matrix",
  "borigin_matrix"
)


##### Data #####

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

start_matrix <- matrix(rep(0, 1), nrow = nstates, ncol = 1)
rownames(start_matrix) <- sim_states
# Start the individual in BON to MCN
start_matrix["mainstem, BON to MCN", 1] <- 1


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


# Copy this matrix for each term in the model;
model_terms <- c("b0", "bflow", "btemp", "brear", "borigin")
for (i in 1:length(parameters)){
  assign(paste0(model_terms[i],"_matrix"), transition_matrix[,1:9])
}

# Get the indices that are 1s (except loss, since that's 1 minus the others)
movements <- which(transition_matrix[,1:9] == 1, arr.ind = TRUE)
nmovements <- dim(movements)[1]
# Now get all of the movements which are fixed to zero
not_movements <- which(transition_matrix[,1:9] == 0, arr.ind = TRUE)
n_notmovements <- dim(not_movements)[1]

data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
             states_mat = states_mat, origin = fish_sim_cat_data[,2], rear = fish_sim_cat_data[,3], 
             movements = movements, not_movements = not_movements, temp_sim = temp_sim, flow_sim = flow_sim,
             nmovements = nmovements, dates = dates, flow_index = flow_index, temp_index = temp_index,
             n_notmovements = n_notmovements)

###### Initial values #####
# Loop through all of the non-zero values in each matrix
# inits <- function(){
#   for (i in 1:length(model_terms)){
#     mat <- eval(parse(text = paste0(model_terms[i], "_matrix")))
#     for (j in 1:dim(movements)[1]){
#       mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
#     }
#     assign(paste0(model_terms[i],"_matrix"), mat)
#   }
# }

inits <- function(){
  mat <- matrix(0, nrow = 10, ncol = 9)
  for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  b0_matrix = mat
  
  mat <- matrix(0, nrow = 10, ncol = 9)
  for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  bflow_matrix = mat
  
  mat <- matrix(0, nrow = 10, ncol = 9)
  for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  btemp_matrix = mat
  
  mat <- matrix(0, nrow = 10, ncol = 9)
  for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  brear_matrix = mat
  
  mat <- matrix(0, nrow = 10, ncol = 9)
  for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  borigin_matrix = mat
}



##### RUN JAGS #####
out.jags = jags(data, inits, parameters, model.file=here::here("simulation", "sim_model.txt"),
                n.chains=3, n.iter=20000, n.burnin=5000, n.thin=1)

# out.jags$summary
