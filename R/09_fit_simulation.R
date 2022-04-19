# 09_fit_simulation

# Here we will be fitting our model to our simulated data

# Load libraries
library(tidyverse)
library(lubridate)
library(janitor)
library(rjags)
library(jagsUI)

# Load data to fit model to
sim_data <- readRDS(here::here("simulation", "sim_600.rds"))

# Read in covariates
temp_sim_zscore_df <- read.csv(here::here("simulation", "temp_600.csv"))
flow_sim_zscore_df <- read.csv(here::here("simulation", "flow_600.csv"))
fish_sim_cat_data <- read.csv(here::here("simulation", "origin_rear_600.csv"))


# Store quantities for loop
# Store the total number of individuals
n.ind <- length(sim_data)

# Store the number of observations per individual
n.obs <- unlist(lapply(sim_data, ncol))

# Get number of possible movements from each site
rownames(sim_data[[1]])
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
  for (j in 1:(n.obs[i]-1)){
    # states_list[[i]][j] <- rownames(as.data.frame(which(sim_data[[i]][,j] == 1)))
    states_list[[i]][j] <- which(sim_data[[i]][,j] == 1) # Get the index of the site instead of the name
  }
}

# Fit the model in JAGS


cat("
model {

#state-space likelihood 
for(i in 1:n.ind){ # Loop through the detection matrices for each individual
  for(j in 1:(n.obs[i]-1)){ # Loop through each of the observations, stopping at the loss column (-1)
    
    # Get the current state
    cur_state <- states_list[[i]][j]
    
    # Get number of possible movements
    n_movements <- possible_movements[cur_state]
    
    # Create a design matrix for number of possible movements
    X <- matrix(0, nrow = (n_movements - 1), ncol = (n_movements - 1) * 5)
    
    # Populate the matrix
    for (m in 1:(n_movements-1)){
      for (k in 1:5){
        X[m, (k-1)*(n_movements-1) + m] <- 1
      }
      
    }
    
    
    # Probability of ascending a dam, for an individual i at location l at time t
    # This way, each ascension probability can be generalized to a[l], rather than having
    # to write out separate ones
    
    # Ascend dam
    logit(a[i,l,t]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) 
    bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
    # Descend dam
    logit(d[i,l,t]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) 
    bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
    # Enter tributary
    # May need one line here for each tributary
    logit([i,l,t]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) 
    bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
    
    # Loss (end - use e to not confuse with location l index)
    logit(e[i,l,t]) <- b0 + bOrigin[natal_origin[i]] + # origin (categorical) 
    bYear[run_year[i]] + # run year (categorical)
    bTemp * temp[l, t] + # temperature (continuous)
    bFlow * flow[l, t] + # flow (continuous)
    bSpill * spill[l, t] # spill (continuous)
  
    
    # Put the movement probabilities into a vector
    # All other movements probabilities are zero
    p <- rep(0, nstates)
    
    # Get the index of the ascend, descend, etc. states
    p_ascend_index <- state_relationships[l,ascend]
    p_descend_index <- state_relationships[l,descend]
    
    p[p_ascend_index] <- a
    p[p_descend_index] <- d
    # put in the loss values
    p[nstates] <- e
    # Evaluate the multinomial likelihood for the counts of detection probabilities
    y[1:nstates] ~ dmulti(p[1:nstates], 1)
  }

}
", fill=TRUE, file=here::here("simulation", "sim_model.txt"))


