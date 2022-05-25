# 12.1 FULL MODEL

# In this initial fit to the full dataset, we will be fitting an intercept-only model.

# This script will write the JAGS code and run the model.

# Load libraries
library(tidyverse)
library(lubridate)
library(janitor)
library(rjags)
library(jagsUI)
library(R2jags)
library(coda)


##### Get model data inputs setup #####


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
  
  # Tributary sites (17)
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
  
  # Loss
  "loss"
)

# 26 states, plus loss
nstates <- length(model_states)


##### Create the transition matrix #####
# Need to update this with the additional states

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
transition_matrix["Asotin Creek", "mainstem, ICH to LGR"] <- 1


# 23: Clearwater River
transition_matrix["Clearwater River", "loss"] <- 1
transition_matrix["Clearwater River", "mainstem, ICH to LGR"] <- 1


# 24: Salmon River
transition_matrix["Salmon River", "loss"] <- 1
transition_matrix["Salmon River", "mainstem, ICH to LGR"] <- 1


# 25: Grande Ronde River
transition_matrix["Grande Ronde River", "loss"] <- 1
transition_matrix["Grande Ronde River", "mainstem, ICH to LGR"] <- 1


# 26: Imnaha River
transition_matrix["Imnaha River", "loss"] <- 1
transition_matrix["Imnaha River", "mainstem, ICH to LGR"] <- 1




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































##### Get states mat #####

##### JAGS model - ORIGIN #####
cat("
model {

#state-space likelihood
for(i in 1:n.ind){ # Loop through the detection matrices for each individual
  for(j in 1:(n.obs[i])){ # Loop through each of the observations, stopping at the loss column (-1)
  
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
  y[,j+1,i] ~ dmulti(c(
  # State 1
  possible_states[states_mat[i,j],1] * exp(b0_matrix[states_mat[i,j], 1])/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
  # State 2
  possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 3
  possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 4
  possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 5
  possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
    
    # State 6
  possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 7
  possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 8
  possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 9
  possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
  # Loss
  (1 - sum(
  # State 1
  possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
  # State 2
  possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 3
  possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 4
  possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 5
  possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
    
    # State 6
  possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 7
  possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 8
  possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))),
  
    # State 9
  possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9]))/
  (1 +   
    possible_states[states_mat[i,j],1] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 1], borigin1_matrix[states_mat[i,j], 1], borigin2_matrix[states_mat[i,j], 1])) + 
    possible_states[states_mat[i,j],2] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 2], borigin1_matrix[states_mat[i,j], 2], borigin2_matrix[states_mat[i,j], 2])) + 
    possible_states[states_mat[i,j],3] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 3], borigin1_matrix[states_mat[i,j], 3], borigin2_matrix[states_mat[i,j], 3])) + 
    possible_states[states_mat[i,j],4] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 4], borigin1_matrix[states_mat[i,j], 4], borigin2_matrix[states_mat[i,j], 4])) + 
    possible_states[states_mat[i,j],5] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 5], borigin1_matrix[states_mat[i,j], 5], borigin2_matrix[states_mat[i,j], 5])) + 
    possible_states[states_mat[i,j],6] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 6], borigin1_matrix[states_mat[i,j], 6], borigin2_matrix[states_mat[i,j], 6])) + 
    possible_states[states_mat[i,j],7] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 7], borigin1_matrix[states_mat[i,j], 7], borigin2_matrix[states_mat[i,j], 7])) + 
    possible_states[states_mat[i,j],8] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 8], borigin1_matrix[states_mat[i,j], 8], borigin2_matrix[states_mat[i,j], 8])) + 
    possible_states[states_mat[i,j],9] * exp(cat_X_mat[i,] %*% c(b0_matrix[states_mat[i,j], 9], borigin1_matrix[states_mat[i,j], 9], borigin2_matrix[states_mat[i,j], 9])))))
  ), 1)
  }

}

    ##### PRIORS #####
    
    ### Set the priors by matrix
    
        for (i in 1:nmovements){
    b0_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)
        }
    
        for (i in 1:nmovements){
    borigin1_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)
    borigin2_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)
        }
    
    
    
    ### Set every other element to zero by using dnorm with a very high precision - this will help with initial values, I believe
    # Needs to be a very high NEGATIVE value, not 0, because exp(0) is 1
            for (i in 1:n_notmovements){
    b0_matrix[not_movements[i,1], not_movements[i,2]] <- -9999
        }
    
        for (i in 1:n_notmovements){
    borigin1_matrix[not_movements[i,1], not_movements[i,2]] <- -9999
    borigin2_matrix[not_movements[i,1], not_movements[i,2]] <- -9999
        }

}", fill=TRUE, file=here::here("simulation", "sim_model_cov_origin.txt")) 



