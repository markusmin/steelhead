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

# Create lists for covariates, to inform which covariates to choose for each transition
cov_index_list <- list()
for (i in 1:(length(possible_movements)-1)){
  vec <- vector(length = (possible_movements[i]-1))
  
  # Store in list
  cov_index_list[[i]] <- vec
}

# Mainstem, mouth to BON
cov_index_list[[1]] <- c(1) # BON
# Mainstem, BON to MCN
cov_index_list[[2]] <- c(1,2,6,5) # BON, MCN, DES, JDR
# Mainstem, MCN to ICH or PRA
cov_index_list[[3]] <- c(2,3,4,7) # MCN, PRA, ICH, YAK
# Mainstem, PRA to RIS
cov_index_list[[4]] <- c(3) # PRA
# Mainstem, ICH to LGR
cov_index_list[[5]] <- c(4, 8) # ICH, TUC
# Deschutes River
cov_index_list[[6]] <- c(6) # DES
# John Day River
cov_index_list[[7]] <- c(5) # JDR
# Yakima River
cov_index_list[[8]] <- c(7) # YAK
# Tucannon
cov_index_list[[9]] <- c(8) # TUC

# Turn this into a matrix for JAGS
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



# Get the beta vectors - will be different for each state

# # Mainstem, mouth to BON
# B_vec_1 <- c(b0_MB_BM, bflow_MB_BM, btemp_MB_BM, brear_MB_BM[rear], borigin_MB_BM[origin])
# 
# # Mainstem, BON to MCN
# B_vec_2 <- c(b0_BM_MB, b0_BM_MIP,  b0_BM_DES,  b0_BM_JDR,
#              bflow_BM_MB, bflow_BM_MIP,  bflow_BM_DES,  bflow_BM_JDR,
#              btemp_BM_MB, btemp_BM_MIP,  btemp_BM_DES,  btemp_BM_JDR,
#              brear_BM_MB[rear], brear_BM_MIP[rear],  brear_BM_DES[rear],  brear_BM_JDR[rear],
#              borigin_BM_MB[origin], borigin_BM_MIP[origin],  borigin_BM_DES[origin], borigin_BM_JDR[origin])
# 
# # Mainstem, MCN to ICH or PRA
# B_vec_3 <- c(b0_MIP_BM, b0_MIP_PR, b0_MIP_IL, b0_MIP_YAK,
#              bflow_MIP_BM, bflow_MIP_PR, bflow_MIP_IL, bflow_MIP_YAK,
#              btemp_MIP_BM, btemp_MIP_PR, btemp_MIP_IL, btemp_MIP_YAK,
#              brear_MIP_BM[rear], brear_MIP_PR[rear], brear_MIP_IL[rear], brear_MIP_YAK[rear],
#              borigin_MIP_BM[origin], borigin_MIP_PR[origin], borigin_MIP_IL[origin], borigin_MIP_YAK[origin])
# 
# # Mainstem, PRA to RIS
# B_vec_4 <- c(b0_PR_MIP, bflow_PR_MIP, btemp_PR_MIP, brear_PR_MIP[rear], borigin_PR_MIP[origin])
# 
# # Mainstem, ICH to LGR
# B_vec_5 <- c(b0_IL_MIP, b0_IL_TUC, 
#              bflow_IL_MIP, bflow_IL_TUC, 
#              btemp_IL_MIP, btemp_IL_TUC, 
#              brear_IL_MIP[rear], brear_IL_TUC[rear], 
#              borigin_IL_MIP[origin], borigin_IL_TUC[origin])
# 
# # Deschutes River
# B_vec_6 <- c(b0_DES_BM, bflow_DES_BM, btemp_DES_BM, brear_DES_BM[rear], borigin_DES_BM[origin])
# 
# # John Day River
# B_vec_7 <- c(b0_JDR_BM, bflow_JDR_BM, btemp_JDR_BM, brear_JDR_BM[rear], borigin_JDR_BM[origin])
# 
# # Yakima River
# B_vec_8 <- c(b0_YAK_MIP, bflow_YAK_MIP, btemp_YAK_MIP, brear_YAK_MIP[rear], borigin_YAK_MIP[origin])
# 
# # Tucannon River
# B_vec_9 <- c(b0_TUC_IL, bflow_TUC_IL, btemp_TUC_IL, brear_TUC_IL[rear], borigin_TUC_IL[origin])
# 
# # Store beta vectors in a list
# B_vec_list <- list(B_vec_1, B_vec_2, B_vec_3, B_vec_4, B_vec_5, 
#                    B_vec_6, B_vec_7, B_vec_8, B_vec_9)









































##### JAGS model #####
cat("
model {

#state-space likelihood
for(i in 1:n.ind){ # Loop through the detection matrices for each individual
  for(j in 1:(n.obs[i]-1)){ # Loop through each of the observations, stopping at the loss column (-1)
  
  # Vectorized state transitions
  exp(b0_matrix[states_mat[i,j],] + bflow_matrix[states_mat[i,j],] + btemp_matrix[states_mat[i,j],] + 
  brear_matrix[states_mat[i,j],] + borigin_matrix[states_mat[i,j],])/ (1 + sum(exp(b0_matrix[states_mat[i,j],] + bflow_matrix[states_mat[i,j],] + btemp_matrix[states_mat[i,j],] + 
  brear_matrix[states_mat[i,j],] + borigin_matrix[states_mat[i,j],])))
  
  y[,,i+1] <- dmulti( ,1)
  
  }

}

    ##### PRIORS #####
    
    # Set the priors by matrix
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

}", fill=TRUE, file=here::here("simulation", "sim_model.txt"))


###### Parameters monitored #####
# Here we have a matrix for each term (b0, bflow, btemp, etc.)
# Where each matrix has the term for each from (row) to (column) movement
parameters <- c(
  "b0",
  "bflow",
  "btemp",
  "brear",
  "borigin"
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

data <- list(y = sim_data,n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
             states_mat = states_mat, origin = fish_sim_cat_data[,2], rear = fish_sim_cat_data[,3], 
             movements = movements, temp_sim = temp_sim, flow_sim = flow_sim)

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
  b0_matrix =    for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  bflow_matrix =    for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  btemp_matrix =    for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  brear_matrix =    for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  borigin_matrix =    for (j in 1:dim(movements)[1]){
    mat[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
}



##### RUN JAGS #####
out.jags = jags(data, inits, parameters, model.file=here::here("simulation", "sim_model.txt"),
                n.chains=3, n.iter=20000, n.burnin=5000, n.thin=1)

# out.jags$summary
