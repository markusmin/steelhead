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

# Create rear and origin vectors, which we can index to


# Get the beta vectors - will be different for each state
# NOTE: every single one of these will require priors

# Mainstem, mouth to BON
B_vec_1 <- c(b0_MB_BM, bflow_MB_BM, btemp_MB_BM, brear_MB_BM[rear], borigin_MB_BM[origin])
temp_vec_1 <- c(temp_BON)
flow_vec_1 <- c(temp_BON)

# Mainstem, BON to MCN
B_vec_2 <- c(b0_BM_MB, b0_BM_MIP,  b0_BM_DES,  b0_BM_JDR,
             bflow_BM_MB, bflow_BM_MIP,  bflow_BM_DES,  bflow_BM_JDR,
             btemp_BM_MB, btemp_BM_MIP,  btemp_BM_DES,  btemp_BM_JDR,
             brear_BM_MB[rear], brear_BM_MIP[rear],  brear_BM_DES[rear],  brear_BM_JDR[rear],
             borigin_BM_MB[origin], borigin_BM_MIP[origin],  borigin_BM_DES[origin], borigin_BM_JDR[origin])

# Mainstem, MCN to ICH or PRA
B_vec_3 <- c(b0_MIP_BM, b0_MIP_PR, b0_MIP_IL, b0_MIP_YAK,
             bflow_MIP_BM, bflow_MIP_PR, bflow_MIP_IL, bflow_MIP_YAK,
             btemp_MIP_BM, btemp_MIP_PR, btemp_MIP_IL, btemp_MIP_YAK,
             brear_MIP_BM[rear], brear_MIP_PR[rear], brear_MIP_IL[rear], brear_MIP_YAK[rear],
             borigin_MIP_BM[origin], borigin_MIP_PR[origin], borigin_MIP_IL[origin], borigin_MIP_YAK[origin])

# Mainstem, PRA to RIS
B_vec_4 <- c(b0_PR_MIP, bflow_PR_MIP, btemp_PR_MIP, brear_PR_MIP[rear], borigin_PR_MIP[origin])

# Mainstem, ICH to LGR
B_vec_5 <- c(b0_IL_MIP, b0_IL_TUC, 
             bflow_IL_MIP, bflow_IL_TUC, 
             btemp_IL_MIP, btemp_IL_TUC, 
             brear_IL_MIP[rear], brear_IL_TUC[rear], 
             borigin_IL_MIP[origin], borigin_IL_TUC[origin])

# Deschutes River
B_vec_6 <- c(b0_DES_BM, bflow_DES_BM, btemp_DES_BM, brear_DES_BM[rear], borigin_DES_BM[origin])

# John Day River
B_vec_7 <- c(b0_JDR_BM, bflow_JDR_BM, btemp_JDR_BM, brear_JDR_BM[rear], borigin_JDR_BM[origin])

# Yakima River
B_vec_8 <- c(b0_YAK_MIP, bflow_YAK_MIP, btemp_YAK_MIP, brear_YAK_MIP[rear], borigin_YAK_MIP[origin])

# Tucannon River
B_vec_9 <- c(b0_TUC_IL, bflow_TUC_IL, btemp_TUC_IL, brear_TUC_IL[rear], borigin_TUC_IL[origin])

B_vec_list <- list(B_vec_1, B_vec_2, B_vec_3, B_vec_4, B_vec_5, 
                   B_vec_6, B_vec_7, B_vec_8, B_vec_9)

# Testing vectors
B_vec_2 <- seq(1,20,1)

B_vec_1 <- seq(1, 5, 1)

# Create a design matrix for number of possible movements
X <- matrix(0, nrow = (n_movements - 1), ncol = (n_movements - 1) * 5)

# Populate the matrix
for (m in 1:(n_movements-1)){
  for (k in 1:5){
    X[m, (k-1)*(n_movements-1) + m] <- 1
  }
}

X %*% B_vec_1

X %*% B_vec_2

# Fit the model in JAGS

cat("
model {

#state-space likelihood 
for(i in 1:n.ind){ # Loop through the detection matrices for each individual
  for(j in 1:(n.obs[i]-1)){ # Loop through each of the observations, stopping at the loss column (-1)
  
    # Get the rear type
    rear <- fish_sim_cat_data[i,3]
    
    # Get the origin
    origin <- fish_sim_cat_data[i,2]
    
    
    # Get the current state
    cur_state <- states_list[[i]][j]
    
    # Get number of possible movements
    n_movements <- possible_movements[cur_state]
    
    # Get the temperature values for this state
    temp <- vector(length = (n_movements - 1))
    for (m in 1:(n_movements-1)){
      temp[m] <- temp_mat[cur_state, m]
    }
    
    
    # Get the flow values for this state
    flow <- vector(length = (n_movements - 1))
    for (m in 1:(n_movements-1)){
      flow[m] <- flow_mat[cur_state, m]
    }
    
    
    # Create a design matrix for number of possible movements
    X <- matrix(0, nrow = (n_movements - 1), ncol = (n_movements - 1) * 5)
    
    # Populate the matrix
    for (m in 1:(n_movements-1)){
      for (k in 1:5){
        X[m, (k-1)*(n_movements-1) + m] <- 1
      }
      
    }
    
    # Populate the temperature and flow elements of the matrix with the covariate values
    
    # Temperature
    for (m in 1:(n_movements-1)){
        X[m, m + (n_movements-1)] <- temp[m]
    }
      
      for (m in 1:(n_movements-1)){
        X[m, m + (n_movements-1)] <- flow[m]
      }
      
      
    
    # Multiply design matrix by beta vector to get movement probabilities
    
    p <- X %*% B_vec_list[[cur_state]]
    
    # Evaluate the multinomial likelihood for the counts of detection probabilities
    y[[i]][,j] ~ dmulti(p, 1)
    
    ##### PRIORS #####
    
    # Mainstem, mouth to BON
    b0_MB_BM ~ dnorm(0,0.001)
    bflow_MB_BM ~ dnorm(0,0.001)
    btemp_MB_BM ~ dnorm(0,0.001)
    for (i in 1:2){
      brear_MB_BM[i] ~ dnorm(0,0.001)
    }
    for (i in 1:3){
      borigin_MB_BM[i] ~ dnorm(0,0.001)
    }

  # Mainstem, BON to MCN
  b0_BM_MB ~ dnorm(0,0.001)
  b0_BM_MIP ~ dnorm(0,0.001)
  b0_BM_DES ~ dnorm(0,0.001)
  b0_BM_JDR ~ dnorm(0,0.001)
  bflow_BM_MB ~ dnorm(0,0.001)
  bflow_BM_MIP ~ dnorm(0,0.001)
  bflow_BM_DES ~ dnorm(0,0.001)
  bflow_BM_JDR ~ dnorm(0,0.001)
  btemp_BM_MB ~ dnorm(0,0.001)
  btemp_BM_MIP ~ dnorm(0,0.001)
  btemp_BM_DES ~ dnorm(0,0.001)
  btemp_BM_JDR ~ dnorm(0,0.001)
  for (i in 1:2){
    brear_BM_MB[i] ~ dnorm(0,0.001)
  {
  for (i in 1:2){
    brear_BM_MIP[i] ~ dnorm(0,0.001)
  }
  for (i in 1:2){
    brear_BM_DES[i] ~ dnorm(0,0.001)
  }
  for (i in 1:2){
    brear_BM_JDR[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_BM_MB[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_BM_MIP[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_BM_DES[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_BM_JDR[i] ~ dnorm(0,0.001)
  }

  # Mainstem, MCN to ICH or PRA
  b0_MIP_BM ~ dnorm(0,0.001)
  b0_MIP_PR ~ dnorm(0,0.001)
  b0_MIP_IL ~ dnorm(0,0.001)
  b0_MIP_YAK ~ dnorm(0,0.001)
  bflow_MIP_BM ~ dnorm(0,0.001)
  bflow_MIP_PR ~ dnorm(0,0.001)
  bflow_MIP_IL ~ dnorm(0,0.001)
  bflow_MIP_YAK ~ dnorm(0,0.001)
  btemp_MIP_BM ~ dnorm(0,0.001)
  btemp_MIP_PR ~ dnorm(0,0.001)
  btemp_MIP_IL ~ dnorm(0,0.001)
  btemp_MIP_YAK
  for (i in 1:2){
    brear_MIP_BM[i] ~ dnorm(0,0.001)
  }
  for (i in 1:2){
    brear_MIP_PR[i] ~ dnorm(0,0.001)
  }
  for (i in 1:2){
    brear_MIP_IL[i] ~ dnorm(0,0.001)
  }
  for (i in 1:2){
    brear_MIP_YAK[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_MIP_BM[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_MIP_PR[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_MIP_IL[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_MIP_YAK[i] ~ dnorm(0,0.001)
  }

  # Mainstem, PRA to RIS
  b0_PR_MIP ~ dnorm(0,0.001)
  bflow_PR_MIP ~ dnorm(0,0.001)
  btemp_PR_MIP ~ dnorm(0,0.001)
  for (i in 1:2){
    brear_PR_MIP[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_PR_MIP[i] ~ dnorm(0,0.001)
  }

  # Mainstem, ICH to LGR
  b0_IL_MIP ~ dnorm(0,0.001)
  b0_IL_TUC ~ dnorm(0,0.001)
  bflow_IL_MIP ~ dnorm(0,0.001)
  bflow_IL_TUC ~ dnorm(0,0.001)
  btemp_IL_MIP ~ dnorm(0,0.001)
  btemp_IL_TUC ~ dnorm(0,0.001)
  for (i in 1:3){
    brear_IL_MIP[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    brear_IL_TUC[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_IL_MIP[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_IL_TUC[i] ~ dnorm(0,0.001)
  }

  # Deschutes River
  b0_DES_BM ~ dnorm(0,0.001)
  bflow_DES_BM ~ dnorm(0,0.001)
  btemp_DES_BM ~ dnorm(0,0.001)
  for (i in 1:3){
    brear_DES_BM[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_DES_BM[i] ~ dnorm(0,0.001)
  }

  # John Day River
  b0_JDR_BM ~ dnorm(0,0.001)
  bflow_JDR_BM ~ dnorm(0,0.001)
  btemp_JDR_BM ~ dnorm(0,0.001)
  for (i in 1:2){
    brear_JDR_BM[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_JDR_BM[i] ~ dnorm(0,0.001)
  }

  # Yakima River
  b0_YAK_MIP ~ dnorm(0,0.001)
  bflow_YAK_MIP ~ dnorm(0,0.001)
  btemp_YAK_MIP ~ dnorm(0,0.001)
  for (i in 1:3){
    brear_YAK_MIP[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_YAK_MIP[i] ~ dnorm(0,0.001)
  }

  # Tucannon River
  b0_TUC_IL ~ dnorm(0,0.001)
  bflow_TUC_IL ~ dnorm(0,0.001)
  btemp_TUC_IL ~ dnorm(0,0.001)
  for (i in 1:3){
    brear_TUC_IL[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_TUC_IL[i] ~ dnorm(0,0.001)
  }
  
  }

}
", fill=TRUE, file=here::here("simulation", "sim_model.txt"))


# Get data

y = sim_states