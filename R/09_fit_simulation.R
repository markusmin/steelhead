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
n.obs <- unlist(lapply(sim_data_list, ncol))

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
  for (j in 1:(n.obs[i]-1)){
    # states_list[[i]][j] <- rownames(as.data.frame(which(sim_data[[i]][,j] == 1)))
    states_list[[i]][j] <- which(sim_data[,j,i] == 1) # Get the index of the site instead of the name
  }
}

# Turn into matrix for JAGS
states_mat <- matrix(nrow = n.ind, ncol = max(n.obs))
for (i in 1:n.ind){
  states_mat[i,1:(n.obs[i]-1)] <- states_list[[i]]
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
    cur_state <- states_mat[i,j]

    # Get number of possible movements
    n_movements <- possible_movements[cur_state]

    # index the covariate data
    temps <- as.numeric(temp_sim[dates[i,j],])
    flows <- as.numeric(flow_sim[dates[i,j],])

    # Get the temperature values for this state
    # temp[] <- vector(length = (n_movements - 1))
    for (m in 1:(n_movements-1)){
      temp[m] <- temps[cov_index_mat[cur_state,m]]
    }


    # Get the flow values for this state
    # flow <- vector(length = (n_movements - 1))
    for (m in 1:(n_movements-1)){
      flow[m] <- flows[temps[cov_index_mat[cur_state,m]]]
    }


    # Create a design matrix for number of possible movements
    # X <- matrix(0, nrow = (n_movements - 1), ncol = (n_movements - 1) * 5)

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

    # Flow
    for (m in 1:(n_movements-1)){
        X[m, m + (n_movements-1)*2] <- flow[m]
    }

    # Multiply design matrix by beta vector to get movement probabilities

    # Get the beta vectors - will be different for each state

# Mainstem, mouth to BON
B_vec_1 <- c(b0_MB_BM, bflow_MB_BM, btemp_MB_BM, brear_MB_BM[rear], borigin_MB_BM[origin])

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

# Store beta vectors in a list
# B_vec_list <- list(B_vec_1, B_vec_2, B_vec_3, B_vec_4, B_vec_5,
#                    B_vec_6, B_vec_7, B_vec_8, B_vec_9)
# We can't use vectors in JAGS

# Store all beta vectors as one big vector
B_vec_all <- c(B_vec_1, B_vec_2, B_vec_3, B_vec_4, B_vec_5, B_vec_6, B_vec_7, B_vec_8, B_vec_9)
# Get lengths of each individual vector and store
B_vec_lengths <- c(length(B_vec_1), length(B_vec_2), length(B_vec_3), length(B_vec_4), length(B_vec_5), length(B_vec_6), length(B_vec_7), length(B_vec_8), length(B_vec_9))

# Now, index to the current vector of betas
# Basically what we're doing here is making sure we start and end with the length of the current vector
B_vec_current <- B_vec_all[(sum(B_vec_lengths[1:cur_state])-B_vec_lengths[1]+1):sum(B_vec_lengths[1:cur_state])]
  

    # p <- X %*% B_vec_list[[cur_state]]
    p <- X %*% B_vec_current

    # Evaluate the multinomial likelihood for the counts of detection probabilities
    y[,j,i] ~ dmulti(p, 1)
  }

}

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
  }
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
  btemp_MIP_YAK ~ dnorm(0,0.001)
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
  for (i in 1:2){
    brear_IL_MIP[i] ~ dnorm(0,0.001)
  }
  for (i in 1:2){
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
  for (i in 1:2){
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
  for (i in 1:2){
    brear_YAK_MIP[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_YAK_MIP[i] ~ dnorm(0,0.001)
  }

  # Tucannon River
  b0_TUC_IL ~ dnorm(0,0.001)
  bflow_TUC_IL ~ dnorm(0,0.001)
  btemp_TUC_IL ~ dnorm(0,0.001)
  for (i in 1:2){
    brear_TUC_IL[i] ~ dnorm(0,0.001)
  }
  for (i in 1:3){
    borigin_TUC_IL[i] ~ dnorm(0,0.001)
  }

}", fill=TRUE, file=here::here("simulation", "sim_model.txt"))


##### Data #####
data <- list(y = sim_data, fish_sim_cat_data = fish_sim_cat_data, temp_sim = temp_sim,
             n.ind = n.ind, n.obs = n.obs, possible_movements = possible_movements,
             flow_sim = flow_sim, states_mat = states_mat, cov_index_mat = cov_index_mat)

###### Initial values #####
inits <- function(){list(
                        # Mainstem, mouth to BON
                        b0_MB_BM = runif(1, -1, 1),
                         bflow_MB_BM = runif(1, -1, 1),
                         btemp_MB_BM = runif(1, -1, 1),
                           brear_MB_BM = runif(2, -1, 1),
                           borigin_MB_BM = runif(3, -1, 1),

                         # Mainstem, BON to MCN
                         b0_BM_MB = runif(1, -1, 1),
                         b0_BM_MIP = runif(1, -1, 1),
                         b0_BM_DES = runif(1, -1, 1),
                         b0_BM_JDR = runif(1, -1, 1),
                         bflow_BM_MB = runif(1, -1, 1),
                         bflow_BM_MIP = runif(1, -1, 1),
                         bflow_BM_DES = runif(1, -1, 1),
                         bflow_BM_JDR = runif(1, -1, 1),
                         btemp_BM_MB = runif(1, -1, 1),
                         btemp_BM_MIP = runif(1, -1, 1),
                         btemp_BM_DES = runif(1, -1, 1),
                         btemp_BM_JDR = runif(1, -1, 1),
                           brear_BM_MB = runif(2, -1, 1),
                               brear_BM_MIP = runif(2, -1, 1),
                               brear_BM_DES = runif(2, -1, 1),
                               brear_BM_JDR = runif(2, -1, 1),
                               borigin_BM_MB = runif(3, -1, 1),
                               borigin_BM_MIP = runif(3, -1, 1),
                               borigin_BM_DES = runif(3, -1, 1),
                               borigin_BM_JDR = runif(3, -1, 1),

                             # Mainstem, MCN to ICH or PRA
                             b0_MIP_BM = runif(1, -1, 1),
                             b0_MIP_PR = runif(1, -1, 1),
                             b0_MIP_IL = runif(1, -1, 1),
                             b0_MIP_YAK = runif(1, -1, 1),
                             bflow_MIP_BM = runif(1, -1, 1),
                             bflow_MIP_PR = runif(1, -1, 1),
                             bflow_MIP_IL = runif(1, -1, 1),
                             bflow_MIP_YAK = runif(1, -1, 1),
                             btemp_MIP_BM = runif(1, -1, 1),
                             btemp_MIP_PR = runif(1, -1, 1),
                             btemp_MIP_IL = runif(1, -1, 1),
                             btemp_MIP_YAK = runif(1, -1, 1),
                               brear_MIP_BM = runif(2, -1, 1),
                               brear_MIP_PR = runif(2, -1, 1),
                               brear_MIP_IL = runif(2, -1, 1),
                               brear_MIP_YAK = runif(2, -1, 1),
                               borigin_MIP_BM = runif(3, -1, 1),
                               borigin_MIP_PR = runif(3, -1, 1),
                               borigin_MIP_IL = runif(3, -1, 1),
                               borigin_MIP_YAK = runif(3, -1, 1),

                             # Mainstem, PRA to RIS
                             b0_PR_MIP = runif(1, -1, 1),
                             bflow_PR_MIP = runif(1, -1, 1),
                             btemp_PR_MIP = runif(1, -1, 1),
                               brear_PR_MIP = runif(2, -1, 1),
                               borigin_PR_MIP = runif(3, -1, 1),

                             # Mainstem, ICH to LGR
                             b0_IL_MIP = runif(1, -1, 1),
                             b0_IL_TUC = runif(1, -1, 1),
                             bflow_IL_MIP = runif(1, -1, 1),
                             bflow_IL_TUC = runif(1, -1, 1),
                             btemp_IL_MIP = runif(1, -1, 1),
                             btemp_IL_TUC = runif(1, -1, 1),
                               brear_IL_MIP = runif(2, -1, 1),
                               brear_IL_TUC = runif(2, -1, 1),
                               borigin_IL_MIP = runif(3, -1, 1),
                               borigin_IL_TUC = runif(3, -1, 1),

                             # Deschutes River
                             b0_DES_BM = runif(1, -1, 1),
                             bflow_DES_BM = runif(1, -1, 1),
                             btemp_DES_BM = runif(1, -1, 1),
                               brear_DES_BM = runif(2, -1, 1),
                               borigin_DES_BM = runif(3, -1, 1),

                             # John Day River
                             b0_JDR_BM = runif(1, -1, 1),
                             bflow_JDR_BM = runif(1, -1, 1),
                             btemp_JDR_BM = runif(1, -1, 1),
                               brear_JDR_BM = runif(2, -1, 1),
                               borigin_JDR_BM = runif(3, -1, 1),

                             # Yakima River
                             b0_YAK_MIP = runif(1, -1, 1),
                             bflow_YAK_MIP = runif(1, -1, 1),
                             btemp_YAK_MIP = runif(1, -1, 1),
                               brear_YAK_MIP = runif(2, -1, 1),
                               borigin_YAK_MIP = runif(3, -1, 1),

                             # Tucannon River
                             b0_TUC_IL = runif(1, -1, 1),
                             bflow_TUC_IL = runif(1, -1, 1),
                             btemp_TUC_IL = runif(1, -1, 1),
                              brear_TUC_IL = runif(2, -1, 1),
                              borigin_TUC_IL = runif(3, -1, 1)
                             )}

###### Parameters monitored #####
parameters <- c(
  # Mouth to BON
  "b0_MB_BM",
  "bflow_MB_BM",
  "btemp_MB_BM",
  "brear_MB_BM",
  "borigin_MB_BM",

  # Mainstem, BON to MCN
  "b0_BM_MB",
  "b0_BM_MIP",
  "b0_BM_DES",
  "b0_BM_JDR",
  "bflow_BM_MB",
  "bflow_BM_MIP",
  "bflow_BM_DES",
  "bflow_BM_JDR",
  "btemp_BM_MB",
  "btemp_BM_MIP",
  "btemp_BM_DES",
  "btemp_BM_JDR",
  "brear_BM_MB",
  "brear_BM_MIP",
  "brear_BM_DES",
  "brear_BM_JDR",
  "borigin_BM_MB",
  "borigin_BM_MIP",
  "borigin_BM_DES",
  "borigin_BM_JDR",

  # Mainstem, MCN to ICH or PRA
  "b0_MIP_BM",
  "b0_MIP_PR",
  "b0_MIP_IL",
  "b0_MIP_YAK",
  "bflow_MIP_BM",
  "bflow_MIP_PR",
  "bflow_MIP_IL",
  "bflow_MIP_YAK",
  "btemp_MIP_BM",
  "btemp_MIP_PR",
  "btemp_MIP_IL",
  "btemp_MIP_YAK",
  "brear_MIP_BM",
  "brear_MIP_PR",
  "brear_MIP_IL",
  "brear_MIP_YAK",
  "borigin_MIP_BM",
  "borigin_MIP_PR",
  "borigin_MIP_IL",
  "borigin_MIP_YAK",

  # Mainstem, PRA to RIS
  "b0_PR_MIP",
  "bflow_PR_MIP",
  "btemp_PR_MIP",
  "brear_PR_MIP",
  "borigin_PR_MIP",

  # Mainstem, ICH to LGR
  "b0_IL_MIP",
  "b0_IL_TUC",
  "bflow_IL_MIP",
  "bflow_IL_TUC",
  "btemp_IL_MIP",
  "btemp_IL_TUC",
  "brear_IL_MIP",
  "brear_IL_TUC",
  "borigin_IL_MIP",
  "borigin_IL_TUC",

  # Deschutes River
  "b0_DES_BM",
  "bflow_DES_BM",
  "btemp_DES_BM",
  "brear_DES_BM",
  "borigin_DES_BM",

  # John Day River
  "b0_JDR_BM",
  "bflow_JDR_BM",
  "btemp_JDR_BM",
  "brear_JDR_BM",
  "borigin_JDR_BM",

  # Yakima River
  "b0_YAK_MIP",
  "bflow_YAK_MIP",
  "btemp_YAK_MIP",
  "brear_YAK_MIP",
  "borigin_YAK_MIP",

  # Tucannon River
  "b0_TUC_IL",
  "bflow_TUC_IL",
  "btemp_TUC_IL",
  "brear_TUC_IL",
  "borigin_TUC_IL"
)

##### RUN JAGS #####
out.jags = jags(data, inits, parameters, model.file=here::here("simulation", "sim_model.txt"),
                n.chains=3, n.iter=20000, n.burnin=5000, n.thin=1)

# out.jags$summary