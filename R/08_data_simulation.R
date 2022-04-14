# 08: Data simulation

# Here we are generating a simulated dataset, in order to test our model

# Load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)

# We will be using a reduced number of states:

# 5 mainstem sites:
# 1: Mouth to BON
# 2: BON to MCN
# 3: MCN to ICH or PRA
# 4: PRA to RIS
# 5: ICH to LGR

# 4 tributary sites:
# 1: John Day River
# 2: Deschutes River
# 3: Yakima River
# 4: Tucannon River

# We will have 3 natal origins:
# John Day River, Yakima River, Tucannon River


# We will have two continuous covariates:
# Temperature, flow

# We will have two categorical covariates:
# Natal origin, rear type
# 3000 fish total: 1000 from each natal origin, half wild and half hatchery

# 3 years: 2017-2019

# For this simulation, we will assume that every movement takes place one month after the previous one


##### Simulate covariate data #####


# Simulate temperature data
# Get dates from 2011-2015

# Simulate as sin wave
x <- seq(-pi/2,(pi*6 - (pi/2)),length.out=365*3)
y <- sin(x) +1
plot(x, y)
dates <- seq(ymd("2017-01-01"), ymd("2019-12-31"), by = "days")

# Offset temperature by 1 degree per dam upstream
# BON = 5-25
# MCN = 4-24
# PRA = 3-23
# ICH = 3-23

BON_temp_sim <- y * 10 + 5
MCN_temp_sim <- y * 10 + 4
PRA_temp_sim <- y * 10 + 3
ICH_temp_sim <- y * 10 + 3

# Create temperature df
temp_sim_df <- data.frame(date = dates, 
                          BON = BON_temp_sim, 
                          MCN = MCN_temp_sim, 
                          PRA = PRA_temp_sim, 
                          ICH = ICH_temp_sim)

# Simulate flow data

# Simulate as sin wave
x <- seq(-pi/2,(pi*6 - (pi/2)),length.out=365*3)
y <- sin(x) + 1
dates <- seq(ymd("2017-01-01"), ymd("2019-12-31"), by = "days")

# Ranges per dam
# BON = 10-300
# MCN = 10-450
# PRA = 10-400
# ICH = 10-200

BON_flow_sim <- y * 145 + 10
MCN_flow_sim <- y * 220 + 10
PRA_flow_sim <- y * 195 + 10
ICH_flow_sim <- y * 95 + 10

# Create flow df
flow_sim_df <- data.frame(date = dates, 
                          BON = BON_flow_sim, 
                          MCN = MCN_flow_sim, 
                          PRA = PRA_flow_sim, 
                          ICH = ICH_flow_sim)

# Transform flow and temperature data by z-scoring
flow_sim_df %>% 
  mutate(BON = (BON - mean(BON))/sd(BON)) %>% 
  mutate(MCN = (MCN - mean(MCN))/sd(MCN)) %>% 
  mutate(PRA = (PRA - mean(PRA))/sd(PRA)) %>% 
  mutate(ICH = (ICH - mean(ICH))/sd(ICH)) -> flow_sim_zscore_df

temp_sim_df %>% 
  mutate(BON = (BON - mean(BON))/sd(BON)) %>% 
  mutate(MCN = (MCN - mean(MCN))/sd(MCN)) %>% 
  mutate(PRA = (PRA - mean(PRA))/sd(PRA)) %>% 
  mutate(ICH = (ICH - mean(ICH))/sd(ICH)) -> temp_sim_zscore_df




# Simulate the categorical covariates (rear type and natal origin)

fish_sim_cat_data <- data.frame(tag_code = seq(1, 3000, 1),
                            natal_origin = c(rep("JDR", 1000),
                                             rep("Yakima", 1000),
                                             rep("Tucannon", 1000)),
                            rear_type = rep(c(rep("H", 500), rep("W", 500)), 3))

# Simulate all of the starting dates

# Take a uniform distribution between July 1 and September 30
BON_arrival_dates_sim <- seq(ymd("2017-07-01"), ymd("2019-09-30"), by = "days")[round(runif(3000, min = 1, max = 92), 0)]



##### Base model: no covariates #####

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

##### State transition code

nfish <- length(fish_sim_cat_data$tag_code)
# Create an array, with dimensions nfish x nstates x noccurrences
# Set noccurrences to a number higher than the number of possible states - use 100
# noccurrences will be udpated in the for loop as the fish visits a new state; start with 2; populate first column with BON to MCN
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

# Create a matrix to store the dates in each state; make ncol 100 to accommodate 100 possible state visits
state_date <- matrix(nrow = nfish, ncol = 100)

for (i in 1:nfish){ # for each fish
  # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
  movement_array["mainstem, BON to MCN", 1, i] <- 1
  # Populate state date matrix with date at arrival at BON
  state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
  
  # populate the rest of the state dates for this fish by adding a month
  for (k in 2:100){
    state_date[i,k] <- format(as.Date(ymd(state_date[i,k-1])+ months(1), origin = "1970-01-01"))
  }
  
  # For the 2nd state and onwards, loop through and evaluate multinomial logit
  # But stop when the fish is lost
  for (j in 2:100){
    while (sum(movement_array[i, "loss",]) == 0){
      # Calculate all movement probabilities via multinomial logit in the transition matrix
      
      
      # Index to row of transition matrix containing the correct "from" state to get movement probabilities
      p <- transition_matrix[rownames(as.data.frame(which(movement_array[,j-1,i] == 1))),]
      
      # Choose next state using rmultinom
      movement_array[i,,j] <- rmultinom(1, size = 1, prob = p)
      
    }
  }

  
  movement_array[i,,]
}




##### Model including covariates #####


###### To understand the dynamics of the multinomial logit, let's loop through some different values and plot them #####

# Start them all with equal probabilities

b0_BM_MB <- 1
bflow_BM_MB <- 0
# btemp_BM_MB <- 5
btemp_BM_MB <- 0
brear_BM_MB <- c(0,0) # c(hatchery, wild)
borigin_BM_MB <- c(0, 0, 0) # c(John Day River, Yakima River, Tucannon River)

# BON to MCN to MCN to ICH or PRA transition
b0_BM_MIP <- 1
bflow_BM_MIP <- 3
btemp_BM_MIP <- 5
# bflow_BM_MIP <- 0
# btemp_BM_MIP <- 0
brear_BM_MIP <- c(0,0)
borigin_BM_MIP <- c(0.5, 5, 5) # JDR fish have lower probability of overshooting MCN than the fish from origins upstream of MCN
borigin_BM_MIP <- c(0, 0, 0)

# BON to MCN to DES R transition
b0_BM_DES <- 1
bflow_BM_DES <- 0
btemp_BM_DES <- 0
brear_BM_DES <- c(0,0)
borigin_BM_DES <- c(1, 0.5, 0.5) # JDR has higher probability of stray to DES because it's closer
borigin_BM_DES <- c(0,0,0)

# BON to MCN to JDR transition
b0_BM_JDR <- 1
bflow_BM_JDR <- 0
btemp_BM_JDR <- 0
brear_BM_JDR <- c(0,0)
borigin_BM_JDR <- c(1, 0.1, 0.1)
borigin_BM_JDR <- c(0,0,0)

# store test covariate values
i <- 1
temp_BON <- 0.5
flow_BON <- -0.5
temp_MCN <- 0.5
flow_MCN <- -0.5
temp_DES <- 0.5
flow_DES <- -0.5
temp_JDR <- 0
flow_JDR <- 0

# Evaluate

# BON to MCN
phi_BM_MB <- exp(b0_BM_MB + btemp_BM_MB*temp_BON + bflow_BM_MB*flow_BON + brear_BM_MB[i] + borigin_BM_MB[i])
phi_BM_MIP <- exp(b0_BM_MIP + btemp_BM_MIP*temp_MCN + bflow_BM_MIP*flow_MCN + brear_BM_MIP[i] + borigin_BM_MIP[i])
phi_BM_DES <- exp(b0_BM_DES + btemp_BM_DES*temp_DES + bflow_BM_DES*flow_DES + brear_BM_DES[i] + borigin_BM_DES[i])
phi_BM_JDR <- exp(b0_BM_JDR + btemp_BM_JDR*temp_JDR + bflow_BM_JDR*flow_JDR + brear_BM_JDR[i] + borigin_BM_JDR[i])

transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- phi_BM_MB/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- phi_BM_MIP/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- phi_BM_DES/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
transition_matrix["mainstem, BON to MCN", "John Day River"] <- phi_BM_JDR/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
transition_matrix["mainstem, BON to MCN", "loss"] <- 1 - sum(transition_matrix["mainstem, BON to MCN", c("mainstem, mouth to BON", "mainstem, MCN to ICH or PRA",
                                                                                                         "Deschutes River", "John Day River")])
transition_matrix["mainstem, BON to MCN",]

# turn this evaluate into a function for one variable
### Let's start with just the intercept term for homing to the JDR
evalute_mlogit_b0_BM_JDR <- function(b0_BM_JDR){
  phi_BM_MB <- exp(b0_BM_MB + btemp_BM_MB*temp_BON + bflow_BM_MB*flow_BON + brear_BM_MB[i] + borigin_BM_MB[i])
  phi_BM_MIP <- exp(b0_BM_MIP + btemp_BM_MIP*temp_MCN + bflow_BM_MIP*flow_MCN + brear_BM_MIP[i] + borigin_BM_MIP[i])
  phi_BM_DES <- exp(b0_BM_DES + btemp_BM_DES*temp_DES + bflow_BM_DES*flow_DES + brear_BM_DES[i] + borigin_BM_DES[i])
  phi_BM_JDR <- exp(b0_BM_JDR + btemp_BM_JDR*temp_JDR + bflow_BM_JDR*flow_JDR + brear_BM_JDR[i] + borigin_BM_JDR[i])
  
  transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- phi_BM_MB/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
  transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- phi_BM_MIP/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
  transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- phi_BM_DES/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
  transition_matrix["mainstem, BON to MCN", "John Day River"] <- phi_BM_JDR/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
  transition_matrix["mainstem, BON to MCN", "loss"] <- 1 - sum(transition_matrix["mainstem, BON to MCN", c("mainstem, mouth to BON", "mainstem, MCN to ICH or PRA",
                                                                                                           "Deschutes River", "John Day River")])
  probs <- transition_matrix["mainstem, BON to MCN",]
  return(probs)
}

probs_df <- as.data.frame(matrix(nrow = 10, ncol = length(seq(0.1 ,10, 0.1))))

for (j in 1:length(seq(0.1,10, 0.1))){
  values <- seq(0.1,10, 0.1)
  print(paste0(j))
  probs <- evalute_mlogit_b0_BM_JDR(b0_BM_JDR = values[j])
  probs_df[,j] <- probs
}

rownames(probs_df) <- sim_states

probs_df %>% 
  rownames_to_column("state") %>% 
  pivot_longer(., colnames(probs_df)) %>% 
  dplyr::rename(b0_BM_JDR = name, prob = value) %>% 
  subset(state %in% c("mainstem, mouth to BON", "mainstem, MCN to ICH or PRA",
                      "Deschutes River", "John Day River", "loss"))-> probs_df_long

probs_df_long$b0_BM_JDR <- rep(seq(0.1, 10, 0.1), 5)

ggplot(probs_df_long, aes(x = b0_BM_JDR, y = prob, color = state)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 10")


### What if you did two intercepts at the same time?
# intercepts for going to JDR and DES
evalute_mlogit_b0_BM_JDR_b0_BM_DES <- function(b0_BM_JDR, b0_BM_DES){
  phi_BM_MB <- exp(b0_BM_MB + btemp_BM_MB*temp_BON + bflow_BM_MB*flow_BON + brear_BM_MB[i] + borigin_BM_MB[i])
  phi_BM_MIP <- exp(b0_BM_MIP + btemp_BM_MIP*temp_MCN + bflow_BM_MIP*flow_MCN + brear_BM_MIP[i] + borigin_BM_MIP[i])
  phi_BM_DES <- exp(b0_BM_DES + btemp_BM_DES*temp_DES + bflow_BM_DES*flow_DES + brear_BM_DES[i] + borigin_BM_DES[i])
  phi_BM_JDR <- exp(b0_BM_JDR + btemp_BM_JDR*temp_JDR + bflow_BM_JDR*flow_JDR + brear_BM_JDR[i] + borigin_BM_JDR[i])
  
  transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- phi_BM_MB/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
  transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- phi_BM_MIP/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
  transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- phi_BM_DES/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
  transition_matrix["mainstem, BON to MCN", "John Day River"] <- phi_BM_JDR/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
  transition_matrix["mainstem, BON to MCN", "loss"] <- 1 - sum(transition_matrix["mainstem, BON to MCN", c("mainstem, mouth to BON", "mainstem, MCN to ICH or PRA",
                                                                                                           "Deschutes River", "John Day River")])
  probs <- transition_matrix["mainstem, BON to MCN",]
  return(probs)
}

probs_df <- as.data.frame(matrix(nrow = 10, ncol = length(seq(0.1 ,10, 0.1))))

for (j in 1:length(seq(0.1,10, 0.1))){
  values <- seq(0.1,10, 0.1)
  print(paste0(j))
  probs <- evalute_mlogit_b0_BM_JDR_b0_BM_DES(b0_BM_JDR = values[j], b0_BM_DES = values[j])
  probs_df[,j] <- probs
}

rownames(probs_df) <- sim_states

probs_df %>% 
  rownames_to_column("state") %>% 
  pivot_longer(., colnames(probs_df)) %>% 
  dplyr::rename(b0_BM_JDR = name, prob = value) %>% 
  subset(state %in% c("mainstem, mouth to BON", "mainstem, MCN to ICH or PRA",
                      "Deschutes River", "John Day River", "loss"))-> probs_df_long

probs_df_long$b0_BM_JDR <- rep(seq(0.1, 10, 0.1), 5)

ggplot(probs_df_long, aes(x = b0_BM_JDR, y = prob, color = state)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 10")

##### ACTUAL SIMULATION #####


# Store simulation parameter values

# BON to MCN to Mouth to BON transition
b0_BM_MB <- 1
bflow_BM_MB <- 0.5 # Make this vary with flow, so higher flow at BON = higher fallback over BON
btemp_BM_MB <- 0 # No relationship with temperature
brear_BM_MB <- c(0,0) # c(hatchery, wild) # No relationship with rear type
borigin_BM_MB <- c(0, 0, 0) # c(John Day River, Yakima River, Tucannon River) # No relationship with origin

# BON to MCN to MCN to ICH or PRA transition
b0_BM_MIP <- 1
bflow_BM_MIP <- 0 # 
btemp_BM_MIP <- 0.5 # When it's hotter, more likely to go upstream
brear_BM_MIP <- c(0,0) # No relationship with rear type
borigin_BM_MIP <- c(0.5, 2, 2) # JDR fish have lower probability of overshooting MCN than the fish from origins upstream of MCN

# BON to MCN to DES R transition
b0_BM_DES <- 1
bflow_BM_DES <- 0 # flow not relevant for tributary entry
btemp_BM_DES <- 0 # Temperature not relevant for tributary entry
brear_BM_DES <- c(0.1,0) # Slightly higher probability of straying to Deschutes if hatchery
borigin_BM_DES <- c(0.25, 0, 0) # JDR has higher probability of stray to DES because it's closer

# BON to MCN to JDR transition
b0_BM_JDR <- 1
bflow_BM_JDR <- 0 # flow not relevant for tribuary entry
btemp_BM_JDR <- 0 # Temperature not relevant for tribuary entry
brear_BM_JDR <- c(0,0) # no effect of rear type
borigin_BM_JDR <- c(2, 0.1, 0.1) # JDR fish have much higher chance of homing

# store test covariate values
i <- 1
temp_BON <- 0.5
flow_BON <- -0.5
temp_MCN <- 0.5
flow_MCN <- -0.5
temp_DES <- 0.5
flow_DES <- -0.5
temp_JDR <- 0
flow_JDR <- 0
# Evaluate

# BON to MCN
phi_BM_MB <- exp(b0_BM_MB + btemp_BM_MB*temp_BON + bflow_BM_MB*flow_BON + brear_BM_MB[i] + borigin_BM_MB[i])
phi_BM_MIP <- exp(b0_BM_MIP + btemp_BM_MIP*temp_MCN + bflow_BM_MIP*flow_MCN + brear_BM_MIP[i] + borigin_BM_MIP[i])
phi_BM_DES <- exp(b0_BM_DES + btemp_BM_DES*temp_DES + bflow_BM_DES*flow_DES + brear_BM_DES[i] + borigin_BM_DES[i])
phi_BM_JDR <- exp(b0_BM_JDR + btemp_BM_JDR*temp_JDR + bflow_BM_JDR*flow_JDR + brear_BM_JDR[i] + borigin_BM_JDR[i])

transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- phi_BM_MB/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- phi_BM_MIP/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- phi_BM_DES/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
transition_matrix["mainstem, BON to MCN", "John Day River"] <- phi_BM_JDR/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
transition_matrix["mainstem, BON to MCN", "loss"] <- 1 - sum(transition_matrix["mainstem, BON to MCN", c("mainstem, mouth to BON", "mainstem, MCN to ICH or PRA",
                                                                                                     "Deschutes River", "John Day River")])
transition_matrix["mainstem, BON to MCN",]

#####

# MCN to ICH or PRA
# transition_matrix["mainstem, MCN to ICH or PRA", "loss"] <- 1
# transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, BON to MCN"] <- 1
# transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS"] <- 1
# transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR"] <- 1
# transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River"] <- 1

# PRA to RIS
# transition_matrix["mainstem, PRA to RIS", "loss"] <- 1
# transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- 1

# ICH to LGR
# transition_matrix["mainstem, ICH to LGR", "loss"] <- 1
# transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- 1

# Deschutes River
# transition_matrix["Deschutes River", "loss"] <- 1
# transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- 1

# John Day River
# transition_matrix["John Day River", "loss"] <- 1
# transition_matrix["John Day River", "mainstem, BON to MCN"] <- 1

# Yakima River
# transition_matrix["Yakima River", "loss"] <- 1
# transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- 1

# Tucannon River
# transition_matrix["Tucannon River", "loss"] <- 1
# transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- 1


# All other relationships will not be a product of covariates














