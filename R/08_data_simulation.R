# 08: Data simulation

# Here we are generating a simulated dataset, in order to test our model

# Load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)
library(ggthemes)

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
# Get dates from 2017-2019

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

plot(dates, BON_temp_sim)


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
  mutate(ICH = (ICH - mean(ICH))/sd(ICH)) %>% 
  # Store dummy values (all 0) for flow in tribs
  mutate(JDR = 0) %>% 
  mutate(DES = 0) %>% 
  mutate(YAK = 0) %>% 
  mutate(TUC = 0)-> flow_sim_zscore_df

# Change dates to numeric (for JAGS)
flow_sim_zscore_df %>% 
  mutate(date_numeric = as.numeric(date - ymd("2016-12-31"))) %>% 
  relocate(date_numeric) %>% 
  dplyr::select(-date) %>% 
  dplyr::rename(date = date_numeric) -> flow_sim_zscore_datenumeric

temp_sim_df %>% 
  mutate(BON = (BON - mean(BON))/sd(BON)) %>% 
  mutate(MCN = (MCN - mean(MCN))/sd(MCN)) %>% 
  mutate(PRA = (PRA - mean(PRA))/sd(PRA)) %>% 
  mutate(ICH = (ICH - mean(ICH))/sd(ICH)) %>%  
  # Store dummy values (all 0) for temp in tribs
  mutate(JDR = 0) %>% 
  mutate(DES = 0) %>% 
  mutate(YAK = 0) %>% 
  mutate(TUC = 0)-> temp_sim_zscore_df

# Change dates to numeric (for JAGS)
temp_sim_zscore_df %>% 
  mutate(date_numeric = as.numeric(date - ymd("2016-12-31"))) %>% 
  relocate(date_numeric) %>% 
  dplyr::select(-date) %>% 
  dplyr::rename(date = date_numeric) -> temp_sim_zscore_datenumeric




# Simulate the categorical covariates (rear type and natal origin)
# 1 = JDR, 2 = Yakima, 3 = Tucannon
# 1 = hatchery, 2 = wild

# Start with just 600 fish (200 from each origin, 100 each hatchery and wild)

fish_sim_cat_data <- data.frame(tag_code = seq(1, 600, 1),
                            natal_origin = c(rep(1, 200),
                                             rep(2, 200),
                                             rep(3, 200)),
                            rear_type = rep(c(rep(1, 100), rep(2, 100)), 3))

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

transition_date <- matrix(nrow = nfish, ncol = 100)

set.seed(123)
for (i in 1:nfish){ # for each fish
# for (i in 1:1){
  print(paste0("i = ", i))
  # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
  movement_array["mainstem, BON to MCN", 1, i] <- 1
  # Populate state date matrix with date at arrival at BON
  state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
  # Store this in the transition date matrix
  transition_date[i,1] <- state_date[i,1]
  
  # populate the rest of the state dates for this fish by adding two weeks
  for (k in 2:100){
    # state_date[i,k] <- format(as.Date(ymd(state_date[i,k-1])+ months(1), origin = "1970-01-01"))
    # state_date[i,k] <- format(as.Date(state_date[i,k-1], origin = "1970-01-01")) + months(1)
    # state_date[i,k] <- ymd(state_date[i,k-1]) + days(14)
    state_date[i,k] <- format(as.Date(ymd(state_date[i,k-1]) + days(14), origin = "1970-01-01"))
  }
  
  # Get the categorical data
  origin <- fish_sim_cat_data$natal_origin[i]
  rear <- fish_sim_cat_data$rear_type[i]
  
  # For the 2nd state and onwards, loop through and evaluate multinomial logit
  # But stop when the fish is lost
  j <- 2
    while (sum(movement_array["loss",,i]) == 0){
      # Calculate all movement probabilities via multinomial logit in the transition matrix
      print(paste0("j = ", j))
      # Get covariates
      # index the covariate data
      temps <- temp_sim_zscore_df[temp_sim_zscore_df$date == state_date[i,j-1],][2:ncol(temp_sim_zscore_df)]
      temp_BON <- temps[1,1]; temp_MCN <- temps[1,2]; temp_PRA <- temps[1,3]; temp_ICH <- temps[1,4]; temp_JDR <- temps[1,5]; temp_DES <- temps[1,6]; temp_YAK <- temps[1,7]; temp_TUC <- temps[1,8]
      flows <- flow_sim_zscore_df[flow_sim_zscore_df$date == state_date[i,j-1],][2:ncol(flow_sim_zscore_df)]
      flow_BON <- flows[1,1]; flow_MCN <- flows[1,2]; flow_PRA <- flows[1,3]; flow_ICH <- flows[1,4]; flow_JDR <- flows[1,5]; flow_DES <- flows[1,6]; flow_YAK <- flows[1,7]; flow_TUC <- flows[1,8]
      
      
      # The function:
      #####
      ##### FROM MOUTH TO BON STATE #####
      
      # Mouth to BON to BON to MCN transition
      b0_MB_BM <- 1
      bflow_MB_BM <- 0 # No relationship with flow
      btemp_MB_BM <- 0 # No relationship with temperature
      brear_MB_BM <- c(0,0) # No effect of rear type
      borigin_MB_BM <- c(2, 2, 2) # Highly positive for each rear type
      
      phi_MB_BM <- exp(b0_MB_BM + btemp_MB_BM*temp_BON + bflow_MB_BM*flow_BON + brear_MB_BM[rear] + borigin_MB_BM[origin])
      
      transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- phi_MB_BM/(1 + phi_MB_BM)
      transition_matrix["mainstem, mouth to BON", "loss"] <- 1 - transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"]
      
      # Check values
      transition_matrix["mainstem, mouth to BON", ]
      
      ##### FROM BON TO MCN STATE #####
      
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
      brear_BM_DES <- c(0.5,0) # higher probability of straying to Deschutes if hatchery
      borigin_BM_DES <- c(0.25, 0, 0) # JDR has higher probability of stray to DES because it's closer
      
      # BON to MCN to JDR transition
      b0_BM_JDR <- 1
      bflow_BM_JDR <- 0 # flow not relevant for tributary entry
      btemp_BM_JDR <- 0 # Temperature not relevant for tributary entry
      brear_BM_JDR <- c(0,0) # no effect of rear type
      borigin_BM_JDR <- c(2, 0.1, 0.1) # JDR fish have much higher chance of homing
      
      # Evaluate
      
      # BON to MCN
      phi_BM_MB <- exp(b0_BM_MB + btemp_BM_MB*temp_BON + bflow_BM_MB*flow_BON + brear_BM_MB[rear] + borigin_BM_MB[origin])
      phi_BM_MIP <- exp(b0_BM_MIP + btemp_BM_MIP*temp_MCN + bflow_BM_MIP*flow_MCN + brear_BM_MIP[rear] + borigin_BM_MIP[origin])
      phi_BM_DES <- exp(b0_BM_DES + btemp_BM_DES*temp_DES + bflow_BM_DES*flow_DES + brear_BM_DES[rear] + borigin_BM_DES[origin])
      phi_BM_JDR <- exp(b0_BM_JDR + btemp_BM_JDR*temp_JDR + bflow_BM_JDR*flow_JDR + brear_BM_JDR[rear] + borigin_BM_JDR[origin])
      
      transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- phi_BM_MB/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
      transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- phi_BM_MIP/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
      transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- phi_BM_DES/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
      transition_matrix["mainstem, BON to MCN", "John Day River"] <- phi_BM_JDR/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
      transition_matrix["mainstem, BON to MCN", "loss"] <- 1 - sum(transition_matrix["mainstem, BON to MCN", c("mainstem, mouth to BON", "mainstem, MCN to ICH or PRA",
                                                                                                               "Deschutes River", "John Day River")])
      # Check values
      # transition_matrix["mainstem, BON to MCN",]
      
      ##### FROM MCN TO ICH OR PRA STATE #####
      
      # MCN TO ICH OR PRA to BON to MCN transition
      b0_MIP_BM <- 1
      bflow_MIP_BM <- 0.05 # Make this vary with flow, so higher flow at MCN = higher fallback - but lower
      btemp_MIP_BM <- -0.5 # Negative relationship with temperature - when it's colder, they fall back more
      brear_MIP_BM <- c(0,0) # c(hatchery, wild) # No relationship with rear type
      borigin_MIP_BM <- c(1, 0, 0) # c(John Day River, Yakima River, Tucannon River) # JDR fish have much higher chance of falling back
      
      # BON to MCN to PRA to RIS transition
      b0_MIP_PR <- 1
      bflow_MIP_PR <- 0 # No relationship with flow
      btemp_MIP_PR <- 0.3 # When it's hotter, more likely to go upstream
      brear_MIP_PR <- c(0,0) # No relationship with rear type
      borigin_MIP_PR <- c(0,0,0) # No relationship with origin - no origins are up here
      
      # BON to MCN to ICH to LGR transition
      b0_MIP_IL <- 1
      bflow_MIP_IL <- -0.3 # Negative relationship - higher flow = lower chance of ascending
      btemp_MIP_IL <- 0.2 # Positive relationship with temperature
      brear_MIP_IL <- c(0,0) # No effect of rear type
      borigin_MIP_IL <- c(-0.5, 0, 1) # Negative for JDR, no effect on YAK, positive for Tucannon
      
      # BON to MCN to Yakima transition
      b0_MIP_YAK <- 1
      bflow_MIP_YAK <- 0 # flow not relevant for tributary entry
      btemp_MIP_YAK <- 0 # Temperature not relevant for tributary entry
      brear_MIP_YAK <- c(0,0.2) # Wild fish more likely to enter the Yakima
      borigin_MIP_YAK <- c(0.1, 2, 0.1) # Yakima fish have a much higher chance of homing
      
      phi_MIP_BM <- exp(b0_MIP_BM + btemp_MIP_BM*temp_MCN + bflow_MIP_BM*flow_MCN + brear_MIP_BM[rear] + borigin_MIP_BM[origin])
      phi_MIP_PR <- exp(b0_MIP_PR + btemp_MIP_PR*temp_PRA + bflow_MIP_PR*flow_PRA + brear_MIP_PR[rear] + borigin_MIP_PR[origin])
      phi_MIP_IL <- exp(b0_MIP_IL + btemp_MIP_IL*temp_ICH + bflow_MIP_IL*flow_ICH + brear_MIP_IL[rear] + borigin_MIP_IL[origin])
      phi_MIP_YAK <- exp(b0_MIP_YAK + btemp_MIP_YAK*temp_YAK + bflow_MIP_YAK*flow_YAK + brear_MIP_YAK[rear] + borigin_MIP_YAK[origin])
      
      transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, BON to MCN"] <- phi_MIP_BM/(1 + phi_MIP_BM + phi_MIP_PR + phi_MIP_IL + phi_MIP_YAK)
      transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS"] <- phi_MIP_PR/(1 + phi_MIP_BM + phi_MIP_PR + phi_MIP_IL + phi_MIP_YAK)
      transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR"] <- phi_MIP_IL/(1 + phi_MIP_BM + phi_MIP_PR + phi_MIP_IL + phi_MIP_YAK)
      transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River"] <- phi_MIP_YAK/(1 + phi_MIP_BM + phi_MIP_PR + phi_MIP_IL + phi_MIP_YAK)
      transition_matrix["mainstem, MCN to ICH or PRA", "loss"]  <- 1 - sum(transition_matrix["mainstem, MCN to ICH or PRA", c("mainstem, BON to MCN", "mainstem, PRA to RIS",
                                                                                                                              "mainstem, ICH to LGR", "Yakima River")])
      
      # Check values
      # transition_matrix["mainstem, MCN to ICH or PRA",]
      
      
      ##### FROM PRA to RIS STATE #####
      
      # PRA to RIS to MCN to ICH or PRA transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(2, 2, 2) # Highly positive for each rear type
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_PRA + bflow_PR_MIP*flow_PRA + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["mainstem, PRA to RIS", "loss"] <- 1 - transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["mainstem, PRA to RIS", ]
      
      
      
      
      ##### FROM ICH TO LGR STATE #####
      
      # ICH to LGR to MCN to ICH or PRA transition #
      b0_IL_MIP <- 1
      bflow_IL_MIP <- 0 # No relationship with flow
      btemp_IL_MIP <- 0.2 # Negative relationship with temperature - when it's colder, tend to go downstream more
      brear_IL_MIP <- c(0,0) # No effect of rear type
      borigin_IL_MIP <- c(1, 1, -0.5) # Positive for JDR and YAK, negative for TUC
      
      # ICH to LGR to Tucannon River transition #
      b0_IL_TUC <- 1
      bflow_IL_TUC <- 0 # No relationship with flow (tributary)
      btemp_IL_TUC <- 0 # No relationship with temperature (tributary)
      brear_IL_TUC <- c(0,0) # No effect of rear type
      borigin_IL_TUC <- c(0, 0, 2) # Highly positive for TUC to home
      
      
      phi_IL_MIP <- exp(b0_IL_MIP + btemp_IL_MIP*temp_ICH + bflow_IL_MIP*flow_ICH + brear_IL_MIP[rear] + borigin_IL_MIP[origin])
      phi_IL_TUC <- exp(b0_IL_TUC + btemp_IL_TUC*temp_TUC + bflow_IL_TUC*flow_TUC + brear_IL_TUC[rear] + borigin_IL_TUC[origin])
      
      transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- phi_IL_MIP/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- phi_IL_TUC/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "loss"] <- 1 - sum(transition_matrix["mainstem, ICH to LGR", c("mainstem, MCN to ICH or PRA", "Tucannon River")])
      
      # Check values
      # transition_matrix["mainstem, ICH to LGR",]
      
      
      ##### FROM DESCHUTES RIVER STATE #####
      # Deschutes River to mainstem, BON to MCN transition
      b0_DES_BM <- 0 # Make all of the return intercepts 1, instead of 0 - increase chance of loss param
      bflow_DES_BM <- 0 # No relationship with flow
      btemp_DES_BM <- 0 # No relationship with temperature
      brear_DES_BM <- c(-0.5,0) # Hatchery fish less likely to return to mainstem
      borigin_DES_BM <- c(0,0,0) # No relationship with rear type
      
      phi_DES_BM <- exp(b0_DES_BM + btemp_DES_BM*temp_DES + bflow_DES_BM*flow_DES + brear_DES_BM[rear] + borigin_DES_BM[origin])
      
      transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- phi_DES_BM/(1 + phi_DES_BM)
      transition_matrix["Deschutes River", "loss"] <- 1 - transition_matrix["Deschutes River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["Deschutes River", ]
      
      
      ##### FROM JOHN DAY RIVER STATE #####
      # JDR to mainstem, BON to MCN transition
      b0_PR_MIP <- 0
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(-2, 0, 0) # JDR fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_JDR + bflow_PR_MIP*flow_JDR + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["John Day River", "mainstem, BON to MCN"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["John Day River", "loss"] <- 1 - transition_matrix["John Day River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["John Day River", ]
      
      
      ##### FROM YAKIMA RIVER STATE #####
      # Yakima River to mainstem, MCN to ICH or PRA transition
      b0_PR_MIP <- 0
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, -2, 0) # Yakima R. fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_YAK + bflow_PR_MIP*flow_YAK + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["Yakima River", "loss"] <- 1 - transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["Yakima River", ]
      
      
      ##### FROM TUCANNON RIVER STATE #####
      # From Tucannon River to ICH to LGR transition
      b0_PR_MIP <- 0
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, -2) # Tucannon River fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_TUC + bflow_PR_MIP*flow_TUC + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["Tucannon River", "loss"] <- 1 - transition_matrix["Tucannon River", "mainstem, ICH to LGR"]
      
      # Check values
      # transition_matrix["Tucannon River", ]
      
      
      #####
      
      

      # Index to row of transition matrix containing the correct "from" state to get movement probabilities
      p <- transition_matrix[rownames(as.data.frame(which(movement_array[,j-1,i] == 1))),]
      
      # Choose next state using rmultinom
      movement_array[,j,i] <- rmultinom(1, size = 1, prob = p)
      
      # Store the date in the transition_date matrix
      transition_date[i,j] <- state_date[i,j]
      
      j <- j +1
    }


}

# Trim all of the all-zero columns, store in a list
# Index value of first all zero row (end of detection history)
# Make an empty list
det_hist_sim <- list()

for (i in 1:nfish){
  hist_end_index <- match(0, colSums(movement_array[,,i]))
  
  det_hist_sim[[i]] <- movement_array[,1:(hist_end_index-1),i]
}

# Trim all of the NAs from transition date matrix, store as a list
transition_date_sim <- list()

for (i in 1:nfish){
  hist_end_index <- match(NA, transition_date[i,])
  
  transition_date_sim[[i]] <- transition_date[i,1:(hist_end_index-1)]
}

# Export the data simulation

# Export the detection history
# Here as a list
saveRDS(det_hist_sim, here::here("simulation", "sim_600_list.rds"))

# Now as an array
# Trim to maximum number of observations
n.obs <- unlist(lapply(det_hist_sim, ncol))
max_obs <- max(n.obs)
det_hist_sim_array <- movement_array[,1:max_obs,]
saveRDS(det_hist_sim_array, here::here("simulation", "sim_600_array.rds"))

# Export the transition dates
saveRDS(transition_date_sim, here::here("simulation", "dates_600.rds"))
# As an array
transition_date_sim_matrix <- transition_date[,1:max_obs]

date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2016-12-31"))
transition_date_sim_matrix %>% 
  as_tibble() %>% 
  mutate_all(date_numeric) -> transition_date_sim_numeric
  
saveRDS(transition_date_sim_numeric, here::here("simulation", "dates_600_matrix.rds"))

# Export the covariates
write.csv(temp_sim_zscore_datenumeric, here::here("simulation", "temp_600.csv"), row.names = FALSE)
write.csv(flow_sim_zscore_datenumeric, here::here("simulation", "flow_600.csv"), row.names = FALSE)
write.csv(fish_sim_cat_data, here::here("simulation", "origin_rear_600.csv"), row.names = FALSE)


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

probs_df <- as.data.frame(matrix(nrow = 10, ncol = length(seq(-10,10, 0.1))))

for (j in 1:length(seq(-10,10, 0.1))){
  values <- seq(-10,10, 0.1)
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

probs_df_long$b0_BM_JDR <- rep(seq(-10, 10, 0.1), 5)

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

probs_df <- as.data.frame(matrix(nrow = 10, ncol = length(seq(-10 ,10, 0.1))))

for (j in 1:length(seq(-10,10, 0.1))){
  values <- seq(-10,10, 0.1)
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

probs_df_long$b0_BM_JDR <- rep(seq(-10, 10, 0.1), 5)

ggplot(probs_df_long, aes(x = b0_BM_JDR, y = prob, color = state)) +
  geom_point() +
  scale_color_tableau(palette = "Tableau 10")

##### ACTUAL SIMULATION #####

# store test covariate values
i <- 1
temp_BON <- 0.5
flow_BON <- -0.5
temp_MCN <- 0.5
flow_MCN <- -0.5
temp_ICH <- 0.5
flow_ICH <- -0.5
temp_PRA <- 0.5
flow_PRA <- -0.5
temp_DES <- 0
flow_DES <- 0
temp_JDR <- 0
flow_JDR <- 0
temp_TUC <- 0
flow_TUC <- 0

# Store simulation parameter values

##### FROM BON TO MCN STATE #####

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
brear_BM_DES <- c(0.5,0) # higher probability of straying to Deschutes if hatchery
borigin_BM_DES <- c(0.25, 0, 0) # JDR has higher probability of stray to DES because it's closer

# BON to MCN to JDR transition
b0_BM_JDR <- 1
bflow_BM_JDR <- 0 # flow not relevant for tributary entry
btemp_BM_JDR <- 0 # Temperature not relevant for tributary entry
brear_BM_JDR <- c(0,0) # no effect of rear type
borigin_BM_JDR <- c(2, 0.1, 0.1) # JDR fish have much higher chance of homing

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
# Check values
# transition_matrix["mainstem, BON to MCN",]

##### FROM MCN TO ICH OR PRA STATE #####

# MCN TO ICH OR PRA to BON to MCN transition
b0_MIP_BM <- 1
bflow_MIP_BM <- 0.05 # Make this vary with flow, so higher flow at MCN = higher fallback - but lower
btemp_MIP_BM <- -0.5 # Negative relationship with temperature - when it's colder, they fall back more
brear_MIP_BM <- c(0,0) # c(hatchery, wild) # No relationship with rear type
borigin_MIP_BM <- c(1, 0, 0) # c(John Day River, Yakima River, Tucannon River) # JDR fish have much higher chance of falling back

# BON to MCN to PRA to RIS transition
b0_MIP_PR <- 1
bflow_MIP_PR <- 0 # No relationship with flow
btemp_MIP_PR <- 0.3 # When it's hotter, more likely to go upstream
brear_MIP_PR <- c(0,0) # No relationship with rear type
borigin_MIP_PR <- c(0,0,0) # No relationship with origin - no origins are up here

# BON to MCN to ICH to LGR transition
b0_MIP_IL <- 1
bflow_MIP_IL <- -0.3 # Negative relationship - higher flow = lower chance of ascending
btemp_MIP_IL <- 0.2 # Positive relationship with temperature
brear_MIP_IL <- c(0,0) # No effect of rear type
borigin_MIP_IL <- c(-0.5, 0, 1) # Negative for JDR, no effect on YAK, positive for Tucannon

# BON to MCN to Yakima transition
b0_MIP_YAK <- 1
bflow_MIP_YAK <- 0 # flow not relevant for tributary entry
btemp_MIP_YAK <- 0 # Temperature not relevant for tributary entry
brear_MIP_YAK <- c(0,0.2) # Wild fish more likely to enter the Yakima
borigin_MIP_YAK <- c(0.1, 2, 0.1) # Yakima fish have a much higher chance of homing

phi_MIP_BM <- exp(b0_MIP_BM + btemp_MIP_BM*temp_MCN + bflow_MIP_BM*flow_MCN + brear_MIP_BM[i] + borigin_MIP_BM[i])
phi_MIP_PR <- exp(b0_MIP_PR + btemp_MIP_PR*temp_PRA + bflow_MIP_PR*flow_PRA + brear_MIP_PR[i] + borigin_MIP_PR[i])
phi_MIP_IL <- exp(b0_MIP_IL + btemp_MIP_IL*temp_ICH + bflow_MIP_IL*flow_ICH + brear_MIP_IL[i] + borigin_MIP_IL[i])
phi_MIP_YAK <- exp(b0_MIP_YAK + btemp_MIP_YAK*temp_YAK + bflow_MIP_YAK*flow_YAK + brear_MIP_YAK[i] + borigin_MIP_YAK[i])

transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, BON to MCN"] <- phi_MIP_BM/(1 + phi_MIP_BM + phi_MIP_PR + phi_MIP_IL + phi_MIP_YAK)
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS"] <- phi_MIP_PR/(1 + phi_MIP_BM + phi_MIP_PR + phi_MIP_IL + phi_MIP_YAK)
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR"] <- phi_MIP_IL/(1 + phi_MIP_BM + phi_MIP_PR + phi_MIP_IL + phi_MIP_YAK)
transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River"] <- phi_MIP_YAK/(1 + phi_MIP_BM + phi_MIP_PR + phi_MIP_IL + phi_MIP_YAK)
transition_matrix["mainstem, MCN to ICH or PRA", "loss"]  <- 1 - sum(transition_matrix["mainstem, MCN to ICH or PRA", c("mainstem, BON to MCN", "mainstem, PRA to RIS",
"mainstem, ICH to LGR", "Yakima River")])

# Check values
# transition_matrix["mainstem, MCN to ICH or PRA",]


##### FROM PRA to RIS STATE #####

# PRA to RIS to MCN to ICH or PRA transition
b0_PR_MIP <- 1
bflow_PR_MIP <- 0 # No relationship with flow
btemp_PR_MIP <- 0 # No relationship with temperature
brear_PR_MIP <- c(0,0) # No effect of rear type
borigin_PR_MIP <- c(2, 2, 2) # Highly positive for each rear type

phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_BON + bflow_PR_MIP*flow_BON + brear_PR_MIP[i] + borigin_PR_MIP[i])

transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
transition_matrix["mainstem, PRA to RIS", "loss"] <- 1 - transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"]

# Check values
# transition_matrix["mainstem, PRA to RIS", ]




##### FROM ICH TO LGR STATE #####

# ICH to LGR to MCN to ICH or PRA transition #
b0_IL_MIP <- 1
bflow_IL_MIP <- 0 # No relationship with flow
btemp_IL_MIP <- 0.2 # Negative relationship with temperature - when it's colder, tend to go downstream more
brear_IL_MIP <- c(0,0) # No effect of rear type
borigin_IL_MIP <- c(1, 1, -0.5) # Positive for JDR and YAK, negative for TUC

# ICH to LGR to Tucannon River transition #
b0_IL_TUC <- 1
bflow_IL_TUC <- 0 # No relationship with flow (tributary)
btemp_IL_TUC <- 0 # No relationship with temperature (tributary)
brear_IL_TUC <- c(0,0) # No effect of rear type
borigin_IL_TUC <- c(0, 0, 2) # Highly positive for TUC to home


phi_IL_MIP <- exp(b0_IL_MIP + btemp_IL_MIP*temp_ICH + bflow_IL_MIP*flow_ICH + brear_IL_MIP[i] + borigin_IL_MIP[i])
phi_IL_TUC <- exp(b0_IL_TUC + btemp_IL_TUC*temp_TUC + bflow_IL_TUC*flow_TUC + brear_IL_TUC[i] + borigin_IL_TUC[i])

transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- phi_IL_MIP/(1 + phi_IL_MIP + phi_IL_TUC)
transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- phi_IL_TUC/(1 + phi_IL_MIP + phi_IL_TUC)
transition_matrix["mainstem, ICH to LGR", "loss"] <- 1 - sum(transition_matrix["mainstem, ICH to LGR", c("mainstem, MCN to ICH or PRA", "Tucannon River")])

# Check values
# transition_matrix["mainstem, ICH to LGR",]


##### FROM DESCHUTES RIVER STATE #####
# Deschutes River to mainstem, BON to MCN transition
b0_DES_BM <- 0 # Make all of the return intercepts 1, instead of 0 - increase chance of loss param
bflow_DES_BM <- 0 # No relationship with flow
btemp_DES_BM <- 0 # No relationship with temperature
brear_DES_BM <- c(-0.5,0) # Hatchery fish less likely to return to mainstem
borigin_DES_BM <- c(0,0,0) # No relationship with rear type

phi_DES_BM <- exp(b0_DES_BM + btemp_DES_BM*temp_DES + bflow_DES_BM*flow_DES + brear_DES_BM[i] + borigin_DES_BM[i])

transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- phi_DES_BM/(1 + phi_DES_BM)
transition_matrix["Deschutes River", "loss"] <- 1 - transition_matrix["Deschutes River", "mainstem, BON to MCN"]

# Check values
# transition_matrix["Deschutes River", ]


##### FROM JOHN DAY RIVER STATE #####
# JDR to mainstem, BON to MCN transition
b0_PR_MIP <- 0
bflow_PR_MIP <- 0 # No relationship with flow
btemp_PR_MIP <- 0 # No relationship with temperature
brear_PR_MIP <- c(0,0) # No effect of rear type
borigin_PR_MIP <- c(-2, 0, 0) # JDR fish less likely to return

phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_JDR + bflow_PR_MIP*flow_JDR + brear_PR_MIP[i] + borigin_PR_MIP[i])

transition_matrix["John Day River", "mainstem, BON to MCN"] <- phi_PR_MIP/(1 + phi_PR_MIP)
transition_matrix["John Day River", "loss"] <- 1 - transition_matrix["John Day River", "mainstem, BON to MCN"]

# Check values
# transition_matrix["John Day River", ]


##### FROM YAKIMA RIVER STATE #####
# Yakima River to mainstem, MCN to ICH or PRA transition
b0_PR_MIP <- 0
bflow_PR_MIP <- 0 # No relationship with flow
btemp_PR_MIP <- 0 # No relationship with temperature
brear_PR_MIP <- c(0,0) # No effect of rear type
borigin_PR_MIP <- c(0, -2, 0) # Yakima R. fish less likely to return

phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_YAK + bflow_PR_MIP*flow_YAK + brear_PR_MIP[i] + borigin_PR_MIP[i])

transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
transition_matrix["Yakima River", "loss"] <- 1 - transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"]

# Check values
# transition_matrix["Yakima River", ]


##### FROM TUCANNON RIVER STATE #####
# From Tucannon River to ICH to LGR transition
b0_PR_MIP <- 0
bflow_PR_MIP <- 0 # No relationship with flow
btemp_PR_MIP <- 0 # No relationship with temperature
brear_PR_MIP <- c(0,0) # No effect of rear type
borigin_PR_MIP <- c(0, 0, -2) # Tucannon River fish less likely to return

phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_TUC + bflow_PR_MIP*flow_TUC + brear_PR_MIP[i] + borigin_PR_MIP[i])

transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- phi_PR_MIP/(1 + phi_PR_MIP)
transition_matrix["Tucannon River", "loss"] <- 1 - transition_matrix["Tucannon River", "mainstem, ICH to LGR"]

# Check values
# transition_matrix["Tucannon River", ]



# NOTE: NEED TO UPDATE INDEXING FOR CATEGORICAL VARIABLES










