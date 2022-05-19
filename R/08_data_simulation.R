# 08: Data simulation

# Here we are generating a simulated dataset, in order to test our model

# Load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)
library(ggthemes)
library(foreach)
library(parallel)
library(doParallel)

# Set up parallel backend
n.cores <- detectCores() - 1

#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)

#check cluster definition (optional)
print(my.cluster)

doParallel::registerDoParallel(cl = my.cluster)

#check if it is registered (optional)
foreach::getDoParRegistered()

#how many workers are available? (optional)
foreach::getDoParWorkers()

# When you're done, turn off the cluster:
parallel::stopCluster(cl = my.cluster)

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




##### Load covariate data #####
# The simulated covariate data is too autocorrelated (as in, it's exactly the same)
# Let's load the actual temperature and flow data

flow_data <- read.csv(here::here("covariate_data", "flow_by_state.csv"), row.names = 1)
temp_data <- read.csv(here::here("covariate_data", "temperature_by_state.csv"), row.names = 1)

flow_data %>% 
  dplyr::rename(BON = mouth.to.BON,
                ICH = MCN.to.ICH.or.PRA..ICH.,
                MCN = BON.to.MCN,
                PRA = MCN.to.ICH.or.PRA..PRA.) %>% 
  dplyr::select(date, BON, ICH, MCN, PRA) %>% 
  mutate(date = ymd(date)) %>% 
  subset(date >= ymd("2017-01-01") & date <= ymd("2019-12-31")) -> flow_data

# First, fill NAs using linear interpolation (TEMPORARY)

# Also, remove strange all NA rows
flow_data %>% 
  subset(., !(is.na(BON) & is.na(MCN) & is.na(ICH) & is.na(PRA))) -> flow_data
library(zoo)
flow_data$BON <- na.approx(flow_data$BON)
flow_data$MCN <- na.approx(flow_data$MCN)
flow_data$PRA <- na.approx(flow_data$PRA)
flow_data$ICH <- na.approx(flow_data$ICH)

flow_data %>% 
  mutate(BON = (BON - mean(BON))/sd(BON)) %>% 
  mutate(MCN = (MCN - mean(MCN))/sd(MCN)) %>% 
  mutate(PRA = (PRA - mean(PRA))/sd(PRA)) %>% 
  mutate(ICH = (ICH - mean(ICH))/sd(ICH)) %>% 
  # Store dummy values (all 0) for flow in tribs
  mutate(JDR = 0) %>% 
  mutate(DES = 0) %>% 
  mutate(YAK = 0) %>% 
  mutate(TUC = 0) -> flow_zscore_df

# Change dates to numeric (for JAGS)
flow_zscore_df %>% 
  mutate(date_numeric = as.numeric(date - ymd("2016-12-31"))) %>% 
  relocate(date_numeric) %>% 
  dplyr::select(-date) %>% 
  dplyr::rename(date = date_numeric) -> flow_sim_zscore_datenumeric

rownames(flow_sim_zscore_datenumeric) <- NULL
flow_sim_zscore_datenumeric %>% 
  column_to_rownames("date") %>% 
  relocate(ICH, .after = PRA) -> flow_sim_zscore_datenumeric

write.csv(flow_sim_zscore_datenumeric, here::here("simulation", "flow_600.csv"))


# temp_data %>% 
#   mutate(BON = (BON - mean(BON))/sd(BON)) %>% 
#   mutate(MCN = (MCN - mean(MCN))/sd(MCN)) %>% 
#   mutate(PRA = (PRA - mean(PRA))/sd(PRA)) %>% 
#   mutate(ICH = (ICH - mean(ICH))/sd(ICH)) %>% 
#   # Store dummy values (all 0) for flow in tribs
#   mutate(JDR = 0) %>% 
#   mutate(DES = 0) %>% 
#   mutate(YAK = 0) %>% 
#   mutate(TUC = 0)-> temp_zscore_df

# For now, let's keep the simulated temperature and take the actual flow data
temp_zscore_df <- temp_sim_zscore_df
flow_sim_zscore_df <- flow_zscore_df


# Let's check correlation between temp and flow
bon_tempflow <- data.frame(temp = temp_zscore_df$BON, flow = flow_zscore_df$BON)
mcn_tempflow <- data.frame(temp = temp_zscore_df$MCN, flow = flow_zscore_df$MCN)
ich_tempflow <- data.frame(temp = temp_zscore_df$ICH, flow = flow_zscore_df$ICH)
pra_tempflow <- data.frame(temp = temp_zscore_df$PRA, flow = flow_zscore_df$PRA)

cor(bon_tempflow$temp, bon_tempflow$flow)
cor(mcn_tempflow$temp, mcn_tempflow$flow)
cor(ich_tempflow$temp, ich_tempflow$flow)
cor(pra_tempflow$temp, pra_tempflow$flow)

ggplot(bon_tempflow, aes(x = temp, y = flow)) +
  geom_point()

ggplot(mcn_tempflow, aes(x = temp, y = flow)) +
  geom_point()

ggplot(ich_tempflow, aes(x = temp, y = flow)) +
  geom_point()

ggplot(pra_tempflow, aes(x = temp, y = flow)) +
  geom_point()


##### Set up Base model: no covariates #####

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


##### SIMULATION FUNCTION 1: NO COVARIATES #####
# This is a function to simulate data with no covariates, with every b0 set to 1.
# the only argument needed is "nfish"

nocov_sim <- function(nfish, fish_sim_cat_data){
  # Take a uniform distribution between July 1 and September 30
  BON_arrival_dates_sim <- seq(ymd("2017-07-01"), ymd("2019-09-30"), by = "days")[round(runif(3000, min = 1, max = 92), 0)]
  
  state_date <- matrix(nrow = nfish, ncol = 100)
  transition_date <- matrix(nrow = nfish, ncol = 100)
  
  for (i in 1:nfish){ # for each fish
    print(paste0("i = ", i))
    # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
    movement_array["mainstem, BON to MCN", 1, i] <- 1
    # Populate state date matrix with date at arrival at BON
    state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
    # Store this in the transition date matrix
    transition_date[i,1] <- state_date[i,1]
    
    # populate the rest of the state dates for this fish by adding two weeks
    for (k in 2:100){
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
      borigin_MB_BM <- c(0,0,0) # No effect of origin
      
      phi_MB_BM <- exp(b0_MB_BM + btemp_MB_BM*temp_BON + bflow_MB_BM*flow_BON + brear_MB_BM[rear] + borigin_MB_BM[origin])
      
      transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- phi_MB_BM/(1 + phi_MB_BM)
      transition_matrix["mainstem, mouth to BON", "loss"] <- 1 - transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"]
      
      # Check values
      transition_matrix["mainstem, mouth to BON", ]
      
      ##### FROM BON TO MCN STATE #####
      
      # BON to MCN to Mouth to BON transition
      b0_BM_MB <- 1
      bflow_BM_MB <- 0 # Make this vary with flow, so higher flow at BON = higher fallback over BON
      btemp_BM_MB <- 0 # No relationship with temperature
      brear_BM_MB <- c(0,0) # c(hatchery, wild) # No relationship with rear type
      borigin_BM_MB <- c(0, 0, 0) # c(John Day River, Yakima River, Tucannon River) # No relationship with origin
      
      # BON to MCN to MCN to ICH or PRA transition
      b0_BM_MIP <- 1
      bflow_BM_MIP <- 0 # No relationship with flows
      btemp_BM_MIP <- 0 # When it's hotter, more likely to go upstream
      brear_BM_MIP <- c(0,0) # No relationship with rear type
      borigin_BM_MIP <- c(0,0,0) # JDR fish have lower probability of overshooting MCN than the fish from origins upstream of MCN
      
      # BON to MCN to DES R transition
      b0_BM_DES <- 1
      bflow_BM_DES <- 0 # flow not relevant for tributary entry
      btemp_BM_DES <- 0 # Temperature not relevant for tributary entry
      brear_BM_DES <- c(0,0) # higher probability of straying to Deschutes if hatchery
      borigin_BM_DES <- c(0, 0, 0) # JDR has higher probability of stray to DES because it's closer
      
      # BON to MCN to JDR transition
      b0_BM_JDR <- 1
      bflow_BM_JDR <- 0 # flow not relevant for tributary entry
      btemp_BM_JDR <- 0 # Temperature not relevant for tributary entry
      brear_BM_JDR <- c(0,0) # no effect of rear type
      borigin_BM_JDR <- c(0,0,0) # JDR fish have much higher chance of homing
      
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
      bflow_MIP_BM <- 0 # Make this vary with flow, so higher flow at MCN = higher fallback - but lower
      btemp_MIP_BM <- 0 # Negative relationship with temperature - when it's colder, they fall back more
      brear_MIP_BM <- c(0,0) # c(hatchery, wild) # No relationship with rear type
      borigin_MIP_BM <- c(0, 0, 0) # c(John Day River, Yakima River, Tucannon River) # JDR fish have much higher chance of falling back
      
      # BON to MCN to PRA to RIS transition
      b0_MIP_PR <- 1
      bflow_MIP_PR <- 0 # No relationship with flow
      btemp_MIP_PR <- 0 # When it's hotter, more likely to go upstream
      brear_MIP_PR <- c(0,0) # No relationship with rear type
      borigin_MIP_PR <- c(0,0,0) # No relationship with origin - no origins are up here
      
      # BON to MCN to ICH to LGR transition
      b0_MIP_IL <- 1
      bflow_MIP_IL <- 0 # Negative relationship - higher flow = lower chance of ascending
      btemp_MIP_IL <- 0 # Positive relationship with temperature
      brear_MIP_IL <- c(0,0) # No effect of rear type
      borigin_MIP_IL <- c(0,0,0) # Negative for JDR, no effect on YAK, positive for Tucannon
      
      # BON to MCN to Yakima transition
      b0_MIP_YAK <- 1
      bflow_MIP_YAK <- 0 # flow not relevant for tributary entry
      btemp_MIP_YAK <- 0 # Temperature not relevant for tributary entry
      brear_MIP_YAK <- c(0,0) # Wild fish more likely to enter the Yakima
      borigin_MIP_YAK <- c(0,0,0) # Yakima fish have a much higher chance of homing
      
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
      borigin_PR_MIP <- c(0,0,0) # Highly positive for each rear type
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_PRA + bflow_PR_MIP*flow_PRA + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["mainstem, PRA to RIS", "loss"] <- 1 - transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["mainstem, PRA to RIS", ]
      
      
      
      
      ##### FROM ICH TO LGR STATE #####
      
      # ICH to LGR to MCN to ICH or PRA transition #
      b0_IL_MIP <- 1
      bflow_IL_MIP <- 0 # No relationship with flow
      btemp_IL_MIP <- 0 # Negative relationship with temperature - when it's colder, tend to go downstream more
      brear_IL_MIP <- c(0,0) # No effect of rear type
      borigin_IL_MIP <- c(0,0,0) # Positive for JDR and YAK, negative for TUC
      
      # ICH to LGR to Tucannon River transition #
      b0_IL_TUC <- 1
      bflow_IL_TUC <- 0 # No relationship with flow (tributary)
      btemp_IL_TUC <- 0 # No relationship with temperature (tributary)
      brear_IL_TUC <- c(0,0) # No effect of rear type
      borigin_IL_TUC <- c(0, 0, 0) # Highly positive for TUC to home
      
      
      phi_IL_MIP <- exp(b0_IL_MIP + btemp_IL_MIP*temp_ICH + bflow_IL_MIP*flow_ICH + brear_IL_MIP[rear] + borigin_IL_MIP[origin])
      phi_IL_TUC <- exp(b0_IL_TUC + btemp_IL_TUC*temp_TUC + bflow_IL_TUC*flow_TUC + brear_IL_TUC[rear] + borigin_IL_TUC[origin])
      
      transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- phi_IL_MIP/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- phi_IL_TUC/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "loss"] <- 1 - sum(transition_matrix["mainstem, ICH to LGR", c("mainstem, MCN to ICH or PRA", "Tucannon River")])
      
      # Check values
      # transition_matrix["mainstem, ICH to LGR",]
      
      
      ##### FROM DESCHUTES RIVER STATE #####
      # Deschutes River to mainstem, BON to MCN transition
      b0_DES_BM <- 1 # Make all of the return intercepts 1, instead of 0 - increase chance of loss param
      bflow_DES_BM <- 0 # No relationship with flow
      btemp_DES_BM <- 0 # No relationship with temperature
      brear_DES_BM <- c(0,0) # Hatchery fish less likely to return to mainstem
      borigin_DES_BM <- c(0,0,0) # No relationship with rear type
      
      phi_DES_BM <- exp(b0_DES_BM + btemp_DES_BM*temp_DES + bflow_DES_BM*flow_DES + brear_DES_BM[rear] + borigin_DES_BM[origin])
      
      transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- phi_DES_BM/(1 + phi_DES_BM)
      transition_matrix["Deschutes River", "loss"] <- 1 - transition_matrix["Deschutes River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["Deschutes River", ]
      
      
      ##### FROM JOHN DAY RIVER STATE #####
      # JDR to mainstem, BON to MCN transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # JDR fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_JDR + bflow_PR_MIP*flow_JDR + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["John Day River", "mainstem, BON to MCN"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["John Day River", "loss"] <- 1 - transition_matrix["John Day River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["John Day River", ]
      
      
      ##### FROM YAKIMA RIVER STATE #####
      # Yakima River to mainstem, MCN to ICH or PRA transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # Yakima R. fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_YAK + bflow_PR_MIP*flow_YAK + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["Yakima River", "loss"] <- 1 - transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["Yakima River", ]
      
      
      ##### FROM TUCANNON RIVER STATE #####
      # From Tucannon River to ICH to LGR transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # Tucannon River fish less likely to return
      
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
  # Trim to maximum number of observations
  n.obs <- unlist(lapply(det_hist_sim, ncol))
  max_obs <- max(n.obs)
  det_hist_sim_array <- movement_array[,1:max_obs,]
  
  # As an array
  transition_date_sim_matrix <- transition_date[,1:max_obs]
  
  date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2016-12-31"))
  transition_date_sim_matrix %>% 
    as_tibble() %>% 
    mutate_all(date_numeric) -> transition_date_sim_numeric
  
  # Return a list of two things: [[1]] = the detection history, [[2]] = the dates
  
  return(list(det_hist_sim_array, transition_date_sim_numeric))
}

fish_sim_cat_data_600 <- data.frame(tag_code = seq(1, 600, 1),
                                    natal_origin = c(rep(1, 200),
                                                     rep(2, 200),
                                                     rep(3, 200)),
                                    rear_type = rep(c(rep(1, 100), rep(2, 100)), 3))

fish_sim_cat_data_1200 <- data.frame(tag_code = seq(1, 1200, 1),
                                natal_origin = c(rep(1, 400),
                                                 rep(2, 400),
                                                 rep(3, 400)),
                                rear_type = rep(c(rep(1, 200), rep(2, 200)), 3))

fish_sim_cat_data_3000 <- data.frame(tag_code = seq(1, 3000, 1),
                                natal_origin = c(rep(1, 1000),
                                                 rep(2, 1000),
                                                 rep(3, 1000)),
                                rear_type = rep(c(rep(1, 500), rep(2, 500)), 3))



nocov_sim(nfish = 600, fish_sim_cat_data = fish_sim_cat_data_600)


##### Loop the simulation

### 600 fish
sim_600_hist_list <- list()
sim_600_dates_list <- list()
nfish = 600
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- nocov_sim(nfish = 600, fish_sim_cat_data = fish_sim_cat_data_600)
  
  sim_600_hist_list[[z]] <- sim_run[[1]]
  sim_600_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_600_hist_list, here::here("simulation", "sim_600_hist_list.rds"))

saveRDS(sim_600_dates_list, here::here("simulation", "sim_600_dates_list.rds"))
write.csv(fish_sim_cat_data_600, here::here("simulation", "origin_rear_600.csv"), row.names = FALSE)


### 1200 fish
sim_1200_hist_list <- list()
sim_1200_dates_list <- list()
nfish = 1200
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- nocov_sim(nfish = 1200, fish_sim_cat_data = fish_sim_cat_data_1200)
  
  sim_1200_hist_list[[z]] <- sim_run[[1]]
  sim_1200_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_1200_hist_list, here::here("simulation", "sim_1200_hist_list.rds"))
saveRDS(sim_1200_dates_list, here::here("simulation", "sim_1200_dates_list.rds"))
write.csv(fish_sim_cat_data_1200, here::here("simulation", "origin_rear_1200.csv"), row.names = FALSE)


### 3000 fish
sim_3000_hist_list <- list()
sim_3000_dates_list <- list()
nfish = 3000
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- nocov_sim(nfish = 3000, fish_sim_cat_data = fish_sim_cat_data_3000)
  
  sim_3000_hist_list[[z]] <- sim_run[[1]]
  sim_3000_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_3000_hist_list, here::here("simulation", "sim_3000_hist_list.rds"))
saveRDS(sim_3000_dates_list, here::here("simulation", "sim_3000_dates_list.rds"))
write.csv(fish_sim_cat_data_3000, here::here("simulation", "origin_rear_3000.csv"), row.names = FALSE)


##### NOT USED: SIMULATION FUNCTION 2: All 4 COVARIATES #####
# This is a function to simulate data with no covariates, with every b0 set to 1.
# the only argument needed is "nfish"

cov_sim <- function(nfish, fish_sim_cat_data){
  # Take a uniform distribution between July 1 and September 30
  BON_arrival_dates_sim <- seq(ymd("2017-07-01"), ymd("2019-09-30"), by = "days")[round(runif(3000, min = 1, max = 92), 0)]
  
  state_date <- matrix(nrow = nfish, ncol = 100)
  transition_date <- matrix(nrow = nfish, ncol = 100)
  
  for (i in 1:nfish){ # for each fish
    print(paste0("i = ", i))
    # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
    movement_array["mainstem, BON to MCN", 1, i] <- 1
    # Populate state date matrix with date at arrival at BON
    state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
    # Store this in the transition date matrix
    transition_date[i,1] <- state_date[i,1]
    
    # populate the rest of the state dates for this fish by adding two weeks
    for (k in 2:100){
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
      bflow_BM_MIP <- 0 # No relationship with flows
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
  # Trim to maximum number of observations
  n.obs <- unlist(lapply(det_hist_sim, ncol))
  max_obs <- max(n.obs)
  det_hist_sim_array <- movement_array[,1:max_obs,]
  
  # As an array
  transition_date_sim_matrix <- transition_date[,1:max_obs]
  
  date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2016-12-31"))
  transition_date_sim_matrix %>% 
    as_tibble() %>% 
    mutate_all(date_numeric) -> transition_date_sim_numeric
  
  # Return a list of two things: [[1]] = the detection history, [[2]] = the dates
  
  return(list(det_hist_sim_array, transition_date_sim_numeric))
}

fish_sim_cat_data_600 <- data.frame(tag_code = seq(1, 600, 1),
                                    natal_origin = c(rep(1, 200),
                                                     rep(2, 200),
                                                     rep(3, 200)),
                                    rear_type = rep(c(rep(1, 100), rep(2, 100)), 3))

fish_sim_cat_data_1200 <- data.frame(tag_code = seq(1, 1200, 1),
                                     natal_origin = c(rep(1, 400),
                                                      rep(2, 400),
                                                      rep(3, 400)),
                                     rear_type = rep(c(rep(1, 200), rep(2, 200)), 3))

fish_sim_cat_data_3000 <- data.frame(tag_code = seq(1, 3000, 1),
                                     natal_origin = c(rep(1, 1000),
                                                      rep(2, 1000),
                                                      rep(3, 1000)),
                                     rear_type = rep(c(rep(1, 500), rep(2, 500)), 3))


##### Loop the simulation

### 600 fish
sim_600_hist_list <- list()
sim_600_dates_list <- list()
nfish = 600
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- cov_sim(nfish = 600, fish_sim_cat_data = fish_sim_cat_data_600)
  
  sim_600_hist_list[[z]] <- sim_run[[1]]
  sim_600_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_600_hist_list, here::here("simulation", "sim_600_cov_hist_list.rds"))

saveRDS(sim_600_dates_list, here::here("simulation", "sim_600_cov_dates_list.rds"))
write.csv(fish_sim_cat_data_600, here::here("simulation", "origin_rear_600.csv"), row.names = FALSE)


### 1200 fish
sim_1200_hist_list <- list()
sim_1200_dates_list <- list()
nfish = 1200
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- cov_sim(nfish = 1200, fish_sim_cat_data = fish_sim_cat_data_1200)
  
  sim_1200_hist_list[[z]] <- sim_run[[1]]
  sim_1200_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_1200_hist_list, here::here("simulation", "sim_1200_cov_hist_list.rds"))
saveRDS(sim_1200_dates_list, here::here("simulation", "sim_1200_cov_dates_list.rds"))
write.csv(fish_sim_cat_data_1200, here::here("simulation", "origin_rear_1200.csv"), row.names = FALSE)


### 3000 fish
sim_3000_hist_list <- list()
sim_3000_dates_list <- list()
nfish = 3000
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- cov_sim(nfish = 3000, fish_sim_cat_data = fish_sim_cat_data_3000)
  
  sim_3000_hist_list[[z]] <- sim_run[[1]]
  sim_3000_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_3000_hist_list, here::here("simulation", "sim_3000_cov_hist_list.rds"))
saveRDS(sim_3000_dates_list, here::here("simulation", "sim_3000_cov_dates_list.rds"))
write.csv(fish_sim_cat_data_3000, here::here("simulation", "origin_rear_3000.csv"), row.names = FALSE)






##### SIMULATION FUNCTION 3: ONLY TEMP AND FLOW #####
# This is a function to simulate data with temperature and flow covariates (no rear or origin), with every b0 set to 1.
# the only argument needed is "nfish"
# I think that the reason we're having some issues with estimation is because b0 and brear borigin all function as intercepts.

cov_continuous_sim <- function(nfish, fish_sim_cat_data, movement_array = movement_array, transition_matrix = transition_matrix){
  # Take a uniform distribution between July 1 and September 30
  BON_arrival_dates_sim <- seq(ymd("2017-07-01"), ymd("2019-09-30"), by = "days")[round(runif(3000, min = 1, max = 92), 0)]
  
  state_date <- matrix(nrow = nfish, ncol = 100)
  transition_date <- matrix(nrow = nfish, ncol = 100)
  
  for (i in 1:nfish){ # for each fish
    print(paste0("i = ", i))
    # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
    movement_array["mainstem, BON to MCN", 1, i] <- 1
    # Populate state date matrix with date at arrival at BON
    state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
    # Store this in the transition date matrix
    transition_date[i,1] <- state_date[i,1]
    
    # populate the rest of the state dates for this fish by adding two weeks
    for (k in 2:100){
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
      temps <- temp_zscore_df[temp_zscore_df$date == state_date[i,j-1],][2:ncol(temp_zscore_df)]
      temp_BON <- temps[1,1]; temp_MCN <- temps[1,2]; temp_PRA <- temps[1,3]; temp_ICH <- temps[1,4]; temp_JDR <- temps[1,5]; temp_DES <- temps[1,6]; temp_YAK <- temps[1,7]; temp_TUC <- temps[1,8]
      flows <- flow_zscore_df[flow_zscore_df$date == state_date[i,j-1],][2:ncol(flow_zscore_df)]
      flow_BON <- flows[1,1]; flow_MCN <- flows[1,2]; flow_PRA <- flows[1,3]; flow_ICH <- flows[1,4]; flow_JDR <- flows[1,5]; flow_DES <- flows[1,6]; flow_YAK <- flows[1,7]; flow_TUC <- flows[1,8]
      
      
      # The function:
      ##### FROM MOUTH TO BON STATE #####
      
      # Mouth to BON to BON to MCN transition
      b0_MB_BM <- 1
      bflow_MB_BM <- 0 # No relationship with flow
      btemp_MB_BM <- 0 # No relationship with temperature
      brear_MB_BM <- c(0,0) # No effect of rear type
      borigin_MB_BM <- c(0,0,0) # Highly positive for each rear type
      
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
      bflow_BM_MIP <- 0 # No relationship with flows
      btemp_BM_MIP <- 0.5 # When it's hotter, more likely to go upstream
      brear_BM_MIP <- c(0,0) # No relationship with rear type
      borigin_BM_MIP <- c(0, 0, 0) # JDR fish have lower probability of overshooting MCN than the fish from origins upstream of MCN
      
      # BON to MCN to DES R transition
      b0_BM_DES <- 1
      bflow_BM_DES <- 0 # flow not relevant for tributary entry
      btemp_BM_DES <- 0 # Temperature not relevant for tributary entry
      brear_BM_DES <- c(0,0) # higher probability of straying to Deschutes if hatchery
      borigin_BM_DES <- c(0, 0, 0) # JDR has higher probability of stray to DES because it's closer
      
      # BON to MCN to JDR transition
      b0_BM_JDR <- 1
      bflow_BM_JDR <- 0 # flow not relevant for tributary entry
      btemp_BM_JDR <- 0 # Temperature not relevant for tributary entry
      brear_BM_JDR <- c(0,0) # no effect of rear type
      borigin_BM_JDR <- c(0, 0, 0) # JDR fish have much higher chance of homing
      
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
      borigin_MIP_BM <- c(0, 0, 0) # c(John Day River, Yakima River, Tucannon River) # JDR fish have much higher chance of falling back
      
      # MCN TO ICH OR PRA to PRA to RIS transition
      b0_MIP_PR <- 1
      bflow_MIP_PR <- 0 # No relationship with flow
      btemp_MIP_PR <- 0.3 # When it's hotter, more likely to go upstream
      brear_MIP_PR <- c(0,0) # No relationship with rear type
      borigin_MIP_PR <- c(0,0,0) # No relationship with origin - no origins are up here
      
      # MCN TO ICH OR PRA to ICH to LGR transition
      b0_MIP_IL <- 1
      bflow_MIP_IL <- -0.3 # Negative relationship - higher flow = lower chance of ascending
      btemp_MIP_IL <- 1 # Positive relationship with temperature
      brear_MIP_IL <- c(0,0) # No effect of rear type
      borigin_MIP_IL <- c(0,0,0) # Negative for JDR, no effect on YAK, positive for Tucannon
      
      # MCN TO ICH OR PRA to Yakima transition
      b0_MIP_YAK <- 1
      bflow_MIP_YAK <- 0 # flow not relevant for tributary entry
      btemp_MIP_YAK <- 0 # Temperature not relevant for tributary entry
      brear_MIP_YAK <- c(0,0) # Wild fish more likely to enter the Yakima
      borigin_MIP_YAK <- c(0,0,0) # Yakima fish have a much higher chance of homing
      
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
      borigin_PR_MIP <- c(0,0,0) # Highly positive for each rear type
      
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
      borigin_IL_MIP <- c(0,0,0) # Positive for JDR and YAK, negative for TUC
      
      # ICH to LGR to Tucannon River transition #
      b0_IL_TUC <- 1
      bflow_IL_TUC <- 0 # No relationship with flow (tributary)
      btemp_IL_TUC <- 0 # No relationship with temperature (tributary)
      brear_IL_TUC <- c(0,0) # No effect of rear type
      borigin_IL_TUC <- c(0,0,0) # Highly positive for TUC to home
      
      
      phi_IL_MIP <- exp(b0_IL_MIP + btemp_IL_MIP*temp_ICH + bflow_IL_MIP*flow_ICH + brear_IL_MIP[rear] + borigin_IL_MIP[origin])
      phi_IL_TUC <- exp(b0_IL_TUC + btemp_IL_TUC*temp_TUC + bflow_IL_TUC*flow_TUC + brear_IL_TUC[rear] + borigin_IL_TUC[origin])
      
      transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- phi_IL_MIP/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- phi_IL_TUC/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "loss"] <- 1 - sum(transition_matrix["mainstem, ICH to LGR", c("mainstem, MCN to ICH or PRA", "Tucannon River")])
      
      # Check values
      # transition_matrix["mainstem, ICH to LGR",]
      
      
      ##### FROM DESCHUTES RIVER STATE #####
      # Deschutes River to mainstem, BON to MCN transition
      b0_DES_BM <- 1 # Make all of the return intercepts 1, instead of 0 - increase chance of loss param
      bflow_DES_BM <- 0 # No relationship with flow
      btemp_DES_BM <- 0 # No relationship with temperature
      brear_DES_BM <- c(0,0) # Hatchery fish less likely to return to mainstem
      borigin_DES_BM <- c(0,0,0) # No relationship with rear type
      
      phi_DES_BM <- exp(b0_DES_BM + btemp_DES_BM*temp_DES + bflow_DES_BM*flow_DES + brear_DES_BM[rear] + borigin_DES_BM[origin])
      
      transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- phi_DES_BM/(1 + phi_DES_BM)
      transition_matrix["Deschutes River", "loss"] <- 1 - transition_matrix["Deschutes River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["Deschutes River", ]
      
      
      ##### FROM JOHN DAY RIVER STATE #####
      # JDR to mainstem, BON to MCN transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # JDR fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_JDR + bflow_PR_MIP*flow_JDR + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["John Day River", "mainstem, BON to MCN"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["John Day River", "loss"] <- 1 - transition_matrix["John Day River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["John Day River", ]
      
      
      ##### FROM YAKIMA RIVER STATE #####
      # Yakima River to mainstem, MCN to ICH or PRA transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # Yakima R. fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_YAK + bflow_PR_MIP*flow_YAK + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["Yakima River", "loss"] <- 1 - transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["Yakima River", ]
      
      
      ##### FROM TUCANNON RIVER STATE #####
      # From Tucannon River to ICH to LGR transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # Tucannon River fish less likely to return
      
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
  # Trim to maximum number of observations
  n.obs <- unlist(lapply(det_hist_sim, ncol))
  max_obs <- max(n.obs)
  det_hist_sim_array <- movement_array[,1:max_obs,]
  
  # As an array
  transition_date_sim_matrix <- transition_date[,1:max_obs]
  
  date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2016-12-31"))
  transition_date_sim_matrix %>% 
    as_tibble() %>% 
    mutate_all(date_numeric) -> transition_date_sim_numeric
  
  # Return a list of two things: [[1]] = the detection history, [[2]] = the dates
  
  return(list(det_hist_sim_array, transition_date_sim_numeric))
}

fish_sim_cat_data_600 <- data.frame(tag_code = seq(1, 600, 1),
                                    natal_origin = c(rep(1, 200),
                                                     rep(2, 200),
                                                     rep(3, 200)),
                                    rear_type = rep(c(rep(1, 100), rep(2, 100)), 3))

fish_sim_cat_data_1200 <- data.frame(tag_code = seq(1, 1200, 1),
                                     natal_origin = c(rep(1, 400),
                                                      rep(2, 400),
                                                      rep(3, 400)),
                                     rear_type = rep(c(rep(1, 200), rep(2, 200)), 3))

fish_sim_cat_data_3000 <- data.frame(tag_code = seq(1, 3000, 1),
                                     natal_origin = c(rep(1, 1000),
                                                      rep(2, 1000),
                                                      rep(3, 1000)),
                                     rear_type = rep(c(rep(1, 500), rep(2, 500)), 3))


##### Loop the simulation

### 600 fish
sim_600_hist_list <- list()
sim_600_dates_list <- list()
nfish = 600
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- cov_continuous_sim(nfish = 600, fish_sim_cat_data = fish_sim_cat_data_600)
  
  sim_600_hist_list[[z]] <- sim_run[[1]]
  sim_600_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_600_hist_list, here::here("simulation", "sim_600_cov_continuous_hist_list.rds"))

saveRDS(sim_600_dates_list, here::here("simulation", "sim_600_cov_continuous_dates_list.rds"))
write.csv(fish_sim_cat_data_600, here::here("simulation", "origin_rear_600.csv"), row.names = FALSE)


### 1200 fish
sim_1200_hist_list <- list()
sim_1200_dates_list <- list()
nfish = 1200
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

# Run for loops in parallel
temp_flow_sim_runs <- foreach(z = 1:10) %dopar% {
  cov_continuous_sim(nfish = 1200, fish_sim_cat_data = fish_sim_cat_data_1200,
               movement_array = movement_array, transition_matrix = transition_matrix)
}

for (i in 1:10){
  sim_1200_hist_list[[i]] <- temp_flow_sim_runs[[i]][[1]]
  sim_1200_dates_list[[i]] <- temp_flow_sim_runs[[i]][[2]]
}

# Export these lists as RDS objects
saveRDS(sim_1200_hist_list, here::here("simulation", "simulated_data", "temperature_flow", "sim_1200_temp_flow_hist_list.rds"))
saveRDS(sim_1200_dates_list, here::here("simulation", "simulated_data", "temperature_flow", "sim_1200_temp_flow_dates_list.rds"))


### 3000 fish
sim_3000_hist_list <- list()
sim_3000_dates_list <- list()
nfish = 3000
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- cov_sim(nfish = 3000, fish_sim_cat_data = fish_sim_cat_data_3000)
  
  sim_3000_hist_list[[z]] <- sim_run[[1]]
  sim_3000_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_3000_hist_list, here::here("simulation", "sim_3000_cov_hist_list.rds"))
saveRDS(sim_3000_dates_list, here::here("simulation", "sim_3000_cov_dates_list.rds"))
write.csv(fish_sim_cat_data_3000, here::here("simulation", "origin_rear_3000.csv"), row.names = FALSE)





##### SIMULATION FUNCTION 4: ONLY TEMPERATURE #####
# This is a function to simulate data with temperature (no rear or origin or flow), with every b0 set to 1.
# Right now the temperature and flow values are the same due to the simulation. That creates an identifiability issue
# where we can't include both.
# the only argument needed is "nfish"
# I think that the reason we're having some issues with estimation is because b0 and brear borigin all function as intercepts.

cov_temp_sim <- function(nfish, fish_sim_cat_data, movement_array = movement_array, transition_matrix = transition_matrix){
  # Take a uniform distribution between July 1 and September 30
  BON_arrival_dates_sim <- seq(ymd("2017-07-01"), ymd("2019-09-30"), by = "days")[round(runif(3000, min = 1, max = 92), 0)]
  
  state_date <- matrix(nrow = nfish, ncol = 100)
  transition_date <- matrix(nrow = nfish, ncol = 100)
  
  for (i in 1:nfish){ # for each fish
    print(paste0("i = ", i))
    # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
    movement_array["mainstem, BON to MCN", 1, i] <- 1
    # Populate state date matrix with date at arrival at BON
    state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
    # Store this in the transition date matrix
    transition_date[i,1] <- state_date[i,1]
    
    # populate the rest of the state dates for this fish by adding two weeks
    for (k in 2:100){
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
      # troubleshooting
      # print(paste0(state_date[i, j-1]))
      # Get covariates
      # index the covariate data
      temps <- temp_sim_zscore_df[temp_sim_zscore_df$date == state_date[i,j-1],][2:ncol(temp_sim_zscore_df)]
      temp_BON <- temps[1,1]; temp_MCN <- temps[1,2]; temp_PRA <- temps[1,3]; temp_ICH <- temps[1,4]; temp_JDR <- temps[1,5]; temp_DES <- temps[1,6]; temp_YAK <- temps[1,7]; temp_TUC <- temps[1,8]
      flows <- flow_sim_zscore_df[flow_sim_zscore_df$date == state_date[i,j-1],][2:ncol(flow_sim_zscore_df)]
      flow_BON <- flows[1,1]; flow_MCN <- flows[1,2]; flow_PRA <- flows[1,3]; flow_ICH <- flows[1,4]; flow_JDR <- flows[1,5]; flow_DES <- flows[1,6]; flow_YAK <- flows[1,7]; flow_TUC <- flows[1,8]
      
      
      # The function:
      ##### FROM MOUTH TO BON STATE #####
      
      # Mouth to BON to BON to MCN transition
      b0_MB_BM <- 1
      bflow_MB_BM <- 0 # No relationship with flow
      btemp_MB_BM <- 0 # No relationship with temperature
      brear_MB_BM <- c(0,0) # No effect of rear type
      borigin_MB_BM <- c(0,0,0) # Highly positive for each rear type
      
      phi_MB_BM <- exp(b0_MB_BM + btemp_MB_BM*temp_BON + bflow_MB_BM*flow_BON + brear_MB_BM[rear] + borigin_MB_BM[origin])
      
      transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- phi_MB_BM/(1 + phi_MB_BM)
      transition_matrix["mainstem, mouth to BON", "loss"] <- 1 - transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"]
      
      # Check values
      transition_matrix["mainstem, mouth to BON", ]
      
      ##### FROM BON TO MCN STATE #####
      
      # BON to MCN to Mouth to BON transition
      b0_BM_MB <- 1
      bflow_BM_MB <- 0 # Make this vary with flow, so higher flow at BON = higher fallback over BON
      btemp_BM_MB <- 0 # No relationship with temperature
      brear_BM_MB <- c(0,0) # c(hatchery, wild) # No relationship with rear type
      borigin_BM_MB <- c(0, 0, 0) # c(John Day River, Yakima River, Tucannon River) # No relationship with origin
      
      # BON to MCN to MCN to ICH or PRA transition
      b0_BM_MIP <- 1
      bflow_BM_MIP <- 0 # No relationship with flows
      btemp_BM_MIP <- 0.5 # When it's hotter, more likely to go upstream
      brear_BM_MIP <- c(0,0) # No relationship with rear type
      borigin_BM_MIP <- c(0, 0, 0) # JDR fish have lower probability of overshooting MCN than the fish from origins upstream of MCN
      
      # BON to MCN to DES R transition
      b0_BM_DES <- 1
      bflow_BM_DES <- 0 # flow not relevant for tributary entry
      btemp_BM_DES <- 0 # Temperature not relevant for tributary entry
      brear_BM_DES <- c(0,0) # higher probability of straying to Deschutes if hatchery
      borigin_BM_DES <- c(0, 0, 0) # JDR has higher probability of stray to DES because it's closer
      
      # BON to MCN to JDR transition
      b0_BM_JDR <- 1
      bflow_BM_JDR <- 0 # flow not relevant for tributary entry
      btemp_BM_JDR <- 0 # Temperature not relevant for tributary entry
      brear_BM_JDR <- c(0,0) # no effect of rear type
      borigin_BM_JDR <- c(0, 0, 0) # JDR fish have much higher chance of homing
      
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
      bflow_MIP_BM <- 0 # Make this vary with flow, so higher flow at MCN = higher fallback - but lower
      btemp_MIP_BM <- -0.5 # Negative relationship with temperature - when it's colder, they fall back more
      brear_MIP_BM <- c(0,0) # c(hatchery, wild) # No relationship with rear type
      borigin_MIP_BM <- c(0, 0, 0) # c(John Day River, Yakima River, Tucannon River) # JDR fish have much higher chance of falling back
      
      # BON to MCN to PRA to RIS transition
      b0_MIP_PR <- 1
      bflow_MIP_PR <- 0 # No relationship with flow
      btemp_MIP_PR <- 0.3 # When it's hotter, more likely to go upstream
      brear_MIP_PR <- c(0,0) # No relationship with rear type
      borigin_MIP_PR <- c(0,0,0) # No relationship with origin - no origins are up here
      
      # BON to MCN to ICH to LGR transition
      b0_MIP_IL <- 1
      bflow_MIP_IL <- 0 # Negative relationship - higher flow = lower chance of ascending
      btemp_MIP_IL <- 1 # Positive relationship with temperature
      brear_MIP_IL <- c(0,0) # No effect of rear type
      borigin_MIP_IL <- c(0,0,0) # Negative for JDR, no effect on YAK, positive for Tucannon
      
      # BON to MCN to Yakima transition
      b0_MIP_YAK <- 1
      bflow_MIP_YAK <- 0 # flow not relevant for tributary entry
      btemp_MIP_YAK <- 0 # Temperature not relevant for tributary entry
      brear_MIP_YAK <- c(0,0) # Wild fish more likely to enter the Yakima
      borigin_MIP_YAK <- c(0,0,0) # Yakima fish have a much higher chance of homing
      
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
      borigin_PR_MIP <- c(0,0,0) # Highly positive for each rear type
      
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
      borigin_IL_MIP <- c(0,0,0) # Positive for JDR and YAK, negative for TUC
      
      # ICH to LGR to Tucannon River transition #
      b0_IL_TUC <- 1
      bflow_IL_TUC <- 0 # No relationship with flow (tributary)
      btemp_IL_TUC <- 0 # No relationship with temperature (tributary)
      brear_IL_TUC <- c(0,0) # No effect of rear type
      borigin_IL_TUC <- c(0,0,0) # Highly positive for TUC to home
      
      
      phi_IL_MIP <- exp(b0_IL_MIP + btemp_IL_MIP*temp_ICH + bflow_IL_MIP*flow_ICH + brear_IL_MIP[rear] + borigin_IL_MIP[origin])
      phi_IL_TUC <- exp(b0_IL_TUC + btemp_IL_TUC*temp_TUC + bflow_IL_TUC*flow_TUC + brear_IL_TUC[rear] + borigin_IL_TUC[origin])
      
      transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- phi_IL_MIP/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- phi_IL_TUC/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "loss"] <- 1 - sum(transition_matrix["mainstem, ICH to LGR", c("mainstem, MCN to ICH or PRA", "Tucannon River")])
      
      # Check values
      # transition_matrix["mainstem, ICH to LGR",]
      
      
      ##### FROM DESCHUTES RIVER STATE #####
      # Deschutes River to mainstem, BON to MCN transition
      b0_DES_BM <- 1 # Make all of the return intercepts 1, instead of 0 - increase chance of loss param
      bflow_DES_BM <- 0 # No relationship with flow
      btemp_DES_BM <- 0 # No relationship with temperature
      brear_DES_BM <- c(0,0) # Hatchery fish less likely to return to mainstem
      borigin_DES_BM <- c(0,0,0) # No relationship with rear type
      
      phi_DES_BM <- exp(b0_DES_BM + btemp_DES_BM*temp_DES + bflow_DES_BM*flow_DES + brear_DES_BM[rear] + borigin_DES_BM[origin])
      
      transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- phi_DES_BM/(1 + phi_DES_BM)
      transition_matrix["Deschutes River", "loss"] <- 1 - transition_matrix["Deschutes River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["Deschutes River", ]
      
      
      ##### FROM JOHN DAY RIVER STATE #####
      # JDR to mainstem, BON to MCN transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # JDR fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_JDR + bflow_PR_MIP*flow_JDR + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["John Day River", "mainstem, BON to MCN"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["John Day River", "loss"] <- 1 - transition_matrix["John Day River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["John Day River", ]
      
      
      ##### FROM YAKIMA RIVER STATE #####
      # Yakima River to mainstem, MCN to ICH or PRA transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # Yakima R. fish less likely to return
      
      phi_PR_MIP <- exp(b0_PR_MIP + btemp_PR_MIP*temp_YAK + bflow_PR_MIP*flow_YAK + brear_PR_MIP[rear] + borigin_PR_MIP[origin])
      
      transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["Yakima River", "loss"] <- 1 - transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["Yakima River", ]
      
      
      ##### FROM TUCANNON RIVER STATE #####
      # From Tucannon River to ICH to LGR transition
      b0_PR_MIP <- 1
      bflow_PR_MIP <- 0 # No relationship with flow
      btemp_PR_MIP <- 0 # No relationship with temperature
      brear_PR_MIP <- c(0,0) # No effect of rear type
      borigin_PR_MIP <- c(0, 0, 0) # Tucannon River fish less likely to return
      
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
  # Trim to maximum number of observations
  n.obs <- unlist(lapply(det_hist_sim, ncol))
  max_obs <- max(n.obs)
  det_hist_sim_array <- movement_array[,1:max_obs,]
  
  # As an array
  transition_date_sim_matrix <- transition_date[,1:max_obs]
  
  date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2016-12-31"))
  transition_date_sim_matrix %>% 
    as_tibble() %>% 
    mutate_all(date_numeric) -> transition_date_sim_numeric
  
  # Return a list of two things: [[1]] = the detection history, [[2]] = the dates
  
  return(list(det_hist_sim_array, transition_date_sim_numeric))
}

fish_sim_cat_data_600 <- data.frame(tag_code = seq(1, 600, 1),
                                    natal_origin = c(rep(1, 200),
                                                     rep(2, 200),
                                                     rep(3, 200)),
                                    rear_type = rep(c(rep(1, 100), rep(2, 100)), 3))

fish_sim_cat_data_1200 <- data.frame(tag_code = seq(1, 1200, 1),
                                     natal_origin = c(rep(1, 400),
                                                      rep(2, 400),
                                                      rep(3, 400)),
                                     rear_type = rep(c(rep(1, 200), rep(2, 200)), 3))

fish_sim_cat_data_3000 <- data.frame(tag_code = seq(1, 3000, 1),
                                     natal_origin = c(rep(1, 1000),
                                                      rep(2, 1000),
                                                      rep(3, 1000)),
                                     rear_type = rep(c(rep(1, 500), rep(2, 500)), 3))


##### Loop the simulation

### 600 fish
sim_600_hist_list <- list()
sim_600_dates_list <- list()
nfish = 600
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:3){
  sim_run <- cov_temp_sim(nfish = 600, fish_sim_cat_data = fish_sim_cat_data_600)
  
  sim_600_hist_list[[z]] <- sim_run[[1]]
  sim_600_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_600_hist_list, here::here("simulation", "sim_600_cov_temp_hist_list.rds"))

saveRDS(sim_600_dates_list, here::here("simulation", "sim_600_cov_temp_dates_list.rds"))
write.csv(fish_sim_cat_data_600, here::here("simulation", "origin_rear_600.csv"), row.names = FALSE)


### 1200 fish
sim_1200_hist_list <- list()
sim_1200_dates_list <- list()
nfish = 1200
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states




# Run for loops in parallel
temp_sim_runs <- foreach(z = 1:10) %dopar% {
  cov_temp_sim(nfish = 1200, fish_sim_cat_data = fish_sim_cat_data_1200,
                      movement_array = movement_array, transition_matrix = transition_matrix)
}

for (i in 1:10){
  sim_1200_hist_list[[i]] <- temp_sim_runs[[i]][[1]]
  sim_1200_dates_list[[i]] <- temp_sim_runs[[i]][[2]]
}

# Export these lists as RDS objects
saveRDS(sim_1200_hist_list, here::here("simulation", "simulated_data", "temperature_only", "sim_1200_cov_temp_hist_list.rds"))
saveRDS(sim_1200_dates_list, here::here("simulation", "simulated_data", "temperature_only", "sim_1200_cov_temp_dates_list.rds"))


### 3000 fish
sim_3000_hist_list <- list()
sim_3000_dates_list <- list()
nfish = 3000
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- cov_sim(nfish = 3000, fish_sim_cat_data = fish_sim_cat_data_3000)
  
  sim_3000_hist_list[[z]] <- sim_run[[1]]
  sim_3000_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_3000_hist_list, here::here("simulation", "sim_3000_cov_hist_list.rds"))
saveRDS(sim_3000_dates_list, here::here("simulation", "sim_3000_cov_dates_list.rds"))
write.csv(fish_sim_cat_data_3000, here::here("simulation", "origin_rear_3000.csv"), row.names = FALSE)






##### SIMULATION FUNCTION 5: ORIGIN AND REAR (CATEGORICAL COV) #####



origin_rear_cov_sim <- function(nfish, fish_sim_cat_data, cat_X_mat, movement_array, transition_matrix){
  # Take a uniform distribution between July 1 and September 30
  BON_arrival_dates_sim <- seq(ymd("2017-07-01"), ymd("2019-09-30"), by = "days")[round(runif(3000, min = 1, max = 92), 0)]
  
  state_date <- matrix(nrow = nfish, ncol = 100)
  transition_date <- matrix(nrow = nfish, ncol = 100)
  
  for (i in 1:nfish){ # for each fish
    print(paste0("i = ", i))
    # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
    movement_array["mainstem, BON to MCN", 1, i] <- 1
    # Populate state date matrix with date at arrival at BON
    state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
    # Store this in the transition date matrix
    transition_date[i,1] <- state_date[i,1]
    
    # populate the rest of the state dates for this fish by adding two weeks
    for (k in 2:100){
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
      ##### FROM MOUTH TO BON STATE #####
      
      # Mouth to BON to BON to MCN transition
      b0_MB_BM <- 1
      bflow_MB_BM <- 0 # No relationship with flow
      btemp_MB_BM <- 0 # No relationship with temperature
      brear_MB_BM <- 0 # No difference for hatchery vs. wild
      borigin_MB_BM <- c(0, 0) # No difference for JDR vs. TUC or YAK vs TUC
      
      # Here we've reformatted it so that b0, brear, and borigin are all in a vector that we multiply by the design matrix for the categorical covariates
      phi_MB_BM <- exp(btemp_MB_BM*temp_BON + bflow_MB_BM*flow_BON + cat_X_mat[i,] %*% as.matrix(c(b0_MB_BM, brear_MB_BM, borigin_MB_BM), ncol = 1))
      
      transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- phi_MB_BM/(1 + phi_MB_BM)
      transition_matrix["mainstem, mouth to BON", "loss"] <- 1 - transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"]
      
      # Check values
      transition_matrix["mainstem, mouth to BON", ]
      
      ##### FROM BON TO MCN STATE #####
      
      # BON to MCN to Mouth to BON transition
      b0_BM_MB <- 1
      bflow_BM_MB <- 0 # No effect of flow
      btemp_BM_MB <- 0 # No relationship with temperature
      brear_BM_MB <- 0 # No difference for hatchery vs. wild
      borigin_BM_MB <- c(0, 0) # No difference for JDR vs. TUC or YAK vs TUC
      
      # BON to MCN to MCN to ICH or PRA transition
      b0_BM_MIP <- 1
      bflow_BM_MIP <- 0 # No relationship with flows
      btemp_BM_MIP <- 0 # No effect of temperature
      brear_BM_MIP <- 0 # No difference for hatchery vs. wild
      borigin_BM_MIP <- c(-0.5, 0.25) # JDR have a lower probability of overshooting MCN than TUC fish, but there is no difference for YAK vs TUC
      
      # BON to MCN to DES R transition
      b0_BM_DES <- 1
      bflow_BM_DES <- 0 # flow not relevant for tributary entry
      btemp_BM_DES <- 0 # Temperature not relevant for tributary entry
      brear_BM_DES <- 0.5 # higher probability of straying to Deschutes if hatchery
      borigin_BM_DES <- c(0.25, -0.125) # JDR have a higher probability of straying to DES than TUC fish (closer), but there is no difference for YAK vs TUC
      
      # BON to MCN to JDR transition
      b0_BM_JDR <- 1
      bflow_BM_JDR <- 0 # flow not relevant for tributary entry
      btemp_BM_JDR <- 0 # Temperature not relevant for tributary entry
      brear_BM_JDR <- 0 # No difference for hatchery vs. wild
      borigin_BM_JDR <- c(0.5, -0.25) # JDR have a much higher probability of homing to JDR than TUC fish , but there is no difference for YAK vs TUC
      
      # Evaluate
      
      # BON to MCN
      phi_BM_MB <- exp(btemp_BM_MB*temp_BON + bflow_BM_MB*flow_BON + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_MB, brear_BM_MB, borigin_BM_MB), ncol = 1))
      phi_BM_MIP <- exp(btemp_BM_MIP*temp_MCN + bflow_BM_MIP*flow_MCN + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_MIP, brear_BM_MIP, borigin_BM_MIP), ncol = 1))
      phi_BM_DES <- exp(btemp_BM_DES*temp_DES + bflow_BM_DES*flow_DES + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_DES, brear_BM_DES, borigin_BM_DES), ncol = 1))
      phi_BM_JDR <- exp(btemp_BM_JDR*temp_JDR + bflow_BM_JDR*flow_JDR + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_JDR, brear_BM_JDR, borigin_BM_JDR), ncol = 1))
      
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
      bflow_MIP_BM <- 0 # No relationship with flow
      btemp_MIP_BM <- 0 # No relationship with temperature
      brear_MIP_BM <- 0 # No difference for hatchery vs. wild
      borigin_MIP_BM <- c(0.5, -0.25) # JDR have much higher chance of falling back than TUC, but no difference between YAK and TUC
      
      # BON to MCN to PRA to RIS transition
      b0_MIP_PR <- 1
      bflow_MIP_PR <- 0 # No relationship with flow
      btemp_MIP_PR <- 0 # 0.3 # When it's hotter, more likely to go upstream
      brear_MIP_PR <- 0 # No difference for hatchery vs. wild
      borigin_MIP_PR <- c(0,0) # No relationship with origin - no origins are up here
      
      # BON to MCN to ICH to LGR transition
      b0_MIP_IL <- 1
      bflow_MIP_IL <- 0 # -0.3 # Negative relationship - higher flow = lower chance of ascending
      btemp_MIP_IL <- 0 # 0.2 # Positive relationship with temperature
      brear_MIP_IL <- c(0) # No effect of rear type
      borigin_MIP_IL <- c(-0.5, 0) # JDR have much lower chance than TUC, YAK have slightly lower chance than TUC
      
      # BON to MCN to Yakima transition
      b0_MIP_YAK <- 1
      bflow_MIP_YAK <- 0 # flow not relevant for tributary entry
      btemp_MIP_YAK <- 0 # Temperature not relevant for tributary entry
      brear_MIP_YAK <- -0.2 # Hatchery fish less likely to enter Yakima
      borigin_MIP_YAK <- c(-0.25, 0.5) # No difference for JDR vs. TUC, but much higher for YAK than JDR
      
      phi_MIP_BM <- exp(btemp_MIP_BM*temp_MCN + bflow_MIP_BM*flow_MCN + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_BM, brear_MIP_BM, borigin_MIP_BM), ncol = 1))
      phi_MIP_PR <- exp(btemp_MIP_PR*temp_PRA + bflow_MIP_PR*flow_PRA + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_PR, brear_MIP_PR, borigin_MIP_PR), ncol = 1))
      phi_MIP_IL <- exp(btemp_MIP_IL*temp_ICH + bflow_MIP_IL*flow_ICH + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_IL, brear_MIP_IL, borigin_MIP_IL), ncol = 1))
      phi_MIP_YAK <- exp(btemp_MIP_YAK*temp_YAK + bflow_MIP_YAK*flow_YAK + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_YAK, brear_MIP_YAK, borigin_MIP_YAK), ncol = 1))
      
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
      brear_PR_MIP <- 0 # No effect of rear type
      borigin_PR_MIP <- c(0,0) # No effect of origins vs. TUC 
      
      phi_PR_MIP <- exp(btemp_PR_MIP*temp_PRA + bflow_PR_MIP*flow_PRA + cat_X_mat[i,] %*% as.matrix(c(b0_PR_MIP, brear_PR_MIP, borigin_PR_MIP), ncol = 1))
      
      transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["mainstem, PRA to RIS", "loss"] <- 1 - transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["mainstem, PRA to RIS", ]
      
      
      
      
      ##### FROM ICH TO LGR STATE #####
      
      # ICH to LGR to MCN to ICH or PRA transition #
      b0_IL_MIP <- 1
      bflow_IL_MIP <- 0 # No relationship with flow
      btemp_IL_MIP <- 0 # 0.2 # Negative relationship with temperature - when it's colder, tend to go downstream more
      brear_IL_MIP <- 0 # No effect of rear type
      borigin_IL_MIP <- c(0.25, 0.25) # JDR and YAK fish are both more likely than TUC fish to return to MCN to ICH or PRA state
      
      # ICH to LGR to Tucannon River transition #
      b0_IL_TUC <- 1
      bflow_IL_TUC <- 0 # No relationship with flow (tributary)
      btemp_IL_TUC <- 0 # No relationship with temperature (tributary)
      brear_IL_TUC <- 0 # No effect of rear type
      borigin_IL_TUC <- c(-0.25, -0.25) # Both other origins are less likely to enter TUC than TUC fish
      
      
      phi_IL_MIP <- exp(btemp_IL_MIP*temp_ICH + bflow_IL_MIP*flow_ICH + cat_X_mat[i,] %*% as.matrix(c(b0_IL_MIP, brear_IL_MIP, borigin_IL_MIP), ncol = 1))
      phi_IL_TUC <- exp(btemp_IL_TUC*temp_TUC + bflow_IL_TUC*flow_TUC + cat_X_mat[i,] %*% as.matrix(c(b0_IL_TUC, brear_IL_TUC, borigin_IL_TUC), ncol = 1))
      
      transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- phi_IL_MIP/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- phi_IL_TUC/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "loss"] <- 1 - sum(transition_matrix["mainstem, ICH to LGR", c("mainstem, MCN to ICH or PRA", "Tucannon River")])
      
      # Check values
      # transition_matrix["mainstem, ICH to LGR",]
      
      
      ##### FROM DESCHUTES RIVER STATE #####
      # Deschutes River to mainstem, BON to MCN transition
      b0_DES_BM <- 1 # Make all of the return intercepts 1, instead of 0 - increase chance of loss param
      bflow_DES_BM <- 0 # No relationship with flow
      btemp_DES_BM <- 0 # No relationship with temperature
      brear_DES_BM <- -0.5 # Hatchery fish less likely to return to mainstem
      borigin_DES_BM <- c(0,0) # No relationship with rear type
      
      phi_DES_BM <- exp(btemp_DES_BM*temp_DES + bflow_DES_BM*flow_DES + cat_X_mat[i,] %*% as.matrix(c(b0_DES_BM, brear_DES_BM, borigin_DES_BM), ncol = 1))
      
      transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- phi_DES_BM/(1 + phi_DES_BM)
      transition_matrix["Deschutes River", "loss"] <- 1 - transition_matrix["Deschutes River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["Deschutes River", ]
      
      
      ##### FROM JOHN DAY RIVER STATE #####
      # JDR to mainstem, BON to MCN transition
      b0_JDR_BM <- 1
      bflow_JDR_BM <- 0 # No relationship with flow
      btemp_JDR_BM <- 0 # No relationship with temperature
      brear_JDR_BM <- 0 # No effect of rear type
      borigin_JDR_BM <- c(-0.5, 0.25) # JDR fish less likely to return than TUC, no difference for YAK vs. TUC
      
      phi_JDR_BM <- exp(btemp_JDR_BM*temp_JDR + bflow_JDR_BM*flow_JDR + cat_X_mat[i,] %*% as.matrix(c(b0_JDR_BM, brear_JDR_BM, borigin_JDR_BM), ncol = 1))
      
      transition_matrix["John Day River", "mainstem, BON to MCN"] <- phi_JDR_BM/(1 + phi_JDR_BM)
      transition_matrix["John Day River", "loss"] <- 1 - transition_matrix["John Day River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["John Day River", ]
      
      
      ##### FROM YAKIMA RIVER STATE #####
      # Yakima River to mainstem, MCN to ICH or PRA transition
      b0_YAK_MIP <- 1
      bflow_YAK_MIP <- 0 # No relationship with flow
      btemp_YAK_MIP <- 0 # No relationship with temperature
      brear_YAK_MIP <- 0 # No effect of rear type
      borigin_YAK_MIP <- c(0.25, -0.5) # Yakima R. fish less likely to return
      
      phi_YAK_MIP <- exp(btemp_YAK_MIP*temp_YAK + bflow_YAK_MIP*flow_YAK + cat_X_mat[i,] %*% as.matrix(c(b0_YAK_MIP, brear_YAK_MIP, borigin_YAK_MIP), ncol = 1))
      
      transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- phi_YAK_MIP/(1 + phi_YAK_MIP)
      transition_matrix["Yakima River", "loss"] <- 1 - transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["Yakima River", ]
      
      
      ##### FROM TUCANNON RIVER STATE #####
      # From Tucannon River to ICH to LGR transition
      b0_TUC_IL <- 1
      bflow_TUC_IL <- 0 # No relationship with flow
      btemp_TUC_IL <- 0 # No relationship with temperature
      brear_TUC_IL <- 0 # No effect of rear type
      borigin_TUC_IL <- c(0.25,0.25) # JDR and YAK fish both more likely to return than TUC fish
      
      phi_TUC_IL <- exp(btemp_TUC_IL*temp_TUC + bflow_TUC_IL*flow_TUC + cat_X_mat[i,] %*% as.matrix(c(b0_TUC_IL, brear_TUC_IL, borigin_TUC_IL), ncol = 1))
      
      transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- phi_TUC_IL/(1 + phi_TUC_IL)
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
  # Trim to maximum number of observations
  n.obs <- unlist(lapply(det_hist_sim, ncol))
  max_obs <- max(n.obs)
  det_hist_sim_array <- movement_array[,1:max_obs,]
  
  # As an array
  transition_date_sim_matrix <- transition_date[,1:max_obs]
  
  date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2016-12-31"))
  transition_date_sim_matrix %>% 
    as_tibble() %>% 
    mutate_all(date_numeric) -> transition_date_sim_numeric
  
  # Return a list of two things: [[1]] = the detection history, [[2]] = the dates
  
  return(list(det_hist_sim_array, transition_date_sim_numeric))
}

fish_sim_cat_data_600 <- data.frame(tag_code = seq(1, 600, 1),
                                    natal_origin = c(rep(1, 200),
                                                     rep(2, 200),
                                                     rep(3, 200)),
                                    rear_type = rep(c(rep(1, 100), rep(2, 100)), 3))

fish_sim_cat_data_1200 <- data.frame(tag_code = seq(1, 1200, 1),
                                     natal_origin = c(rep(1, 400),
                                                      rep(2, 400),
                                                      rep(3, 400)),
                                     rear_type = rep(c(rep(1, 200), rep(2, 200)), 3))

fish_sim_cat_data_3000 <- data.frame(tag_code = seq(1, 3000, 1),
                                     natal_origin = c(rep(1, 1000),
                                                      rep(2, 1000),
                                                      rep(3, 1000)),
                                     rear_type = rep(c(rep(1, 500), rep(2, 500)), 3))


##### Loop the simulation

### 600 fish
sim_600_hist_list <- list()
sim_600_dates_list <- list()
nfish = 600
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

# Create the design matrix for categorical variables
cat_X_mat_600 <- matrix(0, nrow = nfish, ncol = 4)
# The first column everyone gets a 1 (this is beta 0/grand mean mu)
cat_X_mat_600[,1] <- 1

for (i in 1:nfish){
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


# Run for loops in parallel
origin_rear_sim_runs <- foreach(z = 1:10) %dopar% {
  origin_rear_cov_sim(nfish = 600, fish_sim_cat_data = fish_sim_cat_data_600,cat_X_mat = cat_X_mat_600, 
                 movement_array = movement_array, transition_matrix = transition_matrix)
}

for (i in 1:10){
  sim_600_hist_list[[i]] <- origin_rear_sim_runs[[i]][[1]]
  sim_600_dates_list[[i]] <- origin_rear_sim_runs[[i]][[2]]
}

# Export these lists as RDS objects
saveRDS(sim_600_hist_list, here::here("simulation", "sim_600_origin_rear_hist_list.rds"))

saveRDS(sim_600_dates_list, here::here("simulation", "sim_600_origin_rear_dates_list.rds"))
write.csv(fish_sim_cat_data_600, here::here("simulation", "origin_rear_600.csv"), row.names = FALSE)


### 1200 fish
fish_sim_cat_data_1200 <- data.frame(tag_code = seq(1, 1200, 1),
                                     natal_origin = c(rep(1, 400),
                                                      rep(2, 400),
                                                      rep(3, 400)),
                                     rear_type = rep(c(rep(1, 200), rep(2, 200)), 3))

sim_1200_hist_list <- list()
sim_1200_dates_list <- list()
nfish = 1200
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

# Create the design matrix for categorical variables
cat_X_mat_1200 <- matrix(0, nrow = nfish, ncol = 4)
# The first column everyone gets a 1 (this is beta 0/grand mean mu)
cat_X_mat_1200[,1] <- 1

for (i in 1:nfish){
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



# Run for loops in parallel
origin_rear_sim_runs <- foreach(z = 1:10) %dopar% {
  origin_rear_cov_sim(nfish = 1200, fish_sim_cat_data = fish_sim_cat_data_1200,cat_X_mat = cat_X_mat_1200, 
                      movement_array = movement_array, transition_matrix = transition_matrix)
}

for (i in 1:10){
  sim_1200_hist_list[[i]] <- origin_rear_sim_runs[[i]][[1]]
  sim_1200_dates_list[[i]] <- origin_rear_sim_runs[[i]][[2]]
}

# Export these lists as RDS objects
saveRDS(sim_1200_hist_list, here::here("simulation", "simulated_data", "origin_rear", "sim_1200_origin_rear_hist_list.rds"))

saveRDS(sim_1200_dates_list, here::here("simulation", "simulated_data", "origin_rear", "sim_1200_origin_rear_dates_list.rds"))
write.csv(fish_sim_cat_data_1200, here::here("simulation", "simulated_data", "origin_rear", "origin_rear_1200.csv"), row.names = FALSE)

### 3000 fish
sim_3000_hist_list <- list()
sim_3000_dates_list <- list()
nfish = 3000
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- cov_sim(nfish = 3000, fish_sim_cat_data = fish_sim_cat_data_3000)
  
  sim_3000_hist_list[[z]] <- sim_run[[1]]
  sim_3000_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_3000_hist_list, here::here("simulation", "sim_3000_cov_categorical_hist_list.rds"))
saveRDS(sim_3000_dates_list, here::here("simulation", "sim_3000_cov_categorical_dates_list.rds"))
write.csv(fish_sim_cat_data_3000, here::here("simulation", "origin_rear_3000.csv"), row.names = FALSE)







##### SIMULATION FUNCTION 6: ORIGIN ONLY #####



origin_cov_sim <- function(nfish, fish_sim_cat_data, cat_X_mat, movement_array, transition_matrix){
  # Take a uniform distribution between July 1 and September 30
  BON_arrival_dates_sim <- seq(ymd("2017-07-01"), ymd("2019-09-30"), by = "days")[round(runif(3000, min = 1, max = 92), 0)]
  
  state_date <- matrix(nrow = nfish, ncol = 100)
  transition_date <- matrix(nrow = nfish, ncol = 100)
  
  for (i in 1:nfish){ # for each fish
  # foreach(i = 1:nfish, .combine = "c") %dopar% {
    print(paste0("i = ", i))
    # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
    movement_array["mainstem, BON to MCN", 1, i] <- 1
    # Populate state date matrix with date at arrival at BON
    state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
    # Store this in the transition date matrix
    transition_date[i,1] <- state_date[i,1]
    
    # populate the rest of the state dates for this fish by adding two weeks
    for (k in 2:100){
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
      ##### FROM MOUTH TO BON STATE #####
      
      # Mouth to BON to BON to MCN transition
      b0_MB_BM <- 1
      bflow_MB_BM <- 0 # No relationship with flow
      btemp_MB_BM <- 0 # No relationship with temperature
      brear_MB_BM <- 0 # No difference for hatchery vs. wild
      borigin_MB_BM <- c(0, 0) # No difference for JDR vs. TUC or YAK vs TUC
      
      # Here we've reformatted it so that b0, brear, and borigin are all in a vector that we multiply by the design matrix for the categorical covariates
      phi_MB_BM <- exp(btemp_MB_BM*temp_BON + bflow_MB_BM*flow_BON + cat_X_mat[i,] %*% as.matrix(c(b0_MB_BM, brear_MB_BM, borigin_MB_BM), ncol = 1))
      
      transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- phi_MB_BM/(1 + phi_MB_BM)
      transition_matrix["mainstem, mouth to BON", "loss"] <- 1 - transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"]
      
      # Check values
      transition_matrix["mainstem, mouth to BON", ]
      
      ##### FROM BON TO MCN STATE #####
      
      # BON to MCN to Mouth to BON transition
      b0_BM_MB <- 1
      bflow_BM_MB <- 0 # No effect of flow
      btemp_BM_MB <- 0 # No relationship with temperature
      brear_BM_MB <- 0 # No difference for hatchery vs. wild
      borigin_BM_MB <- c(0, 0) # No difference for JDR vs. TUC or YAK vs TUC
      
      # BON to MCN to MCN to ICH or PRA transition
      b0_BM_MIP <- 1
      bflow_BM_MIP <- 0 # No relationship with flows
      btemp_BM_MIP <- 0 # No effect of temperature
      brear_BM_MIP <- 0 # No difference for hatchery vs. wild
      borigin_BM_MIP <- c(-0.5, 0.25) # JDR have a lower probability of overshooting MCN than TUC fish, but there is no difference for YAK vs TUC
      
      # BON to MCN to DES R transition
      b0_BM_DES <- 1
      bflow_BM_DES <- 0 # flow not relevant for tributary entry
      btemp_BM_DES <- 0 # Temperature not relevant for tributary entry
      brear_BM_DES <- 0 # 0.5 # higher probability of straying to Deschutes if hatchery
      borigin_BM_DES <- c(0.25, -0.125) # JDR have a higher probability of straying to DES than TUC fish (closer), but there is no difference for YAK vs TUC
      
      # BON to MCN to JDR transition
      b0_BM_JDR <- 1
      bflow_BM_JDR <- 0 # flow not relevant for tributary entry
      btemp_BM_JDR <- 0 # Temperature not relevant for tributary entry
      brear_BM_JDR <- 0 # No difference for hatchery vs. wild
      borigin_BM_JDR <- c(0.5, -0.25) # JDR have a much higher probability of homing to JDR than TUC fish , but there is no difference for YAK vs TUC
      
      # Evaluate
      
      # BON to MCN
      phi_BM_MB <- exp(btemp_BM_MB*temp_BON + bflow_BM_MB*flow_BON + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_MB, brear_BM_MB, borigin_BM_MB), ncol = 1))
      phi_BM_MIP <- exp(btemp_BM_MIP*temp_MCN + bflow_BM_MIP*flow_MCN + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_MIP, brear_BM_MIP, borigin_BM_MIP), ncol = 1))
      phi_BM_DES <- exp(btemp_BM_DES*temp_DES + bflow_BM_DES*flow_DES + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_DES, brear_BM_DES, borigin_BM_DES), ncol = 1))
      phi_BM_JDR <- exp(btemp_BM_JDR*temp_JDR + bflow_BM_JDR*flow_JDR + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_JDR, brear_BM_JDR, borigin_BM_JDR), ncol = 1))
      
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
      bflow_MIP_BM <- 0 # No relationship with flow
      btemp_MIP_BM <- 0 # No relationship with temperature
      brear_MIP_BM <- 0 # No difference for hatchery vs. wild
      borigin_MIP_BM <- c(0.5, -0.25) # JDR have much higher chance of falling back than TUC, but no difference between YAK and TUC
      
      # MCN TO ICH OR PRA to RIS transition
      b0_MIP_PR <- 1
      bflow_MIP_PR <- 0 # No relationship with flow
      btemp_MIP_PR <- 0 # 0.3 # When it's hotter, more likely to go upstream
      brear_MIP_PR <- 0 # No difference for hatchery vs. wild
      borigin_MIP_PR <- c(0,0) # No relationship with origin - no origins are up here
      
      # MCN TO ICH OR PRA to ICH to LGR transition
      b0_MIP_IL <- 1
      bflow_MIP_IL <- 0 # -0.3 # Negative relationship - higher flow = lower chance of ascending
      btemp_MIP_IL <- 0 # 0.2 # Positive relationship with temperature
      brear_MIP_IL <- 0 # No effect of rear type
      borigin_MIP_IL <- c(-0.5, 0) # JDR have much lower chance than TUC, YAK have slightly lower chance than TUC
      
      # MCN TO ICH OR PRA to Yakima transition
      b0_MIP_YAK <- 1
      bflow_MIP_YAK <- 0 # flow not relevant for tributary entry
      btemp_MIP_YAK <- 0 # Temperature not relevant for tributary entry
      brear_MIP_YAK <- 0 # -0.2 # Hatchery fish less likely to enter Yakima
      borigin_MIP_YAK <- c(-0.25, 0.5) # No difference for JDR vs. TUC, but much higher for YAK than JDR
      
      phi_MIP_BM <- exp(btemp_MIP_BM*temp_MCN + bflow_MIP_BM*flow_MCN + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_BM, brear_MIP_BM, borigin_MIP_BM), ncol = 1))
      phi_MIP_PR <- exp(btemp_MIP_PR*temp_PRA + bflow_MIP_PR*flow_PRA + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_PR, brear_MIP_PR, borigin_MIP_PR), ncol = 1))
      phi_MIP_IL <- exp(btemp_MIP_IL*temp_ICH + bflow_MIP_IL*flow_ICH + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_IL, brear_MIP_IL, borigin_MIP_IL), ncol = 1))
      phi_MIP_YAK <- exp(btemp_MIP_YAK*temp_YAK + bflow_MIP_YAK*flow_YAK + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_YAK, brear_MIP_YAK, borigin_MIP_YAK), ncol = 1))
      
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
      brear_PR_MIP <- 0 # No effect of rear type
      borigin_PR_MIP <- c(0,0) # No effect of origins vs. TUC 
      
      phi_PR_MIP <- exp(btemp_PR_MIP*temp_PRA + bflow_PR_MIP*flow_PRA + cat_X_mat[i,] %*% as.matrix(c(b0_PR_MIP, brear_PR_MIP, borigin_PR_MIP), ncol = 1))
      
      transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["mainstem, PRA to RIS", "loss"] <- 1 - transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["mainstem, PRA to RIS", ]
      
      
      
      
      ##### FROM ICH TO LGR STATE #####
      
      # ICH to LGR to MCN to ICH or PRA transition #
      b0_IL_MIP <- 1
      bflow_IL_MIP <- 0 # No relationship with flow
      btemp_IL_MIP <- 0 # 0.2 # Negative relationship with temperature - when it's colder, tend to go downstream more
      brear_IL_MIP <- 0 # No effect of rear type
      borigin_IL_MIP <- c(0.25, 0.25) # JDR and YAK fish are both more likely than TUC fish to return to MCN to ICH or PRA state
      
      # ICH to LGR to Tucannon River transition #
      b0_IL_TUC <- 1
      bflow_IL_TUC <- 0 # No relationship with flow (tributary)
      btemp_IL_TUC <- 0 # No relationship with temperature (tributary)
      brear_IL_TUC <- 0 # No effect of rear type
      borigin_IL_TUC <- c(-0.25, -0.25) # Both other origins are less likely to enter TUC than TUC fish
      
      
      phi_IL_MIP <- exp(btemp_IL_MIP*temp_ICH + bflow_IL_MIP*flow_ICH + cat_X_mat[i,] %*% as.matrix(c(b0_IL_MIP, brear_IL_MIP, borigin_IL_MIP), ncol = 1))
      phi_IL_TUC <- exp(btemp_IL_TUC*temp_TUC + bflow_IL_TUC*flow_TUC + cat_X_mat[i,] %*% as.matrix(c(b0_IL_TUC, brear_IL_TUC, borigin_IL_TUC), ncol = 1))
      
      transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- phi_IL_MIP/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- phi_IL_TUC/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "loss"] <- 1 - sum(transition_matrix["mainstem, ICH to LGR", c("mainstem, MCN to ICH or PRA", "Tucannon River")])
      
      # Check values
      # transition_matrix["mainstem, ICH to LGR",]
      
      
      ##### FROM DESCHUTES RIVER STATE #####
      # Deschutes River to mainstem, BON to MCN transition
      b0_DES_BM <- 1 # Make all of the return intercepts 1, instead of 0 - increase chance of loss param
      bflow_DES_BM <- 0 # No relationship with flow
      btemp_DES_BM <- 0 # No relationship with temperature
      brear_DES_BM <- 0 # -0.5 # Hatchery fish less likely to return to mainstem
      borigin_DES_BM <- c(0,0) # No relationship with rear type
      
      phi_DES_BM <- exp(btemp_DES_BM*temp_DES + bflow_DES_BM*flow_DES + cat_X_mat[i,] %*% as.matrix(c(b0_DES_BM, brear_DES_BM, borigin_DES_BM), ncol = 1))
      
      transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- phi_DES_BM/(1 + phi_DES_BM)
      transition_matrix["Deschutes River", "loss"] <- 1 - transition_matrix["Deschutes River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["Deschutes River", ]
      
      
      ##### FROM JOHN DAY RIVER STATE #####
      # JDR to mainstem, BON to MCN transition
      b0_JDR_BM <- 1
      bflow_JDR_BM <- 0 # No relationship with flow
      btemp_JDR_BM <- 0 # No relationship with temperature
      brear_JDR_BM <- 0 # No effect of rear type
      borigin_JDR_BM <- c(-0.5, 0.25) # JDR fish less likely to return than TUC, no difference for YAK vs. TUC
      
      phi_JDR_BM <- exp(btemp_JDR_BM*temp_JDR + bflow_JDR_BM*flow_JDR + cat_X_mat[i,] %*% as.matrix(c(b0_JDR_BM, brear_JDR_BM, borigin_JDR_BM), ncol = 1))
      
      transition_matrix["John Day River", "mainstem, BON to MCN"] <- phi_JDR_BM/(1 + phi_JDR_BM)
      transition_matrix["John Day River", "loss"] <- 1 - transition_matrix["John Day River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["John Day River", ]
      
      
      ##### FROM YAKIMA RIVER STATE #####
      # Yakima River to mainstem, MCN to ICH or PRA transition
      b0_YAK_MIP <- 1
      bflow_YAK_MIP <- 0 # No relationship with flow
      btemp_YAK_MIP <- 0 # No relationship with temperature
      brear_YAK_MIP <- 0 # No effect of rear type
      borigin_YAK_MIP <- c(0.25, -0.5) # Yakima R. fish less likely to return
      
      phi_YAK_MIP <- exp(btemp_YAK_MIP*temp_YAK + bflow_YAK_MIP*flow_YAK + cat_X_mat[i,] %*% as.matrix(c(b0_YAK_MIP, brear_YAK_MIP, borigin_YAK_MIP), ncol = 1))
      
      transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- phi_YAK_MIP/(1 + phi_YAK_MIP)
      transition_matrix["Yakima River", "loss"] <- 1 - transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["Yakima River", ]
      
      
      ##### FROM TUCANNON RIVER STATE #####
      # From Tucannon River to ICH to LGR transition
      b0_TUC_IL <- 1
      bflow_TUC_IL <- 0 # No relationship with flow
      btemp_TUC_IL <- 0 # No relationship with temperature
      brear_TUC_IL <- 0 # No effect of rear type
      borigin_TUC_IL <- c(0.25,0.25) # JDR and YAK fish both more likely to return than TUC fish
      
      phi_TUC_IL <- exp(btemp_TUC_IL*temp_TUC + bflow_TUC_IL*flow_TUC + cat_X_mat[i,] %*% as.matrix(c(b0_TUC_IL, brear_TUC_IL, borigin_TUC_IL), ncol = 1))
      
      transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- phi_TUC_IL/(1 + phi_TUC_IL)
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
  # Trim to maximum number of observations
  n.obs <- unlist(lapply(det_hist_sim, ncol))
  max_obs <- max(n.obs)
  det_hist_sim_array <- movement_array[,1:max_obs,]
  
  # As an array
  transition_date_sim_matrix <- transition_date[,1:max_obs]
  
  date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2016-12-31"))
  transition_date_sim_matrix %>% 
    as_tibble() %>% 
    mutate_all(date_numeric) -> transition_date_sim_numeric
  
  # Return a list of two things: [[1]] = the detection history, [[2]] = the dates
  
  return(list(det_hist_sim_array, transition_date_sim_numeric))
}

fish_sim_cat_data_600 <- data.frame(tag_code = seq(1, 600, 1),
                                    natal_origin = c(rep(1, 200),
                                                     rep(2, 200),
                                                     rep(3, 200)),
                                    rear_type = rep(c(rep(1, 100), rep(2, 100)), 3))

fish_sim_cat_data_1200 <- data.frame(tag_code = seq(1, 1200, 1),
                                     natal_origin = c(rep(1, 400),
                                                      rep(2, 400),
                                                      rep(3, 400)),
                                     rear_type = rep(c(rep(1, 200), rep(2, 200)), 3))

fish_sim_cat_data_3000 <- data.frame(tag_code = seq(1, 3000, 1),
                                     natal_origin = c(rep(1, 1000),
                                                      rep(2, 1000),
                                                      rep(3, 1000)),
                                     rear_type = rep(c(rep(1, 500), rep(2, 500)), 3))


##### Loop the simulation

### 600 fish
sim_600_hist_list <- list()
sim_600_dates_list <- list()
nfish = 600
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

# Create the design matrix for categorical variables
cat_X_mat_600 <- matrix(0, nrow = nfish, ncol = 4)
# The first column everyone gets a 1 (this is beta 0/grand mean mu)
cat_X_mat_600[,1] <- 1

for (i in 1:nfish){
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


# Run for loops in parallel
origin_sim_runs <- foreach(z = 1:10) %dopar% {
  origin_cov_sim(nfish = 600, fish_sim_cat_data = fish_sim_cat_data_600,cat_X_mat = cat_X_mat_600, 
                            movement_array = movement_array, transition_matrix = transition_matrix)
}

for (i in 1:10){
  sim_600_hist_list[[i]] <- origin_sim_runs[[i]][[1]]
  sim_600_dates_list[[i]] <- origin_sim_runs[[i]][[2]]
}


# Export these lists as RDS objects
saveRDS(sim_600_hist_list, here::here("simulation", "sim_600_cov_origin_hist_list.rds"))

saveRDS(sim_600_dates_list, here::here("simulation", "sim_600_cov_origin_dates_list.rds"))
write.csv(fish_sim_cat_data_600, here::here("simulation", "origin_rear_600.csv"), row.names = FALSE)


### 1200 fish
sim_1200_hist_list <- list()
sim_1200_dates_list <- list()
nfish = 1200
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

# Create the design matrix for categorical variables
cat_X_mat_1200 <- matrix(0, nrow = nfish, ncol = 4)
# The first column everyone gets a 1 (this is beta 0/grand mean mu)
cat_X_mat_1200[,1] <- 1

for (i in 1:nfish){
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


# Run for loops in parallel
origin_sim_runs <- foreach(z = 1:10) %dopar% {
  origin_cov_sim(nfish = 1200, fish_sim_cat_data = fish_sim_cat_data_1200,cat_X_mat = cat_X_mat_1200, 
                 movement_array = movement_array, transition_matrix = transition_matrix)
}

for (i in 1:10){
  sim_1200_hist_list[[i]] <- origin_sim_runs[[i]][[1]]
  sim_1200_dates_list[[i]] <- origin_sim_runs[[i]][[2]]
}

# Export these lists as RDS objects
saveRDS(sim_1200_hist_list, here::here("simulation", "simulated_data", "origin_only", "sim_1200_cov_origin_hist_list.rds"))

saveRDS(sim_1200_dates_list, here::here("simulation", "simulated_data", "origin_only", "sim_1200_cov_origin_dates_list.rds"))
write.csv(fish_sim_cat_data_1200, here::here("simulation", "simulated_data", "origin_only", "origin_rear_1200.csv"), row.names = FALSE)


### 3000 fish
sim_3000_hist_list <- list()
sim_3000_dates_list <- list()
nfish = 3000
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

for (z in 1:10){
  sim_run <- cov_sim(nfish = 3000, fish_sim_cat_data = fish_sim_cat_data_3000)
  
  sim_3000_hist_list[[z]] <- sim_run[[1]]
  sim_3000_dates_list[[z]] <- sim_run[[2]]
}

# Export these lists as RDS objects
saveRDS(sim_3000_hist_list, here::here("simulation", "sim_3000_cov_categorical_hist_list.rds"))
saveRDS(sim_3000_dates_list, here::here("simulation", "sim_3000_cov_categorical_dates_list.rds"))
write.csv(fish_sim_cat_data_3000, here::here("simulation", "origin_rear_3000.csv"), row.names = FALSE)






##### SIMULATION FUNCTION 7: TEMP, FLOW, REAR, ORIGIN #####

all4_cov_sim <- function(nfish, fish_sim_cat_data, cat_X_mat, movement_array, transition_matrix){
  # Take a uniform distribution between July 1 and September 30
  BON_arrival_dates_sim <- seq(ymd("2017-07-01"), ymd("2019-09-30"), by = "days")[round(runif(3000, min = 1, max = 92), 0)]
  
  state_date <- matrix(nrow = nfish, ncol = 100)
  transition_date <- matrix(nrow = nfish, ncol = 100)
  
  for (i in 1:nfish){ # for each fish
    print(paste0("i = ", i))
    # Populate the first state of each matrix as BON to MCN (since we first see each fish at BON)
    movement_array["mainstem, BON to MCN", 1, i] <- 1
    # Populate state date matrix with date at arrival at BON
    state_date[i,1] <- format(as.Date(BON_arrival_dates_sim[i], origin = "1970-01-01"))
    # Store this in the transition date matrix
    transition_date[i,1] <- state_date[i,1]
    
    # populate the rest of the state dates for this fish by adding two weeks
    for (k in 2:100){
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
      ##### FROM MOUTH TO BON STATE #####
      
      # Mouth to BON to BON to MCN transition
      b0_MB_BM <- 1
      bflow_MB_BM <- 0 # No relationship with flow
      btemp_MB_BM <- 0 # No relationship with temperature
      brear_MB_BM <- 0 # No difference for hatchery vs. wild
      borigin_MB_BM <- c(0, 0) # No difference for JDR vs. TUC or YAK vs TUC
      
      # Here we've reformatted it so that b0, brear, and borigin are all in a vector that we multiply by the design matrix for the categorical covariates
      phi_MB_BM <- exp(btemp_MB_BM*temp_BON + bflow_MB_BM*flow_BON + cat_X_mat[i,] %*% as.matrix(c(b0_MB_BM, brear_MB_BM, borigin_MB_BM), ncol = 1))
      
      transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- phi_MB_BM/(1 + phi_MB_BM)
      transition_matrix["mainstem, mouth to BON", "loss"] <- 1 - transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"]
      
      # Check values
      transition_matrix["mainstem, mouth to BON", ]
      
      ##### FROM BON TO MCN STATE #####
      
      # BON to MCN to Mouth to BON transition
      b0_BM_MB <- 1
      bflow_BM_MB <- 0.5 # 
      btemp_BM_MB <- 0 # No relationship with temperature
      brear_BM_MB <- 0 # No difference for hatchery vs. wild
      borigin_BM_MB <- c(0, 0) # No difference for JDR vs. TUC or YAK vs TUC
      
      # BON to MCN to MCN to ICH or PRA transition
      b0_BM_MIP <- 1
      bflow_BM_MIP <- 0 # No relationship with flows
      btemp_BM_MIP <- 0.5 # 
      brear_BM_MIP <- 0 # No difference for hatchery vs. wild
      borigin_BM_MIP <- c(-0.5, 0.25) # JDR have a lower probability of overshooting MCN than TUC fish, but there is no difference for YAK vs TUC
      
      # BON to MCN to DES R transition
      b0_BM_DES <- 1
      bflow_BM_DES <- 0 # flow not relevant for tributary entry
      btemp_BM_DES <- 0 # Temperature not relevant for tributary entry
      brear_BM_DES <- 0.5 # higher probability of straying to Deschutes if hatchery
      borigin_BM_DES <- c(0.25, -0.125) # JDR have a higher probability of straying to DES than TUC fish (closer), but there is no difference for YAK vs TUC
      
      # BON to MCN to JDR transition
      b0_BM_JDR <- 1
      bflow_BM_JDR <- 0 # flow not relevant for tributary entry
      btemp_BM_JDR <- 0 # Temperature not relevant for tributary entry
      brear_BM_JDR <- 0 # No difference for hatchery vs. wild
      borigin_BM_JDR <- c(0.5, -0.25) # JDR have a much higher probability of homing to JDR than TUC fish , but there is no difference for YAK vs TUC
      
      # Evaluate
      
      # BON to MCN
      phi_BM_MB <- exp(btemp_BM_MB*temp_BON + bflow_BM_MB*flow_BON + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_MB, brear_BM_MB, borigin_BM_MB), ncol = 1))
      phi_BM_MIP <- exp(btemp_BM_MIP*temp_MCN + bflow_BM_MIP*flow_MCN + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_MIP, brear_BM_MIP, borigin_BM_MIP), ncol = 1))
      phi_BM_DES <- exp(btemp_BM_DES*temp_DES + bflow_BM_DES*flow_DES + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_DES, brear_BM_DES, borigin_BM_DES), ncol = 1))
      phi_BM_JDR <- exp(btemp_BM_JDR*temp_JDR + bflow_BM_JDR*flow_JDR + cat_X_mat[i,] %*%  as.matrix(c(b0_BM_JDR, brear_BM_JDR, borigin_BM_JDR), ncol = 1))
      
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
      bflow_MIP_BM <- 0.05 # 
      btemp_MIP_BM <- -0.5 # 
      brear_MIP_BM <- 0 # No difference for hatchery vs. wild
      borigin_MIP_BM <- c(0.5, -0.25) # JDR have much higher chance of falling back than TUC, but no difference between YAK and TUC
      
      # BON to MCN to PRA to RIS transition
      b0_MIP_PR <- 1
      bflow_MIP_PR <- 0 # No relationship with flow
      btemp_MIP_PR <- 0.3 # When it's hotter, more likely to go upstream
      brear_MIP_PR <- 0 # No difference for hatchery vs. wild
      borigin_MIP_PR <- c(0,0) # No relationship with origin - no origins are up here
      
      # BON to MCN to ICH to LGR transition
      b0_MIP_IL <- 1
      bflow_MIP_IL <- -0.3 # -0.3 # Negative relationship - higher flow = lower chance of ascending
      btemp_MIP_IL <- 1 # 0.2 # Positive relationship with temperature
      brear_MIP_IL <- c(0) # No effect of rear type
      borigin_MIP_IL <- c(-0.5, 0) # JDR have much lower chance than TUC, YAK have slightly lower chance than TUC
      
      # BON to MCN to Yakima transition
      b0_MIP_YAK <- 1
      bflow_MIP_YAK <- 0 # flow not relevant for tributary entry
      btemp_MIP_YAK <- 0 # Temperature not relevant for tributary entry
      brear_MIP_YAK <- -0.2 # Hatchery fish less likely to enter Yakima
      borigin_MIP_YAK <- c(-0.25, 0.5) # No difference for JDR vs. TUC, but much higher for YAK than JDR
      
      phi_MIP_BM <- exp(btemp_MIP_BM*temp_MCN + bflow_MIP_BM*flow_MCN + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_BM, brear_MIP_BM, borigin_MIP_BM), ncol = 1))
      phi_MIP_PR <- exp(btemp_MIP_PR*temp_PRA + bflow_MIP_PR*flow_PRA + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_PR, brear_MIP_PR, borigin_MIP_PR), ncol = 1))
      phi_MIP_IL <- exp(btemp_MIP_IL*temp_ICH + bflow_MIP_IL*flow_ICH + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_IL, brear_MIP_IL, borigin_MIP_IL), ncol = 1))
      phi_MIP_YAK <- exp(btemp_MIP_YAK*temp_YAK + bflow_MIP_YAK*flow_YAK + cat_X_mat[i,] %*% as.matrix(c(b0_MIP_YAK, brear_MIP_YAK, borigin_MIP_YAK), ncol = 1))
      
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
      brear_PR_MIP <- 0 # No effect of rear type
      borigin_PR_MIP <- c(0,0) # No effect of origins vs. TUC 
      
      phi_PR_MIP <- exp(btemp_PR_MIP*temp_PRA + bflow_PR_MIP*flow_PRA + cat_X_mat[i,] %*% as.matrix(c(b0_PR_MIP, brear_PR_MIP, borigin_PR_MIP), ncol = 1))
      
      transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- phi_PR_MIP/(1 + phi_PR_MIP)
      transition_matrix["mainstem, PRA to RIS", "loss"] <- 1 - transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["mainstem, PRA to RIS", ]
      
      
      
      
      ##### FROM ICH TO LGR STATE #####
      
      # ICH to LGR to MCN to ICH or PRA transition #
      b0_IL_MIP <- 1
      bflow_IL_MIP <- 0 # No relationship with flow
      btemp_IL_MIP <- 0.2 # 0.2 # Negative relationship with temperature - when it's colder, tend to go downstream more
      brear_IL_MIP <- 0 # No effect of rear type
      borigin_IL_MIP <- c(0.25, 0.25) # JDR and YAK fish are both more likely than TUC fish to return to MCN to ICH or PRA state
      
      # ICH to LGR to Tucannon River transition #
      b0_IL_TUC <- 1
      bflow_IL_TUC <- 0 # No relationship with flow (tributary)
      btemp_IL_TUC <- 0 # No relationship with temperature (tributary)
      brear_IL_TUC <- 0 # No effect of rear type
      borigin_IL_TUC <- c(-0.25, -0.25) # Both other origins are less likely to enter TUC than TUC fish
      
      
      phi_IL_MIP <- exp(btemp_IL_MIP*temp_ICH + bflow_IL_MIP*flow_ICH + cat_X_mat[i,] %*% as.matrix(c(b0_IL_MIP, brear_IL_MIP, borigin_IL_MIP), ncol = 1))
      phi_IL_TUC <- exp(btemp_IL_TUC*temp_TUC + bflow_IL_TUC*flow_TUC + cat_X_mat[i,] %*% as.matrix(c(b0_IL_TUC, brear_IL_TUC, borigin_IL_TUC), ncol = 1))
      
      transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- phi_IL_MIP/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- phi_IL_TUC/(1 + phi_IL_MIP + phi_IL_TUC)
      transition_matrix["mainstem, ICH to LGR", "loss"] <- 1 - sum(transition_matrix["mainstem, ICH to LGR", c("mainstem, MCN to ICH or PRA", "Tucannon River")])
      
      # Check values
      # transition_matrix["mainstem, ICH to LGR",]
      
      
      ##### FROM DESCHUTES RIVER STATE #####
      # Deschutes River to mainstem, BON to MCN transition
      b0_DES_BM <- 1 # Make all of the return intercepts 1, instead of 0 - increase chance of loss param
      bflow_DES_BM <- 0 # No relationship with flow
      btemp_DES_BM <- 0 # No relationship with temperature
      brear_DES_BM <- -0.5 # Hatchery fish less likely to return to mainstem
      borigin_DES_BM <- c(0,0) # No relationship with rear type
      
      phi_DES_BM <- exp(btemp_DES_BM*temp_DES + bflow_DES_BM*flow_DES + cat_X_mat[i,] %*% as.matrix(c(b0_DES_BM, brear_DES_BM, borigin_DES_BM), ncol = 1))
      
      transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- phi_DES_BM/(1 + phi_DES_BM)
      transition_matrix["Deschutes River", "loss"] <- 1 - transition_matrix["Deschutes River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["Deschutes River", ]
      
      
      ##### FROM JOHN DAY RIVER STATE #####
      # JDR to mainstem, BON to MCN transition
      b0_JDR_BM <- 1
      bflow_JDR_BM <- 0 # No relationship with flow
      btemp_JDR_BM <- 0 # No relationship with temperature
      brear_JDR_BM <- 0 # No effect of rear type
      borigin_JDR_BM <- c(-0.5, 0.25) # JDR fish less likely to return than TUC, no difference for YAK vs. TUC
      
      phi_JDR_BM <- exp(btemp_JDR_BM*temp_JDR + bflow_JDR_BM*flow_JDR + cat_X_mat[i,] %*% as.matrix(c(b0_JDR_BM, brear_JDR_BM, borigin_JDR_BM), ncol = 1))
      
      transition_matrix["John Day River", "mainstem, BON to MCN"] <- phi_JDR_BM/(1 + phi_JDR_BM)
      transition_matrix["John Day River", "loss"] <- 1 - transition_matrix["John Day River", "mainstem, BON to MCN"]
      
      # Check values
      # transition_matrix["John Day River", ]
      
      
      ##### FROM YAKIMA RIVER STATE #####
      # Yakima River to mainstem, MCN to ICH or PRA transition
      b0_YAK_MIP <- 1
      bflow_YAK_MIP <- 0 # No relationship with flow
      btemp_YAK_MIP <- 0 # No relationship with temperature
      brear_YAK_MIP <- 0 # No effect of rear type
      borigin_YAK_MIP <- c(0.25, -0.5) # Yakima R. fish less likely to return
      
      phi_YAK_MIP <- exp(btemp_YAK_MIP*temp_YAK + bflow_YAK_MIP*flow_YAK + cat_X_mat[i,] %*% as.matrix(c(b0_YAK_MIP, brear_YAK_MIP, borigin_YAK_MIP), ncol = 1))
      
      transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- phi_YAK_MIP/(1 + phi_YAK_MIP)
      transition_matrix["Yakima River", "loss"] <- 1 - transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"]
      
      # Check values
      # transition_matrix["Yakima River", ]
      
      
      ##### FROM TUCANNON RIVER STATE #####
      # From Tucannon River to ICH to LGR transition
      b0_TUC_IL <- 1
      bflow_TUC_IL <- 0 # No relationship with flow
      btemp_TUC_IL <- 0 # No relationship with temperature
      brear_TUC_IL <- 0 # No effect of rear type
      borigin_TUC_IL <- c(0.25,0.25) # JDR and YAK fish both more likely to return than TUC fish
      
      phi_TUC_IL <- exp(btemp_TUC_IL*temp_TUC + bflow_TUC_IL*flow_TUC + cat_X_mat[i,] %*% as.matrix(c(b0_TUC_IL, brear_TUC_IL, borigin_TUC_IL), ncol = 1))
      
      transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- phi_TUC_IL/(1 + phi_TUC_IL)
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
  # Trim to maximum number of observations
  n.obs <- unlist(lapply(det_hist_sim, ncol))
  max_obs <- max(n.obs)
  det_hist_sim_array <- movement_array[,1:max_obs,]
  
  # As an array
  transition_date_sim_matrix <- transition_date[,1:max_obs]
  
  date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(x) - ymd("2016-12-31"))
  transition_date_sim_matrix %>% 
    as_tibble() %>% 
    mutate_all(date_numeric) -> transition_date_sim_numeric
  
  # Return a list of two things: [[1]] = the detection history, [[2]] = the dates
  
  return(list(det_hist_sim_array, transition_date_sim_numeric))
}

### 1200 fish
fish_sim_cat_data_1200 <- data.frame(tag_code = seq(1, 1200, 1),
                                     natal_origin = c(rep(1, 400),
                                                      rep(2, 400),
                                                      rep(3, 400)),
                                     rear_type = rep(c(rep(1, 200), rep(2, 200)), 3))

sim_1200_hist_list <- list()
sim_1200_dates_list <- list()
nfish = 1200
movement_array <- array(0, dim = c(nstates, 100, nfish))
# Set rownames of each matrix to the states
rownames(movement_array) <- sim_states

# Create the design matrix for categorical variables
cat_X_mat_1200 <- matrix(0, nrow = nfish, ncol = 4)
# The first column everyone gets a 1 (this is beta 0/grand mean mu)
cat_X_mat_1200[,1] <- 1

for (i in 1:nfish){
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



# Run for loops in parallel
all4_sim_runs <- foreach(z = 1:10) %dopar% {
  all4_cov_sim(nfish = 1200, fish_sim_cat_data = fish_sim_cat_data_1200,cat_X_mat = cat_X_mat_1200, 
                      movement_array = movement_array, transition_matrix = transition_matrix)
}

for (i in 1:10){
  sim_1200_hist_list[[i]] <- all4_sim_runs[[i]][[1]]
  sim_1200_dates_list[[i]] <- all4_sim_runs[[i]][[2]]
}

# Export these lists as RDS objects
saveRDS(sim_1200_hist_list, here::here("simulation", "simulated_data", "all4", "sim_1200_all4_hist_list.rds"))

saveRDS(sim_1200_dates_list, here::here("simulation", "simulated_data", "all4", "sim_1200_all4_dates_list.rds"))
write.csv(fish_sim_cat_data_1200, here::here("simulation", "simulated_data", "all4", "all4_1200.csv"), row.names = FALSE)
