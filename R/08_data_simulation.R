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


##### Simulate covariate data #####


# Simulate temperature data
# Get dates from 2011-2015

# Simulate as sin wave plus noise
x <- seq(0,pi*6,length.out=365*3)
y <- sin(x)
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

# Create a temperature df
temp_sim_df <- data.frame(date = dates, 
                          BON = BON_temp_sim, 
                          MCN = MCN_temp_sim, 
                          PRA = PRA_temp_sim, 
                          ICH = ICH_temp_sim)



# Simulate the categorical covariates (rear type and natal origin)

fish_sim_cat_data <- data.frame(tag_code = seq(1, 3000, 1),
                            natal_origin = c(rep("JDR", 1000),
                                             rep("Yakima", 1000),
                                             rep("Tucannon", 1000)),
                            rear_type = rep(c(rep("H", 500), rep("W", 500)), 3))



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






##### Model including covariates #####


# Our movement probabilities will follow the following functional relationships:

# For the mouth to BON state
# transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- 0.9 (no covariates)
# transition_matrix["mainstem, mouth to BON", "loss"] <- 0.1 (no covariates)

# BON to MCN
# phi_BM_MB <- exp(b0_BM_MB + btemp_BM_MB*temp + bflow_BM_MB*flow + brear_BM_MB[i] + borigin_BM_MB[i])
# phi_BM_MIP <- exp(b0_BM_MIP + btemp_BM_MIP*temp + bflow_BM_MIP*flow + brear_BM_MIP[i] + borigin_BM_MIP[i])
# phi_BM_DES <- exp(b0_BM_DES + btemp_BM_DES*temp + bflow_BM_DES*flow + brear_BM_DES[i] + borigin_BM_DES[i])
# phi_BM_JDR <- exp(b0_BM_JDR + btemp_BM_JDR*temp + bflow_BM_JDR*flow + brear_BM_JDR[i] + borigin_BM_JDR[i])


# transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- phi_BM_MB/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
# transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- phi_BM_MIP/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
# transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- phi_BM_DES/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
# transition_matrix["mainstem, BON to MCN", "John Day River"] <- phi_BM_JDR/(1 + phi_BM_MB + phi_BM_MIP + phi_BM_DES + phi_BM_JDR)
# transition_matrix["mainstem, BON to MCN", "loss"] <- 1 - the others


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














