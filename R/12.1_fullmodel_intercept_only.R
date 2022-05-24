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
  "mainstem, upstream of ICH",
  
  # Tributary sites
  "Deschutes River",
  "John Day River",
  "Hood River",
  "Fifteenmile Creek",
  "Umatilla River",
  "Yakima River",
  "Walla Walla River",
  "Wenatchee River",
  "Entiat River",
  "Okangon River",
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

















