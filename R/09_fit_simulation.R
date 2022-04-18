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
