# 11 Full dataset - reformat for JAGS model

# In this R script, we will reformat all of the detection histories so that we can fit the JAGS model

### Load libraries
library(here)
library(tidyverse)
library(lubridate)
library(janitor)

### Load data
states_times <- read.csv(here::here("model_files", "states.csv"))

unique(states_times$tag_code)

sort(unique(states_times$state))
