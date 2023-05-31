### Temperature MAR model

# This script will fit a MAR model to temperature data from different dams
# to estimate temperature for every date in our time series.

# The MAR model that we will fit is a state-space model, with forebay and tailrace
# temperatures as observations of a "Columbia River Basin" temperature, with
# each dam estimated as an offset of it
# We will have sixteen observations of the same process (forebay and tailrace temperatures for eight dams)

# Load libraries
library(tidyverse)
library(MARSS)

# Load data
# note that these data files were produced by script 22, which filtered out bad data
# in addition to some other steps. The only columns that we need from these
# files are date, tailrace_temp, and forebay_temp
BON_temp <- read.csv(here::here("covariate_data", "for_model", "temp", "BON_window_temp.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
MCN_temp <- read.csv(here::here("covariate_data", "for_model", "temp", "MCN_window_temp.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
PRA_temp <- read.csv(here::here("covariate_data", "for_model", "temp", "PRA_window_temp.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
RIS_temp <- read.csv(here::here("covariate_data", "for_model", "temp", "RIS_window_temp.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
RRE_temp <- read.csv(here::here("covariate_data", "for_model", "temp", "RRE_window_temp.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
ICH_temp <- read.csv(here::here("covariate_data", "for_model", "temp", "ICH_window_temp.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
WEL_temp <- read.csv(here::here("covariate_data", "for_model", "temp", "WEL_window_temp.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
LGR_temp <- read.csv(here::here("covariate_data", "for_model", "temp", "LGR_window_temp.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)

# reformat data for model - drop date and transpose, bind all dams together
t(BON_temp[,2:3]) %>% 
  rbind(., t(MCN_temp[,2:3])) %>% 
  rbind(., t(PRA_temp[,2:3])) %>% 
  rbind(., t(RIS_temp[,2:3])) %>% 
  rbind(., t(RRE_temp[,2:3])) %>% 
  rbind(., t(ICH_temp[,2:3])) %>% 
  rbind(., t(WEL_temp[,2:3])) %>% 
  rbind(., t(LGR_temp[,2:3])) -> temp_for_MAR

dat = temp_for_MAR
# fit the model

# no bias term, no errors for observation process, no 
R <- A <- U <- "zero"
B <- Z <- "identity"
Q <- "equalvarcov"
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
kem <- MARSS(dat, model = model.list)










