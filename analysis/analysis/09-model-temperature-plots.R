# 09-model-temperature-plots

# This script takes the output from the stan model runs in 05-stan-runs and
# plots the effects of temperature

#### Load libraries, state information ####
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(dplyr)
library(tidyr)
library(here)
library(ggpubr)
library(stringr)
library(rlang)
library(tibble)
library(forcats)
library(lubridate)

# Load the temperature summary
window_temp_summary <- read.csv(here::here("stan_actual", "window_temps_summary.csv"), row.names = 1)


# get the model states into a df, to help with interpretation
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
  "mainstem, upstream of LGR",
  
  # Tributary states ()
  # With detection efficiencies in the model, we now have more tributary states,
  # since we have an upstream and a river mouth state
  
  # "Deschutes River", 
  "Deschutes River Mouth", "Deschutes River Upstream",
  # "John Day River", 
  "John Day River Mouth", "John Day River Upstream",
  # "Hood River",
  "Hood River Mouth", "Hood River Upstream",
  # "Fifteenmile Creek", 
  "Fifteenmile Creek Mouth", "Fifteenmile Creek Upstream",
  # "Umatilla River",
  "Umatilla River Mouth", "Umatilla River Upstream",
  # "Yakima River",
  "Yakima River Mouth", "Yakima River Upstream",
  # "Walla Walla River",
  "Walla Walla River Mouth", "Walla Walla River Upstream",
  # "Wenatchee River", 
  "Wenatchee River Mouth", "Wenatchee River Upstream",
  # "Entiat River", 
  "Entiat River Mouth", "Entiat River Upstream",
  # "Okanogan River", 
  "Okanogan River Mouth", "Okanogan River Upstream",
  # "Methow River", 
  "Methow River Mouth", "Methow River Upstream",
  # "Tucannon River",
  "Tucannon River Mouth", "Tucannon River Upstream",
  # "Asotin Creek", 
  "Asotin Creek Mouth", "Asotin Creek Upstream",
  "Clearwater River",
  "Salmon River",
  "Grande Ronde River",
  # "Imnaha River",
  "Imnaha River Mouth", "Imnaha River Upstream",
  "BON to MCN other tributaries",
  "Upstream WEL other tributaries",
  
  # Loss
  "loss"
)

# Get info about state names and numbers
from_state_number_names <- data.frame(from = seq(1,43,1), from_name = model_states)
to_state_number_names <- data.frame(to = seq(1,43,1), to_name = model_states)


# get the info on transitions
UCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_wild","UCW_transition_counts.csv"))
UCW_movements <- paste0("_", UCW_transition_counts$from, "_", UCW_transition_counts$to)
UCW_movements <- UCW_movements[!(grepl("NA", UCW_movements))]

UCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_hatchery","UCH_transition_counts.csv"))
UCH_movements <- paste0("_", UCH_transition_counts$from, "_", UCH_transition_counts$to)
UCH_movements <- UCH_movements[!(grepl("NA", UCH_movements))]

MCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_wild","MCW_transition_counts.csv"))
MCW_movements <- paste0("_", MCW_transition_counts$from, "_", MCW_transition_counts$to)
MCW_movements <- MCW_movements[!(grepl("NA", MCW_movements))]

MCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_hatchery","MCH_transition_counts.csv"))
MCH_movements <- paste0("_", MCH_transition_counts$from, "_", MCH_transition_counts$to)
MCH_movements <- MCH_movements[!(grepl("NA", MCH_movements))]

SRW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_wild","SRW_transition_counts.csv"))
SRW_movements <- paste0("_", SRW_transition_counts$from, "_", SRW_transition_counts$to)
SRW_movements <- SRW_movements[!(grepl("NA", SRW_movements))]

SRH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_hatchery","SRH_transition_counts.csv"))
SRH_movements <- paste0("_", SRH_transition_counts$from, "_", SRH_transition_counts$to)
SRH_movements <- SRH_movements[!(grepl("NA", SRH_movements))]


##### Load the model runs #####

# Load the model data associated with each run (necessary to load covariates)
# Store these each in an environment, because most things share names
UCW_envir <- new.env()
UCH_envir <- new.env()
MCW_envir <- new.env()
MCH_envir <- new.env()
SRW_envir <- new.env()
SRH_envir <- new.env()
load(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_wild", "model_data.rda"),
     envir = UCW_envir)
load(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_hatchery", "model_data.rda"),
     envir = UCH_envir)
load(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_wild", "model_data.rda"),
     envir = MCW_envir)
load(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_hatchery", "model_data.rda"),
     envir = MCH_envir)
load(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_wild", "model_data.rda"),
     envir = SRW_envir)
load(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_hatchery", "model_data.rda"),
     envir = SRH_envir)




# Function to bind four chains together
bind4chains <- function(chain1, chain2, chain3, chain4){
  bound_draws <- bind_draws(chain1$draws(),
                            chain2$draws(),
                            chain3$draws(),
                            chain4$draws(), along = "chain")
  
  return(bound_draws)
}

## Upper Columbia, Wild
UCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_wild", "chain1_UCW_reparam_v2_fit.rds"))
UCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_wild", "chain2_UCW_reparam_v2_fit.rds"))
UCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_wild", "chain3_UCW_reparam_v2_fit.rds"))
UCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_wild", "chain4_UCW_reparam_v2_fit.rds"))

# bind chains together
UCW_fit_raw <- bind4chains(UCW_chain1, UCW_chain2, UCW_chain3, UCW_chain4)
# thin2
thin_draws(UCW_fit_raw, thin = 2) -> UCW_fit
# summarise
UCW_fit_summary <- summarise_draws(UCW_fit)

## Upper Columbia, Hatchery
UCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_hatchery", "chain1_UCH_reparam_v2_fit.rds"))
UCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_hatchery", "chain2_UCH_reparam_v2_fit.rds"))
UCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_hatchery", "chain3_UCH_reparam_v2_fit.rds"))
UCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "upper_columbia_hatchery", "chain4_UCH_reparam_v2_fit.rds"))

# bind chains together
UCH_fit_raw <- bind4chains(UCH_chain1, UCH_chain2, UCH_chain3, UCH_chain4)
# thin2
thin_draws(UCH_fit_raw, thin = 2) -> UCH_fit
# summarise
UCH_fit_summary <- summarise_draws(UCH_fit)

## Middle Columbia, Wild
MCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_wild", "chain1_MCW_reparam_v2_fit.rds"))
MCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_wild", "chain2_MCW_reparam_v2_fit.rds"))
MCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_wild", "chain3_MCW_reparam_v2_fit.rds"))
MCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_wild", "chain4_MCW_reparam_v2_fit.rds"))

# bind chains together
MCW_fit_raw <- bind4chains(MCW_chain1, MCW_chain2, MCW_chain3, MCW_chain4)
# thin2
thin_draws(MCW_fit_raw, thin = 2) -> MCW_fit
# summarise
MCW_fit_summary <- summarise_draws(MCW_fit)

## Middle Columbia, Hatchery
MCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_hatchery", "chain1_MCH_reparam_v2_fit.rds"))
MCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_hatchery", "chain2_MCH_reparam_v2_fit.rds"))
MCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_hatchery", "chain3_MCH_reparam_v2_fit.rds"))
MCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "middle_columbia_hatchery", "chain4_MCH_reparam_v2_fit.rds"))

# bind chains together
MCH_fit_raw <- bind4chains(MCH_chain1, MCH_chain2, MCH_chain3, MCH_chain4)
# thin2
thin_draws(MCH_fit_raw, thin = 2) -> MCH_fit
# summarise
MCH_fit_summary <- summarise_draws(MCH_fit)

## Snake River, Wild
SRW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_wild", "chain1_SRW_reparam_v2_fit.rds"))
SRW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_wild", "chain2_SRW_reparam_v2_fit.rds"))
SRW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_wild", "chain3_SRW_reparam_v2_fit.rds"))
SRW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_wild", "chain4_SRW_reparam_v2_fit.rds"))

# bind chains together
SRW_fit_raw <- bind4chains(SRW_chain1, SRW_chain2, SRW_chain3, SRW_chain4)
# thin2
thin_draws(SRW_fit_raw, thin = 2) -> SRW_fit
# summarise
SRW_fit_summary <- summarise_draws(SRW_fit)

## Snake River, Hatchery
SRH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_hatchery", "chain1_SRH_reparam_v2_fit.rds"))
SRH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_hatchery", "chain2_SRH_reparam_v2_fit.rds"))
SRH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_hatchery", "chain3_SRH_reparam_v2_fit.rds"))
SRH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_spillrescaled", "snake_river_hatchery", "chain4_SRH_reparam_v2_fit.rds"))

# bind chains together
SRH_fit_raw <- bind4chains(SRH_chain1, SRH_chain2, SRH_chain3, SRH_chain4)
# thin2
thin_draws(SRH_fit_raw, thin = 2) -> SRH_fit
# summarise
SRH_fit_summary <- summarise_draws(SRH_fit)

#### Extract all parameter values from the model fit objects ####

# Here, you need to make sure to check the origin params, because those change by DPS

# function to take a parameter type (fixed effect) and store all of them in an array
make_parameter_draws_array <- function(parameter_prefix, fit, fit_summary){
  # extract b0 as an array
  parameters <- fit_summary$variable
  parameters[grepl(paste0(parameter_prefix, "_"), parameters)] -> param_subset
  param_subset[!(grepl("vector", param_subset))] -> param_subset
  # drop the NDE parameters
  param_subset <- param_subset[!(grepl("_NDE", param_subset))]
  
  
  param_subset_from = as.numeric(sub("[^_]*_[^_]*_", "", str_extract(param_subset, "[^_]*_[^_]*_[^_]*")))
  param_subset_to = as.numeric(str_extract(sub("[^_]*_[^_]*_[^_]*_", "", param_subset), "\\d+"))
  
  param_subset_indices <- data.frame(parameter = param_subset, from = param_subset_from, to = param_subset_to)
  
  
  # arrange all parameter values into an array
  param_array <- array(data = 0, dim = c(length(model_states), length(model_states),
                                         length(as.matrix(fit[,,1]))))
  
  # 0 is meaningful for loss (this is what is used in the stan code)
  # 0s for all other movements that are not overwritten will not be used
  
  
  for(i in 1:nrow(param_subset_indices)){
    param_array[param_subset_indices[i, "from"], param_subset_indices[i, "to"], ] <- as.matrix(fit[,,param_subset_indices[i, "parameter"]])
  }
  
  return(param_array)
}

### UCW ###
UCW_parameters <- UCW_fit_summary$variable
UCW_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", UCW_parameters)] -> UCW_fixed_effects
UCW_fixed_effects[!(grepl("vector", UCW_fixed_effects))] -> UCW_fixed_effects

b0_array_UCW <- make_parameter_draws_array(parameter_prefix = "b0", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp0_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = UCW_fit, fit_summary = UCW_fit_summary)
bspillwindow_array_UCW <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = UCW_fit, fit_summary = UCW_fit_summary)
bwinterspill_array_UCW <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp0xorigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp1xorigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp0xorigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp1xorigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp0xorigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)
btemp1xorigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)
borigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
borigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
borigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)

### UCH ###
UCH_parameters <- UCH_fit_summary$variable
UCH_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", UCH_parameters)] -> UCH_fixed_effects
UCH_fixed_effects[!(grepl("vector", UCH_fixed_effects))] -> UCH_fixed_effects

b0_array_UCH <- make_parameter_draws_array(parameter_prefix = "b0", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp0_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = UCH_fit, fit_summary = UCH_fit_summary)
bspillwindow_array_UCH <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = UCH_fit, fit_summary = UCH_fit_summary)
bwinterspill_array_UCH <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp0xorigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp1xorigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp0xorigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp1xorigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp0xorigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)
btemp1xorigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)
borigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
borigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
borigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)

### MCW ###
MCW_parameters <- MCW_fit_summary$variable
MCW_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", MCW_parameters)] -> MCW_fixed_effects
MCW_fixed_effects[!(grepl("vector", MCW_fixed_effects))] -> MCW_fixed_effects

b0_array_MCW <- make_parameter_draws_array(parameter_prefix = "b0", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = MCW_fit, fit_summary = MCW_fit_summary)
bspillwindow_array_MCW <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = MCW_fit, fit_summary = MCW_fit_summary)
bwinterspill_array_MCW <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin1_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin1_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin2_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin2_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin3_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin3_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin4_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin4_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin5_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin5_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp0xorigin6_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin6", fit = MCW_fit, fit_summary = MCW_fit_summary)
btemp1xorigin6_array_MCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin6", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin1_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin2_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin3_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin4_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin5_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = MCW_fit, fit_summary = MCW_fit_summary)
borigin6_array_MCW <- make_parameter_draws_array(parameter_prefix = "borigin6", fit = MCW_fit, fit_summary = MCW_fit_summary)

### MCH ###
MCH_parameters <- MCH_fit_summary$variable
MCH_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", MCH_parameters)] -> MCH_fixed_effects
MCH_fixed_effects[!(grepl("vector", MCH_fixed_effects))] -> MCH_fixed_effects

b0_array_MCH <- make_parameter_draws_array(parameter_prefix = "b0", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp0_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp1_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = MCH_fit, fit_summary = MCH_fit_summary)
bspillwindow_array_MCH <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = MCH_fit, fit_summary = MCH_fit_summary)
bwinterspill_array_MCH <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp0xorigin1_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp1xorigin1_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp0xorigin2_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = MCH_fit, fit_summary = MCH_fit_summary)
btemp1xorigin2_array_MCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = MCH_fit, fit_summary = MCH_fit_summary)
borigin1_array_MCH <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = MCH_fit, fit_summary = MCH_fit_summary)
borigin2_array_MCH <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = MCH_fit, fit_summary = MCH_fit_summary)

### SRW ###
SRW_parameters <- SRW_fit_summary$variable
SRW_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRW_parameters)] -> SRW_fixed_effects
SRW_fixed_effects[!(grepl("vector", SRW_fixed_effects))] -> SRW_fixed_effects

b0_array_SRW <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRW_fit, fit_summary = SRW_fit_summary)
bspillwindow_array_SRW <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRW_fit, fit_summary = SRW_fit_summary)
bwinterspill_array_SRW <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp0xorigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)
btemp1xorigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
borigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)

### SRH ###
SRH_parameters <- SRH_fit_summary$variable
SRH_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRH_parameters)] -> SRH_fixed_effects
SRH_fixed_effects[!(grepl("vector", SRH_fixed_effects))] -> SRH_fixed_effects

b0_array_SRH <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRH_fit, fit_summary = SRH_fit_summary)
bspillwindow_array_SRH <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRH_fit, fit_summary = SRH_fit_summary)
bwinterspill_array_SRH <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp0xorigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)
btemp1xorigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
borigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)

#### Temperature plot ####

# Here, we will probably want to look at only certain movements (?)
# For example, overshoot, holding in Deschutes, straying
# But we could plot all of them

# We have DPS-wide temperature effects and origin-specific ones, 
# and we have a winter/spring and summer/fall param (temp0 or temp1)

# to capture posterior correlations, we'll need to select that same iteration
# for all parameters that affect movement


# create a list that maps origin numbers (params) to what they actually are
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(grepl("other tributaries", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

# Use the parameter map to index the right effects
origin_param_map <- data.frame(
  natal_origin = natal_origins,
  hatchery = c(NA, NA, NA, NA, 1, NA, 2, # MC
               1,NA,2,3, # UC
               5,NA,1,4,2,3), # SR,
  wild = c(1,3,NA,2,4,6,5, # MC
           1,2,NA,3, # UC
           6,1,2,5,3,4)) # SR

# extract_covariate_experiences <- function(envir, rear, origin_select){
#   origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
#   for(i in 1:nrow(envir$data$cat_X_mat)){
#     origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
#   }
#   
#   pop_states_dates <- data.frame(state = as.vector(envir$data$y),
#                                  date = as.vector(envir$data$transition_dates),
#                                  year = ceiling(as.vector(envir$data$transition_dates)/365.25)+1,
#                                  origin = rep(origin_vector, ncol(envir$data$y)))
#   # add mainstem dam for each state
#   dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
#                           state = seq(2,9))
#   pop_states_dates %>% 
#     left_join(., dam_index, by = "state") -> pop_states_dates
#   
#   
#   # reformat covariates so that they can be joined
#   as.data.frame(envir$data$spill_window_data) %>% 
#     rownames_to_column("date") %>% 
#     mutate(date = as.numeric(date)) %>% 
#     pivot_longer(cols = -c(date), names_to = "dam", values_to = "spill_window") -> spill_window_long
#   
#   as.data.frame(envir$data$temperature_data) %>% 
#     rownames_to_column("date") %>% 
#     mutate(date = as.numeric(date)) %>% 
#     pivot_longer(cols = -c(date), names_to = "dam", values_to = "temperature") -> temp_long
#   
#   as.data.frame(envir$data$winter_spill_days_data) %>% 
#     rownames_to_column("year") %>% 
#     mutate(year = as.numeric(year)) %>% 
#     pivot_longer(cols = -c(year), names_to = "dam", values_to = "winter_spill") -> spill_days_long
#   
#   
#   # add temperature, spill window, winter spill days
#   pop_states_dates %>% 
#     left_join(., spill_window_long, by = c("date", "dam")) %>% 
#     left_join(., temp_long, by = c("date", "dam")) %>% 
#     left_join(., spill_days_long, by = c("year", "dam")) -> pop_states_dates_covariates
#   
#   # Now, keep only the origin selected
#   if(rear == "wild"){
#     origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
#   } else {
#     origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
#   }
#   
#   subset(pop_states_dates_covariates, origin == origin_numeric) -> origin_states_dates_covariates
#   
#   return(origin_states_dates_covariates)
#   
# }

extract_covariate_experiences <- function(envir, rear, origin_select){
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  
  # for spill days - include the winter post-overshoot vector, which contains
  # info on whether they could have experienced winter spill conditions or not
  pop_states_dates <- data.frame(state = as.vector(t(envir$data$y)),
                                 date = as.vector(t(envir$data$transition_dates)),
                                 year = ceiling(as.vector(t(envir$data$transition_dates))/365.25)+1,
                                 origin = rep(origin_vector, each = ncol(envir$data$y)))
  
  
  # add mainstem dam for each state
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                          state = seq(2,9))
  pop_states_dates %>% 
    left_join(., dam_index, by = "state") -> pop_states_dates
  
  
  # reformat covariates so that they can be joined
  as.data.frame(envir$data$spill_window_data) %>% 
    rownames_to_column("date") %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(cols = -c(date), names_to = "dam", values_to = "spill_window") -> spill_window_long
  
  as.data.frame(envir$data$temperature_data) %>% 
    rownames_to_column("date") %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(cols = -c(date), names_to = "dam", values_to = "temperature") -> temp_long
  
  as.data.frame(envir$data$winter_spill_days_data) %>% 
    rownames_to_column("year") %>% 
    mutate(year = as.numeric(year)) %>% 
    pivot_longer(cols = -c(year), names_to = "dam", values_to = "winter_spill") -> spill_days_long
  
  
  # add temperature, spill window, winter spill days
  pop_states_dates %>% 
    left_join(., spill_window_long, by = c("date", "dam")) %>% 
    left_join(., temp_long, by = c("date", "dam")) %>% 
    left_join(., spill_days_long, by = c("year", "dam")) -> pop_states_dates_covariates
  
  # drop observations in the loss state and with index 0
  pop_states_dates_covariates %>% 
    filter(!(state %in% c(0,43))) -> pop_states_dates_covariates
  
  # now add winter post-overshoot vector
  pop_states_dates_covariates$winter_post_overshoot_vector <- as.vector(envir$data$winter_post_overshoot_vector)
  
  
  # Now, keep only the origin selected
  if(rear == "wild"){
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
  } else {
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
  }
  
  subset(pop_states_dates_covariates, origin == origin_numeric) -> origin_states_dates_covariates
  
  return(origin_states_dates_covariates)
  
}



# Tell it the movements for which you want to estimate temperature effects
# movements are formatted as matrix, with column for from and column for to
estimate_temp_effect_UCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # get median winter spill for this state
    UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                                   date = as.vector(UCW_envir$data$transition_dates))
    spillwindow_data <- UCW_envir$data$spill_window_data
    med_spillwindow <- median(spillwindow_data[subset(UCW_states_dates, state == from)$date,from])
    
    # get median spill window for this state
    winterspill_data <- UCW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter, i] <- exp(b0_array_UCW[from,to,iter] +
                                                       # btemp0_array_UCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_UCW[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_UCW[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_UCW[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_UCW[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_UCW[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_UCW[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_UCW[from,to,iter]*temp_predict[j]*origin2 + 
                                                       # btemp0xorigin3_array_UCW[from,to,iter]*origin3 +
                                                       btemp1xorigin3_array_UCW[from,to,iter]*temp_predict[j]*origin3 +
                                                       borigin1_array_UCW[from,to,iter]*origin1 +
                                                       borigin2_array_UCW[from,to,iter]*origin2 +
                                                       borigin3_array_UCW[from,to,iter]*origin3)/
          sum(exp(b0_array_UCW[from,possible_movements,iter] +
                    # btemp0_array_UCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_UCW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_UCW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_UCW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_UCW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    borigin1_array_UCW[from,possible_movements,iter]*origin1 +
                    borigin2_array_UCW[from,possible_movements,iter]*origin2 +
                    borigin3_array_UCW[from,possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_UCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # get median winter spill for this state
    UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                                   date = as.vector(UCH_envir$data$transition_dates))
    spillwindow_data <- UCH_envir$data$spill_window_data
    med_spillwindow <- median(spillwindow_data[subset(UCH_states_dates, state == from)$date,from])
    
    # get median spill window for this state
    winterspill_data <- UCH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter, i] <- exp(b0_array_UCH[from,to,iter] +
                                                       # btemp0_array_UCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_UCH[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_UCH[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_UCH[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_UCH[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_UCH[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_UCH[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_UCH[from,to,iter]*temp_predict[j]*origin2 + 
                                                       # btemp0xorigin3_array_UCH[from,to,iter]*origin3 +
                                                       btemp1xorigin3_array_UCH[from,to,iter]*temp_predict[j]*origin3 +
                                                       borigin1_array_UCH[from,to,iter]*origin1 +
                                                       borigin2_array_UCH[from,to,iter]*origin2 +
                                                       borigin3_array_UCH[from,to,iter]*origin3)/
          sum(exp(b0_array_UCH[from,possible_movements,iter] +
                    # btemp0_array_UCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_UCH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_UCH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_UCH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_UCH[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCH[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    borigin1_array_UCH[from,possible_movements,iter]*origin1 +
                    borigin2_array_UCH[from,possible_movements,iter]*origin2 +
                    borigin3_array_UCH[from,possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_MCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  origin4 = wild_origin_params[4]
  origin5 = wild_origin_params[5]
  origin6 = wild_origin_params[6]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # get median winter spill for this state
    MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                                   date = as.vector(MCW_envir$data$transition_dates))
    spillwindow_data <- MCW_envir$data$spill_window_data
    med_spillwindow <- median(spillwindow_data[subset(MCW_states_dates, state == from)$date,from])
    
    # get median spill window for this state
    winterspill_data <- MCW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter,i] <- exp(b0_array_MCW[from,to,iter] +
                                                       # btemp0_array_MCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_MCW[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_MCW[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_MCW[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_MCW[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_MCW[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_MCW[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_MCW[from,to,iter]*temp_predict[j]*origin2 + 
                                                       # btemp0xorigin3_array_MCW[from,to,iter]*origin3 +
                                                       btemp1xorigin3_array_MCW[from,to,iter]*temp_predict[j]*origin3 +
                                                       # btemp0xorigin4_array_MCW[from,to,iter]*origin4 +
                                                       btemp1xorigin4_array_MCW[from,to,iter]*temp_predict[j]*origin4 +
                                                       # btemp0xorigin5_array_MCW[from,to,iter]*origin5 +
                                                       btemp1xorigin5_array_MCW[from,to,iter]*temp_predict[j]*origin5 +
                                                       # btemp0xorigin6_array_MCW[from,to,iter]*origin6 +
                                                       btemp1xorigin6_array_MCW[from,to,iter]*temp_predict[j]*origin6 +
                                                       borigin1_array_MCW[from,to,iter]*origin1 +
                                                       borigin2_array_MCW[from,to,iter]*origin2 +
                                                       borigin3_array_MCW[from,to,iter]*origin3 +
                                                       borigin4_array_MCW[from,to,iter]*origin4 +
                                                       borigin5_array_MCW[from,to,iter]*origin5 +
                                                       borigin6_array_MCW[from,to,iter]*origin6)/
          sum(exp(b0_array_MCW[from,possible_movements,iter] +
                    # btemp0_array_MCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_MCW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_MCW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_MCW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_MCW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_MCW[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_MCW[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    # btemp0xorigin6_array_MCW[from,possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_MCW[from,possible_movements,iter]*temp_predict[j]*origin6 +
                    borigin1_array_MCW[from,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[from,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[from,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[from,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[from,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[from,possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_MCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # get median winter spill for this state
    MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                                   date = as.vector(MCH_envir$data$transition_dates))
    spillwindow_data <- MCH_envir$data$spill_window_data
    med_spillwindow <- median(spillwindow_data[subset(MCH_states_dates, state == from)$date,from])
    
    # get median spill window for this state
    winterspill_data <- MCH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j,iter,i] <- exp(b0_array_MCH[from,to,iter] +
                                                       # btemp0_array_MCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_MCH[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_MCH[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_MCH[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_MCH[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_MCH[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_MCH[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_MCH[from,to,iter]*temp_predict[j]*origin2 +
                                                       borigin1_array_MCH[from,to,iter]*origin1 +
                                                       borigin2_array_MCH[from,to,iter]*origin2)/
          sum(exp(b0_array_MCH[from,possible_movements,iter] +
                    # btemp0_array_MCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_MCH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_MCH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_MCH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    borigin1_array_MCH[from,possible_movements,iter]*origin1 +
                    borigin2_array_MCH[from,possible_movements,iter]*origin2))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_SRW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,6)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  origin4 = wild_origin_params[4]
  origin5 = wild_origin_params[5]
  origin6 = wild_origin_params[6]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # get median winter spill for this state
    SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
                                   date = as.vector(SRW_envir$data$transition_dates))
    spillwindow_data <- SRW_envir$data$spill_window_data
    med_spillwindow <- median(spillwindow_data[subset(SRW_states_dates, state == from)$date,from])
    
    # get median spill window for this state
    winterspill_data <- SRW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
  
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j,iter,i] <- exp(b0_array_SRW[from,to,iter] +
              # btemp0_array_SRW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
              btemp1_array_SRW[from,to,iter]*temp_predict[j] + 
              bspillwindow_array_SRW[from,to,iter]*med_spillwindow +
              bwinterspill_array_SRW[from,to,iter]*med_winterspill +
              # btemp0xorigin1_array_SRW[from,to,iter]*origin1 +
              btemp1xorigin1_array_SRW[from,to,iter]*temp_predict[j]*origin1 +
              # btemp0xorigin2_array_SRW[from,to,iter]*origin2 + 
              btemp1xorigin2_array_SRW[from,to,iter]*temp_predict[j]*origin2 + 
              # btemp0xorigin3_array_SRW[from,to,iter]*origin3 +
              btemp1xorigin3_array_SRW[from,to,iter]*temp_predict[j]*origin3 +
              # btemp0xorigin4_array_SRW[from,to,iter]*origin4 +
              btemp1xorigin4_array_SRW[from,to,iter]*temp_predict[j]*origin4 +
              # btemp0xorigin5_array_SRW[from,to,iter]*origin5 +
              btemp1xorigin5_array_SRW[from,to,iter]*temp_predict[j]*origin5 +
              # btemp0xorigin6_array_SRW[from,to,iter]*origin6 +
              btemp1xorigin6_array_SRW[from,to,iter]*temp_predict[j]*origin6 +
              borigin1_array_SRW[from,to,iter]*origin1 +
              borigin2_array_SRW[from,to,iter]*origin2 +
              borigin3_array_SRW[from,to,iter]*origin3 +
              borigin4_array_SRW[from,to,iter]*origin4 +
              borigin5_array_SRW[from,to,iter]*origin5 +
              borigin6_array_SRW[from,to,iter]*origin6)/
          sum(exp(b0_array_SRW[from,possible_movements,iter] +
                    # btemp0_array_SRW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRW[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_SRW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_SRW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_SRW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_SRW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_SRW[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_SRW[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    # btemp0xorigin6_array_SRW[from,possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_SRW[from,possible_movements,iter]*temp_predict[j]*origin6 +
                    borigin1_array_SRW[from,possible_movements,iter]*origin1 +
                    borigin2_array_SRW[from,possible_movements,iter]*origin2 +
                    borigin3_array_SRW[from,possible_movements,iter]*origin3 +
                    borigin4_array_SRW[from,possible_movements,iter]*origin4 +
                    borigin5_array_SRW[from,possible_movements,iter]*origin5 +
                    borigin6_array_SRW[from,possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

estimate_temp_effect_SRH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,5)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  origin4 = hatchery_origin_params[4]
  origin5 = hatchery_origin_params[5]
  
  # set up temperature values to predict across
  temp_predict <- seq(-2,2,length = 100)
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different temperatures
  # temp_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(temp_predict), niter))
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  temp_move_prob_array <- array(dim = c(length(temp_predict), niter, nrow(movements)))
  dimnames(temp_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # get median winter spill for this state
    SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
                                   date = as.vector(SRH_envir$data$transition_dates))
    spillwindow_data <- SRH_envir$data$spill_window_data
    med_spillwindow <- median(spillwindow_data[subset(SRH_states_dates, state == from)$date,from])
    
    # get median spill window for this state
    winterspill_data <- SRH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(temp_predict)){
        # evaluation movement 
        temp_move_prob_array[j, iter,i] <- exp(b0_array_SRH[from,to,iter] +
                                                       # btemp0_array_SRH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                       btemp1_array_SRH[from,to,iter]*temp_predict[j] + 
                                                       bspillwindow_array_SRH[from,to,iter]*med_spillwindow +
                                                       bwinterspill_array_SRH[from,to,iter]*med_winterspill +
                                                       # btemp0xorigin1_array_SRH[from,to,iter]*origin1 +
                                                       btemp1xorigin1_array_SRH[from,to,iter]*temp_predict[j]*origin1 +
                                                       # btemp0xorigin2_array_SRH[from,to,iter]*origin2 + 
                                                       btemp1xorigin2_array_SRH[from,to,iter]*temp_predict[j]*origin2 + 
                                                       # btemp0xorigin3_array_SRH[from,to,iter]*origin3 +
                                                       btemp1xorigin3_array_SRH[from,to,iter]*temp_predict[j]*origin3 +
                                                       # btemp0xorigin4_array_SRH[from,to,iter]*origin4 +
                                                       btemp1xorigin4_array_SRH[from,to,iter]*temp_predict[j]*origin4 +
                                                       # btemp0xorigin5_array_SRH[from,to,iter]*origin5 +
                                                       btemp1xorigin5_array_SRH[from,to,iter]*temp_predict[j]*origin5 +
                                                       borigin1_array_SRH[from,to,iter]*origin1 +
                                                       borigin2_array_SRH[from,to,iter]*origin2 +
                                                       borigin3_array_SRH[from,to,iter]*origin3 +
                                                       borigin4_array_SRH[from,to,iter]*origin4 +
                                                       borigin5_array_SRH[from,to,iter]*origin5)/
          sum(exp(b0_array_SRH[from,possible_movements,iter] +
                    # btemp0_array_SRH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRH[from,possible_movements,iter]*temp_predict[j] + 
                    bspillwindow_array_SRH[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_SRH[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRH[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin1 +
                    # btemp0xorigin2_array_SRH[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin2 + 
                    # btemp0xorigin3_array_SRH[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin3 +
                    # btemp0xorigin4_array_SRH[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin4 +
                    # btemp0xorigin5_array_SRH[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRH[from,possible_movements,iter]*temp_predict[j]*origin5 +
                    borigin1_array_SRH[from,possible_movements,iter]*origin1 +
                    borigin2_array_SRH[from,possible_movements,iter]*origin2 +
                    borigin3_array_SRH[from,possible_movements,iter]*origin3 +
                    borigin4_array_SRH[from,possible_movements,iter]*origin4 +
                    borigin5_array_SRH[from,possible_movements,iter]*origin5))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(temp_move_prob_array)
  
}

# plot temp effect function

# create df to index to right dam temp for plotting
dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                        state = seq(2,9))


plot_temp_effect <- function(move_prob_array, from, to, plot_title = NULL){
  temp_move_prob <- as.data.frame(move_prob_array[from, to,,])
  
  colnames(temp_move_prob) <- paste0("iter", 1:niter) 
  temp_predict <- seq(-2,2,length = 100)
  temp_move_prob$temp <- temp_predict
  
  # Add a column with the actual temperatures
  temp_move_prob %>% 
    mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
             window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> temp_move_prob
  
  temp_move_prob %>% 
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
    group_by(temp_actual) %>% 
    summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> temp_move_prob_quantiles
  
  temp_move_prob_plot <- ggplot(temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
    xlab(expression(~"Temperature" ~ ("°C"))) +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(temp_move_prob_plot)
}

# another plot option to compare between hatchery and wild
rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
plot_compare_rear_temp_effect <- function(origin_select,
                                          wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                          wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                          movements_evaluated,
                                          from, to, plot_title = NULL){
  
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  niter <- 4000 # this is the number of draws we have
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])

    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    temp_predict <- seq(-2,2,length = 100)
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    wild_temp_move_prob %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    wild_temp_move_prob_quantiles  -> rear_temp_move_prob_quantiles

  } 
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    hatchery_temp_move_prob %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    # combine wild and hatchery
    hatchery_temp_move_prob_quantiles -> rear_temp_move_prob_quantiles
    
  } 
  # else run both
  else {
    wild_temp_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    colnames(wild_temp_move_prob) <- paste0("iter", 1:niter) 
    temp_predict <- seq(-2,2,length = 100)
    wild_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    wild_temp_move_prob %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> wild_temp_move_prob
    
    wild_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_temp_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> wild_covariate_experiences
    
    hatchery_temp_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    colnames(hatchery_temp_move_prob) <- paste0("iter", 1:niter) 
    hatchery_temp_move_prob$temp <- temp_predict
    
    # Add a column with the actual temperatures
    hatchery_temp_move_prob %>% 
      mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
               window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temp) -> hatchery_temp_move_prob
    
    hatchery_temp_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(temp_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_temp_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> hatchery_covariate_experiences
    
    # combine wild and hatchery
    wild_temp_move_prob_quantiles %>% 
      bind_rows(., hatchery_temp_move_prob_quantiles) -> rear_temp_move_prob_quantiles
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) -> covariate_experiences
    
  }
  
  # convert temp to temp_actual in covariate experiences
  covariate_experiences %>% 
    mutate(temp_actual = window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_mean")] + 
             window_temp_summary[, paste0(subset(dam_index, state == from)$dam, "_sd")]*temperature) -> covariate_experiences
  
  # Remove any instances for rug plot of fish that would have received temp0 covariate
  covariate_experiences %>% 
    dplyr::rename(date_numeric = date) %>% 
    # keep only jan/feb/mar 
    mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
    mutate(month = month(date)) %>% 
    filter(!(month %in% c(1,2,3,4,5))) -> covariate_experiences
  
  
  rear_temp_move_prob_plot <- ggplot(rear_temp_move_prob_quantiles, aes(x = temp_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                        color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_rug(data = subset(covariate_experiences, rear == "wild"), aes(x = temp_actual, color = rear), inherit.aes = FALSE,
             sides = "t", length = unit(0.3, "cm")) +
    geom_rug(data = subset(covariate_experiences, rear == "hatchery"), aes(x = temp_actual, color = rear), inherit.aes = FALSE,
             sides = "b", length = unit(0.3, "cm")) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab(expression(~"Temperature" ~ ("°C"))) +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(rear_temp_move_prob_plot)
}


#### Plot rear type comparison overshoot for all ####

### Upper Columbia ###
# Wenatchee
WEN_movements <- data.frame(from = c(5), to = c(6))
WEN_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Wenatchee River", movements = WEN_movements)
WEN_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Wenatchee River", movements = WEN_movements)
# Wenatchee - get covariate experiences
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")


WEN_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Wenatchee River",
                                                            wild_move_prob_array = WEN_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = WEN_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = WEN_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
                                                            movements_evaluated = WEN_movements,
                                                            from = WEN_movements$from, to = WEN_movements$to, plot_title = "Wenatchee - overshoot RRE")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WEN_compare_overshoot_temp.png"), WEN_compare_overshoot_temp, height = 8, width = 8)

# Entiat
ENT_movements <- data.frame(from = c(6), to = c(7))
ENT_hatchery_temp_move_prob_array <- estimate_temp_effect_UCH(origin_select = "Eniat River", movements = ENT_movements)
ENT_wild_temp_move_prob_array <- estimate_temp_effect_UCW(origin_select = "Eniat River", movements = ENT_movements)

# Entiat - get covariate experiences
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
ENT_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Entiat River")


ENT_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Entiat River",
                                                            wild_move_prob_array = ENT_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = ENT_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = ENT_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = ENT_hatchery_covariate_experiences,
                                                            movements_evaluated = ENT_movements,
                                                            from = ENT_movements$from, to = ENT_movements$to, plot_title = "Entiat - overshoot WEL")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "ENT_compare_overshoot_temp.png"), ENT_compare_overshoot_temp, height = 8, width = 8)






### Middle Columbia ###

# Deschutes River
DES_movements <- data.frame(from = c(2), to = c(3))
DES_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Deschutes River", movements = DES_movements)
DES_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Deschutes River", movements = DES_movements)

# Deschutes - get covariate experiences
DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")
DES_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Deschutes River")


DES_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Deschutes River",
                                                            wild_move_prob_array = DES_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = DES_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = DES_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = DES_hatchery_covariate_experiences,
                                                            movements_evaluated = DES_movements,
                                                            from = DES_movements$from, to = DES_movements$to, plot_title = "Deschutes - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "DES_compare_overshoot_temp.png"), DES_compare_overshoot_temp, height = 8, width = 8)


# John Day River
JDR_movements <- data.frame(from = c(2), to = c(3))
JDR_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "John Day River", movements = JDR_movements)
JDR_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "John Day River", movements = JDR_movements)

# John Day - get covariate experiences
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")
JDR_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "John Day River")


JDR_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "John Day River",
                                                            wild_move_prob_array = JDR_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = JDR_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = JDR_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = JDR_hatchery_covariate_experiences,
                                                            movements_evaluated = JDR_movements,
                                                            from = JDR_movements$from, to = JDR_movements$to, plot_title = "John Day - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "JDR_compare_overshoot_temp.png"), JDR_compare_overshoot_temp, height = 8, width = 8)


# Fifteenmile Creek
FIF_movements <- data.frame(from = c(2), to = c(3))
FIF_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Fifteenmile Creek", movements = FIF_movements)
FIF_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Fifteenmile Creek", movements = FIF_movements)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
FIF_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Fifteenmile Creek")


FIF_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Fifteenmile Creek",
                                                            wild_move_prob_array = FIF_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = FIF_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = FIF_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = FIF_hatchery_covariate_experiences,
                                                            movements_evaluated = FIF_movements,
                                                            from = FIF_movements$from, to = FIF_movements$to, plot_title = "Fifteenmile - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "FIF_compare_overshoot_temp.png"), FIF_compare_overshoot_temp, height = 8, width = 8)


# Umatilla River
UMA_movements <- data.frame(from = c(2), to = c(3))
UMA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Umatilla River", movements = UMA_movements)
UMA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Umatilla River", movements = UMA_movements)

# Umatilla - get covariate experiences
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")


UMA_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Umatilla River",
                                                            wild_move_prob_array = UMA_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = UMA_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = UMA_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                            movements_evaluated = UMA_movements,
                                                            from = UMA_movements$from, to = UMA_movements$to, plot_title = "Umatilla - overshoot MCN")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "UMA_compare_overshoot_temp.png"), UMA_compare_overshoot_temp, height = 8, width = 8)


# Yakima River
YAK_movements <- data.frame(from = c(3), to = c(4))
YAK_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Yakima River", movements = YAK_movements)
YAK_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = YAK_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
YAK_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Yakima River")


YAK_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Yakima River",
                                                            wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = YAK_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = YAK_hatchery_covariate_experiences,
                                                            movements_evaluated = YAK_movements,
                                                            from = YAK_movements$from, to = YAK_movements$to, plot_title = "Yakima - overshoot PRA")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "YAK_compare_overshoot_PRA_temp.png"), YAK_compare_overshoot_temp, height = 8, width = 8)

YAK_movements <- data.frame(from = c(3), to = c(8))
YAK_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Yakima River", movements = YAK_movements)
YAK_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Yakima River", movements = YAK_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
YAK_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Yakima River")


YAK_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Yakima River",
                                                            wild_move_prob_array = YAK_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = YAK_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = YAK_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = YAK_hatchery_covariate_experiences,
                                                            movements_evaluated = YAK_movements,
                                                            from = YAK_movements$from, to = YAK_movements$to, plot_title = "Yakima - overshoot ICH")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "YAK_compare_overshoot_ICH_temp.png"), YAK_compare_overshoot_temp, height = 8, width = 8)


# Walla Walla River
WAWA_movements <- data.frame(from = c(3), to = c(4))
WAWA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")


WAWA_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Walla Walla River",
                                                            wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                            movements_evaluated = WAWA_movements,
                                                            from = WAWA_movements$from, to = WAWA_movements$to, plot_title = "Walla Walla - overshoot PRA")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WAWA_compare_overshoot_PRA_temp.png"), WAWA_compare_overshoot_temp, height = 8, width = 8)


WAWA_movements <- data.frame(from = c(3), to = c(8))
WAWA_hatchery_temp_move_prob_array <- estimate_temp_effect_MCH(origin_select = "Walla Walla River", movements = WAWA_movements)
WAWA_wild_temp_move_prob_array <- estimate_temp_effect_MCW(origin_select = "Walla Walla River", movements = WAWA_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")


WAWA_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Walla Walla River",
                                                            wild_move_prob_array = WAWA_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = WAWA_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = WAWA_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                            movements_evaluated = WAWA_movements,
                                                             from = WAWA_movements$from, to = WAWA_movements$to, plot_title = "Walla Walla - overshoot ICH")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "WAWA_compare_overshoot_ICH_temp.png"), WAWA_compare_overshoot_temp, height = 8, width = 8)








### Snake River ###

## Tucannon River
TUC_movements <- data.frame(from = c(8), to = c(9))
TUC_hatchery_temp_move_prob_array <- estimate_temp_effect_SRH(origin_select = "Tucannon River", movements = TUC_movements)
TUC_wild_temp_move_prob_array <- estimate_temp_effect_SRW(origin_select = "Tucannon River", movements = TUC_movements)

# Tucannon - get covariate experiences
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")


TUC_compare_overshoot_temp <- plot_compare_rear_temp_effect(origin_select = "Tucannon River",
                                                            wild_move_prob_array = TUC_wild_temp_move_prob_array,
                                                            hatchery_move_prob_array = TUC_hatchery_temp_move_prob_array,
                                                            wild_covariate_experiences = TUC_wild_covariate_experiences,
                                                            hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
                                                            movements_evaluated = TUC_movements,
                                                            from = TUC_movements$from, to = TUC_movements$to, plot_title = "Tucannon - overshoot LGR")

ggsave(here::here("stan_actual", "output", "covariate_effects", "temperature", "TUC_compare_overshoot_temp.png"), TUC_compare_overshoot_temp, height = 8, width = 8)



# Snake only has one overshoot - for Tucannon
