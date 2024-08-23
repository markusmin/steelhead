# 14-model-final-fates-covariates

# This script takes the output from the stan model runs in 05-stan-runs and
# generates estimates of final fate distributions under different scenarios of covariates

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
UCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2", "upper_columbia_wild","UCW_transition_counts.csv"))
UCW_movements <- paste0("_", UCW_transition_counts$from, "_", UCW_transition_counts$to)
UCW_movements <- UCW_movements[!(grepl("NA", UCW_movements))]

UCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2", "upper_columbia_hatchery","UCH_transition_counts.csv"))
UCH_movements <- paste0("_", UCH_transition_counts$from, "_", UCH_transition_counts$to)
UCH_movements <- UCH_movements[!(grepl("NA", UCH_movements))]

MCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2", "middle_columbia_wild","MCW_transition_counts.csv"))
MCW_movements <- paste0("_", MCW_transition_counts$from, "_", MCW_transition_counts$to)
MCW_movements <- MCW_movements[!(grepl("NA", MCW_movements))]

MCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2", "middle_columbia_hatchery","MCH_transition_counts.csv"))
MCH_movements <- paste0("_", MCH_transition_counts$from, "_", MCH_transition_counts$to)
MCH_movements <- MCH_movements[!(grepl("NA", MCH_movements))]

SRW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2", "snake_river_wild","SRW_transition_counts.csv"))
SRW_movements <- paste0("_", SRW_transition_counts$from, "_", SRW_transition_counts$to)
SRW_movements <- SRW_movements[!(grepl("NA", SRW_movements))]

SRH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery","SRH_transition_counts.csv"))
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
load(here::here("stan_actual", "reparameterization_v2", "upper_columbia_wild", "model_data.rda"),
                       envir = UCW_envir)
load(here::here("stan_actual", "reparameterization_v2", "upper_columbia_hatchery", "model_data.rda"),
                       envir = UCH_envir)
load(here::here("stan_actual", "reparameterization_v2", "middle_columbia_wild", "model_data.rda"),
     envir = MCW_envir)
load(here::here("stan_actual", "reparameterization_v2", "middle_columbia_hatchery", "model_data.rda"),
     envir = MCH_envir)
load(here::here("stan_actual", "reparameterization_v2", "snake_river_wild", "model_data.rda"),
     envir = SRW_envir)
load(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery", "model_data.rda"),
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
UCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2", "upper_columbia_wild", "chain1_UCW_reparam_v2_fit.rds"))
UCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2", "upper_columbia_wild", "chain2_UCW_reparam_v2_fit.rds"))
UCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2", "upper_columbia_wild", "chain3_UCW_reparam_v2_fit.rds"))
UCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2", "upper_columbia_wild", "chain4_UCW_reparam_v2_fit.rds"))

# bind chains together
UCW_fit_raw <- bind4chains(UCW_chain1, UCW_chain2, UCW_chain3, UCW_chain4)
# thin2
thin_draws(UCW_fit_raw, thin = 2) -> UCW_fit
# summarise
UCW_fit_summary <- summarise_draws(UCW_fit)

## Upper Columbia, Hatchery
UCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2", "upper_columbia_hatchery", "chain1_UCH_reparam_v2_fit.rds"))
UCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2", "upper_columbia_hatchery", "chain2_UCH_reparam_v2_fit.rds"))
UCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2", "upper_columbia_hatchery", "chain3_UCH_reparam_v2_fit.rds"))
UCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2", "upper_columbia_hatchery", "chain4_UCH_reparam_v2_fit.rds"))

# bind chains together
UCH_fit_raw <- bind4chains(UCH_chain1, UCH_chain2, UCH_chain3, UCH_chain4)
# thin2
thin_draws(UCH_fit_raw, thin = 2) -> UCH_fit
# summarise
UCH_fit_summary <- summarise_draws(UCH_fit)

## Middle Columbia, Wild
MCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2", "middle_columbia_wild", "chain1_MCW_reparam_v2_fit.rds"))
MCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2", "middle_columbia_wild", "chain2_MCW_reparam_v2_fit.rds"))
MCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2", "middle_columbia_wild", "chain3_MCW_reparam_v2_fit.rds"))
MCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2", "middle_columbia_wild", "chain4_MCW_reparam_v2_fit.rds"))

# bind chains together
MCW_fit_raw <- bind4chains(MCW_chain1, MCW_chain2, MCW_chain3, MCW_chain4)
# thin2
thin_draws(MCW_fit_raw, thin = 2) -> MCW_fit
# summarise
MCW_fit_summary <- summarise_draws(MCW_fit)

## Middle Columbia, Hatchery
MCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2", "middle_columbia_hatchery", "chain1_MCH_reparam_v2_fit.rds"))
MCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2", "middle_columbia_hatchery", "chain2_MCH_reparam_v2_fit.rds"))
MCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2", "middle_columbia_hatchery", "chain3_MCH_reparam_v2_fit.rds"))
MCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2", "middle_columbia_hatchery", "chain4_MCH_reparam_v2_fit.rds"))

# bind chains together
MCH_fit_raw <- bind4chains(MCH_chain1, MCH_chain2, MCH_chain3, MCH_chain4)
# thin2
thin_draws(MCH_fit_raw, thin = 2) -> MCH_fit
# summarise
MCH_fit_summary <- summarise_draws(MCH_fit)

## Snake River, Wild
SRW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_wild", "chain1_SRW_reparam_v2_fit.rds"))
SRW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_wild", "chain2_SRW_reparam_v2_fit.rds"))
SRW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_wild", "chain3_SRW_reparam_v2_fit.rds"))
SRW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_wild", "chain4_SRW_reparam_v2_fit.rds"))

# bind chains together
SRW_fit_raw <- bind4chains(SRW_chain1, SRW_chain2, SRW_chain3, SRW_chain4)
# thin2
thin_draws(SRW_fit_raw, thin = 2) -> SRW_fit
# summarise
SRW_fit_summary <- summarise_draws(SRW_fit)

## Snake River, Hatchery
SRH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery", "chain1_SRH_reparam_v2_fit.rds"))
SRH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery", "chain2_SRH_reparam_v2_fit.rds"))
SRH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery", "chain3_SRH_reparam_v2_fit.rds"))
SRH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery", "chain4_SRH_reparam_v2_fit.rds"))

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

#### Select representative warm vs. cold years ####

run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
run_year_numeric = seq(4, 22, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)

temp_ts <- as.data.frame(UCW_envir$data$temperature_data)
dates <- seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days")
years <- year(dates)
temp_ts$date <- dates
temp_ts$year <- years

temp_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  # drop 22/23 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "22/23")) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> annual_temp_medians

spill_ts <- as.data.frame(UCW_envir$data$spill_window_data)
dates <- seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days")
years <- year(dates)
spill_ts$date <- dates
spill_ts$year <- years

spill_ts %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  # drop 22/23 because it's not a complete year (and it isn't modeled)
  filter(!(run_year == "22/23")) %>% 
  group_by(run_year) %>% 
  summarise_all(median) %>% 
  arrange(BON) -> annual_spill_medians

median_temp_by_run_year_BON <- ggplot(annual_temp_medians, aes(x = year, y = BON, label = run_year)) +
  geom_point() +
  ylab("Z-scored basin temperature") +
  geom_text(aes(y = BON + 0.01))

ggsave(here::here("stan_actual", "output", "data_plots", "median_temp_by_run_year_BON.png"), median_temp_by_run_year_BON, height = 8, width = 8)

ggplot(annual_spill_medians, aes(x = year, y = BON)) +
  geom_point()

plot(x = annual_temp_medians$BON, y = annual_spill_medians$BON)
# spill volume and temperature are clearly correlated: lower temperature = lower spill
# So you have to plot spill and temp together


# pick out median, 25%, 75%
summary(1:17)
annual_temp_medians[c(1,5,9,13,17),]

# cool year: 08/09
# average year: 19/20
# warm year: 17/18

# coldest year on record: 07/08
# hottest year on record: 15/16

rep_years <- data.frame(fix_run_year = c("07/08", "08/09", "19/20", "17/18", "15/16"),
                        conditions = c("coldest", "cool", "average", "warm", "warmest"))





#### Final fates functions under different covariate values ####

# This function takes a model fit object, a number of simulated fish,
# the identities of those fish (origin, rear), the conditions throughout the
# basin (spill, temperature, year)
# and then simulates the final fates of those fish
# By default, the conditions (spill and temperature) will be taken from the data
# itself, but note that new values of these conditions could be simulated in order

# Note that another option here would be to use average (median) conditions experienced
# in each state - and that would remove the variability associated with taking random 
# conditions. But note that you'll have to decide what the data is that you're taking
# the median of - is it the whole DPS? Just the population? Currently, we have it set up
# to run by DPS, not population
# Another option would be to make the number of fish per simulation smaller and
# run it more times - that would give you a wider spread of covariate values
# and therefore reduce the chance that a handful of outlier covariate values
# lead to crazy distributions


# Final states simulation under different covariate values

# this function will be very similar to the other final fates simulations, but
# with additional arguments that allow you to fix conditions in certain reaches




# This function will have to be re-written for each DPS, because each DPS has different params
# currently we are not including a random effect of year, but we could easily
# amend the code to simulate final fates in a particular year
# origin1/origin2/origin3 are indicator variables, so set the origin you want to be 1
final_fates_simulation_fixcov_UCW <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       fix_state, fix_temp_value = NA, 
                                       fix_spillwindow_value = NA, fix_winterspill_value = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(UCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  
  
  temp_upstream <- rep(0, length(model_states))
  temp_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "upstream")) == 0){
          temp_upstream[i] <- 0
        } else {
          temp_upstream[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "upstream")$date,i])
        }
      } else {
        temp_upstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "downstream")) == 0){
          temp_downstream[i] <- 0
        } else {
          temp_downstream[i] <- median(temp_data[subset(UCW_states_dates, state == i & direction == "downstream")$date,i])
        }
      } else {
        temp_downstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
      }
  }
  
  
  spillwindow_upstream <- rep(0, length(model_states))
  spillwindow_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "upstream")) == 0){
          spillwindow_upstream[i] <- 0
        } else {
          spillwindow_upstream[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "upstream")$date,i])
        }
      } else {
        spillwindow_upstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(UCW_states_dates, state == i & direction == "downstream")) == 0){
          spillwindow_downstream[i] <- 0
        } else {
          spillwindow_downstream[i] <- median(spillwindow_data[subset(UCW_states_dates, state == i & direction == "downstream")$date,i])
        }
      } else {
        spillwindow_downstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
      }

  }
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_value))){
    temp_upstream[fix_state] <- fix_temp_value
    temp_downstream[fix_state] <- fix_temp_value
  }
  
  if(!(is.na(fix_spillwindow_value))){
    spillwindow_upstream[fix_state] <- fix_spillwindow_value
    spillwindow_downstream[fix_state] <- fix_spillwindow_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                  # btemp0_array_UCW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                  btemp1_array_UCW[i,possible_movements[j],iter]*temp_upstream[i] + 
                                                                  bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_upstream[i] + 
                                                                  bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                  # btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                  btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp_upstream[i]*origin1 +
                                                                  # btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*origin2 + 
                                                                  btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp_upstream[i]*origin2 + 
                                                                  # btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*origin3 +
                                                                  btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp_upstream[i]*origin3 +
                                                                  borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                  borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                  borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
        sum(exp(b0_array_UCW[i,possible_movements,iter] +
                  # btemp0_array_UCW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_UCW[i,possible_movements,iter]*temp_upstream[i] + 
                  bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_upstream[i] + 
                  bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_UCW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp_upstream[i]*origin1 +
                  # btemp0xorigin2_array_UCW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp_upstream[i]*origin2 + 
                  # btemp0xorigin3_array_UCW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp_upstream[i]*origin3 +
                  borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                  borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                  borigin3_array_UCW[i,possible_movements,iter]*origin3))
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep("upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                                    # btemp0_array_UCW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                    btemp1_array_UCW[i,possible_movements[j],iter]*temp_downstream[i] + 
                                                                    bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow_downstream[i] + 
                                                                    bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    # btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                    btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp_downstream[i]*origin1 +
                                                                    # btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*origin2 + 
                                                                    btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp_downstream[i]*origin2 + 
                                                                    # btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*origin3 +
                                                                    btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp_downstream[i]*origin3 +
                                                                    borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
        sum(exp(b0_array_UCW[i,possible_movements,iter] +
                  # btemp0_array_UCW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_UCW[i,possible_movements,iter]*temp_downstream[i] + 
                  bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow_downstream[i] + 
                  bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_UCW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp_downstream[i]*origin1 +
                  # btemp0xorigin2_array_UCW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp_downstream[i]*origin2 + 
                  # btemp0xorigin3_array_UCW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp_downstream[i]*origin3 +
                  borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                  borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                  borigin3_array_UCW[i,possible_movements,iter]*origin3))
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_UCH <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(UCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  
  
  temp_upstream <- rep(0, length(model_states))
  temp_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp_upstream[i] <- temp_data[sample_date[i],i]
      temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "upstream")) == 0){
          temp_upstream[i] <- 0
        } else {
          temp_upstream[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "upstream")$date,i])
        }
      } else {
        temp_upstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "downstream")) == 0){
          temp_downstream[i] <- 0
        } else {
          temp_downstream[i] <- median(temp_data[subset(UCH_states_dates, state == i & direction == "downstream")$date,i])
        }
      } else {
        temp_downstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream <- rep(0, length(model_states))
  spillwindow_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "upstream")) == 0){
          spillwindow_upstream[i] <- 0
        } else {
          spillwindow_upstream[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "upstream")$date,i])
        }
      } else {
        spillwindow_upstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(UCH_states_dates, state == i & direction == "downstream")) == 0){
          spillwindow_downstream[i] <- 0
        } else {
          spillwindow_downstream[i] <- median(spillwindow_data[subset(UCH_states_dates, state == i & direction == "downstream")$date,i])
        }
      } else {
        spillwindow_downstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                  # btemp0_array_UCH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                  btemp1_array_UCH[i,possible_movements[j],iter]*temp_upstream[i] + 
                                                                  bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_upstream[i] + 
                                                                  bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                  # btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                  btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp_upstream[i]*origin1 +
                                                                  # btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*origin2 + 
                                                                  btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp_upstream[i]*origin2 + 
                                                                  # btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*origin3 +
                                                                  btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp_upstream[i]*origin3 +
                                                                  borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                  borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                  borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
        sum(exp(b0_array_UCH[i,possible_movements,iter] +
                  # btemp0_array_UCH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_UCH[i,possible_movements,iter]*temp_upstream[i] + 
                  bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_upstream[i] + 
                  bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_UCH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp_upstream[i]*origin1 +
                  # btemp0xorigin2_array_UCH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp_upstream[i]*origin2 + 
                  # btemp0xorigin3_array_UCH[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp_upstream[i]*origin3 +
                  borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                  borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                  borigin3_array_UCH[i,possible_movements,iter]*origin3))
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep("upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                                    # btemp0_array_UCH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                    btemp1_array_UCH[i,possible_movements[j],iter]*temp_downstream[i] + 
                                                                    bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow_downstream[i] + 
                                                                    bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    # btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                    btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp_downstream[i]*origin1 +
                                                                    # btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*origin2 + 
                                                                    btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp_downstream[i]*origin2 + 
                                                                    # btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*origin3 +
                                                                    btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp_downstream[i]*origin3 +
                                                                    borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
        sum(exp(b0_array_UCH[i,possible_movements,iter] +
                  # btemp0_array_UCH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_UCH[i,possible_movements,iter]*temp_downstream[i] + 
                  bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow_downstream[i] + 
                  bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_UCH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp_downstream[i]*origin1 +
                  # btemp0xorigin2_array_UCH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp_downstream[i]*origin2 + 
                  # btemp0xorigin3_array_UCH[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp_downstream[i]*origin3 +
                  borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                  borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                  borigin3_array_UCH[i,possible_movements,iter]*origin3))
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_MCW <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       origin4 = 0, origin5 = 0, origin6 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       fix_state, fix_temp_value = NA, 
                                       fix_spillwindow_value = NA, fix_winterspill_value = NA,
                                       fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(MCW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  
  
  temp_upstream <- rep(0, length(model_states))
  temp_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream")) == 0){
        temp_upstream[i] <- 0
      } else {
        temp_upstream[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "upstream")$date,i])
      }
    } else {
      temp_upstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream")) == 0){
        temp_downstream[i] <- 0
      } else {
        temp_downstream[i] <- median(temp_data[subset(MCW_states_dates, state == i & direction == "downstream")$date,i])
      }
    } else {
      temp_downstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
    }
  }
  
  
  spillwindow_upstream <- rep(0, length(model_states))
  spillwindow_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "upstream")) == 0){
        spillwindow_upstream[i] <- 0
      } else {
        spillwindow_upstream[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "upstream")$date,i])
      }
    } else {
      spillwindow_upstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCW_states_dates, state == i & direction == "downstream")) == 0){
        spillwindow_downstream[i] <- 0
      } else {
        spillwindow_downstream[i] <- median(spillwindow_data[subset(MCW_states_dates, state == i & direction == "downstream")$date,i])
      }
    } else {
      spillwindow_downstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> MCW_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    
    
    
    MCW_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> MCW_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    temp_upstream <- rep(0, length(model_states))
    temp_downstream <- rep(0, length(model_states))
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_month_day, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream")) == 0){
          temp_upstream[i] <- 0
        } else {
          temp_upstream[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "upstream")$date_run_year,i])
        }
      } else {
        temp_upstream[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream")$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream")) == 0){
          temp_downstream[i] <- 0
        } else {
          temp_downstream[i] <- median(temp_data[subset(MCW_states_month_day, state == i & direction == "downstream")$date_run_year,i])
        }
      } else {
        temp_downstream[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream")$date_run_year,i])
      }
    }
    
    spill_window_upstream <- rep(0, length(model_states))
    spill_window_downstream <- rep(0, length(model_states))
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_month_day, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream")) == 0){
          spill_window_upstream[i] <- 0
        } else {
          spill_window_upstream[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "upstream")$date_run_year,i])
        }
      } else {
        spill_window_upstream[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream")$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream")) == 0){
          spill_window_downstream[i] <- 0
        } else {
          spill_window_downstream[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "downstream")$date_run_year,i])
        }
      } else {
        spill_window_downstream[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream")$date_run_year,i])
      }
    }
    
    
  }
  
  
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_value))){
    temp_upstream[fix_state] <- fix_temp_value
    temp_downstream[fix_state] <- fix_temp_value
  }
  
  if(!(is.na(fix_spillwindow_value))){
    spillwindow_upstream[fix_state] <- fix_spillwindow_value
    spillwindow_downstream[fix_state] <- fix_spillwindow_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  

  
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                             nrow = length(model_states),
                             ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                         # btemp0_array_MCW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                         btemp1_array_MCW[i,possible_movements[j],iter]*temp_upstream[i] + 
                                                         bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_upstream[i] + 
                                                         bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                         # btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                         btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp_upstream[i]*origin1 +
                                                         # btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*origin2 + 
                                                         btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp_upstream[i]*origin2 + 
                                                         # btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                         btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp_upstream[i]*origin3 +
                                                         # btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                         btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp_upstream[i]*origin4 +
                                                         # btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                         btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp_upstream[i]*origin5 +
                                                         # btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*origin6 +
                                                         btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp_upstream[i]*origin6 +
                                                         borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                         borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                         borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                         borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                         borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                         borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
        sum(exp(b0_array_MCW[i,possible_movements,iter] +
                  # btemp0_array_MCW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_MCW[i,possible_movements,iter]*temp_upstream[i] + 
                  bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_upstream[i] + 
                  bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_MCW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp_upstream[i]*origin1 +
                  # btemp0xorigin2_array_MCW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp_upstream[i]*origin2 + 
                  # btemp0xorigin3_array_MCW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp_upstream[i]*origin3 +
                  # btemp0xorigin4_array_MCW[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp_upstream[i]*origin4 +
                  # btemp0xorigin5_array_MCW[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp_upstream[i]*origin5 +
                  # btemp0xorigin6_array_MCW[i,possible_movements,iter]*origin6 +
                  btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp_upstream[i]*origin6 +
                  borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                  borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                  borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                  borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                  borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                  borigin6_array_MCW[i,possible_movements,iter]*origin6))
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep("upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                                  # btemp0_array_MCW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                  btemp1_array_MCW[i,possible_movements[j],iter]*temp_downstream[i] + 
                                                                  bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow_downstream[i] + 
                                                                  bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                                  # btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                  btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp_downstream[i]*origin1 +
                                                                  # btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*origin2 + 
                                                                  btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp_downstream[i]*origin2 + 
                                                                  # btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                  btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp_downstream[i]*origin3 +
                                                                  # btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                  btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp_downstream[i]*origin4 +
                                                                  # btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                  btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp_downstream[i]*origin5 +
                                                                  # btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*origin6 +
                                                                  btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp_downstream[i]*origin6 +
                                                                  borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                                  borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                                  borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                                  borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                                  borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                                  borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
        sum(exp(b0_array_MCW[i,possible_movements,iter] +
                  # btemp0_array_MCW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_MCW[i,possible_movements,iter]*temp_downstream[i] + 
                  bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow_downstream[i] + 
                  bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_MCW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp_downstream[i]*origin1 +
                  # btemp0xorigin2_array_MCW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp_downstream[i]*origin2 + 
                  # btemp0xorigin3_array_MCW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp_downstream[i]*origin3 +
                  # btemp0xorigin4_array_MCW[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp_downstream[i]*origin4 +
                  # btemp0xorigin5_array_MCW[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp_downstream[i]*origin5 +
                  # btemp0xorigin6_array_MCW[i,possible_movements,iter]*origin6 +
                  btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp_downstream[i]*origin6 +
                  borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                  borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                  borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                  borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                  borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                  borigin6_array_MCW[i,possible_movements,iter]*origin6))
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  

  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_fixcov_MCH <- function(nsim,
                                              start_state = 2, states_dates,
                                              origin1 = 0, origin2 = 0, 
                                              temp_data, spillwindow_data, winterspill_data,
                                              fix_state, fix_temp_value = NA, 
                                              fix_spillwindow_value = NA, fix_winterspill_value = NA,
                                              fix_run_year = NA){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(MCH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state across all years
    winterspill[i] <- median(winterspill_data[,i])
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  
  
  temp_upstream <- rep(0, length(model_states))
  temp_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream")) == 0){
        temp_upstream[i] <- 0
      } else {
        temp_upstream[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "upstream")$date,i])
      }
    } else {
      temp_upstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream")) == 0){
        temp_downstream[i] <- 0
      } else {
        temp_downstream[i] <- median(temp_data[subset(MCH_states_dates, state == i & direction == "downstream")$date,i])
      }
    } else {
      temp_downstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
    }
  }
  
  
  spillwindow_upstream <- rep(0, length(model_states))
  spillwindow_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    # approach 2: Take the median conditions for each state
    if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "upstream")) == 0){
        spillwindow_upstream[i] <- 0
      } else {
        spillwindow_upstream[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "upstream")$date,i])
      }
    } else {
      spillwindow_upstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
    }
    if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
      if(nrow(subset(MCH_states_dates, state == i & direction == "downstream")) == 0){
        spillwindow_downstream[i] <- 0
      } else {
        spillwindow_downstream[i] <- median(spillwindow_data[subset(MCH_states_dates, state == i & direction == "downstream")$date,i])
      }
    } else {
      spillwindow_downstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
    }
    
  }
  
  # get temperature and spill values for a specific run year
  if (!(is.na(fix_run_year))){
    # get the date numeric for each day/month of that run year
    run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
    run_year_start <- seq(ymd("2004-06-01"), ymd("2022-06-01"), by = "years")
    run_year_end <- seq(ymd("2005-05-31"), ymd("2023-05-31"), by = "years")
    run_year_numeric = seq(4, 22, 1)
    
    run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)
    
    
    
    # rowwise() %>%
    # dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_actual & run_year_end >= date_actual)$run_year) %>% 
    # filter(run_year == fix_run_year)-> MCH_states_dates_run_year
    
    # subset states_dates for origin and population to pull out just month day, and then
    # use those months and days to index temps from that year
    dates_actual <- data.frame(date_actual = seq(ymd("2005-06-01"), ymd("2022-12-31"), by = "days"),
                               date = 1:6423)
    dates_actual %>% 
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) -> dates_actual
    
    
    subset(run_year_df, run_year == fix_run_year) -> fix_run_year_df
    subset(dates_actual, date_actual <= fix_run_year_df$run_year_end & date_actual >= fix_run_year_df$run_year_start) %>% 
      dplyr::rename(date_run_year = date) -> dates_fix_run_year
    
    # Deal with leap years by using 2/28 twice
    if(!("02-29" %in% dates_fix_run_year$month_day)){
      subset(dates_fix_run_year, month_day == "02-28") %>% 
        mutate(month_day = "02-29") -> leap_year
      
      dates_fix_run_year %>% 
        bind_rows(., leap_year) -> dates_fix_run_year
      
    }
    
    
    MCH_states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> MCH_states_month_day
    
    states_dates %>% 
      left_join(dates_actual, by = "date") %>% 
      # extract just the day of year for transitions
      mutate(month_day = format(as.Date(date_actual), "%m-%d")) %>% 
      # turn this into month/day for that run year
      left_join(., dates_fix_run_year, by = "month_day") -> states_month_day
    
    temp_upstream <- rep(0, length(model_states))
    temp_downstream <- rep(0, length(model_states))
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_month_day, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "upstream")) == 0){
          temp_upstream[i] <- 0
        } else {
          temp_upstream[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "upstream")$date_run_year,i])
        }
      } else {
        temp_upstream[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "upstream")$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(MCH_states_month_day, state == i & direction == "downstream")) == 0){
          temp_downstream[i] <- 0
        } else {
          temp_downstream[i] <- median(temp_data[subset(MCH_states_month_day, state == i & direction == "downstream")$date_run_year,i])
        }
      } else {
        temp_downstream[i] <- median(temp_data[subset(states_month_day, state == i  & direction == "downstream")$date_run_year,i])
      }
    }
    
    spill_window_upstream <- rep(0, length(model_states))
    spill_window_downstream <- rep(0, length(model_states))
    for (i in 1:9){ #1:9 because 9 mainstem states
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_month_day, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "upstream")) == 0){
          spill_window_upstream[i] <- 0
        } else {
          spill_window_upstream[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "upstream")$date_run_year,i])
        }
      } else {
        spill_window_upstream[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "upstream")$date_run_year,i])
      }
      if(nrow(subset(states_month_day, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(MCW_states_month_day, state == i & direction == "downstream")) == 0){
          spill_window_downstream[i] <- 0
        } else {
          spill_window_downstream[i] <- median(spillwindow_data[subset(MCW_states_month_day, state == i & direction == "downstream")$date_run_year,i])
        }
      } else {
        spill_window_downstream[i] <- median(spillwindow_data[subset(states_month_day, state == i  & direction == "downstream")$date_run_year,i])
      }
    }
    
    
  }
  
  
  
  
  # Fix covariate values to the ones that you are setting
  if(!(is.na(fix_temp_value))){
    temp_upstream[fix_state] <- fix_temp_value
    temp_downstream[fix_state] <- fix_temp_value
  }
  
  if(!(is.na(fix_spillwindow_value))){
    spillwindow_upstream[fix_state] <- fix_spillwindow_value
    spillwindow_downstream[fix_state] <- fix_spillwindow_value
  }
  
  if(!(is.na(fix_winterspill_value))){
    winterspill[fix_state] <- fix_winterspill_value
  }
  
  
  
  
  
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                  # btemp0_array_MCH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                  btemp1_array_MCH[i,possible_movements[j],iter]*temp_upstream[i] + 
                                                                  bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_upstream[i] + 
                                                                  bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                  # btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                  btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp_upstream[i]*origin1 +
                                                                  # btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*origin2 + 
                                                                  btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp_upstream[i]*origin2 + 
                                                                  borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                  borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
        sum(exp(b0_array_MCH[i,possible_movements,iter] +
                  # btemp0_array_MCH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_MCH[i,possible_movements,iter]*temp_upstream[i] + 
                  bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_upstream[i] + 
                  bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_MCH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp_upstream[i]*origin1 +
                  # btemp0xorigin2_array_MCH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp_upstream[i]*origin2 + 
                  borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                  borigin2_array_MCH[i,possible_movements,iter]*origin2))
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep("upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                                    # btemp0_array_MCH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                    btemp1_array_MCH[i,possible_movements[j],iter]*temp_downstream[i] + 
                                                                    bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow_downstream[i] + 
                                                                    bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    # btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                    btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp_downstream[i]*origin1 +
                                                                    # btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*origin2 + 
                                                                    btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp_downstream[i]*origin2 + 
                                                                    borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
        sum(exp(b0_array_MCH[i,possible_movements,iter] +
                  # btemp0_array_MCH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_MCH[i,possible_movements,iter]*temp_downstream[i] + 
                  bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow_downstream[i] + 
                  bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_MCH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp_downstream[i]*origin1 +
                  # btemp0xorigin2_array_MCH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp_downstream[i]*origin2 + 
                  borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                  borigin2_array_MCH[i,possible_movements,iter]*origin2))
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_SRW <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       origin4 = 0, origin5 = 0, origin6 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(SRW_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  
  
  temp_upstream <- rep(0, length(model_states))
  temp_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp_upstream[i] <- temp_data[sample_date[i],i]
      temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "upstream")) == 0){
          temp_upstream[i] <- 0
        } else {
          temp_upstream[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "upstream")$date,i])
        }
      } else {
        temp_upstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "downstream")) == 0){
          temp_downstream[i] <- 0
        } else {
          temp_downstream[i] <- median(temp_data[subset(SRW_states_dates, state == i & direction == "downstream")$date,i])
        }
      } else {
        temp_downstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream <- rep(0, length(model_states))
  spillwindow_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "upstream")) == 0){
          spillwindow_upstream[i] <- 0
        } else {
          spillwindow_upstream[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "upstream")$date,i])
        }
      } else {
        spillwindow_upstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(SRW_states_dates, state == i & direction == "downstream")) == 0){
          spillwindow_downstream[i] <- 0
        } else {
          spillwindow_downstream[i] <- median(spillwindow_data[subset(SRW_states_dates, state == i & direction == "downstream")$date,i])
        }
      } else {
        spillwindow_downstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                  # btemp0_array_SRW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                  btemp1_array_SRW[i,possible_movements[j],iter]*temp_upstream[i] + 
                                                                  bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_upstream[i] + 
                                                                  bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                  # btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                  btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp_upstream[i]*origin1 +
                                                                  # btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*origin2 + 
                                                                  btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp_upstream[i]*origin2 + 
                                                                  # btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                  btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp_upstream[i]*origin3 +
                                                                  # btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                  btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp_upstream[i]*origin4 +
                                                                  # btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                  btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp_upstream[i]*origin5 +
                                                                  # btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*origin6 +
                                                                  btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp_upstream[i]*origin6 +
                                                                  borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                  borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                  borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                  borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                  borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                  borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
        sum(exp(b0_array_SRW[i,possible_movements,iter] +
                  # btemp0_array_SRW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_SRW[i,possible_movements,iter]*temp_upstream[i] + 
                  bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_upstream[i] + 
                  bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_SRW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp_upstream[i]*origin1 +
                  # btemp0xorigin2_array_SRW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp_upstream[i]*origin2 + 
                  # btemp0xorigin3_array_SRW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp_upstream[i]*origin3 +
                  # btemp0xorigin4_array_SRW[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp_upstream[i]*origin4 +
                  # btemp0xorigin5_array_SRW[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp_upstream[i]*origin5 +
                  # btemp0xorigin6_array_SRW[i,possible_movements,iter]*origin6 +
                  btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp_upstream[i]*origin6 +
                  borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                  borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                  borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                  borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                  borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                  borigin6_array_SRW[i,possible_movements,iter]*origin6))
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep("upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                                    # btemp0_array_SRW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                    btemp1_array_SRW[i,possible_movements[j],iter]*temp_downstream[i] + 
                                                                    bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow_downstream[i] + 
                                                                    bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                                    # btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                    btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp_downstream[i]*origin1 +
                                                                    # btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*origin2 + 
                                                                    btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp_downstream[i]*origin2 + 
                                                                    # btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                    btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp_downstream[i]*origin3 +
                                                                    # btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                    btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp_downstream[i]*origin4 +
                                                                    # btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                    btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp_downstream[i]*origin5 +
                                                                    # btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*origin6 +
                                                                    btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp_downstream[i]*origin6 +
                                                                    borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                                    borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
        sum(exp(b0_array_SRW[i,possible_movements,iter] +
                  # btemp0_array_SRW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_SRW[i,possible_movements,iter]*temp_downstream[i] + 
                  bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow_downstream[i] + 
                  bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_SRW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp_downstream[i]*origin1 +
                  # btemp0xorigin2_array_SRW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp_downstream[i]*origin2 + 
                  # btemp0xorigin3_array_SRW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp_downstream[i]*origin3 +
                  # btemp0xorigin4_array_SRW[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp_downstream[i]*origin4 +
                  # btemp0xorigin5_array_SRW[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp_downstream[i]*origin5 +
                  # btemp0xorigin6_array_SRW[i,possible_movements,iter]*origin6 +
                  btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp_downstream[i]*origin6 +
                  borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                  borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                  borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                  borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                  borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                  borigin6_array_SRW[i,possible_movements,iter]*origin6))
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_SRH <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       origin4 = 0, origin5 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # for all conditions:
  # if there are ever any states that aren't visited by this particular origin/rear
  # combo, then borrow states_dates from the DPS
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    if(nrow(subset(states_dates, state == i)) == 0){
      sample_date[i] <- sample(subset(SRH_states_dates, state == i)$date, 1)
    } else{
      sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
    }
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)+1
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      winterspill[i] <- winterspill_data[sample_year[i],i] 
    } else {
      # approach 2: Take the median conditions for each state across all years
      winterspill[i] <- median(winterspill_data[,i])
    }
    
    
  }
  
  
  # get temp and spill window data
  # 2024-07-17 edit: We are now going to have temperature and spill window data change in
  # the simulation depending on which direction the fish is coming from. Otherwise,
  # the median will in effect default to the upstream movement conditions (higher temp)
  # which was resulting in overestimated overshoot final fates
  
  
  temp_upstream <- rep(0, length(model_states))
  temp_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp_upstream[i] <- temp_data[sample_date[i],i]
      temp_downstream[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "upstream")) == 0){
          temp_upstream[i] <- 0
        } else {
          temp_upstream[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "upstream")$date,i])
        }
      } else {
        temp_upstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "downstream")) == 0){
          temp_downstream[i] <- 0
        } else {
          temp_downstream[i] <- median(temp_data[subset(SRH_states_dates, state == i & direction == "downstream")$date,i])
        }
      } else {
        temp_downstream[i] <- median(temp_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
      }
      
    }
  }
  
  
  spillwindow_upstream <- rep(0, length(model_states))
  spillwindow_downstream <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow_upstream[i] <- spillwindow_data[sample_date[i],i]
      spillwindow_downstream[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      if(nrow(subset(states_dates, state == i  & direction == "upstream")) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "upstream")) == 0){
          spillwindow_upstream[i] <- 0
        } else {
          spillwindow_upstream[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "upstream")$date,i])
        }
      } else {
        spillwindow_upstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "upstream")$date,i])
      }
      if(nrow(subset(states_dates, state == i  & direction == "downstream")) == 0){
        if(nrow(subset(SRH_states_dates, state == i & direction == "downstream")) == 0){
          spillwindow_downstream[i] <- 0
        } else {
          spillwindow_downstream[i] <- median(spillwindow_data[subset(SRH_states_dates, state == i & direction == "downstream")$date,i])
        }
      } else {
        spillwindow_downstream[i] <- median(spillwindow_data[subset(states_dates, state == i  & direction == "downstream")$date,i])
      }
      
    }
  }
  
  
  # make an upstream and a downstream prob matrix
  
  upstream_move_prob_matrix <- matrix(data = 0,
                                      nrow = length(model_states),
                                      ncol = length(model_states))
  
  for (i in 1:nrow(upstream_move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      upstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                  # btemp0_array_SRH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                  btemp1_array_SRH[i,possible_movements[j],iter]*temp_upstream[i] + 
                                                                  bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_upstream[i] + 
                                                                  bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                  # btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                  btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp_upstream[i]*origin1 +
                                                                  # btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*origin2 + 
                                                                  btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp_upstream[i]*origin2 + 
                                                                  # btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                  btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp_upstream[i]*origin3 +
                                                                  # btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                  btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp_upstream[i]*origin4 +
                                                                  # btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*origin5 +
                                                                  btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp_upstream[i]*origin5 +
                                                                  borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                  borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                  borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                  borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                  borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
        sum(exp(b0_array_SRH[i,possible_movements,iter] +
                  # btemp0_array_SRH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_SRH[i,possible_movements,iter]*temp_upstream[i] + 
                  bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_upstream[i] + 
                  bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_SRH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp_upstream[i]*origin1 +
                  # btemp0xorigin2_array_SRH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp_upstream[i]*origin2 + 
                  # btemp0xorigin3_array_SRH[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp_upstream[i]*origin3 +
                  # btemp0xorigin4_array_SRH[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp_upstream[i]*origin4 +
                  # btemp0xorigin5_array_SRH[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp_upstream[i]*origin5 +
                  borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                  borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                  borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                  borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                  borigin5_array_SRH[i,possible_movements,iter]*origin5))
    }
    
  }
  
  downstream_move_prob_matrix <- matrix(data = 0,
                                        nrow = length(model_states),
                                        ncol = length(model_states))
  
  for (i in 1:nrow(downstream_move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep("upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    for (j in 1:length(possible_movements)){
      downstream_move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                                    # btemp0_array_SRH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                                    btemp1_array_SRH[i,possible_movements[j],iter]*temp_downstream[i] + 
                                                                    bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow_downstream[i] + 
                                                                    bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                                    # btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                    btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp_downstream[i]*origin1 +
                                                                    # btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*origin2 + 
                                                                    btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp_downstream[i]*origin2 + 
                                                                    # btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                    btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp_downstream[i]*origin3 +
                                                                    # btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                    btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp_downstream[i]*origin4 +
                                                                    # btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*origin5 +
                                                                    btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp_downstream[i]*origin5 +
                                                                    borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                                    borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                                    borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                                    borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                                    borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
        sum(exp(b0_array_SRH[i,possible_movements,iter] +
                  # btemp0_array_SRH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_SRH[i,possible_movements,iter]*temp_downstream[i] + 
                  bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow_downstream[i] + 
                  bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_SRH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp_downstream[i]*origin1 +
                  # btemp0xorigin2_array_SRH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp_downstream[i]*origin2 + 
                  # btemp0xorigin3_array_SRH[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp_downstream[i]*origin3 +
                  # btemp0xorigin4_array_SRH[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp_downstream[i]*origin4 +
                  # btemp0xorigin5_array_SRH[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp_downstream[i]*origin5 +
                  borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                  borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                  borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                  borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                  borigin5_array_SRH[i,possible_movements,iter]*origin5))
    }
    
  }
  
  # version 2 of while loop - tracking individual fish
  # Now, use the movement probability matrices to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  # differentiate between fish coming from upstream and fish coming from downstream
  upstream_state_matrix <- state_matrix
  downstream_state_matrix <- state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  # they all are going upstream to start
  upstream_state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  
  
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states to both versions of the state matrix
    upstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> upstream_state_matrix
    
    downstream_state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> downstream_state_matrix
    
    # Now, calculate the probabilities from each starting state
    # do this once for fish moving upstream and once for fish moving downstream
    upstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(upstream_movement_matrix) <- model_states
    colnames(upstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      upstream_movement_matrix[,j] <- rmultinom(n = 1, size = upstream_state_matrix[j,i], prob = upstream_move_prob_matrix[j,])
    }
    
    downstream_movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(downstream_movement_matrix) <- model_states
    colnames(downstream_movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      downstream_movement_matrix[,j] <- rmultinom(n = 1, size = downstream_state_matrix[j,i], prob = downstream_move_prob_matrix[j,])
    }
    
    for (k in 1:ncol(upstream_movement_matrix)){
      for (j in 1:nrow(upstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + upstream_movement_matrix[j,k]
        }
      }
    }
    
    for (k in 1:ncol(downstream_movement_matrix)){
      for (j in 1:nrow(downstream_movement_matrix)){
        if (j < k) { # if the row (to) is less than the column (from), then the fish is moving downstream
          downstream_state_matrix[j, i+1] <- downstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        } else { # if the row (to) is NOT less than the column (from), then the fish is moving upstream
          upstream_state_matrix[j, i+1] <- upstream_state_matrix[j, i+1] + downstream_movement_matrix[j,k]
        }
      }
    }
    
    # add loss row to final fate from each
    final_fate_matrix %>% 
      rbind(., upstream_movement_matrix['loss',] + downstream_movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    # if 
    if (upstream_state_matrix[length(model_states), ncol(upstream_state_matrix)] + 
        downstream_state_matrix[length(model_states), ncol(downstream_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, upstream_move_prob_matrix, downstream_move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

#### Run final fates simulation - H v W comparisons ####


# In order to simulate covariate values, we are going to determine the dates
# where the fish were in each state, and then sample from those dates to 
# get a representative value for temperature/spill when fish are in those states
get_states_dates_direction <- function(envir){
  transitions <- envir$data$y
  transition_dates <- envir$data$transition_dates
  
  state_history <- as.vector(t(transitions))
  transition_date_history <- as.vector(t(transition_dates))
  
  # drop the non-data
  nostate <- which(state_history == 0)
  state_history <- state_history[-nostate]
  transition_date_history <- transition_date_history[-nostate]
  
  transitions_df <- data.frame(state = rep(NA, length(state_history)-1),
                               previous_state = rep(NA, length(state_history)-1),
                               date = rep(NA, length(state_history)-1))
  
  
  for (i in 1:nrow(transitions_df)){
    if (i == 1){ # first observation has to come from the ether
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- 0
      transitions_df$date[i] <- transition_date_history[i]
    } else {
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- ifelse(state_history[i-1]==43, 0, state_history[i-1]) # all first detections come from state 0 (fish purgatory)
      transitions_df$date[i] <- transition_date_history[i]
    }
  }
  
  # drop the non-data
  transitions_df %>% 
    filter(!(state %in% c(0,43))) -> transitions_df
  
  # determine if the fish is going upstream or downstream
  transitions_df %>% 
    mutate(direction = ifelse(state > previous_state, "upstream", "downstream")) -> transitions_df
  
  return(transitions_df)
}


UCW_states_dates <- get_states_dates_direction(envir = UCW_envir)
UCH_states_dates <- get_states_dates_direction(envir = UCH_envir)
MCW_states_dates <- get_states_dates_direction(envir = MCW_envir)
MCH_states_dates <- get_states_dates_direction(envir = MCH_envir)
SRW_states_dates <- get_states_dates_direction(envir = SRW_envir)
SRH_states_dates <- get_states_dates_direction(envir = SRH_envir)

# UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
#                                date = as.vector(UCW_envir$data$transition_dates))
# UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
#                                date = as.vector(UCH_envir$data$transition_dates))
# MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
#                                date = as.vector(MCW_envir$data$transition_dates))
# MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
#                                date = as.vector(MCH_envir$data$transition_dates))
# SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
#                                date = as.vector(SRW_envir$data$transition_dates))
# SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
#                                date = as.vector(SRH_envir$data$transition_dates))

# Get states/dates by origin as well
get_origin_states_dates <- function(envir, origin_select, rear){
  
  transitions <- envir$data$y
  transition_dates <- envir$data$transition_dates
  
  state_history <- as.vector(t(transitions))
  transition_date_history <- as.vector(t(transition_dates))
  
  # get origin info
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  origin_vector_matched <- rep(origin_vector, each = ncol(envir$data$y))
  
  # Now, keep only the origin selected
  if(rear == "wild"){
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
  } else {
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
  }
  
  origin_to_keep <- which(origin_vector_matched == origin_numeric)
  
  state_history <- state_history[origin_to_keep]
  transition_date_history <- transition_date_history[origin_to_keep]
  
  # drop the non-data
  nostate <- which(state_history == 0)
  state_history <- state_history[-nostate]
  transition_date_history <- transition_date_history[-nostate]
  
  transitions_df <- data.frame(state = rep(NA, length(state_history)-1),
                               previous_state = rep(NA, length(state_history)-1),
                               date = rep(NA, length(state_history)-1))
  
  for (i in 1:nrow(transitions_df)){
    if (i == 1){ # first observation has to come from the ether
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- 0
      transitions_df$date[i] <- transition_date_history[i]
    } else {
      transitions_df$state[i] <- state_history[i]
      transitions_df$previous_state[i] <- ifelse(state_history[i-1]==43, 0, state_history[i-1]) # all first detections come from state 0 (fish purgatory)
      transitions_df$date[i] <- transition_date_history[i]
    }
  }
  
  # drop the non-data
  transitions_df %>% 
    filter(!(state %in% c(0,43))) -> transitions_df
  
  # determine if the fish is going upstream or downstream
  transitions_df %>% 
    mutate(direction = ifelse(state > previous_state, "upstream", "downstream")) -> transitions_df
  
  return(transitions_df)
  
}

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


# wen_wild_states_dates <- get_origin_states_dates(envir = UCW_envir, origin_select = "Wenatchee River", rear = "wild")


compare_final_fate_rear_type_UC <- function(niter, nsim, condition_jitter,
                                            origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){

    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- c(0,0,0)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_UCW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                             winterspill_data = UCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
  
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- c(0,0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_UCH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3],
                                                 temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                                 winterspill_data = UCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
  
  
  
  ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
  ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- c(0,0,0)
  hatchery_origin_params <- c(0,0,0)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  for(i in 1:niter) {
    # Run final fates simulation for wild
    sim_wild <- final_fates_simulation_UCW(nsim = nsim,
                                      start_state = 2, states_dates = get_origin_states_dates(envir = UCW_envir, origin_select = origin_select, rear = "wild"),
                                      origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                      temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                      winterspill_data = UCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
    ff_wild %>% 
      bind_cols(sim_wild[[2]]) -> ff_wild
    
    # Run final fates simulation for hatchery
    sim_hatchery <- final_fates_simulation_UCH(nsim = nsim,
                                      start_state = 2, states_dates = get_origin_states_dates(envir = UCH_envir, origin_select = origin_select, rear = "hatchery"),
                                      origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], origin3 = hatchery_origin_params[3],
                                      temp_data = UCH_envir$data$temperature_data, spillwindow_data = UCH_envir$data$spill_window_data, 
                                      winterspill_data = UCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
    ff_hatchery %>% 
      bind_cols(sim_hatchery[[2]]) -> ff_hatchery
  }
  
  
  # Reformat final fates simulation for wild
  rownames(ff_wild) <- NULL
  column_to_rownames(ff_wild, "state") -> ff_wild
  colnames(ff_wild) <- paste0("iter", 1:niter)
  
  rownames_to_column(ff_wild, "state") -> ff_wild
  
  ff_wild %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> ff_wild
  
  ff_wild %>% 
    group_by(state) %>% 
    summarise(across(where(is.numeric), sum)) -> ff_wild
  
  ff_wild %>%
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
  
  ff_wild_long %>% 
    mutate(prop = count/nsim) %>% 
    group_by(state) %>%
    summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) %>% 
    mutate(rear_type = "wild") -> ff_wild_quantiles
  
  # Reformat final fates simulation for hatchery
  rownames(ff_hatchery) <- NULL
  column_to_rownames(ff_hatchery, "state") -> ff_hatchery
  colnames(ff_hatchery) <- paste0("iter", 1:niter)
  
  rownames_to_column(ff_hatchery, "state") -> ff_hatchery
  
  ff_hatchery %>% 
    mutate(state = gsub(" Upstream", "", state)) %>% 
    mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
  
  ff_hatchery %>% 
    group_by(state) %>% 
    summarise(across(where(is.numeric), sum)) -> ff_hatchery
  
  ff_hatchery %>%
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
  
  ff_hatchery_long %>% 
    mutate(prop = count/nsim) %>% 
    group_by(state) %>%
    summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) %>% 
    mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
  
  ff_wild_quantiles %>% 
    bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_fixcov_rear_type_MC <- function(niter, nsim,
                                            origin_select,
                                            fix_state, fix_temp_value = NA, 
                                            fix_spillwindow_value = NA, fix_winterspill_value = NA,
                                            fix_run_year = NA){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_MCW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                             winterspill_data = MCW_envir$data$winter_spill_days_data,
                                             fix_state = fix_state, fix_temp_value = fix_temp_value, 
                                             fix_spillwindow_value = fix_spillwindow_value, fix_winterspill_value = fix_winterspill_value,
                                             fix_run_year = fix_run_year)
      

        
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_MCH(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                    origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                    temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                    winterspill_data = MCH_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_value = fix_temp_value, 
                                                    fix_spillwindow_value = fix_spillwindow_value, fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- c(0,0)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_fixcov_MCW(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = MCW_envir, origin_select = origin_select, rear = "wild"),
                                                    origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                                    origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                                    temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                                    winterspill_data = MCW_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_value = fix_temp_value, 
                                                    fix_spillwindow_value = fix_spillwindow_value, fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_fixcov_MCH(nsim = nsim,
                                                    start_state = 2, states_dates = get_origin_states_dates(envir = MCH_envir, origin_select = origin_select, rear = "hatchery"),
                                                    origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                    temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                    winterspill_data = MCH_envir$data$winter_spill_days_data,
                                                    fix_state = fix_state, fix_temp_value = fix_temp_value, 
                                                    fix_spillwindow_value = fix_spillwindow_value, fix_winterspill_value = fix_winterspill_value,
                                                    fix_run_year = fix_run_year)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}

compare_final_fate_rear_type_SR <- function(niter, nsim, condition_jitter,
                                            origin_select){
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    
    # index to the right origin param to turn it on
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_SRW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                             winterspill_data = SRW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    ff_wild_quantiles -> ff_rear_quantiles
    
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_SRH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                 origin5 = hatchery_origin_params[5],
                                                 temp_data = SRH_envir$data$temperature_data, spillwindow_data = SRH_envir$data$spill_window_data, 
                                                 winterspill_data = SRH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_hatchery_quantiles -> ff_rear_quantiles
  } 
  # else run both
  else {
    
    
    
    ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
    ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
    
    # set up vectors to control which origin parameters are turned on
    wild_origin_params <- rep(0,6)
    hatchery_origin_params <- rep(0,5)
    
    # index to the right origin param to turn it on
    hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
    wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
    
    for(i in 1:niter) {
      # Run final fates simulation for wild
      sim_wild <- final_fates_simulation_SRW(nsim = nsim,
                                             start_state = 2, states_dates = get_origin_states_dates(envir = SRW_envir, origin_select = origin_select, rear = "wild"),
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                             winterspill_data = SRW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_SRH(nsim = nsim,
                                                 start_state = 2, states_dates = get_origin_states_dates(envir = SRH_envir, origin_select = origin_select, rear = "hatchery"),
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 origin3 = hatchery_origin_params[3], origin4 = hatchery_origin_params[4], 
                                                 origin5 = hatchery_origin_params[5],
                                                 temp_data = SRH_envir$data$temperature_data, spillwindow_data = SRH_envir$data$spill_window_data, 
                                                 winterspill_data = SRH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_hatchery %>% 
        bind_cols(sim_hatchery[[2]]) -> ff_hatchery
    }
    
    
    # Reformat final fates simulation for wild
    rownames(ff_wild) <- NULL
    column_to_rownames(ff_wild, "state") -> ff_wild
    colnames(ff_wild) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_wild, "state") -> ff_wild
    
    ff_wild %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_wild
    
    ff_wild %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_wild
    
    ff_wild %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_wild_long
    
    ff_wild_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "wild") -> ff_wild_quantiles
    
    # Reformat final fates simulation for hatchery
    rownames(ff_hatchery) <- NULL
    column_to_rownames(ff_hatchery, "state") -> ff_hatchery
    colnames(ff_hatchery) <- paste0("iter", 1:niter)
    
    rownames_to_column(ff_hatchery, "state") -> ff_hatchery
    
    ff_hatchery %>% 
      mutate(state = gsub(" Upstream", "", state)) %>% 
      mutate(state = gsub(" Mouth", "", state)) -> ff_hatchery
    
    ff_hatchery %>% 
      group_by(state) %>% 
      summarise(across(where(is.numeric), sum)) -> ff_hatchery
    
    ff_hatchery %>%
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "count") -> ff_hatchery_long
    
    ff_hatchery_long %>% 
      mutate(prop = count/nsim) %>% 
      group_by(state) %>%
      summarise(prob = quantile(prop, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear_type = "hatchery") -> ff_hatchery_quantiles
    
    ff_wild_quantiles %>% 
      bind_rows(ff_hatchery_quantiles) -> ff_rear_quantiles
  }
  
  
  return(ff_rear_quantiles)
  
}


# ORDER THE STATES FOR PLOTTING
states_order_for_plot <- gsub(" Mouth| Upstream", "", model_states)
states_order_for_plot <- states_order_for_plot[!(duplicated(states_order_for_plot))]
# Make a couple of changes to make them be in the order from most downstream to most upstream
states_order_for_plot[10] <- "Hood River"
states_order_for_plot[11] <- "Fifteenmile Creek"
states_order_for_plot[12] <- "Deschutes River"
states_order_for_plot[13] <- "John Day River"
states_order_for_plot[15] <- "Walla Walla River"
states_order_for_plot[16] <- "Yakima River"
states_order_for_plot[19] <- "Methow River"
states_order_for_plot[20] <- "Okanogan River"

states_order_for_plot[16:29] <- states_order_for_plot[15:28]
states_order_for_plot[15] <- "BON to MCN other tributaries"
states_order_for_plot[23:29] <- states_order_for_plot[22:28]
states_order_for_plot[22] <- "Upstream WEL other tributaries"
states_order_for_plot[29] <- "loss"

states_order_for_plot[24] <- "Clearwater River"
states_order_for_plot[25] <- "Asotin Creek"
states_order_for_plot[26] <- "Grande Ronde River"
states_order_for_plot[27] <- "Salmon River"


plot_final_fate_rear_type <- function(ff_comp, natal_origin){
  rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
  
  ff_comp$state <- fct_rev(factor(ff_comp$state, levels = states_order_for_plot))
  
  ff_comp_plot <- ggplot(ff_comp, aes(x = state, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = rear_type)) +
    geom_point(size = 3.5, shape = 18, position=position_dodge(width=0.5)) +
    geom_linerange(position=position_dodge(width=0.5)) +
    coord_flip() +
    ylab("Final Distribution Probability") +
    xlab("Model State") +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = rear_colors) +
    theme(plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    ggtitle(natal_origin)
  
  return(ff_comp_plot)
  
}


ff_iter <- 1000
ff_nsim <- 1000

#### Winter spill days comparison ####


### Walla Walla River ###

# Walla Walla River Steelhead: compare how homing changes based on winter spill days at Ice Harbor Dam
# compare homing at 0, 20, 30, 50 winter spill days
ICH_winter_spill_values <- c(0, 0.10, 0.20,0.30, 0.40, 0.50)

# cool year: 08/09
# average year: 19/20
# warm year: 17/18

# coldest year on record: 07/08
# hottest year on record: 15/16
fix_run_years <- c("08/09", "19/20", "17/18")

WAWA_ICH_winterspill_homing <- data.frame()

for (i in 1:length(ICH_winter_spill_values)){
  for (j in 1:length(fix_run_years)){
  WAWA_FF <- compare_final_fate_fixcov_rear_type_MC(niter = ff_iter, nsim = ff_nsim,
                                         origin_select = "Walla Walla River",
                                         fix_state = 8, fix_temp_value = NA, 
                                         fix_spillwindow_value = NA, fix_winterspill_value = ICH_winter_spill_values[i],
                                         fix_run_year = fix_run_years[j])
  WAWA_FF$ICH_winterspill <-  ICH_winter_spill_values[i]
  WAWA_FF$fix_run_year <-  fix_run_years[j]
  
  WAWA_homing <- subset(WAWA_FF, state == "Walla Walla River")
  
  WAWA_ICH_winterspill_homing %>% 
    bind_rows(., WAWA_homing) -> WAWA_ICH_winterspill_homing
  
  }
}

WAWA_ICH_winterspill_homing %>% 
  mutate(ICH_winterspill_actual = ICH_winterspill*100) -> WAWA_ICH_winterspill_homing

WAWA_ICH_winterspill_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> WAWA_ICH_winterspill_homing

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
WAWA_ICH_winterspill_homing_plot <- ggplot(WAWA_ICH_winterspill_homing, aes(x = ICH_winterspill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = rear_type,
                                                                            shape = conditions)) +
         geom_point(size = 3.5, position=position_dodge(width=3)) +
         geom_linerange(position=position_dodge(width=3)) +
  scale_shape_manual(values = c(15,16,17)) +
  ylab("Homing Probability") +
  xlab("Days of Winter Spill at Ice Harbor Dam") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = rear_colors) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  ggtitle("Homing by Walla Walla River Steelhead under different winter spill conditions at Ice Harbor Dam")

ggsave(here::here("stan_actual", "output", "final_fates_covariates", "WAWA_ICH_winterspill_homing_plot_temps.png"), WAWA_ICH_winterspill_homing_plot, height = 8, width = 8)

### John Day River ###

# John Day River Steelhead: compare how homing changes based on winter spill days at McNary Dam
# compare homing at 0, 20, 30, 50 winter spill days
MCN_winter_spill_values <- c(0, 0.10, 0.20,0.30, 0.40, 0.50)

# cool year: 08/09
# average year: 19/20
# warm year: 17/18

# coldest year on record: 07/08
# hottest year on record: 15/16
fix_run_years <- c("08/09", "19/20", "17/18")

JDR_MCN_winterspill_homing <- data.frame()

for (i in 1:length(MCN_winter_spill_values)){
  for (j in 1:length(fix_run_years)){
    JDR_FF <- compare_final_fate_fixcov_rear_type_MC(niter = ff_iter, nsim = ff_nsim,
                                                     origin_select = "John Day River",
                                                     fix_state = 3, fix_temp_value = NA, 
                                                     fix_spillwindow_value = NA, fix_winterspill_value = MCN_winter_spill_values[i],
                                                     fix_run_year = fix_run_years[j])
    JDR_FF$MCN_winterspill <-  MCN_winter_spill_values[i]
    JDR_FF$fix_run_year <-  fix_run_years[j]
    
    JDR_homing <- subset(JDR_FF, state == "John Day River")
    
    JDR_MCN_winterspill_homing %>% 
      bind_rows(., JDR_homing) -> JDR_MCN_winterspill_homing
    
  }
}

JDR_MCN_winterspill_homing %>% 
  mutate(MCN_winterspill_actual = MCN_winterspill*100) -> JDR_MCN_winterspill_homing

JDR_MCN_winterspill_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> JDR_MCN_winterspill_homing

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
JDR_MCN_winterspill_homing_plot <- ggplot(JDR_MCN_winterspill_homing, aes(x = MCN_winterspill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = rear_type,
                                                                          shape = conditions)) +
  geom_point(size = 3.5, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3)) +
  scale_shape_manual(values = c(15,16,17)) +
  ylab("Homing Probability") +
  xlab("Days of Winter Spill at McNary Dam") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = rear_colors) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  ggtitle("Homing by John Day River Steelhead under different winter spill conditions at McNary Dam")

ggsave(here::here("stan_actual", "output", "final_fates_covariates", "JDR_MCN_winterspill_homing_plot_temps.png"), JDR_MCN_winterspill_homing_plot, height = 8, width = 8)