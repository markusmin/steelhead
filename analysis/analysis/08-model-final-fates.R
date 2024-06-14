# 08-model-final-fates

# This script takes the output from the stan model runs in 05-stan-runs and
# generates estimates of final fate distributions

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
  parameters[grepl(parameter_prefix, parameters)] -> param_subset
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

#### Final fates functions ####

# This function takes a model fit object, a number of simulated fish,
# the identities of those fish (origin, rear), the conditions throughout the
# basin (spill, temperature, year)
# and then simulates the final fates of those fish
# By default, the conditions (spill and temperature) will be taken from the data
# itself, but note that new values of these conditions could be simulated in order

# Final states simulation
# This function will have to be re-written for each DPS, because each DPS has different params
# currently we are not including a random effect of year, but we could easily
# amend the code to simulate final fates in a particular year
# origin1/origin2/origin3 are indicator variables, so set the origin you want to be 1
final_fates_simulation_UCW <- function(nsim,
                                   start_state = 2, states_dates,
                                   origin1 = 0, origin2 = 0, origin3 = 0,
                                   temp_data, spillwindow_data, winterspill_data){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)

  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    winterspill[i] <- winterspill_data[sample_year[i],i]
  }

  
  # get temp and spill window data
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    temp[i] <- temp_data[sample_date[i],i]
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    spillwindow[i] <- spillwindow_data[sample_date[i],i]
  }
    
  
  move_prob_matrix <- matrix(data = 0,
                             nrow = length(model_states),
                             ncol = length(model_states))
  
  for (i in 1:nrow(move_prob_matrix)){
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    for (j in 1:length(possible_movements)){
      move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCW[i,possible_movements[j],iter] +
                                                         # btemp0_array_UCW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                         btemp1_array_UCW[i,possible_movements[j],iter]*temp[i] + 
                                                         bspillwindow_array_UCW[i,possible_movements[j],iter]*spillwindow[i] + 
                                                         # bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
                                                         # btemp0xorigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                         btemp1xorigin1_array_UCW[i,possible_movements[j],iter]*temp[i]*origin1 +
                                                         # btemp0xorigin2_array_UCW[i,possible_movements[j],iter]*origin2 + 
                                                         btemp1xorigin2_array_UCW[i,possible_movements[j],iter]*temp[i]*origin2 + 
                                                         # btemp0xorigin3_array_UCW[i,possible_movements[j],iter]*origin3 +
                                                         btemp1xorigin3_array_UCW[i,possible_movements[j],iter]*temp[i]*origin3 +
                                                         borigin1_array_UCW[i,possible_movements[j],iter]*origin1 +
                                                         borigin2_array_UCW[i,possible_movements[j],iter]*origin2 +
                                                         borigin3_array_UCW[i,possible_movements[j],iter]*origin3)/
        sum(exp(b0_array_UCW[i,possible_movements,iter] +
                  # btemp0_array_UCW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_UCW[i,possible_movements,iter]*temp[i] + 
                  bspillwindow_array_UCW[i,possible_movements,iter]*spillwindow[i] + 
                  # bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_UCW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_UCW[i,possible_movements,iter]*temp[i]*origin1 +
                  # btemp0xorigin2_array_UCW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_UCW[i,possible_movements,iter]*temp[i]*origin2 + 
                  # btemp0xorigin3_array_UCW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_UCW[i,possible_movements,iter]*temp[i]*origin3 +
                  borigin1_array_UCW[i,possible_movements,iter]*origin1 +
                  borigin2_array_UCW[i,possible_movements,iter]*origin2 +
                  borigin3_array_UCW[i,possible_movements,iter]*origin3))
    }
    
  }
  
  # Now, use the movement probability matrix to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states

  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  # run while loop until all fish have entered loss state
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(movement_matrix) <- model_states
    colnames(movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = state_matrix[j,i], prob = move_prob_matrix[j,])
    }
    # sum across rows, store in the state matrix
    state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (state_matrix[length(model_states), ncol(state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_UCH <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    winterspill[i] <- winterspill_data[sample_year[i],i]
  }
  
  
  # get temp and spill window data
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    temp[i] <- temp_data[sample_date[i],i]
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    spillwindow[i] <- spillwindow_data[sample_date[i],i]
  }
  
  
  move_prob_matrix <- matrix(data = 0,
                             nrow = length(model_states),
                             ncol = length(model_states))
  
  for (i in 1:nrow(move_prob_matrix)){
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    for (j in 1:length(possible_movements)){
      move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_UCH[i,possible_movements[j],iter] +
                                                         # btemp0_array_UCH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                         btemp1_array_UCH[i,possible_movements[j],iter]*temp[i] + 
                                                         bspillwindow_array_UCH[i,possible_movements[j],iter]*spillwindow[i] + 
                                                         # bwinterspill_array_UCH[i,possible_movements[j],iter]*winterspill[i] +
                                                         # btemp0xorigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                         btemp1xorigin1_array_UCH[i,possible_movements[j],iter]*temp[i]*origin1 +
                                                         # btemp0xorigin2_array_UCH[i,possible_movements[j],iter]*origin2 + 
                                                         btemp1xorigin2_array_UCH[i,possible_movements[j],iter]*temp[i]*origin2 + 
                                                         # btemp0xorigin3_array_UCH[i,possible_movements[j],iter]*origin3 +
                                                         btemp1xorigin3_array_UCH[i,possible_movements[j],iter]*temp[i]*origin3 +
                                                         borigin1_array_UCH[i,possible_movements[j],iter]*origin1 +
                                                         borigin2_array_UCH[i,possible_movements[j],iter]*origin2 +
                                                         borigin3_array_UCH[i,possible_movements[j],iter]*origin3)/
        sum(exp(b0_array_UCH[i,possible_movements,iter] +
                  # btemp0_array_UCH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_UCH[i,possible_movements,iter]*temp[i] + 
                  bspillwindow_array_UCH[i,possible_movements,iter]*spillwindow[i] + 
                  # bwinterspill_array_UCH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_UCH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_UCH[i,possible_movements,iter]*temp[i]*origin1 +
                  # btemp0xorigin2_array_UCH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_UCH[i,possible_movements,iter]*temp[i]*origin2 + 
                  # btemp0xorigin3_array_UCH[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_UCH[i,possible_movements,iter]*temp[i]*origin3 +
                  borigin1_array_UCH[i,possible_movements,iter]*origin1 +
                  borigin2_array_UCH[i,possible_movements,iter]*origin2 +
                  borigin3_array_UCH[i,possible_movements,iter]*origin3))
    }
    
  }
  
  # Now, use the movement probability matrix to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  # run while loop until all fish have entered loss state
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(movement_matrix) <- model_states
    colnames(movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = state_matrix[j,i], prob = move_prob_matrix[j,])
    }
    # sum across rows, store in the state matrix
    state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (state_matrix[length(model_states), ncol(state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_MCW <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    winterspill[i] <- winterspill_data[sample_year[i],i]
  }
  
  
  # get temp and spill window data
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    temp[i] <- temp_data[sample_date[i],i]
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    spillwindow[i] <- spillwindow_data[sample_date[i],i]
  }
  
  
  move_prob_matrix <- matrix(data = 0,
                             nrow = length(model_states),
                             ncol = length(model_states))
  
  for (i in 1:nrow(move_prob_matrix)){
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    for (j in 1:length(possible_movements)){
      move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCW[i,possible_movements[j],iter] +
                                                         # btemp0_array_MCW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                         btemp1_array_MCW[i,possible_movements[j],iter]*temp[i] + 
                                                         bspillwindow_array_MCW[i,possible_movements[j],iter]*spillwindow[i] + 
                                                         # bwinterspill_array_MCW[i,possible_movements[j],iter]*winterspill[i] +
                                                         # btemp0xorigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                         btemp1xorigin1_array_MCW[i,possible_movements[j],iter]*temp[i]*origin1 +
                                                         # btemp0xorigin2_array_MCW[i,possible_movements[j],iter]*origin2 + 
                                                         btemp1xorigin2_array_MCW[i,possible_movements[j],iter]*temp[i]*origin2 + 
                                                         # btemp0xorigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                         btemp1xorigin3_array_MCW[i,possible_movements[j],iter]*temp[i]*origin3 +
                                                         # btemp0xorigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                         btemp1xorigin4_array_MCW[i,possible_movements[j],iter]*temp[i]*origin4 +
                                                         # btemp0xorigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                         btemp1xorigin5_array_MCW[i,possible_movements[j],iter]*temp[i]*origin5 +
                                                         # btemp0xorigin6_array_MCW[i,possible_movements[j],iter]*origin6 +
                                                         btemp1xorigin6_array_MCW[i,possible_movements[j],iter]*temp[i]*origin6 +
                                                         borigin1_array_MCW[i,possible_movements[j],iter]*origin1 +
                                                         borigin2_array_MCW[i,possible_movements[j],iter]*origin2 +
                                                         borigin3_array_MCW[i,possible_movements[j],iter]*origin3 +
                                                         borigin4_array_MCW[i,possible_movements[j],iter]*origin4 +
                                                         borigin5_array_MCW[i,possible_movements[j],iter]*origin5 +
                                                         borigin6_array_MCW[i,possible_movements[j],iter]*origin6)/
        sum(exp(b0_array_MCW[i,possible_movements,iter] +
                  # btemp0_array_MCW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_MCW[i,possible_movements,iter]*temp[i] + 
                  bspillwindow_array_MCW[i,possible_movements,iter]*spillwindow[i] + 
                  # bwinterspill_array_MCW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_MCW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_MCW[i,possible_movements,iter]*temp[i]*origin1 +
                  # btemp0xorigin2_array_MCW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_MCW[i,possible_movements,iter]*temp[i]*origin2 + 
                  # btemp0xorigin3_array_MCW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_MCW[i,possible_movements,iter]*temp[i]*origin3 +
                  # btemp0xorigin4_array_MCW[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_MCW[i,possible_movements,iter]*temp[i]*origin4 +
                  # btemp0xorigin5_array_MCW[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_MCW[i,possible_movements,iter]*temp[i]*origin5 +
                  # btemp0xorigin6_array_MCW[i,possible_movements,iter]*origin6 +
                  btemp1xorigin6_array_MCW[i,possible_movements,iter]*temp[i]*origin6 +
                  borigin1_array_MCW[i,possible_movements,iter]*origin1 +
                  borigin2_array_MCW[i,possible_movements,iter]*origin2 +
                  borigin3_array_MCW[i,possible_movements,iter]*origin3 +
                  borigin4_array_MCW[i,possible_movements,iter]*origin4 +
                  borigin5_array_MCW[i,possible_movements,iter]*origin5 +
                  borigin6_array_MCW[i,possible_movements,iter]*origin6))
    }
    
  }
  
  # Now, use the movement probability matrix to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  # run while loop until all fish have entered loss state
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(movement_matrix) <- model_states
    colnames(movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = state_matrix[j,i], prob = move_prob_matrix[j,])
    }
    # sum across rows, store in the state matrix
    state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (state_matrix[length(model_states), ncol(state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_MCH <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    winterspill[i] <- winterspill_data[sample_year[i],i]
  }
  
  
  # get temp and spill window data
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    temp[i] <- temp_data[sample_date[i],i]
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    spillwindow[i] <- spillwindow_data[sample_date[i],i]
  }
  
  
  move_prob_matrix <- matrix(data = 0,
                             nrow = length(model_states),
                             ncol = length(model_states))
  
  for (i in 1:nrow(move_prob_matrix)){
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    for (j in 1:length(possible_movements)){
      move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_MCH[i,possible_movements[j],iter] +
                                                         # btemp0_array_MCH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                         btemp1_array_MCH[i,possible_movements[j],iter]*temp[i] + 
                                                         bspillwindow_array_MCH[i,possible_movements[j],iter]*spillwindow[i] + 
                                                         # bwinterspill_array_MCH[i,possible_movements[j],iter]*winterspill[i] +
                                                         # btemp0xorigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                         btemp1xorigin1_array_MCH[i,possible_movements[j],iter]*temp[i]*origin1 +
                                                         # btemp0xorigin2_array_MCH[i,possible_movements[j],iter]*origin2 + 
                                                         btemp1xorigin2_array_MCH[i,possible_movements[j],iter]*temp[i]*origin2 + 
                                                         
                                                         borigin1_array_MCH[i,possible_movements[j],iter]*origin1 +
                                                         borigin2_array_MCH[i,possible_movements[j],iter]*origin2)/
        sum(exp(b0_array_MCH[i,possible_movements,iter] +
                  # btemp0_array_MCH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_MCH[i,possible_movements,iter]*temp[i] + 
                  bspillwindow_array_MCH[i,possible_movements,iter]*spillwindow[i] + 
                  # bwinterspill_array_MCH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_MCH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_MCH[i,possible_movements,iter]*temp[i]*origin1 +
                  # btemp0xorigin2_array_MCH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_MCH[i,possible_movements,iter]*temp[i]*origin2 + 
                  borigin1_array_MCH[i,possible_movements,iter]*origin1 +
                  borigin2_array_MCH[i,possible_movements,iter]*origin2))
    }
    
  }
  
  # Now, use the movement probability matrix to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  # run while loop until all fish have entered loss state
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(movement_matrix) <- model_states
    colnames(movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = state_matrix[j,i], prob = move_prob_matrix[j,])
    }
    # sum across rows, store in the state matrix
    state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (state_matrix[length(model_states), ncol(state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_SRW <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    winterspill[i] <- winterspill_data[sample_year[i],i]
  }
  
  
  # get temp and spill window data
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    temp[i] <- temp_data[sample_date[i],i]
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    spillwindow[i] <- spillwindow_data[sample_date[i],i]
  }
  
  
  move_prob_matrix <- matrix(data = 0,
                             nrow = length(model_states),
                             ncol = length(model_states))
  
  for (i in 1:nrow(move_prob_matrix)){
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    for (j in 1:length(possible_movements)){
      move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRW[i,possible_movements[j],iter] +
                                                         # btemp0_array_SRW[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                         btemp1_array_SRW[i,possible_movements[j],iter]*temp[i] + 
                                                         bspillwindow_array_SRW[i,possible_movements[j],iter]*spillwindow[i] + 
                                                         # bwinterspill_array_SRW[i,possible_movements[j],iter]*winterspill[i] +
                                                         # btemp0xorigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                         btemp1xorigin1_array_SRW[i,possible_movements[j],iter]*temp[i]*origin1 +
                                                         # btemp0xorigin2_array_SRW[i,possible_movements[j],iter]*origin2 + 
                                                         btemp1xorigin2_array_SRW[i,possible_movements[j],iter]*temp[i]*origin2 + 
                                                         # btemp0xorigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                         btemp1xorigin3_array_SRW[i,possible_movements[j],iter]*temp[i]*origin3 +
                                                         # btemp0xorigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                         btemp1xorigin4_array_SRW[i,possible_movements[j],iter]*temp[i]*origin4 +
                                                         # btemp0xorigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                         btemp1xorigin5_array_SRW[i,possible_movements[j],iter]*temp[i]*origin5 +
                                                         # btemp0xorigin6_array_SRW[i,possible_movements[j],iter]*origin6 +
                                                         btemp1xorigin6_array_SRW[i,possible_movements[j],iter]*temp[i]*origin6 +
                                                         borigin1_array_SRW[i,possible_movements[j],iter]*origin1 +
                                                         borigin2_array_SRW[i,possible_movements[j],iter]*origin2 +
                                                         borigin3_array_SRW[i,possible_movements[j],iter]*origin3 +
                                                         borigin4_array_SRW[i,possible_movements[j],iter]*origin4 +
                                                         borigin5_array_SRW[i,possible_movements[j],iter]*origin5 +
                                                         borigin6_array_SRW[i,possible_movements[j],iter]*origin6)/
        sum(exp(b0_array_SRW[i,possible_movements,iter] +
                  # btemp0_array_SRW[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_SRW[i,possible_movements,iter]*temp[i] + 
                  bspillwindow_array_SRW[i,possible_movements,iter]*spillwindow[i] + 
                  # bwinterspill_array_SRW[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_SRW[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_SRW[i,possible_movements,iter]*temp[i]*origin1 +
                  # btemp0xorigin2_array_SRW[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_SRW[i,possible_movements,iter]*temp[i]*origin2 + 
                  # btemp0xorigin3_array_SRW[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_SRW[i,possible_movements,iter]*temp[i]*origin3 +
                  # btemp0xorigin4_array_SRW[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_SRW[i,possible_movements,iter]*temp[i]*origin4 +
                  # btemp0xorigin5_array_SRW[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_SRW[i,possible_movements,iter]*temp[i]*origin5 +
                  # btemp0xorigin6_array_SRW[i,possible_movements,iter]*origin6 +
                  btemp1xorigin6_array_SRW[i,possible_movements,iter]*temp[i]*origin6 +
                  borigin1_array_SRW[i,possible_movements,iter]*origin1 +
                  borigin2_array_SRW[i,possible_movements,iter]*origin2 +
                  borigin3_array_SRW[i,possible_movements,iter]*origin3 +
                  borigin4_array_SRW[i,possible_movements,iter]*origin4 +
                  borigin5_array_SRW[i,possible_movements,iter]*origin5 +
                  borigin6_array_SRW[i,possible_movements,iter]*origin6))
    }
    
  }
  
  # Now, use the movement probability matrix to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  # run while loop until all fish have entered loss state
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(movement_matrix) <- model_states
    colnames(movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = state_matrix[j,i], prob = move_prob_matrix[j,])
    }
    # sum across rows, store in the state matrix
    state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (state_matrix[length(model_states), ncol(state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

final_fates_simulation_SRH <- function(nsim,
                                       start_state = 2, states_dates,
                                       origin1 = 0, origin2 = 0, origin3 = 0,
                                       temp_data, spillwindow_data, winterspill_data){
  # select the iteration you will use
  iter <- sample(1:4000, 1)
  
  # randomly select dates from actual dates at states
  # create a vector of length 9 (for 9 mainstem sites)
  sample_date <- vector(length = 9)
  
  # populate by selecting one date randomly for each state
  for(i in 1:length(sample_date)){
    sample_date[i] <- sample(subset(states_dates, state == i)$date, 1)
  }
  
  # convert these to random years for winter spill
  # From modeling code: We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
  sample_year <- ceiling(sample_date/365.25)
  
  
  winterspill <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    winterspill[i] <- winterspill_data[sample_year[i],i]
  }
  
  
  # get temp and spill window data
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    temp[i] <- temp_data[sample_date[i],i]
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    spillwindow[i] <- spillwindow_data[sample_date[i],i]
  }
  
  
  move_prob_matrix <- matrix(data = 0,
                             nrow = length(model_states),
                             ncol = length(model_states))
  
  for (i in 1:nrow(move_prob_matrix)){
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == i]
    possible_movements <- c(possible_movements, 43)
    
    for (j in 1:length(possible_movements)){
      move_prob_matrix[i,possible_movements[j]] <- exp(b0_array_SRH[i,possible_movements[j],iter] +
                                                         # btemp0_array_SRH[i,possible_movements[j],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                         btemp1_array_SRH[i,possible_movements[j],iter]*temp[i] + 
                                                         bspillwindow_array_SRH[i,possible_movements[j],iter]*spillwindow[i] + 
                                                         # bwinterspill_array_SRH[i,possible_movements[j],iter]*winterspill[i] +
                                                         # btemp0xorigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                         btemp1xorigin1_array_SRH[i,possible_movements[j],iter]*temp[i]*origin1 +
                                                         # btemp0xorigin2_array_SRH[i,possible_movements[j],iter]*origin2 + 
                                                         btemp1xorigin2_array_SRH[i,possible_movements[j],iter]*temp[i]*origin2 + 
                                                         # btemp0xorigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                         btemp1xorigin3_array_SRH[i,possible_movements[j],iter]*temp[i]*origin3 +
                                                         # btemp0xorigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                         btemp1xorigin4_array_SRH[i,possible_movements[j],iter]*temp[i]*origin4 +
                                                         # btemp0xorigin5_array_SRH[i,possible_movements[j],iter]*origin5 +
                                                         btemp1xorigin5_array_SRH[i,possible_movements[j],iter]*temp[i]*origin5 +
                                                         borigin1_array_SRH[i,possible_movements[j],iter]*origin1 +
                                                         borigin2_array_SRH[i,possible_movements[j],iter]*origin2 +
                                                         borigin3_array_SRH[i,possible_movements[j],iter]*origin3 +
                                                         borigin4_array_SRH[i,possible_movements[j],iter]*origin4 +
                                                         borigin5_array_SRH[i,possible_movements[j],iter]*origin5)/
        sum(exp(b0_array_SRH[i,possible_movements,iter] +
                  # btemp0_array_SRH[i,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_SRH[i,possible_movements,iter]*temp[i] + 
                  bspillwindow_array_SRH[i,possible_movements,iter]*spillwindow[i] + 
                  # bwinterspill_array_SRH[i,possible_movements,iter]*winterspill[i] +
                  # btemp0xorigin1_array_SRH[i,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_SRH[i,possible_movements,iter]*temp[i]*origin1 +
                  # btemp0xorigin2_array_SRH[i,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_SRH[i,possible_movements,iter]*temp[i]*origin2 + 
                  # btemp0xorigin3_array_SRH[i,possible_movements,iter]*origin3 +
                  btemp1xorigin3_array_SRH[i,possible_movements,iter]*temp[i]*origin3 +
                  # btemp0xorigin4_array_SRH[i,possible_movements,iter]*origin4 +
                  btemp1xorigin4_array_SRH[i,possible_movements,iter]*temp[i]*origin4 +
                  # btemp0xorigin5_array_SRH[i,possible_movements,iter]*origin5 +
                  btemp1xorigin5_array_SRH[i,possible_movements,iter]*temp[i]*origin5 +
                  borigin1_array_SRH[i,possible_movements,iter]*origin1 +
                  borigin2_array_SRH[i,possible_movements,iter]*origin2 +
                  borigin3_array_SRH[i,possible_movements,iter]*origin3 +
                  borigin4_array_SRH[i,possible_movements,iter]*origin4 +
                  borigin5_array_SRH[i,possible_movements,iter]*origin5))
    }
    
  }
  
  # Now, use the movement probability matrix to simulate movements
  
  
  # create a state matrix to track movements through the model
  state_matrix <- matrix(rep(0, length(model_states)), nrow = length(model_states), ncol = 1)
  rownames(state_matrix) <- model_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, length(model_states)), nrow = 1, ncol = length(model_states))
  colnames(final_fate_matrix) <- model_states
  
  # run while loop until all fish have entered loss state
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    state_matrix %>% 
      bind_cols(., t = rep(0, length(model_states))) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, length(model_states)*length(model_states)), ncol = length(model_states), nrow = length(model_states))
    rownames(movement_matrix) <- model_states
    colnames(movement_matrix) <- model_states
    for (j in 1:nrow(state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = state_matrix[j,i], prob = move_prob_matrix[j,])
    }
    # sum across rows, store in the state matrix
    state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (state_matrix[length(model_states), ncol(state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  # count all of the final fates
  data.frame(colSums(final_fate_matrix[,1:(length(model_states)-1)])) -> final_fates
  colnames(final_fates) <- "count"
  
  outputs <- list(state_matrix, final_fates, move_prob_matrix)
  # return(final_fates)
  return(outputs)
  
}

#### Run final fates simulation - H v W comparisons for UC ####


# In order to simulate covariate values, we are going to determine the dates
# where the fish were in each state, and then sample from those dates to 
# get a representative value for temperature/spill when fish are in those states
UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                               date = as.vector(UCW_envir$data$transition_dates))
UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                               date = as.vector(UCH_envir$data$transition_dates))


# Create a function to compare hatchery and wild final fates for each tributary
# select one origin to compare
compare_final_fate_rear_type_UC <- function(niter, nsim,
                                            origin1, origin2, origin3){
  ff_wild <- data.frame(state = model_states[1:(length(model_states)-1)])
  ff_hatchery <- data.frame(state = model_states[1:(length(model_states)-1)])
  
  for(i in 1:niter) {
    # Run final fates simulation for wild
    sim_wild <- final_fates_simulation_UCW(nsim = nsim,
                                      start_state = 2, states_dates = UCW_states_dates,
                                      origin1 = origin1, origin2 = origin2, origin3 = origin3,
                                      temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, winterspill_data = UCW_envir$data$winter_spill_days_data)
    ff_wild %>% 
      bind_cols(sim_wild[[2]]) -> ff_wild
    
    # Run final fates simulation for hatchery
    sim_hatchery <- final_fates_simulation_UCH(nsim = nsim,
                                      start_state = 2, states_dates = UCW_states_dates,
                                      origin1 = origin1, origin2 = origin2, origin3 = origin3,
                                      temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, winterspill_data = UCW_envir$data$winter_spill_days_data)
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
  
  return(ff_rear_quantiles)
  
}

plot_final_fate_rear_type <- function(ff_comp, natal_origin){
  ff_comp_plot <- ggplot(ff_comp, aes(x = state, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = rear_type)) +
    geom_point(size = 3.5, shape = 18, position=position_dodge(width=0.5)) +
    geom_linerange(position=position_dodge(width=0.5)) +
    coord_flip() +
    ylab("Final Distribution Probability") +
    xlab("Model State") +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    theme(plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    ggtitle(natal_origin)
  
  return(ff_comp_plot)
  
}

# Wenatchee comparison
wen_ff_comp <- compare_final_fate_rear_type_UC(niter = 100, nsim = 100, origin1 = 1, origin2 = 0, origin3 = 0)
plot_final_fate_rear_type(wen_ff_comp, natal_origin = "Wenatchee River")

# Entiat comparison
ent_ff_comp <- compare_final_fate_rear_type_UC(niter = 100, nsim = 100, origin1 = 0, origin2 = 1, origin3 = 0)
plot_final_fate_rear_type(ent_ff_comp, natal_origin = "Entiat River")

# Methow comparison
met_ff_comp <- compare_final_fate_rear_type_UC(niter = 100, nsim = 100, origin1 = 0, origin2 = 0, origin3 = 3)
plot_final_fate_rear_type(met_ff_comp, natal_origin = "Methow River")








