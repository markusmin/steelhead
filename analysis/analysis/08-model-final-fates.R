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
library(forcats)


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

#### Final fates functions ####

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


# Final states simulation
# This function will have to be re-written for each DPS, because each DPS has different params
# currently we are not including a random effect of year, but we could easily
# amend the code to simulate final fates in a particular year
# origin1/origin2/origin3 are indicator variables, so set the origin you want to be 1
final_fates_simulation_UCW <- function(nsim,
                                   start_state = 2, states_dates,
                                   origin1 = 0, origin2 = 0, origin3 = 0,
                                   temp_data, spillwindow_data, winterspill_data,
                                   condition_jitter = FALSE){
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
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      temp[i] <- median(temp_data[subset(states_dates, state == i)$date,i])
    }
    
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      spillwindow[i] <- median(spillwindow_data[subset(states_dates, state == i)$date,i])
    }
    
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
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
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
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      temp[i] <- median(temp_data[subset(states_dates, state == i)$date,i])
    }
    
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      spillwindow[i] <- median(spillwindow_data[subset(states_dates, state == i)$date,i])
    }
    
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
                                       origin4 = 0, origin5 = 0, origin6 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
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
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      temp[i] <- median(temp_data[subset(states_dates, state == i)$date,i])
    }
    
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      spillwindow[i] <- median(spillwindow_data[subset(states_dates, state == i)$date,i])
    }
    
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
                                       origin1 = 0, origin2 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
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
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      temp[i] <- median(temp_data[subset(states_dates, state == i)$date,i])
    }
    
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      spillwindow[i] <- median(spillwindow_data[subset(states_dates, state == i)$date,i])
    }
    
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
                                       origin4 = 0, origin5 = 0, origin6 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
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
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      temp[i] <- median(temp_data[subset(states_dates, state == i)$date,i])
    }
    
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      spillwindow[i] <- median(spillwindow_data[subset(states_dates, state == i)$date,i])
    }
    
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
                                       origin4 = 0, origin5 = 0,
                                       temp_data, spillwindow_data, winterspill_data,
                                       condition_jitter = FALSE){
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
  temp <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      temp[i] <- temp_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      temp[i] <- median(temp_data[subset(states_dates, state == i)$date,i])
    }
    
  }
  
  spillwindow <- rep(0, length(model_states))
  for (i in 1:9){ #1:9 because 9 mainstem states
    if(condition_jitter == TRUE){
      # approach 1: jitter conditions
      spillwindow[i] <- spillwindow_data[sample_date[i],i]
    } else {
      # approach 2: Take the median conditions for each state
      spillwindow[i] <- median(spillwindow_data[subset(states_dates, state == i)$date,i])
    }
    
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

#### Run final fates simulation - H v W comparisons ####


# In order to simulate covariate values, we are going to determine the dates
# where the fish were in each state, and then sample from those dates to 
# get a representative value for temperature/spill when fish are in those states
UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                               date = as.vector(UCW_envir$data$transition_dates))
UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                               date = as.vector(UCH_envir$data$transition_dates))
MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                               date = as.vector(MCW_envir$data$transition_dates))
MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                               date = as.vector(MCH_envir$data$transition_dates))
SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
                               date = as.vector(SRW_envir$data$transition_dates))
SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
                               date = as.vector(SRH_envir$data$transition_dates))


# create a list that maps origin numbers (params) to what they actually are
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

origin_param_map <- data.frame(
  natal_origin = natal_origins,
  hatchery = c(NA, NA, NA, NA, 1, NA, 2, # MC
               1,NA,2,3, # UC
               5,NA,1,4,2,3, # SR
               NA, NA),
  wild = c(1,3,NA,2,4,6,5, # MC
           1,2,NA,3, # UC
           6,1,2,5,3,4, # SR
           NA, NA))


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
                                             start_state = 2, states_dates = UCW_states_dates,
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
                                                 start_state = 2, states_dates = UCH_states_dates,
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
                                      start_state = 2, states_dates = UCW_states_dates,
                                      origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                      temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, 
                                      winterspill_data = UCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
    ff_wild %>% 
      bind_cols(sim_wild[[2]]) -> ff_wild
    
    # Run final fates simulation for hatchery
    sim_hatchery <- final_fates_simulation_UCH(nsim = nsim,
                                      start_state = 2, states_dates = UCH_states_dates,
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

compare_final_fate_rear_type_MC <- function(niter, nsim, condition_jitter,
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
      sim_wild <- final_fates_simulation_MCW(nsim = nsim,
                                             start_state = 2, states_dates = MCW_states_dates,
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                             winterspill_data = MCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
      sim_hatchery <- final_fates_simulation_MCH(nsim = nsim,
                                                 start_state = 2, states_dates = MCH_states_dates,
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2], 
                                                 temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                 winterspill_data = MCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
      sim_wild <- final_fates_simulation_MCW(nsim = nsim,
                                             start_state = 2, states_dates = MCW_states_dates,
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = MCW_envir$data$temperature_data, spillwindow_data = MCW_envir$data$spill_window_data, 
                                             winterspill_data = MCW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_MCH(nsim = nsim,
                                                 start_state = 2, states_dates = MCH_states_dates,
                                                 origin1 = hatchery_origin_params[1], origin2 = hatchery_origin_params[2],
                                                 temp_data = MCH_envir$data$temperature_data, spillwindow_data = MCH_envir$data$spill_window_data, 
                                                 winterspill_data = MCH_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
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
                                             start_state = 2, states_dates = SRW_states_dates,
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
                                                 start_state = 2, states_dates = SRH_states_dates,
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
                                             start_state = 2, states_dates = SRW_states_dates,
                                             origin1 = wild_origin_params[1], origin2 = wild_origin_params[2], origin3 = wild_origin_params[3],
                                             origin4 = wild_origin_params[4], origin5 = wild_origin_params[5], origin6 = wild_origin_params[6],
                                             temp_data = SRW_envir$data$temperature_data, spillwindow_data = SRW_envir$data$spill_window_data, 
                                             winterspill_data = SRW_envir$data$winter_spill_days_data, condition_jitter = condition_jitter)
      ff_wild %>% 
        bind_cols(sim_wild[[2]]) -> ff_wild
      
      # Run final fates simulation for hatchery
      sim_hatchery <- final_fates_simulation_SRH(nsim = nsim,
                                                 start_state = 2, states_dates = SRH_states_dates,
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

## Upper Columbia

# Wenatchee comparison
# wen_ff_comp_jitter <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Wenatchee River", condition_jitter = TRUE)
# wen_ff_comp_jitter_plot <- plot_final_fate_rear_type(wen_ff_comp, natal_origin = "Wenatchee River")
# ggsave(here::here("stan_actual", "output", "final_fates", "wen_ff_comp_jitter_plot.png"), wen_ff_comp_jitter_plot, height = 8, width = 8)
wen_ff_comp_median <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Wenatchee River", condition_jitter = FALSE)
wen_ff_comp_median_plot <- plot_final_fate_rear_type(wen_ff_comp_median, natal_origin = "Wenatchee River")
ggsave(here::here("stan_actual", "output", "final_fates", "wen_ff_comp_median_plot.png"), wen_ff_comp_median_plot, height = 8, width = 8)

# Entiat comparison
ent_ff_comp_median <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Entiat River", condition_jitter = FALSE)
ent_ff_comp_median_plot <- plot_final_fate_rear_type(ent_ff_comp_median, natal_origin = "Entiat River")
ggsave(here::here("stan_actual", "output", "final_fates", "ent_ff_comp_median_plot.png"), ent_ff_comp_median_plot, height = 8, width = 8)

# Okanogan comparison
oka_ff_comp_median <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Okanogan River", condition_jitter = FALSE)
oka_ff_comp_median_plot <- plot_final_fate_rear_type(oka_ff_comp_median, natal_origin = "Okanogan River")
ggsave(here::here("stan_actual", "output", "final_fates", "oka_ff_comp_median_plot.png"), oka_ff_comp_median_plot, height = 8, width = 8)

# Methow comparison
# this used to crash with jitter - let's see if it works now
met_ff_comp_median <- compare_final_fate_rear_type_UC(niter = ff_iter, nsim = ff_nsim, origin_select = "Methow River", condition_jitter = FALSE)
met_ff_comp_median_plot <- plot_final_fate_rear_type(met_ff_comp_median, natal_origin = "Methow River")
ggsave(here::here("stan_actual", "output", "final_fates", "met_ff_comp_median_plot.png"), met_ff_comp_median_plot, height = 8, width = 8)



## Middle Columbia
# Deschutes comparison
des_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Deschutes River", condition_jitter = FALSE)
des_ff_comp_median_plot <- plot_final_fate_rear_type(des_ff_comp_median, natal_origin = "Deschutes River")
ggsave(here::here("stan_actual", "output", "final_fates", "des_ff_comp_median_plot_median_conditions.png"), des_ff_comp_median_plot, height = 8, width = 8)

# John Day comparison
jdr_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "John Day River", condition_jitter = FALSE)
jdr_ff_comp_median_plot <- plot_final_fate_rear_type(jdr_ff_comp_median, natal_origin = "John Day River")
ggsave(here::here("stan_actual", "output", "final_fates", "jdr_ff_comp_median_plot_median_conditions_v2.png"), jdr_ff_comp_median_plot, height = 8, width = 8)

# Fifteenmile Creek comparison
fif_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Fifteenmile Creek", condition_jitter = FALSE)
fif_ff_comp_median_plot <- plot_final_fate_rear_type(fif_ff_comp_median, natal_origin = "Fifteenmile Creek")
ggsave(here::here("stan_actual", "output", "final_fates", "fif_ff_comp_median_plot_median_conditions.png"), fif_ff_comp_median_plot, height = 8, width = 8)

# Umatilla comparison
uma_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Umatilla River", condition_jitter = FALSE)
uma_ff_comp_median_plot <- plot_final_fate_rear_type(uma_ff_comp_median, natal_origin = "Umatilla River")
ggsave(here::here("stan_actual", "output", "final_fates", "uma_ff_comp_median_plot_median_conditions.png"), uma_ff_comp_median_plot, height = 8, width = 8)

# Yakima comparison
yak_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Yakima River", condition_jitter = FALSE)
yak_ff_comp_median_plot <- plot_final_fate_rear_type(yak_ff_comp_median, natal_origin = "Yakima River")
ggsave(here::here("stan_actual", "output", "final_fates", "yak_ff_comp_median_plot_median_conditions.png"), yak_ff_comp_median_plot, height = 8, width = 8)

# Walla Walla comparison
wawa_ff_comp_median <- compare_final_fate_rear_type_MC(niter = ff_iter, nsim = ff_nsim, origin_select = "Walla Walla River", condition_jitter = FALSE)
wawa_ff_comp_median_plot <- plot_final_fate_rear_type(wawa_ff_comp_median, natal_origin = "Walla Walla River")
ggsave(here::here("stan_actual", "output", "final_fates", "wawa_ff_comp_median_plot_median_conditions.png"), wawa_ff_comp_median_plot, height = 8, width = 8)


## Snake River
# Asotin Creek comparison
aso_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Asotin Creek", condition_jitter = FALSE)
aso_ff_comp_median_plot <- plot_final_fate_rear_type(aso_ff_comp_median, natal_origin = "Asotin Creek")
ggsave(here::here("stan_actual", "output", "final_fates", "aso_ff_comp_median_plot.png"), aso_ff_comp_median_plot, height = 8, width = 8)

# Clearwater comparison
cle_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Clearwater River", condition_jitter = FALSE)
cle_ff_comp_median_plot <- plot_final_fate_rear_type(cle_ff_comp_median, natal_origin = "Clearwater River")
ggsave(here::here("stan_actual", "output", "final_fates", "cle_ff_comp_median_plot.png"), cle_ff_comp_median_plot, height = 8, width = 8)

# Salmon comparison
sal_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Salmon River", condition_jitter = FALSE)
sal_ff_comp_median_plot <- plot_final_fate_rear_type(sal_ff_comp_median, natal_origin = "Salmon River")
ggsave(here::here("stan_actual", "output", "final_fates", "sal_ff_comp_median_plot.png"), sal_ff_comp_median_plot, height = 8, width = 8)

# Grande Ronde comparison
gr_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Grande Ronde River", condition_jitter = FALSE)
gr_ff_comp_median_plot <- plot_final_fate_rear_type(gr_ff_comp_median, natal_origin = "Grande Ronde River")
ggsave(here::here("stan_actual", "output", "final_fates", "gr_ff_comp_median_plot.png"), gr_ff_comp_median_plot, height = 8, width = 8)

# Imnaha comparison
imn_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Imnaha River", condition_jitter = FALSE)
imn_ff_comp_median_plot <- plot_final_fate_rear_type(imn_ff_comp_median, natal_origin = "Imnaha River")
ggsave(here::here("stan_actual", "output", "final_fates", "imn_ff_comp_median_plot.png"), imn_ff_comp_median_plot, height = 8, width = 8)

# Tucannon comparison
tuc_ff_comp_median <- compare_final_fate_rear_type_SR(niter = ff_iter, nsim = ff_nsim, origin_select = "Tucannon River", condition_jitter = FALSE)
tuc_ff_comp_median_plot <- plot_final_fate_rear_type(tuc_ff_comp_median, natal_origin = "Tucannon River")
ggsave(here::here("stan_actual", "output", "final_fates", "tuc_ff_comp_median_plot.png"), tuc_ff_comp_median_plot, height = 8, width = 8)


