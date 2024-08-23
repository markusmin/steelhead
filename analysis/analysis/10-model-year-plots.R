# 10-model-year-plot

# This script takes the output from the stan model runs in 05-stan-runs and
# plots the random effect of year

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

#### Extract fixed parameter values from model fit objects ####

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


#### Extract year parameter values from the model fit objects ####

# Here it's a bit weird with NDE vs DE - because those change by year
# We're only going to look at DE, because we trust those more than NDE

# function to store all year effects in an array
make_year_origin_parameter_draws_array <- function(fit, fit_summary, origin_select_numeric, envir){
  
  parameters <- fit_summary$variable
  # Figure out which year effects we're actually estimating
  parameters[grepl(paste0("raw_vector"), parameters)] -> year_raw_params
  year_raw_params_from = as.numeric(sub("[^_]*_[^_]*_[^_]*_", "", str_extract(year_raw_params, "[^_]*_[^_]*_[^_]*_[^_]*")))
  year_raw_params_to = as.numeric(str_extract(sub("[^_]*_[^_]*_[^_]*_[^_]*_", "", year_raw_params), "\\d+"))
  year_raw_params_indices <- data.frame(parameter = year_raw_params, from = year_raw_params_from, to = year_raw_params_to)
  year_raw_params_indices %>% 
    dplyr::select(-parameter) %>% 
    mutate(from_to = paste0(from, "_", to)) %>% 
    filter(!(duplicated(from_to)))-> year_raw_params_indices
    
  
  # get the parameter indexing
  param_indices_matrix <- envir$data$parameter_indices_matrix
  param_indices <- data.frame(row = vector(length = max(param_indices_matrix)-1),
                              col = vector(length = max(param_indices_matrix)-1))
  for (i in 1:(max(param_indices_matrix)-1)){
    param_indices[i,] <- which(param_indices_matrix == i, arr.ind = TRUE)
  }
  param_indices$index <- as.numeric(rownames(param_indices))
  param_indices$from_to <- paste0(param_indices$row, "_", param_indices$col)
  
  # Now, include only the indices that actually have year effects
  year_movement_indices <- filter(param_indices, from_to %in% year_raw_params_indices$from_to)
  year_movement_indices %>% 
    dplyr::rename(from = row, to = col) -> year_movement_indices
  # Loop through this (?) to create the param array below

  
  
  
  # extract year effect as an 4-d array with rows=from, columns=to, years = slices, and iter=4th dimension
  parameters <- fit_summary$variable
  parameters[grepl(paste0("actual"), parameters)] -> year_params
  # drop the NDE parameters
  year_params <- year_params[!(grepl("_NDE", year_params))]
  # select only one origin
  year_origin_params <- year_params[grepl(origin_select_numeric, year_params)]
  
  
  year_origin_params_index <- as.numeric(str_split_i(str_extract(year_origin_params, "(?<=\\[).*(?=\\])"), "[,]", 1))
  year_origin_params_year <-  as.numeric(str_split_i(str_extract(year_origin_params, "(?<=\\[).*(?=\\])"), "[,]", 2))
  
  year_origin_params_indices <- data.frame(parameter = year_origin_params, index = year_origin_params_index, year = year_origin_params_year)
  year_origin_params_indices %>% 
    left_join(., year_movement_indices, by = "index") %>% 
    filter(!(is.na(from))) -> year_origin_params_indices
  
  # get a repeating index for indexing the 3d array below
  year_origin_params_indices$array_index <- rep(1:length(unique(year_origin_params_indices$index)), length(unique(year_origin_params_indices$year)))
  
  
  # update: extract year effect as 3-d array, with rows = movements, columns = years, iter = slices
  # include movements that don't have a year effect, to make our indexing easier later on
  # these will all have zero (for no year effect)
  year_origin_param_array <- array(data = 0, dim = c(nrow(param_indices),
                                                     length(unique(year_origin_params_indices$year)),
                                                     length(as.matrix(fit[,,1]))))
  
  # give them names that indicate movements
  rownames(year_origin_param_array) <- param_indices$from_to

  # now populate the array with the iteration draws
  for(i in 1:nrow(year_origin_params_indices)){
    year_origin_param_array[year_origin_params_indices[i, "index"], year_origin_params_indices[i, "year"], ] <- as.matrix(fit[,,year_origin_params_indices[i, "parameter"]])
  }
  

  
  
  return(year_origin_param_array)
}



### UCW ###
origin1_year_param_array_UCW <- make_year_origin_parameter_draws_array(fit = UCW_fit, fit_summary = UCW_fit_summary, 
                                                                origin_select_numeric = "origin1", envir = UCW_envir)

origin2_year_param_array_UCW <- make_year_origin_parameter_draws_array(fit = UCW_fit, fit_summary = UCW_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = UCW_envir)

origin3_year_param_array_UCW <- make_year_origin_parameter_draws_array(fit = UCW_fit, fit_summary = UCW_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = UCW_envir)

### UCH ###
origin1_year_param_array_UCH <- make_year_origin_parameter_draws_array(fit = UCH_fit, fit_summary = UCH_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = UCH_envir)

origin2_year_param_array_UCH <- make_year_origin_parameter_draws_array(fit = UCH_fit, fit_summary = UCH_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = UCH_envir)

origin3_year_param_array_UCH <- make_year_origin_parameter_draws_array(fit = UCH_fit, fit_summary = UCH_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = UCH_envir)

### MCW ###
origin1_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = MCW_envir)

origin2_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = MCW_envir)

origin3_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = MCW_envir)

origin4_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin4", envir = MCW_envir)

origin5_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin5", envir = MCW_envir)

origin6_year_param_array_MCW <- make_year_origin_parameter_draws_array(fit = MCW_fit, fit_summary = MCW_fit_summary, 
                                                                       origin_select_numeric = "origin6", envir = MCW_envir)

### MCH ###
origin1_year_param_array_MCH <- make_year_origin_parameter_draws_array(fit = MCH_fit, fit_summary = MCH_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = MCH_envir)

origin2_year_param_array_MCH <- make_year_origin_parameter_draws_array(fit = MCH_fit, fit_summary = MCH_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = MCH_envir)

### SRW ###
origin1_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = SRW_envir)

origin2_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = SRW_envir)

origin3_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = SRW_envir)

origin4_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin4", envir = SRW_envir)

origin5_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin5", envir = SRW_envir)

origin6_year_param_array_SRW <- make_year_origin_parameter_draws_array(fit = SRW_fit, fit_summary = SRW_fit_summary, 
                                                                       origin_select_numeric = "origin6", envir = SRW_envir)

### SRH ###
origin1_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin1", envir = SRH_envir)

origin2_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin2", envir = SRH_envir)

origin3_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin3", envir = SRH_envir)

origin4_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin4", envir = SRH_envir)

origin5_year_param_array_SRH <- make_year_origin_parameter_draws_array(fit = SRH_fit, fit_summary = SRH_fit_summary, 
                                                                       origin_select_numeric = "origin5", envir = SRH_envir)





#### Functions to estimate movement probability by year ####

# This function takes an origin and a movement and plots the probability of
# that movement in different years. This uses the median covariate values from
# that year to show the overall probability of movement


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

estimate_year_move_prob_UCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  
  # set up years to predict across
  year_predict <- seq(1,18)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCW_states_dates_years <- data.frame(state = as.vector(UCW_envir$data$y),
                                   date = as.vector(UCW_envir$data$transition_dates))
    
    UCW_states_dates_years$year <- ceiling(UCW_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- UCW_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 18)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(UCW_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- UCW_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 18)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- UCW_envir$data$temperature_data
    
    med_temp <- vector(length = 18)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(UCW_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_UCW[from,to,iter] +
                                                         # btemp0_array_UCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                         btemp1_array_UCW[from,to,iter]*med_temp[j] + 
                                                         bspillwindow_array_UCW[from,to,iter]*med_spillwindow[j] +
                                                         bwinterspill_array_UCW[from,to,iter]*med_winterspill[j] +
                                                         # btemp0xorigin1_array_UCW[from,to,iter]*origin1 +
                                                         btemp1xorigin1_array_UCW[from,to,iter]*med_temp[j]*origin1 +
                                                         # btemp0xorigin2_array_UCW[from,to,iter]*origin2 + 
                                                         btemp1xorigin2_array_UCW[from,to,iter]*med_temp[j]*origin2 + 
                                                         # btemp0xorigin3_array_UCW[from,to,iter]*origin3 +
                                                         btemp1xorigin3_array_UCW[from,to,iter]*med_temp[j]*origin3 +
                                                         borigin1_array_UCW[from,to,iter]*origin1 +
                                                         borigin2_array_UCW[from,to,iter]*origin2 +
                                                         borigin3_array_UCW[from,to,iter]*origin3 +
                                                         
                                                         # year effects
                                                         origin1_year_param_array_UCW[paste0(from, "_", to), j, iter]*origin1 +
                                                         origin2_year_param_array_UCW[paste0(from, "_", to), j, iter]*origin2 +
                                                         origin3_year_param_array_UCW[paste0(from, "_", to), j, iter]*origin3)/
            sum(exp(b0_array_UCW[from,possible_movements,iter] +
                      # btemp0_array_UCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_UCW[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_UCW[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_UCW[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_UCW[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_UCW[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_UCW[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_UCW[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_UCW[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_UCW[from,possible_movements,iter]*med_temp[j]*origin3 +
                      borigin1_array_UCW[from,possible_movements,iter]*origin1 +
                      borigin2_array_UCW[from,possible_movements,iter]*origin2 +
                      borigin3_array_UCW[from,possible_movements,iter]*origin3+
                  # year effects
                  c(origin1_year_param_array_UCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                  c(origin2_year_param_array_UCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                  c(origin3_year_param_array_UCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_UCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  
  # set up years to predict across
  year_predict <- seq(1,18)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCH_states_dates_years <- data.frame(state = as.vector(UCH_envir$data$y),
                                         date = as.vector(UCH_envir$data$transition_dates))
    
    UCH_states_dates_years$year <- ceiling(UCH_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- UCH_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 18)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(UCH_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- UCH_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 18)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- UCH_envir$data$temperature_data
    
    med_temp <- vector(length = 18)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(UCH_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_UCH[from,to,iter] +
                                                    # btemp0_array_UCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_UCH[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_UCH[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_UCH[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_UCH[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_UCH[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_UCH[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_UCH[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_UCH[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_UCH[from,to,iter]*med_temp[j]*origin3 +
                                                    borigin1_array_UCH[from,to,iter]*origin1 +
                                                    borigin2_array_UCH[from,to,iter]*origin2 +
                                                    borigin3_array_UCH[from,to,iter]*origin3 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_UCH[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_UCH[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_UCH[paste0(from, "_", to), j, iter]*origin3)/
            sum(exp(b0_array_UCH[from,possible_movements,iter] +
                      # btemp0_array_UCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_UCH[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_UCH[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_UCH[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_UCH[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_UCH[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_UCH[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_UCH[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_UCH[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_UCH[from,possible_movements,iter]*med_temp[j]*origin3 +
                      borigin1_array_UCH[from,possible_movements,iter]*origin1 +
                      borigin2_array_UCH[from,possible_movements,iter]*origin2 +
                      borigin3_array_UCH[from,possible_movements,iter]*origin3+
                      # year effects
                      c(origin1_year_param_array_UCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_UCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_UCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_MCW <- function(origin_select, movements){
  
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
  
  # set up years to predict across
  year_predict <- seq(1,18)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCW_states_dates_years <- data.frame(state = as.vector(MCW_envir$data$y),
                                         date = as.vector(MCW_envir$data$transition_dates))
    
    MCW_states_dates_years$year <- ceiling(MCW_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- MCW_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 18)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(MCW_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- MCW_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 18)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- MCW_envir$data$temperature_data
    
    med_temp <- vector(length = 18)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(MCW_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_MCW[from,to,iter] +
                                                    # btemp0_array_MCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_MCW[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_MCW[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_MCW[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_MCW[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_MCW[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_MCW[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_MCW[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_MCW[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_MCW[from,to,iter]*med_temp[j]*origin3 +
                                                    # btemp0xorigin4_array_MCW[from,to,iter]*origin4 +
                                                    btemp1xorigin4_array_MCW[from,to,iter]*med_temp[j]*origin4 +
                                                    # btemp0xorigin5_array_MCW[from,to,iter]*origin5 +
                                                    btemp1xorigin5_array_MCW[from,to,iter]*med_temp[j]*origin5 +
                                                    # btemp0xorigin6_array_MCW[from,to,iter]*origin6 +
                                                    btemp1xorigin6_array_MCW[from,to,iter]*med_temp[j]*origin6 +
                                                    borigin1_array_MCW[from,to,iter]*origin1 +
                                                    borigin2_array_MCW[from,to,iter]*origin2 +
                                                    borigin3_array_MCW[from,to,iter]*origin3 +
                                                    borigin4_array_MCW[from,to,iter]*origin4 +
                                                    borigin5_array_MCW[from,to,iter]*origin5 +
                                                    borigin6_array_MCW[from,to,iter]*origin6 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin3 +
                                                    origin4_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin4 +
                                                    origin5_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin5 +
                                                    origin6_year_param_array_MCW[paste0(from, "_", to), j, iter]*origin6)/
            sum(exp(b0_array_MCW[from,possible_movements,iter] +
                      # btemp0_array_MCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_MCW[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_MCW[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_MCW[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_MCW[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_MCW[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_MCW[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_MCW[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_MCW[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_MCW[from,possible_movements,iter]*med_temp[j]*origin3 +
                      # btemp0xorigin4_array_MCW[from,possible_movements,iter]*origin4 +
                      btemp1xorigin4_array_MCW[from,possible_movements,iter]*med_temp[j]*origin4 +
                      # btemp0xorigin5_array_MCW[from,possible_movements,iter]*origin5 +
                      btemp1xorigin5_array_MCW[from,possible_movements,iter]*med_temp[j]*origin5 +
                      # btemp0xorigin6_array_MCW[from,possible_movements,iter]*origin6 +
                      btemp1xorigin6_array_MCW[from,possible_movements,iter]*med_temp[j]*origin6 +
                      borigin1_array_MCW[from,possible_movements,iter]*origin1 +
                      borigin2_array_MCW[from,possible_movements,iter]*origin2 +
                      borigin3_array_MCW[from,possible_movements,iter]*origin3 +
                      borigin4_array_MCW[from,possible_movements,iter]*origin4 +
                      borigin5_array_MCW[from,possible_movements,iter]*origin5 +
                      borigin6_array_MCW[from,possible_movements,iter]*origin6 +
                      # year effects
                      c(origin1_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3 +
                      c(origin4_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin4 +
                      c(origin5_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin5 +
                      c(origin6_year_param_array_MCW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin6))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_MCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  # set up years to predict across
  year_predict <- seq(1,18)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCH_states_dates_years <- data.frame(state = as.vector(MCH_envir$data$y),
                                         date = as.vector(MCH_envir$data$transition_dates))
    
    MCH_states_dates_years$year <- ceiling(MCH_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- MCH_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 18)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(MCH_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- MCH_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 18)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- MCH_envir$data$temperature_data
    
    med_temp <- vector(length = 18)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(MCH_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_MCH[from,to,iter] +
                                                    # btemp0_array_MCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_MCH[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_MCH[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_MCH[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_MCH[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_MCH[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_MCH[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_MCH[from,to,iter]*med_temp[j]*origin2 + 
                                                    borigin1_array_MCH[from,to,iter]*origin1 +
                                                    borigin2_array_MCH[from,to,iter]*origin2 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_MCH[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_MCH[paste0(from, "_", to), j, iter]*origin2)/
            sum(exp(b0_array_MCH[from,possible_movements,iter] +
                      # btemp0_array_MCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_MCH[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_MCH[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_MCH[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_MCH[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_MCH[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_MCH[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_MCH[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      borigin1_array_MCH[from,possible_movements,iter]*origin1 +
                      borigin2_array_MCH[from,possible_movements,iter]*origin2 +
                      # year effects
                      c(origin1_year_param_array_MCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_MCH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_SRW <- function(origin_select, movements){
  
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
  
  # set up years to predict across
  year_predict <- seq(1,18)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRW_states_dates_years <- data.frame(state = as.vector(SRW_envir$data$y),
                                         date = as.vector(SRW_envir$data$transition_dates))
    
    SRW_states_dates_years$year <- ceiling(SRW_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- SRW_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 18)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(SRW_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- SRW_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 18)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- SRW_envir$data$temperature_data
    
    med_temp <- vector(length = 18)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(SRW_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_SRW[from,to,iter] +
                                                    # btemp0_array_SRW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_SRW[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_SRW[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_SRW[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_SRW[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_SRW[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_SRW[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_SRW[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_SRW[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_SRW[from,to,iter]*med_temp[j]*origin3 +
                                                    # btemp0xorigin4_array_SRW[from,to,iter]*origin4 +
                                                    btemp1xorigin4_array_SRW[from,to,iter]*med_temp[j]*origin4 +
                                                    # btemp0xorigin5_array_SRW[from,to,iter]*origin5 +
                                                    btemp1xorigin5_array_SRW[from,to,iter]*med_temp[j]*origin5 +
                                                    # btemp0xorigin6_array_SRW[from,to,iter]*origin6 +
                                                    btemp1xorigin6_array_SRW[from,to,iter]*med_temp[j]*origin6 +
                                                    borigin1_array_SRW[from,to,iter]*origin1 +
                                                    borigin2_array_SRW[from,to,iter]*origin2 +
                                                    borigin3_array_SRW[from,to,iter]*origin3 +
                                                    borigin4_array_SRW[from,to,iter]*origin4 +
                                                    borigin5_array_SRW[from,to,iter]*origin5 +
                                                    borigin6_array_SRW[from,to,iter]*origin6 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin3 +
                                                    origin4_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin4 +
                                                    origin5_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin5 +
                                                    origin6_year_param_array_SRW[paste0(from, "_", to), j, iter]*origin6)/
            sum(exp(b0_array_SRW[from,possible_movements,iter] +
                      # btemp0_array_SRW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_SRW[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_SRW[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_SRW[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_SRW[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_SRW[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_SRW[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_SRW[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_SRW[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_SRW[from,possible_movements,iter]*med_temp[j]*origin3 +
                      # btemp0xorigin4_array_SRW[from,possible_movements,iter]*origin4 +
                      btemp1xorigin4_array_SRW[from,possible_movements,iter]*med_temp[j]*origin4 +
                      # btemp0xorigin5_array_SRW[from,possible_movements,iter]*origin5 +
                      btemp1xorigin5_array_SRW[from,possible_movements,iter]*med_temp[j]*origin5 +
                      # btemp0xorigin6_array_SRW[from,possible_movements,iter]*origin6 +
                      btemp1xorigin6_array_SRW[from,possible_movements,iter]*med_temp[j]*origin6 +
                      borigin1_array_SRW[from,possible_movements,iter]*origin1 +
                      borigin2_array_SRW[from,possible_movements,iter]*origin2 +
                      borigin3_array_SRW[from,possible_movements,iter]*origin3 +
                      borigin4_array_SRW[from,possible_movements,iter]*origin4 +
                      borigin5_array_SRW[from,possible_movements,iter]*origin5 +
                      borigin6_array_SRW[from,possible_movements,iter]*origin6 +
                      # year effects
                      c(origin1_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3 +
                      c(origin4_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin4 +
                      c(origin5_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin5 +
                      c(origin6_year_param_array_SRW[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin6))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}

estimate_year_move_prob_SRH <- function(origin_select, movements){
  
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
  
  # set up years to predict across
  year_predict <- seq(1,18)
  
  # get an array to store probabilities of movements in different years
  niter <- 4000 # this is the number of draws we have
  
  # set this up where rownames = years, columns = iter, slices = movements
  year_move_prob_array <- array(dim = c(length(year_predict), niter, nrow(movements)))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRH_states_dates_years <- data.frame(state = as.vector(SRH_envir$data$y),
                                         date = as.vector(SRH_envir$data$transition_dates))
    
    SRH_states_dates_years$year <- ceiling(SRH_states_dates_years$date/365.25)+1
    
    # get median spill window for this state, by year - this will be indexed the same as year, using j
    spillwindow_data <- SRH_envir$data$spill_window_data
    
    med_spillwindow <- vector(length = 18)
    
    # for spill and temp, if there were no fish in year/state combination,
    # there will be NAs in the vector - and that's ok
    
    for (x in 1:length(med_spillwindow)){
      med_spillwindow[x] <- median(spillwindow_data[subset(SRH_states_dates_years, state == from & year == x)$date,from])
    }
    
    # get median winter spill for this state (note that this is just one value, since each year already only has one winter spill value), 
    # by year - this will be indexed the same as year, using j
    winterspill_data <- SRH_envir$data$winter_spill_days_data
    med_winterspill <- vector(length = 18)
    
    for (y in 1:length(med_winterspill)){
      med_winterspill[y] <- winterspill_data[y,from]
    }
    
    # get median temperature for this state, by year - this will be indexed the same as year, using j
    temp_data <- SRH_envir$data$temperature_data
    
    med_temp <- vector(length = 18)
    
    for (z in 1:length(med_temp)){
      med_temp[z] <- median(temp_data[subset(SRH_states_dates_years, state == from & year == z)$date,from])
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(year_predict)){
        
        # if there aren't any movements in that year, skip it
        if (!(is.na(med_temp[j]))){
          # evaluation movement 
          year_move_prob_array[j, iter, i] <- exp(b0_array_SRH[from,to,iter] +
                                                    # btemp0_array_SRH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                    btemp1_array_SRH[from,to,iter]*med_temp[j] + 
                                                    bspillwindow_array_SRH[from,to,iter]*med_spillwindow[j] +
                                                    bwinterspill_array_SRH[from,to,iter]*med_winterspill[j] +
                                                    # btemp0xorigin1_array_SRH[from,to,iter]*origin1 +
                                                    btemp1xorigin1_array_SRH[from,to,iter]*med_temp[j]*origin1 +
                                                    # btemp0xorigin2_array_SRH[from,to,iter]*origin2 + 
                                                    btemp1xorigin2_array_SRH[from,to,iter]*med_temp[j]*origin2 + 
                                                    # btemp0xorigin3_array_SRH[from,to,iter]*origin3 +
                                                    btemp1xorigin3_array_SRH[from,to,iter]*med_temp[j]*origin3 +
                                                    # btemp0xorigin4_array_SRH[from,to,iter]*origin4 +
                                                    btemp1xorigin4_array_SRH[from,to,iter]*med_temp[j]*origin4 +
                                                    # btemp0xorigin5_array_SRH[from,to,iter]*origin5 +
                                                    btemp1xorigin5_array_SRH[from,to,iter]*med_temp[j]*origin5 +
                                                    borigin1_array_SRH[from,to,iter]*origin1 +
                                                    borigin2_array_SRH[from,to,iter]*origin2 +
                                                    borigin3_array_SRH[from,to,iter]*origin3 +
                                                    borigin4_array_SRH[from,to,iter]*origin4 +
                                                    borigin5_array_SRH[from,to,iter]*origin5 +
                                                    
                                                    # year effects
                                                    origin1_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin1 +
                                                    origin2_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin2 +
                                                    origin3_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin3 +
                                                    origin4_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin4 +
                                                    origin5_year_param_array_SRH[paste0(from, "_", to), j, iter]*origin5)/
            sum(exp(b0_array_SRH[from,possible_movements,iter] +
                      # btemp0_array_SRH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                      btemp1_array_SRH[from,possible_movements,iter]*med_temp[j] + 
                      bspillwindow_array_SRH[from,possible_movements,iter]*med_spillwindow[j] +
                      bwinterspill_array_SRH[from,possible_movements,iter]*med_winterspill[j] +
                      # btemp0xorigin1_array_SRH[from,possible_movements,iter]*origin1 +
                      btemp1xorigin1_array_SRH[from,possible_movements,iter]*med_temp[j]*origin1 +
                      # btemp0xorigin2_array_SRH[from,possible_movements,iter]*origin2 + 
                      btemp1xorigin2_array_SRH[from,possible_movements,iter]*med_temp[j]*origin2 + 
                      # btemp0xorigin3_array_SRH[from,possible_movements,iter]*origin3 +
                      btemp1xorigin3_array_SRH[from,possible_movements,iter]*med_temp[j]*origin3 +
                      # btemp0xorigin4_array_SRH[from,possible_movements,iter]*origin4 +
                      btemp1xorigin4_array_SRH[from,possible_movements,iter]*med_temp[j]*origin4 +
                      # btemp0xorigin5_array_SRH[from,possible_movements,iter]*origin5 +
                      btemp1xorigin5_array_SRH[from,possible_movements,iter]*med_temp[j]*origin5 +
                      borigin1_array_SRH[from,possible_movements,iter]*origin1 +
                      borigin2_array_SRH[from,possible_movements,iter]*origin2 +
                      borigin3_array_SRH[from,possible_movements,iter]*origin3 +
                      borigin4_array_SRH[from,possible_movements,iter]*origin4 +
                      borigin5_array_SRH[from,possible_movements,iter]*origin5 +
                      # year effects
                      c(origin1_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin1 +
                      c(origin2_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin2 +
                      c(origin3_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin3 +
                      c(origin4_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin4 +
                      c(origin5_year_param_array_SRH[paste0(from, "_", possible_movements[1:(length(possible_movements)-1)]), j, iter],0)*origin5))
          
        } else {
          
        }
        
        
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(year_move_prob_array)
  
}


#### Function to plot probability of movements in different years ####

plot_prob_by_year <- function(year_move_prob_array, from, to, plot_title = NULL){
  year_move_prob <- as.data.frame(year_move_prob_array[from, to,,])
  
  niter <- 4000 # for the number of draws
  
  colnames(year_move_prob) <- paste0("iter", 1:niter) 
  year_predict <- 1:18
  year_move_prob$year <- year_predict
  
  # Add a column with the actual temperatures
  year_move_prob$year_actual <- 2004:2021
  
  # drop years without observations
  year_move_prob <- na.omit(year_move_prob)
  
  
  year_move_prob %>% 
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
    group_by(year_actual) %>% 
    summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> year_move_prob_quantiles
  
  
  # start here!
  year_move_prob_plot <- ggplot(year_move_prob_quantiles, aes(x = year_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`)) +
    geom_line() +
    geom_ribbon(alpha = 0.2) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(2003, 2022), breaks = seq(2005, 2020, 5), expand = c(0,0)) +
    xlab("Run Year") +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(year_move_prob_plot)
}
rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
plot_compare_prob_by_year <- function(origin_select,
                                      wild_year_move_prob_array = NULL, hatchery_year_move_prob_array = NULL,
                                      movements_evaluated = NULL,
                                      from, to, plot_title = NULL){
  if(is.null(movements_evaluated)){
    movements_evaluated <- data.frame(from = from, to = to)
  }
  
  niter <- 4000 # for the number of draws
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    wild_year_move_prob <- as.data.frame(wild_year_move_prob_array[,,array_index])
    
    colnames(wild_year_move_prob) <- paste0("iter", 1:niter) 
    year_predict <- 1:18
    wild_year_move_prob$year <- year_predict
    
    # Add a column with the actual years
    wild_year_move_prob$year_actual <- 2004:2021
    
    # drop years without observations
    wild_year_move_prob <- na.omit(wild_year_move_prob)
    
    
    wild_year_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_year_move_prob_quantiles
    
    wild_year_move_prob_quantiles -> rear_year_move_prob_quantiles
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    
    hatchery_year_move_prob <- as.data.frame(hatchery_year_move_prob_array[,,array_index])
    
    colnames(hatchery_year_move_prob) <- paste0("iter", 1:niter) 
    year_predict <- 1:18
    hatchery_year_move_prob$year <- year_predict
    
    # Add a column with the actual temperatures
    hatchery_year_move_prob$year_actual <- 2004:2021
    
    # drop years without observations
    hatchery_year_move_prob <- na.omit(hatchery_year_move_prob)
    
    
    hatchery_year_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_year_move_prob_quantiles
    
    # join hatchery and wild
    hatchery_year_move_prob_quantiles -> rear_year_move_prob_quantiles
  }
  # else run both
  else {
    wild_year_move_prob <- as.data.frame(wild_year_move_prob_array[,,array_index])
    
    colnames(wild_year_move_prob) <- paste0("iter", 1:niter) 
    year_predict <- 1:18
    wild_year_move_prob$year <- year_predict
    
    # Add a column with the actual temperatures
    wild_year_move_prob$year_actual <- 2004:2021
    
    # drop years without observations
    wild_year_move_prob <- na.omit(wild_year_move_prob)
    
    
    wild_year_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_year_move_prob_quantiles
    
    hatchery_year_move_prob <- as.data.frame(hatchery_year_move_prob_array[,,array_index])
    
    colnames(hatchery_year_move_prob) <- paste0("iter", 1:niter) 
    year_predict <- 1:18
    hatchery_year_move_prob$year <- year_predict
    
    # Add a column with the actual temperatures
    hatchery_year_move_prob$year_actual <- 2004:2021
    
    # drop years without observations
    hatchery_year_move_prob <- na.omit(hatchery_year_move_prob)
    
    
    hatchery_year_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_year_move_prob_quantiles
    
    # join hatchery and wild
    wild_year_move_prob_quantiles %>% 
      bind_rows(., hatchery_year_move_prob_quantiles) -> rear_year_move_prob_quantiles
    
  }
  
  
  rear_year_move_prob_plot <- ggplot(rear_year_move_prob_quantiles, aes(x = year_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`,
                                                                        color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(2003, 2022), breaks = seq(2005, 2020, 5), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab("Run Year") +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(rear_year_move_prob_plot)
}

#### Function to plot random year effects on their own ####
# separately, this function plots only the parameter value for the random effect of year

plot_RE_year <- function(wild_year_param_array = NULL,
                         hatchery_year_param_array = NULL,
                         from, to, plot_title = NULL){
  rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
  year_predict <- seq(1,18)
  
  if(is.null(hatchery_year_param_array)){
    wild_RE_df <- as.data.frame(wild_year_param_array[paste0(from, "_", to), ,])
    
    niter <- 4000 # for the number of draws
    
    colnames(wild_RE_df) <- paste0("iter", 1:niter) 
    wild_year_predict <- 1:18
    wild_RE_df$year <- year_predict
    
    # Add a column with the actual temperatures
    wild_RE_df$year_actual <- 2004:2021
    
    # drop years without observations
    wild_RE_df <- na.omit(wild_RE_df)
    
    
    wild_RE_df %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) -> wild_RE_df_quantiles
    
    wild_RE_df_quantiles %>% 
      mutate(rear = "wild") -> RE_df_quantiles
    
  } else if(is.null(wild_year_param_array)){
    hatchery_RE_df <- as.data.frame(hatchery_year_param_array[paste0(from, "_", to), ,])
    
    niter <- 4000 # for the number of draws
    
    colnames(hatchery_RE_df) <- paste0("iter", 1:niter) 
    hatchery_year_predict <- 1:18
    hatchery_RE_df$year <- year_predict
    
    # Add a column with the actual temperatures
    hatchery_RE_df$year_actual <- 2004:2021
    
    # drop years without observations
    hatchery_RE_df <- na.omit(hatchery_RE_df)
    
    
    hatchery_RE_df %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) -> hatchery_RE_df_quantiles
    
    hatchery_RE_df_quantiles %>% 
      mutate(rear = "hatchery") -> RE_df_quantiles
    
  } else {
    wild_RE_df <- as.data.frame(wild_year_param_array[paste0(from, "_", to), ,])
    
    niter <- 4000 # for the number of draws
    
    colnames(wild_RE_df) <- paste0("iter", 1:niter) 
    wild_year_predict <- 1:18
    wild_RE_df$year <- year_predict
    
    # Add a column with the actual temperatures
    wild_RE_df$year_actual <- 2004:2021
    
    # drop years without observations
    wild_RE_df <- na.omit(wild_RE_df)
    
    
    wild_RE_df %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) -> wild_RE_df_quantiles
    
    wild_RE_df_quantiles %>% 
      mutate(rear = "wild") -> wild_RE_df_quantiles
    
    hatchery_RE_df <- as.data.frame(hatchery_year_param_array[paste0(from, "_", to), ,])
    
    niter <- 4000 # for the number of draws
    
    colnames(hatchery_RE_df) <- paste0("iter", 1:niter) 
    hatchery_year_predict <- 1:18
    hatchery_RE_df$year <- year_predict
    
    # Add a column with the actual temperatures
    hatchery_RE_df$year_actual <- 2004:2021
    
    # drop years without observations
    hatchery_RE_df <- na.omit(hatchery_RE_df)
    
    
    hatchery_RE_df %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(year_actual) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) -> hatchery_RE_df_quantiles
    
    hatchery_RE_df_quantiles %>% 
      mutate(rear = "hatchery") -> hatchery_RE_df_quantiles
    
    wild_RE_df_quantiles %>% 
      bind_rows(., hatchery_RE_df_quantiles) -> RE_df_quantiles
    
  }
  
  
  
  # start here!
  RE_year_plot <- ggplot(RE_df_quantiles, aes(x = year_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`,
                                              color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    # scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(2003, 2022), breaks = seq(2005, 2020, 5), expand = c(0,0)) +
    xlab("Run Year") +
    ylab("Parameter estimate") +
    ggtitle(plot_title)
  
  return(RE_year_plot)
}

#### Plot movement prob by year for different origins ####

#### Upper Columbia ####

## Wenatchee River
wen_wild_year_movement_probs <- estimate_year_move_prob_UCW(origin_select = "Wenatchee River", movements = data.frame(from = c(5), to = c(24)))

wen_hatchery_year_movement_probs <- estimate_year_move_prob_UCH(origin_select = "Wenatchee River", movements = data.frame(from = c(5), to = c(24)))

# comparison plot
wen_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Wenatchee River",
                                                                   wild_year_move_prob_array = wen_wild_year_movement_probs, 
                                                           hatchery_year_move_prob_array = wen_hatchery_year_movement_probs,
                                                               from = 5, to = 24, plot_title = "Wenatchee - movement into Wenatchee")

ggsave(here::here("stan_actual", "output", "annual", "wen_rear_homing_movement_by_year_plot.png"), wen_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Entiat River
ent_wild_year_movement_probs <- estimate_year_move_prob_UCW(origin_select = "Entiat River", movements = data.frame(from = c(6), to = c(26)))

# comparison plot
ent_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Entiat River", 
                                                                   wild_year_move_prob_array = ent_wild_year_movement_probs, 
                                                                   from = 6, to = 26, plot_title = "Entiat - movement into Entiat")

ggsave(here::here("stan_actual", "output", "annual", "ent_rear_homing_movement_by_year_plot.png"), ent_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Okanogan River

oka_hatchery_year_movement_probs <- estimate_year_move_prob_UCH(origin_select = "Okanogan River", movements = data.frame(from = c(7), to = c(28)))

# comparison plot
oka_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Okanogan River",
                                                                   hatchery_year_move_prob_array = oka_hatchery_year_movement_probs,
                                                                   from = 7, to = 28, plot_title = "Okanogan - movement into Okanogan")

ggsave(here::here("stan_actual", "output", "annual", "oka_rear_homing_movement_by_year_plot.png"), oka_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Methow River
met_wild_year_movement_probs <- estimate_year_move_prob_UCW(origin_select = "Methow River", movements = data.frame(from = c(7), to = c(30)))

met_hatchery_year_movement_probs <- estimate_year_move_prob_UCH(origin_select = "Methow River", movements = data.frame(from = c(7), to = c(30)))

# comparison plot
met_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Methow River",
                                                                   wild_year_move_prob_array = met_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = met_hatchery_year_movement_probs,
                                                                   from = 7, to = 30, plot_title = "Methow - movement into Methow")

ggsave(here::here("stan_actual", "output", "annual", "met_rear_homing_movement_by_year_plot.png"), met_rear_homing_movement_by_year_plot, height = 5, width = 8)

#### Middle Columbia ####
## Deschutes River
des_wild_year_movement_probs <- estimate_year_move_prob_MCW(origin_select = "Deschutes River", movements = data.frame(from = c(2), to = c(10)))

# comparison plot
des_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Deschutes River",
                                                                   wild_year_move_prob_array = des_wild_year_movement_probs, 
                                                                   from = 2, to = 10, plot_title = "Deschutes - movement into Deschutes")

ggsave(here::here("stan_actual", "output", "annual", "des_rear_homing_movement_by_year_plot.png"), des_rear_homing_movement_by_year_plot, height = 5, width = 8)

## John Day River
jdr_wild_year_movement_probs <- estimate_year_move_prob_MCW(origin_select = "John Day River", movements = data.frame(from = c(2), to = c(12)))

# comparison plot
jdr_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "John Day River",
                                                                   wild_year_move_prob_array = jdr_wild_year_movement_probs, 
                                                                   from = 2, to = 12, plot_title = "John Day - movement into John Day")

ggsave(here::here("stan_actual", "output", "annual", "jdr_rear_homing_movement_by_year_plot.png"), jdr_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Fifteenmile Creek
fif_wild_year_movement_probs <- estimate_year_move_prob_MCW(origin_select = "Fifteenmile Creek", movements = data.frame(from = c(2), to = c(16)))

# comparison plot
fif_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Fifteenmile Creek",
                                                                   wild_year_move_prob_array = fif_wild_year_movement_probs, 
                                                                   from = 2, to = 16, plot_title = "Fifteenmile Creek - movement into Fifteenmile Creek")

ggsave(here::here("stan_actual", "output", "annual", "fif_rear_homing_movement_by_year_plot.png"), fif_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Umatilla River
uma_wild_year_movement_probs <- estimate_year_move_prob_MCW(origin_select = "Umatilla River", movements = data.frame(from = c(2), to = c(18)))

uma_hatchery_year_movement_probs <- estimate_year_move_prob_MCH(origin_select = "Umatilla River", movements = data.frame(from = c(2), to = c(18)))

# comparison plot
uma_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Umatilla River",
                                                                   wild_year_move_prob_array = uma_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = uma_hatchery_year_movement_probs,
                                                                   from = 2, to = 18, plot_title = "Umatilla - movement into Umatilla")

ggsave(here::here("stan_actual", "output", "annual", "uma_rear_homing_movement_by_year_plot.png"), uma_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Yakima River
yak_wild_year_movement_probs <- estimate_year_move_prob_MCW(origin_select = "Yakima River", movements = data.frame(from = c(3), to = c(20)))

yak_hatchery_year_movement_probs <- estimate_year_move_prob_MCH(origin_select = "Yakima River", movements = data.frame(from = c(3), to = c(20)))

# comparison plot
yak_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Yakima River",
                                                                   wild_year_move_prob_array = yak_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = yak_hatchery_year_movement_probs,
                                                                   from = 3, to = 20, plot_title = "Yakima - movement into Yakima")

ggsave(here::here("stan_actual", "output", "annual", "yak_rear_homing_movement_by_year_plot.png"), yak_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Walla Walla River
wawa_wild_year_movement_probs <- estimate_year_move_prob_MCW(origin_select = "Walla Walla River", movements = data.frame(from = c(3), to = c(22)))

wawa_hatchery_year_movement_probs <- estimate_year_move_prob_MCH(origin_select = "Walla Walla River", movements = data.frame(from = c(3), to = c(22)))

# comparison plot
wawa_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Walla Walla River",
                                                                   wild_year_move_prob_array = wawa_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = wawa_hatchery_year_movement_probs,
                                                                   from = 3, to = 22, plot_title = "Walla Walla - movement into Walla Walla")

ggsave(here::here("stan_actual", "output", "annual", "wawa_rear_homing_movement_by_year_plot.png"), wawa_rear_homing_movement_by_year_plot, height = 5, width = 8)


#### Snake River ####
## Asotin Creek
aso_wild_year_movement_probs <- estimate_year_move_prob_SRW(origin_select = "Asotin Creek", movements = data.frame(from = c(9), to = c(34)))

# comparison plot
aso_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Asotin Creek",
                                                                    wild_year_move_prob_array = aso_wild_year_movement_probs, 
                                                                    from = 9, to = 34, plot_title = "Asotin Creek - movement into Asotin Creek")

ggsave(here::here("stan_actual", "output", "annual", "aso_rear_homing_movement_by_year_plot.png"), aso_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Clearwater River
cle_wild_year_movement_probs <- estimate_year_move_prob_SRW(origin_select = "Clearwater River", movements = data.frame(from = c(9), to = c(36)))

cle_hatchery_year_movement_probs <- estimate_year_move_prob_SRH(origin_select = "Clearwater River", movements = data.frame(from = c(9), to = c(36)))

# comparison plot
cle_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Clearwater River",
                                                                   wild_year_move_prob_array = cle_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = cle_hatchery_year_movement_probs,
                                                                   from = 9, to = 36, plot_title = "Clearwater - movement into Clearwater")

ggsave(here::here("stan_actual", "output", "annual", "cle_rear_homing_movement_by_year_plot.png"), cle_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Imnaha River
imn_wild_year_movement_probs <- estimate_year_move_prob_SRW(origin_select = "Imnaha River", movements = data.frame(from = c(9), to = c(39)))
imn_hatchery_year_movement_probs <- estimate_year_move_prob_SRH(origin_select = "Imnaha River", movements = data.frame(from = c(9), to = c(39)))

# comparison plot
imn_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Imnaha River",
                                                                   wild_year_move_prob_array = imn_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = imn_hatchery_year_movement_probs,
                                                                   from = 9, to = 39, plot_title = "Imnaha - movement into Imnaha")

ggsave(here::here("stan_actual", "output", "annual", "imn_rear_homing_movement_by_year_plot.png"), imn_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Grande Ronde River
gr_wild_year_movement_probs <- estimate_year_move_prob_SRW(origin_select = "Grande Ronde River", movements = data.frame(from = c(9), to = c(38)))
gr_hatchery_year_movement_probs <- estimate_year_move_prob_SRH(origin_select = "Grande Ronde River", movements = data.frame(from = c(9), to = c(38)))

# comparison plot
gr_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Grande Ronde River",
                                                                   wild_year_move_prob_array = gr_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = gr_hatchery_year_movement_probs,
                                                                   from = 9, to = 38, plot_title = "Grande Ronde - movement into Grande Ronde")

ggsave(here::here("stan_actual", "output", "annual", "gr_rear_homing_movement_by_year_plot.png"), gr_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Salmon River
sal_wild_year_movement_probs <- estimate_year_move_prob_SRW(origin_select = "Salmon River", movements = data.frame(from = c(9), to = c(37)))
sal_hatchery_year_movement_probs <- estimate_year_move_prob_SRH(origin_select = "Salmon River", movements = data.frame(from = c(9), to = c(37)))

# comparison plot
sal_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Salmon River",
                                                                   wild_year_move_prob_array = sal_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = sal_hatchery_year_movement_probs,
                                                                   from = 9, to = 37, plot_title = "Salmon - movement into Salmon")

ggsave(here::here("stan_actual", "output", "annual", "sal_rear_homing_movement_by_year_plot.png"), sal_rear_homing_movement_by_year_plot, height = 5, width = 8)

## Tucannon River
tuc_wild_year_movement_probs <- estimate_year_move_prob_SRW(origin_select = "Tucannon River", movements = data.frame(from = c(8), to = c(32)))
tuc_hatchery_year_movement_probs <- estimate_year_move_prob_SRH(origin_select = "Tucannon River", movements = data.frame(from = c(8), to = c(32)))

# comparison plot
tuc_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(origin_select = "Tucannon River",
                                                                   wild_year_move_prob_array = tuc_wild_year_movement_probs, 
                                                                   hatchery_year_move_prob_array = tuc_hatchery_year_movement_probs,
                                                                   from = 8, to = 32, plot_title = "Tucannon - movement into Tucannon")

ggsave(here::here("stan_actual", "output", "annual", "tuc_rear_homing_movement_by_year_plot.png"), tuc_rear_homing_movement_by_year_plot, height = 5, width = 8)








#### Plot year effect parameter estimates for different origins ####
# use origin_param_map to determine which parameters to select when plotting individual REs

#### Upper Columbia ####
# Wenatchee River
subset(origin_param_map, natal_origin == "Wenatchee River")

mainstem_WEN_RE_year_plot <- plot_RE_year(wild_year_param_array = origin1_year_param_array_UCW, 
                                          hatchery_year_param_array = origin1_year_param_array_UCH,
                                          from = 5, to = 24, plot_title = "Year effect for Wenatchee: mainstem into Wenatchee")


ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_WEN_RE_year_plot.png"), mainstem_WEN_RE_year_plot, height = 5, width = 8)

# this is interesting, but what seems to be happening is that maybe there are correlations with
# other parameters elsewhere

# Entiat River
subset(origin_param_map, natal_origin == "Entiat River")

mainstem_ENT_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin2_year_param_array_UCW, from = 6, to = 26, plot_title = "Year effect for wild Entiat: mainstem into Entiat")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_ENT_wild_RE_year_plot.png"), mainstem_ENT_wild_RE_year_plot, height = 5, width = 8)

# Okanogan River
subset(origin_param_map, natal_origin == "Okanogan River")

mainstem_OKA_hatchery_RE_year_plot <- plot_RE_year(hatchery_year_param_array = origin2_year_param_array_UCH, from = 7, to = 28, plot_title = "Year effect for hatchery Okanogan: mainstem into Okanogan")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_OKA_hatchery_RE_year_plot.png"), mainstem_OKA_hatchery_RE_year_plot, height = 5, width = 8)

# Methow River
subset(origin_param_map, natal_origin == "Methow River")

mainstem_MET_RE_year_plot <- plot_RE_year(wild_year_param_array = origin3_year_param_array_UCW,
                                          hatchery_year_param_array = origin3_year_param_array_UCH,
                                          from = 7, to = 30, plot_title = "Year effect for Methow: mainstem into Methow")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_MET_RE_year_plot.png"), mainstem_MET_RE_year_plot, height = 5, width = 8)


#### Middle Columbia ####
# Deschutes River
subset(origin_param_map, natal_origin == "Deschutes River")

mainstem_DES_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin1_year_param_array_MCW, from = 2, to = 10, plot_title = "Year effect for wild Deschutes: mainstem into Deschutes")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_DES_wild_RE_year_plot.png"), mainstem_DES_wild_RE_year_plot, height = 5, width = 8)

# John Day River
subset(origin_param_map, natal_origin == "John Day River")

mainstem_JDR_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin3_year_param_array_MCW, from = 2, to = 12, plot_title = "Year effect for wild John Day: mainstem into John Day")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_JDR_wild_RE_year_plot.png"), mainstem_JDR_wild_RE_year_plot, height = 5, width = 8)


# Fifteenmile Creek
subset(origin_param_map, natal_origin == "Fifteenmile Creek")

mainstem_FIF_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin2_year_param_array_MCW, from = 2, to = 16, plot_title = "Year effect for wild Fifteenmile: mainstem into Fifteenmile")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_FIF_wild_RE_year_plot.png"), mainstem_FIF_wild_RE_year_plot, height = 5, width = 8)

# Umatilla River
subset(origin_param_map, natal_origin == "Umatilla River")

mainstem_UMA_RE_year_plot <- plot_RE_year(wild_year_param_array = origin4_year_param_array_MCW, 
                                          hatchery_year_param_array = origin1_year_param_array_MCH,
                                          from = 2, to = 18, plot_title = "Year effect for Umatilla: mainstem into Umatilla")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_UMA_RE_year_plot.png"), mainstem_UMA_RE_year_plot, height = 5, width = 8)

# Yakima River
subset(origin_param_map, natal_origin == "Yakima River")

mainstem_YAK_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin6_year_param_array_MCW, from = 3, to = 20, plot_title = "Year effect for wild Yakima: mainstem into Yakima")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_YAK_wild_RE_year_plot.png"), mainstem_YAK_wild_RE_year_plot, height = 5, width = 8)

# Walla Walla River
subset(origin_param_map, natal_origin == "Walla Walla River")

mainstem_WAWA_RE_year_plot <- plot_RE_year(wild_year_param_array = origin5_year_param_array_MCW, 
                                           hatchery_year_param_array = origin2_year_param_array_MCH,
                                           from = 3, to = 22, plot_title = "Year effect for Walla Walla: mainstem into Walla Walla")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_WAWA_RE_year_plot.png"), mainstem_WAWA_RE_year_plot, height = 5, width = 8)


#### Snake River ####

# Tucannon River
subset(origin_param_map, natal_origin == "Tucannon River")

mainstem_TUC_RE_year_plot <- plot_RE_year(wild_year_param_array = origin6_year_param_array_SRW, 
                                          hatchery_year_param_array = origin5_year_param_array_SRH,
                                          from = 8, to = 32, plot_title = "Year effect for Tucannon: mainstem into Tucannon")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_TUC_RE_year_plot.png"), mainstem_TUC_RE_year_plot, height = 5, width = 8)

# Asotin Creek
subset(origin_param_map, natal_origin == "Asotin Creek")

mainstem_ASO_wild_RE_year_plot <- plot_RE_year(wild_year_param_array = origin1_year_param_array_SRW, from = 9, to = 34, plot_title = "Year effect for wild Asotin: mainstem into Asotin")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_ASO_wild_RE_year_plot.png"), mainstem_ASO_wild_RE_year_plot, height = 5, width = 8)


# Clearwater River
subset(origin_param_map, natal_origin == "Clearwater River")

mainstem_CLE_RE_year_plot <- plot_RE_year(wild_year_param_array = origin2_year_param_array_SRW,
                                          hatchery_year_param_array = origin1_year_param_array_SRH,
                                          from = 9, to = 36, plot_title = "Year effect for Clearwater: mainstem into Clearwater")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_CLE_RE_year_plot.png"), mainstem_CLE_RE_year_plot, height = 5, width = 8)


# Imnaha River
subset(origin_param_map, natal_origin == "Imnaha River")

mainstem_IMN_RE_year_plot <- plot_RE_year(wild_year_param_array = origin4_year_param_array_SRW, 
                                          hatchery_year_param_array = origin3_year_param_array_SRH,
                                          from = 9, to = 39, plot_title = "Year effect for Imnaha: mainstem into Imnaha")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_IMN_RE_year_plot.png"), mainstem_IMN_RE_year_plot, height = 5, width = 8)



# Grande Ronde River
subset(origin_param_map, natal_origin == "Grande Ronde River")

mainstem_GR_RE_year_plot <- plot_RE_year(wild_year_param_array = origin3_year_param_array_SRW, 
                                              hatchery_year_param_array = origin2_year_param_array_SRH,
                                              from = 9, to = 38, plot_title = "Year effect for Grande Ronde: mainstem into Grande Ronde")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_GR_RE_year_plot.png"), mainstem_GR_RE_year_plot, height = 5, width = 8)


# Salmon River
subset(origin_param_map, natal_origin == "Salmon River")

mainstem_SAL_RE_year_plot <- plot_RE_year(wild_year_param_array = origin5_year_param_array_SRW,
                                               hatchery_year_param_array = origin4_year_param_array_SRH,
                                               from = 9, to = 37, plot_title = "Year effect for Salmon: mainstem into Salmon")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_SAL_RE_year_plot.png"), mainstem_SAL_RE_year_plot, height = 5, width = 8)








