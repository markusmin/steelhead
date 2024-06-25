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
  movement_indices <- filter(param_indices, from_to %in% year_raw_params_indices$from_to)
  movement_indices %>% 
    dplyr::rename(from = row, to = col) -> movement_indices
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
    left_join(., movement_indices, by = "index") %>% 
    filter(!(is.na(from))) -> year_origin_params_indices
  
  
  # arrange all parameter values into an array
  # 0 is meaningful for loss (this is what is used in the stan code)
  # 0s for all other movements that are not overwritten will not be used
  year_origin_param_array <- array(data = 0, dim = c(length(model_states), length(model_states),
                                         length(unique(year_origin_params_indices$year)),
                                         length(as.matrix(fit[,,1]))))
  

  # now populate the array with the iteration draws
  for(i in 1:nrow(year_origin_params_indices)){
    year_origin_param_array[year_origin_params_indices[i, "from"], year_origin_params_indices[i, "to"], year_origin_params_indices[i, "year"], ] <- as.matrix(fit[,,year_origin_params_indices[i, "parameter"]])
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



#### Functions to estimate movement probability by year ####

# This function takes an origin and a movement and plots the probability of
# that movement in different years. This uses the median covariate values from
# that year to show the overall probability of movement


# create a list that maps origin numbers (params) to what they actually are
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

# Use the parameter map to index the right effects
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
  
  # get an array to store probabilities of movements at different temperatures
  niter <- 4000 # this is the number of draws we have
  year_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(year_predict), niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
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
          year_move_prob_array[from,to,j, iter] <- exp(b0_array_UCW[from,to,iter] +
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
                                                         origin1_year_param_array_UCW[from, to, j, iter]*origin1 +
                                                         origin2_year_param_array_UCW[from, to, j, iter]*origin2 +
                                                         origin3_year_param_array_UCW[from, to, j, iter]*origin3)/
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
                  origin1_year_param_array_UCW[from, possible_movements, j, iter]*origin1 +
                  origin2_year_param_array_UCW[from, possible_movements, j, iter]*origin2 +
                  origin3_year_param_array_UCW[from, possible_movements, j, iter]*origin3))
          
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
  
  # get an array to store probabilities of movements at different temperatures
  niter <- 4000 # this is the number of draws we have
  year_move_prob_array <- array(dim = c(length(model_states), length(model_states), length(year_predict), niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
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
          year_move_prob_array[from,to,j, iter] <- exp(b0_array_UCH[from,to,iter] +
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
                                                         origin1_year_param_array_UCH[from, to, j, iter]*origin1 +
                                                         origin2_year_param_array_UCH[from, to, j, iter]*origin2 +
                                                         origin3_year_param_array_UCH[from, to, j, iter]*origin3)/
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
                      origin1_year_param_array_UCH[from, possible_movements, j, iter]*origin1 +
                      origin2_year_param_array_UCH[from, possible_movements, j, iter]*origin2 +
                      origin3_year_param_array_UCH[from, possible_movements, j, iter]*origin3))
          
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
plot_compare_prob_by_year <- function(wild_year_move_prob_array, hatchery_year_move_prob_array,
                                      from, to, plot_title = NULL){
  niter <- 4000 # for the number of draws
  
  wild_year_move_prob <- as.data.frame(wild_year_move_prob_array[from, to,,])
  
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
  
  hatchery_year_move_prob <- as.data.frame(hatchery_year_move_prob_array[from, to,,])
  
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

plot_RE_year <- function(year_param_array, from, to, plot_title = NULL){
  RE_df <- as.data.frame(year_param_array[from, to, ,])
  
  niter <- 4000 # for the number of draws
  
  colnames(RE_df) <- paste0("iter", 1:niter) 
  year_predict <- 1:18
  RE_df$year <- year_predict
  
  # Add a column with the actual temperatures
  RE_df$year_actual <- 2004:2021
  
  # drop years without observations
  RE_df <- na.omit(RE_df)
  
  
  RE_df %>% 
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
    group_by(year_actual) %>% 
    summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> RE_df_quantiles
  
  
  # start here!
  RE_year_plot <- ggplot(RE_df_quantiles, aes(x = year_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`)) +
    geom_line() +
    geom_ribbon(alpha = 0.2) +
    # scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(2003, 2022), breaks = seq(2005, 2020, 5), expand = c(0,0)) +
    xlab("Run Year") +
    ylab("Parameter estimate") +
    ggtitle(plot_title)
  
  return(RE_year_plot)
}

#### Plot movement prob by year for different origins ####

### UCW ###

## Wenatchee River
# Wild
wen_wild_year_movement_probs <- estimate_year_move_prob_UCW(origin_select = "Wenatchee River", movements = data.frame(from = c(5), to = c(24)))


wen_wild_homing_movement_by_year_plot <- plot_prob_by_year(year_move_prob_array = wen_wild_year_movement_probs, 
                                                      from = 5, to = 24, plot_title = "Wenatchee Wild - movement into Wenatchee")

ggsave(here::here("stan_actual", "output", "year_effects", "wen_wild_homing_movement_by_year_plot.png"), wen_wild_homing_movement_by_year_plot, height = 8, width = 8)

# Hatchery
wen_hatchery_year_movement_probs <- estimate_year_move_prob_UCH(origin_select = "Wenatchee River", movements = data.frame(from = c(5), to = c(24)))


wen_hatchery_homing_movement_by_year_plot <- plot_prob_by_year(year_move_prob_array = wen_hatchery_year_movement_probs, 
                                                               from = 5, to = 24, plot_title = "Wenatchee hatchery - movement into Wenatchee")

ggsave(here::here("stan_actual", "output", "year_effects", "wen_hatchery_homing_movement_by_year_plot.png"), wen_hatchery_homing_movement_by_year_plot, height = 8, width = 8)

# comparison plot
wen_rear_homing_movement_by_year_plot <- plot_compare_prob_by_year(wild_year_move_prob_array = wen_wild_year_movement_probs, 
                                                           hatchery_year_move_prob_array = wen_hatchery_year_movement_probs,
                                                               from = 5, to = 24, plot_title = "Wenatchee - movement into Wenatchee")

ggsave(here::here("stan_actual", "output", "year_effects", "wen_rear_homing_movement_by_year_plot.png"), wen_rear_homing_movement_by_year_plot, height = 8, width = 8)














#### Plot year effect parameter estimates for different origins ####
# use origin_param_map to determine which parameters to select when plotting individual REs

### UCW
# Wenatchee River
subset(origin_param_map, natal_origin == "Wenatchee River")

mainstem_WEN_wild_RE_year_plot <- plot_RE_year(year_param_array = origin1_year_param_array_UCW, from = 5, to = 24, plot_title = "Year effect for wild: moving from mainstem into Wenatchee")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_WEN_wild_RE_year_plot.png"), mainstem_WEN_wild_RE_year_plot, height = 8, width = 8)

mainstem_WEN_hatchery_RE_year_plot <- plot_RE_year(year_param_array = origin1_year_param_array_UCH, from = 5, to = 24, plot_title = "Year effect for hatchery: moving from mainstem into Wenatchee")
ggsave(here::here("stan_actual", "output", "year_effects", "mainstem_WEN_hatchery_RE_year_plot.png"), mainstem_WEN_hatchery_RE_year_plot, height = 8, width = 8)

# this is interesting, but what seems to be happening is that maybe there are correlations with
# other parameters elsewhere




