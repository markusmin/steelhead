# 07-model-results

# This script takes the output from the stan model runs in 05-stan-runs and
# generates some results figures

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
library(ggthemes)

window_temp_summary <- read.csv(here::here("stan_actual", "window_temps_summary.csv"), row.names = 1)
actual_temps <- seq(-2,2,0.01)*window_temp_summary$BON_sd + window_temp_summary$BON_mean

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

# make the origin param map to help with interpretation of parameters
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(grepl("other tributaries", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

origin_param_map <- data.frame(
  natal_origin = natal_origins,
  hatchery = c(NA, NA, NA, NA, 1, NA, 2, # MC
               1,NA,2,3, # UC
               5,NA,1,4,2,3), # SR,
  wild = c(1,3,NA,2,4,6,5, # MC
           1,2,NA,3, # UC
           6,1,2,5,3,4)) # SR

##### Figure 1: Hatchery vs. wild parameter estimates #####

# For each DPS, the exact same parameters (with the exception of the origin parameters)
# are being estimated for both hatchery and wild.


# write a function that plots posteriors for one parameter, estimated for
# two populations

plot_posterior_comparion <- function(parameter, hatchery_fit, wild_fit){
  as.matrix(hatchery_fit[,,parameter]) %>% 
    as.data.frame() %>% 
    dplyr::rename(hatchery = V1) -> hatchery_draws
  
  as.matrix(wild_fit[,,parameter]) %>% 
    as.data.frame() %>% 
    dplyr::rename(wild = V1) -> wild_draws
  
  hatchery_draws %>% 
    bind_cols(wild_draws) %>% 
    pivot_longer(cols = c("hatchery", "wild"), names_to = "rear_type", values_to = parameter) -> draws_comparison
    
  posterior_comparison_plot <- ggplot(draws_comparison, aes(x = eval(parse(text = parameter)), color = rear_type, fill = rear_type)) +
    geom_density(alpha = 0.1) +
    xlab(parameter) +
    scale_color_manual(values = c("hatchery" = "#1f78b4", wild = "#33a02c")) +
    scale_fill_manual(values = c("hatchery" = "#1f78b4", wild = "#33a02c"))
  
  return(posterior_comparison_plot)
  
}

# second function: for comparing origin-specific movements
plot_posterior_comparion_origin <- function(hatchery_parameter, wild_parameter, hatchery_fit, wild_fit){
  as.matrix(hatchery_fit[,,hatchery_parameter]) %>% 
    as.data.frame() %>% 
    dplyr::rename(hatchery = V1) -> hatchery_draws
  
  as.matrix(wild_fit[,,wild_parameter]) %>% 
    as.data.frame() %>% 
    dplyr::rename(wild = V1) -> wild_draws
  
  hatchery_draws %>% 
    bind_cols(wild_draws) %>% 
    pivot_longer(cols = c("hatchery", "wild"), names_to = "rear_type", values_to = "draws") -> draws_comparison
  
  posterior_comparison_plot <- ggplot(draws_comparison, aes(x = draws, color = rear_type, fill = rear_type)) +
    geom_density(alpha = 0.1) +
    xlab(paste0(hatchery_parameter, "/", wild_parameter)) +
    scale_color_manual(values = c("hatchery" = "#1f78b4", wild = "#33a02c")) +
    scale_fill_manual(values = c("hatchery" = "#1f78b4", wild = "#33a02c"))
  
  return(posterior_comparison_plot)
  
}

# third function: plot single posterior
plot_single_posterior <- function(parameter, fit){
  as.matrix(fit[,,parameter]) %>% 
    as.data.frame() %>% 
    dplyr::rename(draws = V1) -> parameter_draws
  
  posterior_comparison_plot <- ggplot(parameter_draws, aes(x = draws)) +
    geom_density(alpha = 0.1) +
    xlab(parameter)
  
  return(posterior_comparison_plot)
  
}


# Posterior comparisons: Upper Columbia
UCW_parameters <- UCW_fit_summary$variable
UCW_parameters[grepl("b0|btemp1_|btemp0_|bspillwindow|bwinterspill", UCW_parameters)] -> UCW_parameters_2
UCW_parameters_2[!(grepl("vector", UCW_parameters_2))] -> UCW_parameters_to_compare

UCH_parameters <- UCH_fit_summary$variable
UCH_parameters[grepl("b0|btemp1_|btemp0_|bspillwindow|bwinterspill", UCH_parameters)] -> UCH_parameters_2
UCH_parameters_2[!(grepl("vector", UCH_parameters_2))] -> UCH_parameters_to_compare

UC_parameters_to_compare <- intersect(UCH_parameters_to_compare, UCW_parameters_to_compare)

UC_parameters_comparison_plotlist <- vector(mode = "list", length = length(UC_parameters_to_compare))

for (i in 1:length(UC_parameters_to_compare)){
  UC_parameters_comparison_plotlist[[i]] <- plot_posterior_comparion(parameter = UC_parameters_to_compare[i], 
                                                                      hatchery_fit = UCH_fit, wild_fit = UCW_fit)
  print(paste0("i = ", i))
}

for (i in 1:length(UC_parameters_to_compare)){
  ggsave(file = paste0("stan_actual/output/posterior_comparison_plots/UC/", UC_parameters_to_compare[[i]], 
                       "_comp_plot.png"), plot = UC_parameters_comparison_plotlist[[i]],
         height = 8, width = 8)
  
  print(paste0("i = ", i))
}

# Posterior comparisons: Middle Columbia

MCW_parameters <- MCW_fit_summary$variable
MCW_parameters[grepl("b0|btemp1_|btemp0_|bspillwindow|bwinterspill", MCW_parameters)] -> MCW_parameters_2
MCW_parameters_2[!(grepl("vector", MCW_parameters_2))] -> MCW_parameters_to_compare

MCH_parameters <- MCH_fit_summary$variable
MCH_parameters[grepl("b0|btemp1_|btemp0_|bspillwindow|bwinterspill", MCH_parameters)] -> MCH_parameters_2
MCH_parameters_2[!(grepl("vector", MCH_parameters_2))] -> MCH_parameters_to_compare

MC_parameters_to_compare <- intersect(MCH_parameters_to_compare, MCW_parameters_to_compare)

MC_parameters_comparison_plotlist <- vector(mode = "list", length = length(MC_parameters_to_compare))

for (i in 1:length(MC_parameters_to_compare)){
  MC_parameters_comparison_plotlist[[i]] <- plot_posterior_comparion(parameter = MC_parameters_to_compare[i], 
                                                                     hatchery_fit = MCH_fit, wild_fit = MCW_fit)
  print(paste0("i = ", i))
}

for (i in 1:length(MC_parameters_to_compare)){
  ggsave(file = paste0("stan_actual/output/posterior_comparison_plots/MC/", MC_parameters_to_compare[[i]], 
                       "_comp_plot.png"), plot = MC_parameters_comparison_plotlist[[i]],
         height = 8, width = 8)
  
  print(paste0("i = ", i))
}

#### Posterior comparisons: Deschutes movement investigation ####
MCW_parameters <- MCW_fit_summary$variable
MCW_parameters[grepl("2_10", MCW_parameters)] -> MCW_DES_movement_parameters
MCW_DES_movement_parameters[!(grepl("vector|NDE", MCW_DES_movement_parameters))] -> MCW_DES_movement_parameters

MCH_parameters <- MCH_fit_summary$variable
MCH_parameters[grepl("2_10", MCH_parameters)] -> MCH_DES_movement_parameters
MCH_DES_movement_parameters[!(grepl("vector|NDE", MCH_DES_movement_parameters))] -> MCH_DES_movement_parameters

hatchery_origins_for_comp <- c("origin1", "origin2")
wild_origins_for_comp <- c("origin4", "origin5")


MCW_DES_movement_parameters[(grepl(paste0(wild_origins_for_comp[1], "|", wild_origins_for_comp[2]), MCW_DES_movement_parameters))] -> MCW_DES_movements_comp
MCH_DES_movement_parameters[(grepl(paste0(hatchery_origins_for_comp[1], "|", hatchery_origins_for_comp[2]), MCH_DES_movement_parameters))] -> MCH_DES_movements_comp

# compare the parameters that are not shared between origins

MCW_DES_movements_solo <- MCW_DES_movement_parameters[!(MCW_DES_movement_parameters %in% c(MCW_DES_movements_comp))]

MCW_DES_movements_solo_plotlist <- vector(mode = "list", length = length(MCW_DES_movements_solo))

for (i in 1:length(MCW_DES_movements_solo)){
  MCW_DES_movements_solo_plotlist[[i]] <- plot_single_posterior(parameter = MCW_DES_movements_solo[i], 
                                                                         fit = MCW_fit)
  print(paste0("i = ", i))
}

for (i in 1:length(MCW_DES_movements_solo)){
  ggsave(file = paste0("stan_actual/output/posterior_comparison_plots/MC/MCW_", MCW_DES_movements_solo[[i]], 
                       "plot.png"), plot = MCW_DES_movements_solo_plotlist[[i]],
         height = 8, width = 8)
  
  print(paste0("i = ", i))
}

# now, compare origin parameters

MCW_DES_movements_comp_plotlist <- vector(mode = "list", length = length(MCH_DES_movements_comp))

for (i in 1:length(MCW_DES_movements_comp)){
  MCW_DES_movements_comp_plotlist[[i]] <- plot_posterior_comparion_origin(hatchery_parameter = MCH_DES_movements_comp[i], 
                                                                          wild_parameter = MCW_DES_movements_comp[i], 
                                                                          hatchery_fit = MCH_fit, 
                                                                          wild_fit = MCW_fit)
  
  ggsave(file = paste0("stan_actual/output/posterior_comparison_plots/MC/comp_", MCW_DES_movements_comp[[i]], "_",
                       MCH_DES_movements_comp[[i]],
                       "plot.png"), plot = MCW_DES_movements_comp_plotlist[[i]],
         height = 8, width = 8)

  print(paste0("i = ", i))
}


# Without plotting, look into individual parameters
MCW_parameters[grepl("_2_", MCW_parameters)] -> MCW_out2_params
MCW_out2_params[!(grepl("vector|NDE|actual", MCW_out2_params))] -> MCW_out2_params
MCH_parameters[grepl("_2_", MCH_parameters)] -> MCH_out2_params
MCH_out2_params[!(grepl("vector|NDE|actual", MCH_out2_params))] -> MCH_out2_params
MCH_out2_params[grepl("origin1", MCH_out2_params)] -> UMA_H_out2_params

subset(MCH_fit_summary, variable %in% UMA_H_out2_params)


##### Posterior comparisons: Snake River ####

SRW_parameters <- SRW_fit_summary$variable
SRW_parameters[grepl("b0|btemp1_|btemp0_|bspillwindow|bwinterspill", SRW_parameters)] -> SRW_parameters_2
SRW_parameters_2[!(grepl("vector", SRW_parameters_2))] -> SRW_parameters_to_compare

SRH_parameters <- SRH_fit_summary$variable
SRH_parameters[grepl("b0|btemp1_|btemp0_|bspillwindow|bwinterspill", SRH_parameters)] -> SRH_parameters_2
SRH_parameters_2[!(grepl("vector", SRH_parameters_2))] -> SRH_parameters_to_compare

SR_parameters_to_compare <- intersect(SRH_parameters_to_compare, SRW_parameters_to_compare)

SR_parameters_comparison_plotlist <- vector(mode = "list", length = length(SR_parameters_to_compare))

for (i in 1:length(SR_parameters_to_compare)){
  SR_parameters_comparison_plotlist[[i]] <- plot_posterior_comparion(parameter = SR_parameters_to_compare[i], 
                                                                     hatchery_fit = SRH_fit, wild_fit = SRW_fit)
  print(paste0("i = ", i))
}

for (i in 1:length(SR_parameters_to_compare)){
  ggsave(file = paste0("stan_actual/output/posterior_comparison_plots/SR/", SR_parameters_to_compare[[i]], 
                       "_comp_plot.png"), plot = SR_parameters_comparison_plotlist[[i]],
         height = 8, width = 8)
  
  print(paste0("i = ", i))
}

#### Calculate movement probabilities ####

# Here we are using draws of parameters to estimate overall movement probabilities

calculate_movement_probabilities <- function(fit, fit_summary, state){
  parameters <- unique(fit_summary$variable)
  
  one_state_params <- parameters[grep(paste0("b0_matrix_", state, "_"), parameters)]
  # remove NDE parameters
  one_state_params <- one_state_params[!(grepl("NDE", one_state_params))]
  # remove year parameters
  one_state_params <- one_state_params[!(grepl("year", one_state_params))]
  
  one_state_param_draws <- data.frame(iter = 1:4000)
  for (i in 1:length(one_state_params)){
    as.matrix(fit[,,one_state_params[i]]) %>% 
        as.data.frame() -> param_draws
    colnames(param_draws) <- one_state_params[i]
    
    one_state_param_draws %>% 
      bind_cols(., param_draws) -> one_state_param_draws
  }
    
    
  
  # as.data.frame(fit[,,one_state_params]) %>% 
  #   as.data.frame() %>% 
  #   dplyr::rename(hatchery = V1) -> hatchery_draws
  
  
}

#### Investigate transition counts from individual states ####

# Load detection efficiency information
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")

deschutes_river_trib_det_eff_capability <- data.frame(natal_origin = "Deschutes River",
                                                      run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                      DE = c(rep(0,9), rep(1,6), rep(0,3)))

john_day_river_trib_det_eff_capability <- data.frame(natal_origin = "John Day River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,8), rep(1,10)))

hood_river_trib_det_eff_capability <- data.frame(natal_origin = "Hood River",
                                                 run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                 DE = c(rep(0,9), rep(1,9)))

fifteenmile_creek_trib_det_eff_capability <- data.frame(natal_origin = "Fifteenmile Creek",
                                                        run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                        DE = c(rep(0,8), rep(1,10)))

umatilla_river_trib_det_eff_capability <- data.frame(natal_origin = "Umatilla River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,3), rep(1,15)))

yakima_river_trib_det_eff_capability <- data.frame(natal_origin = "Yakima River",     
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   # DE = c(rep(1,18)))
                                                   # For now - remove 04/05 run year
                                                   DE = c(0, rep(1,17)))

walla_walla_river_trib_det_eff_capability <- data.frame(natal_origin = "Walla Walla River",
                                                        run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                        DE = c(rep(0,1), rep(1,17)))

wenatchee_river_trib_det_eff_capability <- data.frame(natal_origin = "Wenatchee River",
                                                      run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                      DE = c(rep(0,7), rep(1,11)))

entiat_river_trib_det_eff_capability <- data.frame(natal_origin = "Entiat River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,4), rep(1,14)))

okanogan_river_trib_det_eff_capability <- data.frame(natal_origin = "Okanogan River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,9), rep(1,9)))

methow_river_trib_det_eff_capability <- data.frame(natal_origin = "Methow River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,5), rep(1,13)))

tucannon_river_trib_det_eff_capability <- data.frame(natal_origin = "Tucannon River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,7), rep(1,11)))

asotin_creek_trib_det_eff_capability <- data.frame(natal_origin = "Asotin Creek",      
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,7), rep(1,11)))

imnaha_river_trib_det_eff_capability <- data.frame(natal_origin = "Imnaha River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,7), rep(1,11)))


deschutes_river_trib_det_eff_capability %>% 
  bind_rows(., john_day_river_trib_det_eff_capability) %>% 
  bind_rows(., hood_river_trib_det_eff_capability) %>% 
  bind_rows(., fifteenmile_creek_trib_det_eff_capability) %>% 
  bind_rows(., umatilla_river_trib_det_eff_capability) %>% 
  bind_rows(., yakima_river_trib_det_eff_capability) %>% 
  bind_rows(., walla_walla_river_trib_det_eff_capability) %>% 
  bind_rows(., wenatchee_river_trib_det_eff_capability) %>% 
  bind_rows(., entiat_river_trib_det_eff_capability) %>% 
  bind_rows(., okanogan_river_trib_det_eff_capability) %>% 
  bind_rows(., methow_river_trib_det_eff_capability) %>% 
  bind_rows(., tucannon_river_trib_det_eff_capability) %>% 
  bind_rows(., asotin_creek_trib_det_eff_capability) %>% 
  bind_rows(., imnaha_river_trib_det_eff_capability) -> trib_det_eff_capability



# Load the model data associated with each run (necessary to load covariates)
# Store these each in an environment, because most things share names
UCW_envir <- new.env()
UCH_envir <- new.env()
MCW_envir <- new.env()
MCH_envir <- new.env()
SRW_envir <- new.env()
SRH_envir <- new.env()
# load(here::here("stan_actual", "reparameterization_v2", "upper_columbia_wild", "model_data.rda"),
#      envir = UCW_envir)
# load(here::here("stan_actual", "reparameterization_v2", "upper_columbia_hatchery", "model_data.rda"),
#      envir = UCH_envir)
load(here::here("stan_actual", "reparameterization_v2", "middle_columbia_wild", "model_data.rda"),
     envir = MCW_envir)
load(here::here("stan_actual", "reparameterization_v2", "middle_columbia_hatchery", "model_data.rda"),
     envir = MCH_envir)
# load(here::here("stan_actual", "reparameterization_v2", "snake_river_wild", "model_data.rda"),
#      envir = SRW_envir)
# load(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery", "model_data.rda"),
#      envir = SRH_envir)


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

# function to extract states and dates by origin
get_origin_states_transitions <- function(envir, origin_select, rear){
  
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  # Now, keep only the origin selected
  if(rear == "wild"){
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
  } else {
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
  }
  
  # keep only the transitions for fish from this origin
  transitions <- envir$data$y
  transitions[which(origin_vector == origin_numeric),] -> origin_transitions
  
  state_history <- as.vector(t(origin_transitions))
  
  transitions_df <- data.frame(from = rep(NA, length(state_history)-1),
                               to = rep(NA, length(state_history)-1))
  
  for (i in 1:nrow(transitions_df)){
    transitions_df$from[i] <- state_history[i]
    transitions_df$to[i] <- state_history[i+1]
  }
  # drop the non-data
  transitions_df %>% 
    filter(!(from %in% c(0,43))) -> transitions_df
  
  # count transitions out of each state
  transitions_df %>% 
    group_by(from, to) %>% 
    summarise(count = n()) -> transition_counts
  
  transition_counts %>% 
    ungroup() %>% 
    group_by(from) %>% 
    mutate(total = sum(count)) %>% 
    mutate(prop = count/total) -> transition_counts
  
  # return(transition_counts)

  return(list(transitions_df = transitions_df, transition_counts = transition_counts))
}

# Let's compare movements out of state 2 across different populations
# Deschutes
DES_wild_transitions_ind <- get_origin_states_transitions(envir = MCW_envir, origin_select = "Deschutes River", rear = "wild")[[1]]
DES_wild_transitions <- get_origin_states_transitions(envir = MCW_envir, origin_select = "Deschutes River", rear = "wild")[[2]]
DES_out2 <- subset(DES_wild_transitions, from == 2)
DES_out2$origin <- "Deschutes River"

# John Day
JDR_wild_transitions <- get_origin_states_transitions(envir = MCW_envir, origin_select = "John Day River", rear = "wild")[[2]]
JDR_out2 <- subset(JDR_wild_transitions, from == 2)
JDR_out2$origin <- "John Day River"

# Umatilla
UMA_wild_transitions <- get_origin_states_transitions(envir = MCW_envir, origin_select = "Umatilla River", rear = "wild")[[2]]
UMA_out2 <- subset(UMA_wild_transitions, from == 2)
UMA_out2$origin <- "Umatilla River"

# Walla Walla
WAWA_wild_transitions <- get_origin_states_transitions(envir = MCW_envir, origin_select = "Walla Walla River", rear = "wild")[[2]]
WAWA_out2 <- subset(WAWA_wild_transitions, from == 2)
WAWA_out2$origin <- "Walla Walla River"

# Yakima
YAK_wild_transitions <- get_origin_states_transitions(envir = MCW_envir, origin_select = "Yakima River", rear = "wild")[[2]]
YAK_out2 <- subset(YAK_wild_transitions, from == 2)
YAK_out2$origin <- "Yakima River"

# Fifteenmile
FIF_wild_transitions <- get_origin_states_transitions(envir = MCW_envir, origin_select = "Fifteenmile Creek", rear = "wild")[[2]]
FIF_out2 <- subset(FIF_wild_transitions, from == 2)
FIF_out2$origin <- "Fifteenmile Creek"

DES_out2 %>% 
  bind_rows(., JDR_out2) %>% 
  bind_rows(., UMA_out2) %>% 
  bind_rows(., WAWA_out2) %>% 
  bind_rows(., YAK_out2) %>% 
  bind_rows(., FIF_out2) -> MCW_out2
  



### Inspect some populations that don't seem quite right

# Fifteenmile
FIF_wild_transitions <- get_origin_states_transitions(envir = MCW_envir, origin_select = "Fifteenmile Creek", rear = "wild")[[2]]

# movements out of the Deschutes?
subset(FIF_wild_transitions, from == 10)

# movements out of the mainstem?
subset(FIF_wild_transitions, from == 2)

ggplot(FIF_wild_transitions, aes(x = from, y = prop, fill = as.character(to))) +
  geom_bar(stat = "identity") +
  scale_fill_tableau(palette = "Tableau 20")

# Umatilla
UMA_wild_transitions <- get_origin_states_transitions(envir = MCW_envir, origin_select = "Umatilla River", rear = "wild")[[2]]

# movements out of the Deschutes?
subset(UMA_wild_transitions, from == 10)
subset(UMA_wild_transitions, from == 11)

# movements out of the mainstem?
subset(UMA_wild_transitions, from == 2)

ggplot(UMA_wild_transitions, aes(x = from, y = prop, fill = as.character(to))) +
  geom_bar(stat = "identity") +
  scale_fill_tableau(palette = "Tableau 20")


UMA_hatchery_transitions <- get_origin_states_transitions(envir = MCH_envir, origin_select = "Umatilla River", rear = "hatchery")[[2]]

# movements out of the Deschutes?
subset(UMA_hatchery_transitions, from == 10)
subset(UMA_hatchery_transitions, from == 11)

# Every single Umatilla hatchery fish that went into the Deschutes never came out?

# movements out of the mainstem?
subset(UMA_hatchery_transitions, from == 2)

ggplot(UMA_hatchery_transitions, aes(x = from, y = prop, fill = as.character(to))) +
  geom_bar(stat = "identity") +
  scale_fill_tableau(palette = "Tableau 20")
  


# Plot state transitions by temperature at which they occurred
get_origin_states_transitions_temperature <- function(envir, origin_select, rear){
  
  temp_data <- envir$data$temperature
  
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  # Now, keep only the origin selected
  if(rear == "wild"){
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
  } else {
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
  }
  
  # keep only the transitions for fish from this origin
  transitions <- envir$data$y
  transitions[which(origin_vector == origin_numeric),] -> origin_transitions
  transition_dates <- envir$data$transition_dates
  transition_dates[which(origin_vector == origin_numeric),] -> origin_transition_dates
  
  state_history <- as.vector(t(origin_transitions))
  transition_date_history <- as.vector(t(origin_transition_dates))
  
  transitions_df <- data.frame(from = rep(NA, length(state_history)-1),
                               to = rep(NA, length(state_history)-1),
                               entry_date = rep(NA, length(state_history)-1),
                               from_temp = rep(NA, length(state_history)-1))
  
  
  for (i in 1:nrow(transitions_df)){
    transitions_df$from[i] <- state_history[i]
    transitions_df$to[i] <- state_history[i+1]
    transitions_df$entry_date[i] <- transition_date_history[i]
    # add the temperature data, if transition is from a mainstem state
    if (transitions_df$from[i] %in% c(1:9)){
      transitions_df$from_temp[i] <- temp_data[transitions_df$entry_date[i], transitions_df$from[i]]
    }
    
  }
  # drop the non-data
  transitions_df %>% 
    filter(!(from %in% c(0,43))) -> transitions_df
  
  return(transitions_df)
}


FIF_transitions_temps <- get_origin_states_transitions_temperature(envir = MCW_envir, origin_select = "Fifteenmile Creek", rear = "wild")

ggplot(subset(FIF_transitions_temps, from == 2), aes(x = from_temp, color = as.character(to))) +
  geom_boxplot() +
  scale_color_tableau(palette = "Tableau 10")

#### Investigate movement probabilities out of one state ####
# this function is based largely on the estimation of movement probabilities from covariate values
estimate_move_prob_MCW <- function(origin_select, from_state){
  
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
  
  # get all of the possible movements from that from state
  
  possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from_state]
  possible_movements <- c(possible_movements, 43)
  
  # make sure to drop all movements to upstream states (in DE years, these aren't possible)
  grep(" Upstream", model_states) -> upstream_indices
  possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements to different states
  # store as matrix, with columns = possible movements, rows = iterations/draws
  move_prob_array <- matrix(nrow = niter, ncol = length(possible_movements))
  dimnames(move_prob_array)[[2]] <- possible_movements
  

    
    # If it's a mainstem state, get the spill data. If it's not, ignore it.
    
    if (from_state %in% c(1:9)){
      # get median winter spill for this state
      MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                                     date = as.vector(MCW_envir$data$transition_dates))
      spillwindow_data <- MCW_envir$data$spill_window_data
      med_spillwindow <- median(spillwindow_data[subset(MCW_states_dates, state == from_state)$date,from_state])
      
      # get median spill window for this state
      winterspill_data <- MCW_envir$data$winter_spill_days_data
      med_winterspill <- median(winterspill_data[,from_state])     
      
      # get median temperature for this state
      temp_data <- MCW_envir$data$temperature
      med_temp <- median(temp_data[subset(MCW_states_dates, state == from_state)$date,from_state])
      
    } else {
      med_spillwindow <- 0
      med_winterspill <- 0
    }
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(possible_movements)){
          from <- from_state
          to <- possible_movements[j]
        
        # evaluation movement 
        move_prob_array[iter,j] <- exp(b0_array_MCW[from,to,iter] +
                                                 # btemp0_array_MCW[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_MCW[from,to,iter]*med_temp + 
                                                 bspillwindow_array_MCW[from,to,iter]*med_spillwindow +
                                                 bwinterspill_array_MCW[from,to,iter]*med_winterspill +
                                                 # btemp0xorigin1_array_MCW[from,to,iter]*origin1 +
                                                 btemp1xorigin1_array_MCW[from,to,iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_MCW[from,to,iter]*origin2 + 
                                                 btemp1xorigin2_array_MCW[from,to,iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_MCW[from,to,iter]*origin3 +
                                                 btemp1xorigin3_array_MCW[from,to,iter]*med_temp*origin3 +
                                                 # btemp0xorigin4_array_MCW[from,to,iter]*origin4 +
                                                 btemp1xorigin4_array_MCW[from,to,iter]*med_temp*origin4 +
                                                 # btemp0xorigin5_array_MCW[from,to,iter]*origin5 +
                                                 btemp1xorigin5_array_MCW[from,to,iter]*med_temp*origin5 +
                                                 # btemp0xorigin6_array_MCW[from,to,iter]*origin6 +
                                                 btemp1xorigin6_array_MCW[from,to,iter]*med_temp*origin6 +
                                                 borigin1_array_MCW[from,to,iter]*origin1 +
                                                 borigin2_array_MCW[from,to,iter]*origin2 +
                                                 borigin3_array_MCW[from,to,iter]*origin3 +
                                                 borigin4_array_MCW[from,to,iter]*origin4 +
                                                 borigin5_array_MCW[from,to,iter]*origin5 +
                                                 borigin6_array_MCW[from,to,iter]*origin6)/
          sum(exp(b0_array_MCW[from,possible_movements,iter] +
                    # btemp0_array_MCW[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCW[from,possible_movements,iter]*med_temp + 
                    bspillwindow_array_MCW[from,possible_movements,iter]*med_spillwindow +
                    bwinterspill_array_MCW[from,possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCW[from,possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCW[from,possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_MCW[from,possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCW[from,possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_MCW[from,possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_MCW[from,possible_movements,iter]*med_temp*origin3 +
                    # btemp0xorigin4_array_MCW[from,possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_MCW[from,possible_movements,iter]*med_temp*origin4 +
                    # btemp0xorigin5_array_MCW[from,possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_MCW[from,possible_movements,iter]*med_temp*origin5 +
                    # btemp0xorigin6_array_MCW[from,possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_MCW[from,possible_movements,iter]*med_temp*origin6 +
                    borigin1_array_MCW[from,possible_movements,iter]*origin1 +
                    borigin2_array_MCW[from,possible_movements,iter]*origin2 +
                    borigin3_array_MCW[from,possible_movements,iter]*origin3 +
                    borigin4_array_MCW[from,possible_movements,iter]*origin4 +
                    borigin5_array_MCW[from,possible_movements,iter]*origin5 +
                    borigin6_array_MCW[from,possible_movements,iter]*origin6))
        
      }
      
      
      
    }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(move_prob_array)
  
}

estimate_move_prob_MCH <- function(origin_select, from_state){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  # get all of the possible movements from that from state
  
  possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from_state]
  possible_movements <- c(possible_movements, 43)
  
  # make sure to drop all movements to upstream states (in DE years, these aren't possible)
  grep(" Upstream", model_states) -> upstream_indices
  possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements to different states
  # store as matrix, with columns = possible movements, rows = iterations/draws
  move_prob_array <- matrix(nrow = niter, ncol = length(possible_movements))
  dimnames(move_prob_array)[[2]] <- possible_movements
  
  
  
  # If it's a mainstem state, get the spill data. If it's not, ignore it.
  
  if (from_state %in% c(1:9)){
    # get median winter spill for this state
    MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                                   date = as.vector(MCH_envir$data$transition_dates))
    spillwindow_data <- MCH_envir$data$spill_window_data
    med_spillwindow <- median(spillwindow_data[subset(MCH_states_dates, state == from_state)$date,from_state])
    
    # get median spill window for this state
    winterspill_data <- MCH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from_state])     
    
    # get median temperature for this state
    temp_data <- MCH_envir$data$temperature
    med_temp <- median(temp_data[subset(MCH_states_dates, state == from_state)$date,from_state])
    
  } else {
    med_spillwindow <- 0
    med_winterspill <- 0
  }
  
  # Loop through all of the iterations
  for (iter in 1:niter){
    
    # Loop through a sequence of temperature values to get predicted response
    # temperature was z-scored so we can plot 2 standard deviations
    for (j in 1:length(possible_movements)){
      from <- from_state
      to <- possible_movements[j]
      
      # evaluation movement 
      move_prob_array[iter,j] <- exp(b0_array_MCH[from,to,iter] +
                                       # btemp0_array_MCH[from,to,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                       btemp1_array_MCH[from,to,iter]*med_temp + 
                                       bspillwindow_array_MCH[from,to,iter]*med_spillwindow +
                                       bwinterspill_array_MCH[from,to,iter]*med_winterspill +
                                       # btemp0xorigin1_array_MCH[from,to,iter]*origin1 +
                                       btemp1xorigin1_array_MCH[from,to,iter]*med_temp*origin1 +
                                       # btemp0xorigin2_array_MCH[from,to,iter]*origin2 + 
                                       btemp1xorigin2_array_MCH[from,to,iter]*med_temp*origin2 + 
                                       borigin1_array_MCH[from,to,iter]*origin1 +
                                       borigin2_array_MCH[from,to,iter]*origin2)/
        sum(exp(b0_array_MCH[from,possible_movements,iter] +
                  # btemp0_array_MCH[from,possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                  btemp1_array_MCH[from,possible_movements,iter]*med_temp + 
                  bspillwindow_array_MCH[from,possible_movements,iter]*med_spillwindow +
                  bwinterspill_array_MCH[from,possible_movements,iter]*med_winterspill +
                  # btemp0xorigin1_array_MCH[from,possible_movements,iter]*origin1 +
                  btemp1xorigin1_array_MCH[from,possible_movements,iter]*med_temp*origin1 +
                  # btemp0xorigin2_array_MCH[from,possible_movements,iter]*origin2 + 
                  btemp1xorigin2_array_MCH[from,possible_movements,iter]*med_temp*origin2 + 
                  borigin1_array_MCH[from,possible_movements,iter]*origin1 +
                  borigin2_array_MCH[from,possible_movements,iter]*origin2))
      
    }
    
    
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(move_prob_array)
  
}

# Plot individual movement probabilities
plot_move_probs <- function(move_prob, plot_title){
  
  move_prob %>% 
    as.data.frame() %>% 
    pivot_longer(cols = everything(), names_to = "to", values_to = "prob") %>% 
    group_by(to) %>% 
    summarise(prob = quantile(prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) -> move_prob_quantiles
  
  move_prob %>% 
    as.data.frame() %>% 
    pivot_longer(cols = everything(), names_to = "to", values_to = "prob") -> move_prob_long
  
  # move_prob_plot <- ggplot(subset(move_prob_long, to =="1"), aes(x = prob, fill = to, color = to, group = to)) +
  move_prob_plot <- ggplot(move_prob_long, aes(x = prob, y = after_stat(scaled), fill = to, color = to, group = to)) +
    geom_density(alpha = 0.1) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(0,NA), expand = c(0,0)) +
    scale_fill_tableau(palette = "Tableau 10") +
    scale_color_tableau(palette = "Tableau 10") +
    xlab("Movement probability") +
    ylab("Density") +
    ggtitle(plot_title)
  
  return(move_prob_plot)
}

FIF_from2_move_prob <- estimate_move_prob_MCW(origin_select = "Fifteenmile Creek", from_state = 2)

FIF_from2_move_prob_plot <- plot_move_probs(move_prob = FIF_from2_move_prob, plot_title = "Fifteenmile Creek: Movements out of BON to MCN")

UMA_W_from2_move_prob <- estimate_move_prob_MCW(origin_select = "Umatilla River", from_state = 2)

UMA_W_from2_move_prob_plot <- plot_move_probs(move_prob = UMA_W_from2_move_prob, plot_title = "Umatilla River, Wild: Movements out of BON to MCN")

UMA_H_from2_move_prob <- estimate_move_prob_MCH(origin_select = "Umatilla River", from_state = 2)

UMA_H_from2_move_prob_plot <- plot_move_probs(move_prob = UMA_H_from2_move_prob, plot_title = "Umatilla River, Hatchery: Movements out of BON to MCN")


#### Troubleshooting ####

# Why is the model estimated effect not matching the data?
# We have far too many fish moving into the Deschutes, across origins from the Middle Columbia

# We need to create a df with the following columns: fish, state (2), next state, temp
# Then, we can try fitting a logistic regression to the movement, based on temp

# We can use the get_origin_states_transitions_temperature() function to get this

UMA_H_transitions_temps <- get_origin_states_transitions_temperature(envir = MCH_envir, origin_select = "Umatilla River", rear = "hatchery")

subset(UMA_H_transitions_temps, from == 2) -> UMA_H_out2

# Do a very rough boxplot

model_states_for_join <- data.frame(to = 1:length(model_states),
                                    state = model_states)

UMA_H_out2 %>% 
  left_join(., model_states_for_join, by = "to") -> UMA_H_out2

ggplot(UMA_H_out2, aes(x = from_temp, color = state)) +
  geom_boxplot() +
  scale_color_tableau(palette = "Tableau 10")

# now, change response so that 1 = Deschutes, 0 = not Deschutes

UMA_H_out2 %>% 
  mutate(DES = ifelse(state == "Deschutes River Mouth", 1, 0)) %>% 
  mutate(temp_actual = from_temp*window_temp_summary$BON_sd + window_temp_summary$BON_mean) -> UMA_H_out2

ggplot(UMA_H_out2, aes(x = from_temp, fill = DES)) +
  # geom_bar(position = "stack", stat = "identity") +
  geom_bar(stat = "count") +
  scale_x_binned()

UMA_H_temp_move_plot <- ggplot(UMA_H_out2, aes(x = temp_actual, fill = state)) + 
  # geom_bar(position = "stack", stat = "identity") +
  geom_bar(stat = "count") +
  scale_x_binned(n.breaks = 20) +
  ggtitle("Umatilla Hatchery: Movements out of BON to MCN by temperature")

ggsave(here::here("stan_actual", "output", "data_plots", "UMA_H_temp_move_plot.png"), UMA_H_temp_move_plot, height = 8, width = 8)






# try umatilla again - but now, only keep the years where we have DE ability in the Deschutes
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_numeric = seq(4, 22, 1)
run_year_index = 1:19
run_year_df <- data.frame(run_year,run_year_numeric, run_year_index)

UMA_H_transitions_temps %>% 
  mutate(run_year_index = ceiling(entry_date/365.25)+1) %>% 
  left_join(., run_year_df, by = "run_year_index") %>% 
  filter(run_year %in% c("13/14", "14/15", "15/16", "16/17", "17/18", "18/19")) -> UMA_H_transitions_temps_Deschutes_DE_years

subset(UMA_H_transitions_temps_Deschutes_DE_years, from == 2) -> UMA_H_out2_Deschutes_DE

# Do a very rough boxplot

model_states_for_join <- data.frame(to = 1:length(model_states),
                                    state = model_states)

UMA_H_out2_Deschutes_DE %>% 
  left_join(., model_states_for_join, by = "to") -> UMA_H_out2_Deschutes_DE

ggplot(UMA_H_out2_Deschutes_DE, aes(x = from_temp, color = state)) +
  geom_boxplot() +
  scale_color_tableau(palette = "Tableau 10")

# now, change response so that 1 = Deschutes, 0 = not Deschutes

UMA_H_out2_Deschutes_DE %>% 
  mutate(DES = ifelse(state == "Deschutes River Mouth", 1, 0)) %>% 
  mutate(temp_actual = from_temp*window_temp_summary$BON_sd + window_temp_summary$BON_mean) -> UMA_H_out2_Deschutes_DE

UMA_H_temp_move_plot_Deschutes_DE <- ggplot(UMA_H_out2_Deschutes_DE, aes(x = temp_actual, fill = state)) + 
  # geom_bar(position = "stack", stat = "identity") +
  geom_bar(stat = "count") +
  scale_x_binned(breaks = 0.5:22.5, lim = c(0,25)) +
  ggtitle("Umatilla Hatchery: Movements out of BON to MCN by temperature") +
  theme(axis.text.x = element_text(angle = 45))

ggsave(here::here("stan_actual", "output", "data_plots", "UMA_H_temp_move_plot_Deschutes_DE.png"), UMA_H_temp_move_plot_Deschutes_DE, height = 8, width = 8)





FIF_W_transitions_temps <- get_origin_states_transitions_temperature(envir = MCW_envir, origin_select = "Fifteenmile Creek", rear = "wild")

subset(FIF_W_transitions_temps, from == 2) -> FIF_W_out2

# Do a very rough boxplot

model_states_for_join <- data.frame(to = 1:length(model_states),
                                    state = model_states)

FIF_W_out2 %>% 
  left_join(., model_states_for_join, by = "to") -> FIF_W_out2

ggplot(FIF_W_out2, aes(x = from_temp, color = state)) +
  geom_boxplot() +
  scale_color_tableau(palette = "Tableau 10")

# now, change response so that 1 = Deschutes, 0 = not Deschutes

FIF_W_out2 %>% 
  mutate(DES = ifelse(state == "Deschutes River Mouth", 1, 0)) %>% 
  mutate(temp_actual = from_temp*window_temp_summary$BON_sd + window_temp_summary$BON_mean) -> FIF_W_out2

ggplot(FIF_W_out2, aes(x = from_temp, fill = DES)) +
  # geom_bar(position = "stack", stat = "identity") +
  geom_bar(stat = "count") +
  scale_x_binned()

FIF_W_temp_move_plot <- ggplot(FIF_W_out2, aes(x = temp_actual, fill = state)) + 
  # geom_bar(position = "stack", stat = "identity") +
  geom_bar(stat = "count") +
  scale_x_binned(n.breaks = 20) +
  ggtitle("Fifteenmile Wild: Movements out of BON to MCN by temperature")

ggsave(here::here("stan_actual", "output", "data_plots", "FIF_W_temp_move_plot.png"), FIF_W_temp_move_plot, height = 8, width = 8)

# quick logistic regression
glm(DES~from_temp, data = UMA_H_out2, family = "binomial") -> UMA_H_DES_glm

predict.glm(UMA_H_DES_glm,
            type = "response",
            newdata = data.frame(from_temp = seq(-2,2,0.01))) -> UMA_H_DES_glm_predict

# Plot actual temperatures
window_temp_summary <- read.csv(here::here("stan_actual", "window_temps_summary.csv"), row.names = 1)
actual_temps <- seq(-2,2,0.01)*window_temp_summary$BON_sd + window_temp_summary$BON_mean


plot(x = actual_temps, y = UMA_H_DES_glm_predict)

# Another thing to look into: What about temp0 vs. temp1?
















#### Plot parameter estimates ####


hist(btemp1xorigin2_array_MCW[2,16,])
summary(btemp1xorigin2_array_MCW[2,16,])
hist(btemp1xorigin2_array_MCW[2,10,])
summary(btemp1xorigin2_array_MCW[2,10,])
hist(borigin2_array_MCW[2,16,])
summary(borigin2_array_MCW[2,16,])
hist(borigin2_array_MCW[2,10,])
summary(borigin2_array_MCW[2,10,])


#### Compare parameter estimates for movements out of state 2 - MC populations ####

MCW_parameters <- unique(MCW_fit_summary$variable)
MCH_parameters <- unique(MCW_fit_summary$variable)

MCW_parameters_state2 <- MCW_parameters[grep("_2_", MCW_parameters)]
MCH_parameters_state2 <- MCH_parameters[grep("_2_", MCH_parameters)]
# remove NDE parameters
MCW_parameters_state2 <- MCW_parameters_state2[!(grepl("NDE", MCW_parameters_state2))]
MCH_parameters_state2 <- MCH_parameters_state2[!(grepl("NDE", MCH_parameters_state2))]
# remove year parameters 
MCW_parameters_state2 <- MCW_parameters_state2[!(grepl("year", MCW_parameters_state2))]
MCH_parameters_state2 <- MCH_parameters_state2[!(grepl("year", MCH_parameters_state2))]

# Create DFs of these parameter estimates
subset(MCW_fit_summary, variable %in% MCW_parameters_state2) -> MCW_fit_summary_state2
subset(MCH_fit_summary, variable %in% MCH_parameters_state2) -> MCH_fit_summary_state2

model_states_df <- data.frame(to_state_numeric = 1:length(model_states),
                              state = model_states)

# add origin information and to state information
MCW_fit_summary_state2 %>%
  mutate(to_state_numeric = as.numeric(str_extract(sub("[^_]*_[^_]*_[^_]*_", "", variable), "\\d+"))) %>% 
  left_join(model_states_df, by = "to_state_numeric") %>% 
  mutate(state = gsub(" Mouth", "", state)) %>% 
  bind_cols(str_locate(MCW_fit_summary_state2$variable, "origin")) %>% 
  mutate(origin_numeric = as.numeric(str_sub(variable, end+1, end+1))) %>% 
  left_join(origin_param_map[1:7,], join_by(origin_numeric == wild)) %>% 
  mutate(natal_origin = ifelse(is.na(origin_numeric), NA, natal_origin)) %>% 
  dplyr::select(-c(mean, sd, mad, q5, q95, rhat, ess_bulk, ess_tail, start, end, origin_numeric, hatchery)) %>% 
  mutate(rear = "wild") -> MCW_fit_summary_state2

MCH_fit_summary_state2 %>%
  mutate(to_state_numeric = as.numeric(str_extract(sub("[^_]*_[^_]*_[^_]*_", "", variable), "\\d+"))) %>% 
  left_join(model_states_df, by = "to_state_numeric") %>% 
  mutate(state = gsub(" Mouth", "", state)) %>% 
  bind_cols(str_locate(MCH_fit_summary_state2$variable, "origin")) %>% 
  mutate(origin_numeric = as.numeric(str_sub(variable, end+1, end+1))) %>% 
  left_join(origin_param_map[c(5,7),], join_by(origin_numeric == hatchery)) %>% 
  mutate(natal_origin = ifelse(is.na(origin_numeric), NA, natal_origin)) %>% 
  dplyr::select(-c(mean, sd, mad, q5, q95, rhat, ess_bulk, ess_tail, start, end, origin_numeric, wild)) %>% 
  mutate(rear = "hatchery") -> MCH_fit_summary_state2

MCW_fit_summary_state2 %>% 
  bind_rows(., MCH_fit_summary_state2) %>% 
  mutate(natal_origin = ifelse(is.na(natal_origin), "DPS", natal_origin)) %>% 
  mutate(parameter_type = ifelse(str_detect(variable, "temp0"), "temp0", 
                                 ifelse(str_detect(variable, "temp1"), "temp1", 
                                        ifelse(str_detect(variable, "spillwindow"), "spillwindow", 
                                               ifelse(str_detect(variable, "borigin"), "borigin", "big uh oh!"))))) %>% 
  mutate(population = paste0(natal_origin, " - ", rear)) %>% 
  # drop spillwindow, it's not origin specific so it won't help us diagnose things here
  filter(parameter_type != "spillwindow") -> MC_fit_summary_state_2

MC_state_2_param_plot <- ggplot(MC_fit_summary_state_2, aes(x = state, y = median, color = population)) +
  geom_point(position =position_dodge(width=0.25)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_color_tableau(palette = "Tableau 10") +
  facet_wrap(~parameter_type, nrow = 1) +
  ylab("Parameter Estimate")

ggsave(here::here("stan_actual", "output", "MC_state_2_param_plot.png"), MC_state_2_param_plot, height = 8, width = 8)









