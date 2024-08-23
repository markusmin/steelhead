# 11-model-spill-plots-march

# This script takes the output from the stan model runs in 05-stan-runs and
# plots the effects of spill (window and winter days)
# this looks at the effect of winter spill if only march is included

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
# UCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_wild","UCW_transition_counts.csv"))
# UCW_movements <- paste0("_", UCW_transition_counts$from, "_", UCW_transition_counts$to)
# UCW_movements <- UCW_movements[!(grepl("NA", UCW_movements))]
# 
# UCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_hatchery","UCH_transition_counts.csv"))
# UCH_movements <- paste0("_", UCH_transition_counts$from, "_", UCH_transition_counts$to)
# UCH_movements <- UCH_movements[!(grepl("NA", UCH_movements))]

MCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_wild","MCW_transition_counts.csv"))
MCW_movements <- paste0("_", MCW_transition_counts$from, "_", MCW_transition_counts$to)
MCW_movements <- MCW_movements[!(grepl("NA", MCW_movements))]

MCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_hatchery","MCH_transition_counts.csv"))
MCH_movements <- paste0("_", MCH_transition_counts$from, "_", MCH_transition_counts$to)
MCH_movements <- MCH_movements[!(grepl("NA", MCH_movements))]

# SRW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_wild","SRW_transition_counts.csv"))
# SRW_movements <- paste0("_", SRW_transition_counts$from, "_", SRW_transition_counts$to)
# SRW_movements <- SRW_movements[!(grepl("NA", SRW_movements))]
# 
# SRH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_hatchery","SRH_transition_counts.csv"))
# SRH_movements <- paste0("_", SRH_transition_counts$from, "_", SRH_transition_counts$to)
# SRH_movements <- SRH_movements[!(grepl("NA", SRH_movements))]


##### Load the model runs #####

# Load the model data associated with each run (necessary to load covariates)
# Store these each in an environment, because most things share names
# UCW_envir <- new.env()
# UCH_envir <- new.env()
MCW_envir <- new.env()
MCH_envir <- new.env()
# SRW_envir <- new.env()
# SRH_envir <- new.env()
# load(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_wild", "model_data.rda"),
#      envir = UCW_envir)
# load(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_hatchery", "model_data.rda"),
#      envir = UCH_envir)
load(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_wild", "model_data.rda"),
     envir = MCW_envir)
load(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_hatchery", "model_data.rda"),
     envir = MCH_envir)
# load(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_wild", "model_data.rda"),
#      envir = SRW_envir)
# load(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_hatchery", "model_data.rda"),
#      envir = SRH_envir)




# Function to bind four chains together
bind4chains <- function(chain1, chain2, chain3, chain4){
  bound_draws <- bind_draws(chain1$draws(),
                            chain2$draws(),
                            chain3$draws(),
                            chain4$draws(), along = "chain")
  
  return(bound_draws)
}

# ## Upper Columbia, Wild
# UCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_wild", "chain1_UCW_reparam_v2_fit.rds"))
# UCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_wild", "chain2_UCW_reparam_v2_fit.rds"))
# UCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_wild", "chain3_UCW_reparam_v2_fit.rds"))
# UCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_wild", "chain4_UCW_reparam_v2_fit.rds"))
# 
# # bind chains together
# UCW_fit_raw <- bind4chains(UCW_chain1, UCW_chain2, UCW_chain3, UCW_chain4)
# # thin2
# thin_draws(UCW_fit_raw, thin = 2) -> UCW_fit
# # summarise
# UCW_fit_summary <- summarise_draws(UCW_fit)
# 
# ## Upper Columbia, Hatchery
# UCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_hatchery", "chain1_UCH_reparam_v2_fit.rds"))
# UCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_hatchery", "chain2_UCH_reparam_v2_fit.rds"))
# UCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_hatchery", "chain3_UCH_reparam_v2_fit.rds"))
# UCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "upper_columbia_hatchery", "chain4_UCH_reparam_v2_fit.rds"))
# 
# # bind chains together
# UCH_fit_raw <- bind4chains(UCH_chain1, UCH_chain2, UCH_chain3, UCH_chain4)
# # thin2
# thin_draws(UCH_fit_raw, thin = 2) -> UCH_fit
# # summarise
# UCH_fit_summary <- summarise_draws(UCH_fit)

## Middle Columbia, Wild
MCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_wild", "chain1_MCW_reparam_v2_fit.rds"))
MCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_wild", "chain2_MCW_reparam_v2_fit.rds"))
MCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_wild", "chain3_MCW_reparam_v2_fit.rds"))
MCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_wild", "chain4_MCW_reparam_v2_fit.rds"))

# bind chains together
MCW_fit_raw <- bind4chains(MCW_chain1, MCW_chain2, MCW_chain3, MCW_chain4)
# thin2
thin_draws(MCW_fit_raw, thin = 2) -> MCW_fit
# summarise
MCW_fit_summary <- summarise_draws(MCW_fit)

## Middle Columbia, Hatchery
MCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_hatchery", "chain1_MCH_reparam_v2_fit.rds"))
MCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_hatchery", "chain2_MCH_reparam_v2_fit.rds"))
MCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_hatchery", "chain3_MCH_reparam_v2_fit.rds"))
MCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "middle_columbia_hatchery", "chain4_MCH_reparam_v2_fit.rds"))

# bind chains together
MCH_fit_raw <- bind4chains(MCH_chain1, MCH_chain2, MCH_chain3, MCH_chain4)
# thin2
thin_draws(MCH_fit_raw, thin = 2) -> MCH_fit
# summarise
MCH_fit_summary <- summarise_draws(MCH_fit)

# ## Snake River, Wild
# SRW_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_wild", "chain1_SRW_reparam_v2_fit.rds"))
# SRW_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_wild", "chain2_SRW_reparam_v2_fit.rds"))
# SRW_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_wild", "chain3_SRW_reparam_v2_fit.rds"))
# SRW_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_wild", "chain4_SRW_reparam_v2_fit.rds"))
# 
# # bind chains together
# SRW_fit_raw <- bind4chains(SRW_chain1, SRW_chain2, SRW_chain3, SRW_chain4)
# # thin2
# thin_draws(SRW_fit_raw, thin = 2) -> SRW_fit
# # summarise
# SRW_fit_summary <- summarise_draws(SRW_fit)
# 
# ## Snake River, Hatchery
# SRH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_hatchery", "chain1_SRH_reparam_v2_fit.rds"))
# SRH_chain2 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_hatchery", "chain2_SRH_reparam_v2_fit.rds"))
# SRH_chain3 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_hatchery", "chain3_SRH_reparam_v2_fit.rds"))
# SRH_chain4 <- readRDS(here::here("stan_actual", "reparameterization_v2_marchspilldays", "snake_river_hatchery", "chain4_SRH_reparam_v2_fit.rds"))
# 
# # bind chains together
# SRH_fit_raw <- bind4chains(SRH_chain1, SRH_chain2, SRH_chain3, SRH_chain4)
# # thin2
# thin_draws(SRH_fit_raw, thin = 2) -> SRH_fit
# # summarise
# SRH_fit_summary <- summarise_draws(SRH_fit)

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

# ### UCW ###
# UCW_parameters <- UCW_fit_summary$variable
# UCW_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", UCW_parameters)] -> UCW_fixed_effects
# UCW_fixed_effects[!(grepl("vector", UCW_fixed_effects))] -> UCW_fixed_effects
# 
# b0_array_UCW <- make_parameter_draws_array(parameter_prefix = "b0", fit = UCW_fit, fit_summary = UCW_fit_summary)
# btemp0_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = UCW_fit, fit_summary = UCW_fit_summary)
# btemp1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = UCW_fit, fit_summary = UCW_fit_summary)
# bspillwindow_array_UCW <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = UCW_fit, fit_summary = UCW_fit_summary)
# bwinterspill_array_UCW <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = UCW_fit, fit_summary = UCW_fit_summary)
# btemp0xorigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
# btemp1xorigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
# btemp0xorigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
# btemp1xorigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
# btemp0xorigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)
# btemp1xorigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)
# borigin1_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = UCW_fit, fit_summary = UCW_fit_summary)
# borigin2_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = UCW_fit, fit_summary = UCW_fit_summary)
# borigin3_array_UCW <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = UCW_fit, fit_summary = UCW_fit_summary)
# 
# ### UCH ###
# UCH_parameters <- UCH_fit_summary$variable
# UCH_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", UCH_parameters)] -> UCH_fixed_effects
# UCH_fixed_effects[!(grepl("vector", UCH_fixed_effects))] -> UCH_fixed_effects
# 
# b0_array_UCH <- make_parameter_draws_array(parameter_prefix = "b0", fit = UCH_fit, fit_summary = UCH_fit_summary)
# btemp0_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = UCH_fit, fit_summary = UCH_fit_summary)
# btemp1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = UCH_fit, fit_summary = UCH_fit_summary)
# bspillwindow_array_UCH <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = UCH_fit, fit_summary = UCH_fit_summary)
# bwinterspill_array_UCH <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = UCH_fit, fit_summary = UCH_fit_summary)
# btemp0xorigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
# btemp1xorigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
# btemp0xorigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
# btemp1xorigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
# btemp0xorigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)
# btemp1xorigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)
# borigin1_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = UCH_fit, fit_summary = UCH_fit_summary)
# borigin2_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = UCH_fit, fit_summary = UCH_fit_summary)
# borigin3_array_UCH <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = UCH_fit, fit_summary = UCH_fit_summary)

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
# SRW_parameters <- SRW_fit_summary$variable
# SRW_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRW_parameters)] -> SRW_fixed_effects
# SRW_fixed_effects[!(grepl("vector", SRW_fixed_effects))] -> SRW_fixed_effects
# 
# b0_array_SRW <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp0_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRW_fit, fit_summary = SRW_fit_summary)
# bspillwindow_array_SRW <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRW_fit, fit_summary = SRW_fit_summary)
# bwinterspill_array_SRW <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp0xorigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp1xorigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp0xorigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp1xorigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp0xorigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp1xorigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp0xorigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp1xorigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp0xorigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp1xorigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp0xorigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)
# btemp1xorigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)
# borigin1_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRW_fit, fit_summary = SRW_fit_summary)
# borigin2_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRW_fit, fit_summary = SRW_fit_summary)
# borigin3_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRW_fit, fit_summary = SRW_fit_summary)
# borigin4_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRW_fit, fit_summary = SRW_fit_summary)
# borigin5_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRW_fit, fit_summary = SRW_fit_summary)
# borigin6_array_SRW <- make_parameter_draws_array(parameter_prefix = "borigin6", fit = SRW_fit, fit_summary = SRW_fit_summary)
# 
# ### SRH ###
# SRH_parameters <- SRH_fit_summary$variable
# SRH_parameters[grepl("b0|btemp1|btemp0|bspillwindow|bwinterspill|borigin", SRH_parameters)] -> SRH_fixed_effects
# SRH_fixed_effects[!(grepl("vector", SRH_fixed_effects))] -> SRH_fixed_effects
# 
# b0_array_SRH <- make_parameter_draws_array(parameter_prefix = "b0", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp0_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1", fit = SRH_fit, fit_summary = SRH_fit_summary)
# bspillwindow_array_SRH <- make_parameter_draws_array(parameter_prefix = "bspillwindow", fit = SRH_fit, fit_summary = SRH_fit_summary)
# bwinterspill_array_SRH <- make_parameter_draws_array(parameter_prefix = "bwinterspill", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp0xorigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp1xorigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp0xorigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp1xorigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp0xorigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp1xorigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp0xorigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp1xorigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp0xorigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp0xorigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)
# btemp1xorigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "btemp1xorigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)
# borigin1_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin1", fit = SRH_fit, fit_summary = SRH_fit_summary)
# borigin2_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin2", fit = SRH_fit, fit_summary = SRH_fit_summary)
# borigin3_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin3", fit = SRH_fit, fit_summary = SRH_fit_summary)
# borigin4_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin4", fit = SRH_fit, fit_summary = SRH_fit_summary)
# borigin5_array_SRH <- make_parameter_draws_array(parameter_prefix = "borigin5", fit = SRH_fit, fit_summary = SRH_fit_summary)

#### Extract transition data by spill for rug plot ####


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





#### Spill window plot ####

# We will plot all of these - they affect fallback through mainstem states
# Note that spill window parameter, like spill days, is shared for the DPS. So we 
# will only notice differences in overall movement probability in from
# states that have origin effects for other movements (e.g., fallback over BON
# will have same prob for different UC origins, fallback over RIS will not)

# to capture posterior correlations, we'll need to select that same iteration
# for all parameters that affect movement

# Tell it the movements for which you want to estimate temperature effects
# movements are formatted as matrix, with column for from and column for to
estimate_spillwindow_effect_UCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- UCW_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                                   date = as.vector(UCW_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- UCW_envir$data$temperature
    med_temp <- median(temp_data[subset(UCW_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- UCW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_UCW[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_UCW[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_UCW[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_UCW[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_UCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_UCW[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_UCW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 borigin1_array_UCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_UCW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_UCW[movements$from[i],movements$to[i],iter]*origin3)/
          sum(exp(b0_array_UCW[movements$from[i],possible_movements,iter] +
                    # btemp0_array_UCW[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCW[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_UCW[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_UCW[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCW[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_UCW[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_UCW[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    borigin1_array_UCW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_UCW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_UCW[movements$from[i],possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}
estimate_spillwindow_effect_UCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- UCH_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                                   date = as.vector(UCH_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- UCH_envir$data$temperature
    med_temp <- median(temp_data[subset(UCH_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- UCH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_UCH[movements$from[i],movements$to[i],iter] +
                                                        # btemp0_array_UCH[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                        btemp1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp + 
                                                        bspillwindow_array_UCH[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                        # bwinterspill_array_UCH[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                        # btemp0xorigin1_array_UCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                        btemp1xorigin1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                        # btemp0xorigin2_array_UCH[movements$from[i],movements$to[i],iter]*origin2 + 
                                                        btemp1xorigin2_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                        # btemp0xorigin3_array_UCH[movements$from[i],movements$to[i],iter]*origin3 +
                                                        btemp1xorigin3_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                        borigin1_array_UCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                        borigin2_array_UCH[movements$from[i],movements$to[i],iter]*origin2 +
                                                        borigin3_array_UCH[movements$from[i],movements$to[i],iter]*origin3)/
          sum(exp(b0_array_UCH[movements$from[i],possible_movements,iter] +
                    # btemp0_array_UCH[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_UCH[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_UCH[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_UCH[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_UCH[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_UCH[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_UCH[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    borigin1_array_UCH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_UCH[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_UCH[movements$from[i],possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}
estimate_spillwindow_effect_MCW <- function(origin_select, movements){
  
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
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- MCW_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                                   date = as.vector(MCW_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- MCW_envir$data$temperature
    med_temp <- median(temp_data[subset(MCW_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- MCW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_MCW[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_MCW[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_MCW[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_MCW[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_MCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_MCW[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_MCW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp0xorigin4_array_MCW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 btemp1xorigin4_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp0xorigin5_array_MCW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 btemp1xorigin5_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp0xorigin6_array_MCW[movements$from[i],movements$to[i],iter]*origin6 +
                                                 btemp1xorigin6_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 borigin1_array_MCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_MCW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_MCW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_MCW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_MCW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 borigin6_array_MCW[movements$from[i],movements$to[i],iter]*origin6)/
          sum(exp(b0_array_MCW[movements$from[i],possible_movements,iter] +
                    # btemp0_array_MCW[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCW[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_MCW[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_MCW[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCW[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_MCW[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_MCW[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp0xorigin4_array_MCW[movements$from[i],possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp0xorigin5_array_MCW[movements$from[i],possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp0xorigin6_array_MCW[movements$from[i],possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    borigin1_array_MCW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_MCW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_MCW[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_MCW[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_MCW[movements$from[i],possible_movements,iter]*origin5 +
                    borigin6_array_MCW[movements$from[i],possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spillwindow_effect_MCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- MCH_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                                   date = as.vector(MCH_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- MCH_envir$data$temperature
    med_temp <- median(temp_data[subset(MCH_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- MCH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_MCH[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_MCH[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_MCH[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_MCH[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_MCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_MCH[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 borigin1_array_MCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_MCH[movements$from[i],movements$to[i],iter]*origin2)/
          sum(exp(b0_array_MCH[movements$from[i],possible_movements,iter] +
                    # btemp0_array_MCH[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_MCH[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_MCH[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_MCH[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_MCH[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_MCH[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    borigin1_array_MCH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_MCH[movements$from[i],possible_movements,iter]*origin2))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spillwindow_effect_SRW <- function(origin_select, movements){
  
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
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- SRW_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
                                   date = as.vector(SRW_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- SRW_envir$data$temperature
    med_temp <- median(temp_data[subset(SRW_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- SRW_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_SRW[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_SRW[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_SRW[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_SRW[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_SRW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_SRW[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_SRW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp0xorigin4_array_SRW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 btemp1xorigin4_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp0xorigin5_array_SRW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 btemp1xorigin5_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp0xorigin6_array_SRW[movements$from[i],movements$to[i],iter]*origin6 +
                                                 btemp1xorigin6_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 borigin1_array_SRW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_SRW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_SRW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_SRW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_SRW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 borigin6_array_SRW[movements$from[i],movements$to[i],iter]*origin6)/
          sum(exp(b0_array_SRW[movements$from[i],possible_movements,iter] +
                    # btemp0_array_SRW[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRW[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_SRW[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_SRW[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRW[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_SRW[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_SRW[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp0xorigin4_array_SRW[movements$from[i],possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp0xorigin5_array_SRW[movements$from[i],possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp0xorigin6_array_SRW[movements$from[i],possible_movements,iter]*origin6 +
                    btemp1xorigin6_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    borigin1_array_SRW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_SRW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_SRW[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_SRW[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_SRW[movements$from[i],possible_movements,iter]*origin5 +
                    borigin6_array_SRW[movements$from[i],possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spillwindow_effect_SRH <- function(origin_select, movements){
  
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
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill values to predict across
    # extract the data
    spill_window_data <- SRH_envir$data$spill_window_data
    # create a sequence of values from 1 to max
    spill_predict <- seq(0, max(spill_window_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
                                   date = as.vector(SRH_envir$data$transition_dates))
    
    # get median temperature for this state
    temp_data <- SRH_envir$data$temperature
    med_temp <- median(temp_data[subset(SRH_states_dates, state == from)$date,from])
    
    # get median winter spill for this state
    winterspill_data <- SRH_envir$data$winter_spill_days_data
    med_winterspill <- median(winterspill_data[,from])
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_SRH[movements$from[i],movements$to[i],iter] +
                                                 # btemp0_array_SRH[movements$from[i],movements$to[i],iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                                                 btemp1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp + 
                                                 bspillwindow_array_SRH[movements$from[i],movements$to[i],iter]*spill_predict[j] +
                                                 # bwinterspill_array_SRH[movements$from[i],movements$to[i],iter]*med_winterspill + # indexing is set up so that if spill window is selected, spill days is ignored
                                                 # btemp0xorigin1_array_SRH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 btemp1xorigin1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp0xorigin2_array_SRH[movements$from[i],movements$to[i],iter]*origin2 + 
                                                 btemp1xorigin2_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 # btemp0xorigin3_array_SRH[movements$from[i],movements$to[i],iter]*origin3 +
                                                 btemp1xorigin3_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp0xorigin4_array_SRH[movements$from[i],movements$to[i],iter]*origin4 +
                                                 btemp1xorigin4_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp0xorigin5_array_SRH[movements$from[i],movements$to[i],iter]*origin5 +
                                                 btemp1xorigin5_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 borigin1_array_SRH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_SRH[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_SRH[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_SRH[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_SRH[movements$from[i],movements$to[i],iter]*origin5)/
          sum(exp(b0_array_SRH[movements$from[i],possible_movements,iter] +
                    # btemp0_array_SRH[movements$from[i],possible_movements,iter] + # ignore this parameter - because the vast majority of transitions occur in summer/fall, which corresponds to temp1 param
                    btemp1_array_SRH[movements$from[i],possible_movements,iter]*med_temp + 
                    bspillwindow_array_SRH[movements$from[i],possible_movements,iter]*spill_predict[j] +
                    # bwinterspill_array_SRH[movements$from[i],possible_movements,iter]*med_winterspill +
                    # btemp0xorigin1_array_SRH[movements$from[i],possible_movements,iter]*origin1 +
                    btemp1xorigin1_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp0xorigin2_array_SRH[movements$from[i],possible_movements,iter]*origin2 + 
                    btemp1xorigin2_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    # btemp0xorigin3_array_SRH[movements$from[i],possible_movements,iter]*origin3 +
                    btemp1xorigin3_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp0xorigin4_array_SRH[movements$from[i],possible_movements,iter]*origin4 +
                    btemp1xorigin4_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp0xorigin5_array_SRH[movements$from[i],possible_movements,iter]*origin5 +
                    btemp1xorigin5_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    borigin1_array_SRH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_SRH[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_SRH[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_SRH[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_SRH[movements$from[i],possible_movements,iter]*origin5))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}




#### plot spill window function ####

# create df to index to right dam temp for plotting
dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                        state = seq(2,9))

# another plot option to compare between hatchery and wild
rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")

# set up post-overshoot info so that we don't select fish that got the spill days covariate instead
post_overshoot_combos <- list(c(NA), # Asotin Creek
                              c(NA), # Clearwater River
                              c(3,4,5,6,7,8,9), # Deschutes River
                              c(7), # Entiat River
                              c(3,4,5,6,7,8,9), # Fifteenmile Creek
                              c(NA), # Grande Ronde River
                              c(3,4,5,6,7,8,9), # Hood River
                              c(NA), # Imnaha River
                              c(3,4,5,6,7,8,9), # John Day River
                              c(NA), # Methow River
                              c(NA), # Okanogan River
                              c(NA), # Salmon River
                              c(9), # Tucannon River
                              c(3,4,5,6,7,8,9), # Umatilla River
                              c(4,5,6,7,8,9), # Walla Walla River
                              c(6,7), # Wenatchee River
                              c(4,5,6,7,8,9)) # Yakima River

names(post_overshoot_combos) <- natal_origins[order(natal_origins)]


plot_compare_rear_spill_effect <- function(origin_select,
                                           wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                           wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                          movements_evaluated, spill_predict,
                                          from, to, plot_title = NULL){
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences

    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    hatchery_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
  } 
  # else run both
  else {
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles %>% 
      bind_rows(., hatchery_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    
  }
  
  # Remove any instances for rug plot of fish that would have received winter spill days covariate
  # if there aren't any for this combination, don't try to manipulate it
  if(nrow(covariate_experiences) != 0) {
    covariate_experiences %>% 
      dplyr::rename(date_numeric = date) %>% 
      # keep only jan/feb/mar 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) %>% 
      filter(!(state %in% post_overshoot_combos[[origin_select]] & month %in% c(1,2,3))) -> covariate_experiences
    
  } else {
    
    
  }

  
  rear_spill_move_prob_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_rug(data = subset(covariate_experiences, rear == "wild"), aes(x = spill_window*100, color = rear), inherit.aes = FALSE,
             sides = "t", length = unit(0.3, "cm")) +
    geom_rug(data = subset(covariate_experiences, rear == "hatchery"), aes(x = spill_window*100, color = rear), inherit.aes = FALSE,
             sides = "b", length = unit(0.3, "cm")) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(-5,NA), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab("Spill (kcfs)") +
    ylab("Movement probability") +
    ggtitle(plot_title)
    
  return(rear_spill_move_prob_plot)
}


#### Plot rear type comparison fallback spill window for all ####
mainstem_fallback_movements <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
  from = c(2,3,4,5,6,7,8,9), to = c(1,2,3,4,5,6,3,8))

# #### Upper Columbia ####
# UCH_spill_window_data <- UCH_envir$data$spill_window_data
# UCW_spill_window_data <- UCW_envir$data$spill_window_data
# UC_spill_window_data <- bind_rows(as.data.frame(UCH_spill_window_data), as.data.frame(UCW_spill_window_data))
# 
# # Wenatchee - evaluate all fallback move probs
# WEN_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Wenatchee River", movements = mainstem_fallback_movements)
# WEN_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Wenatchee River", movements = mainstem_fallback_movements)
# 
# # Wenatchee - get covariate experiences
# WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
# WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")
# 
# # Wenatchee - use a loop to plot all of them
# WEN_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   WEN_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Wenatchee River",
#                                                                         wild_move_prob_array = WEN_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = WEN_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = WEN_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Wenatchee - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("WEN_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          WEN_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Entiat - evaluate all fallback move probs
# # Only wild for entiat
# ENT_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Entiat River", movements = mainstem_fallback_movements)
# 
# # Entiat - get covariate experiences
# ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
# 
# # Entiat - use a loop to plot all of them
# ENT_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   ENT_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Entiat River",
#                                                                         wild_move_prob_array = ENT_wild_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = ENT_wild_covariate_experiences, 
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Entiat - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("ENT_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          ENT_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Methow - evaluate all fallback move probs
# MET_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCW(origin_select= "Methow River", movements = mainstem_fallback_movements)
# MET_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Methow River", movements = mainstem_fallback_movements)
# 
# # Methow - get covariate experiences
# MET_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Methow River")
# MET_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Methow River")
# 
# # Methow - use a loop to plot all of them
# MET_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   MET_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Methow River",
#                                                                         wild_move_prob_array = MET_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = MET_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = MET_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = MET_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Methow - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("MET_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          MET_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Okanoagn - evaluate all fallback move probs
# # Only hatchery for Okanogan
# OKA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_UCH(origin_select= "Okanogan River", movements = mainstem_fallback_movements)
# 
# # Okanogan - get covariate experiences
# OKA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Okanogan River")
# 
# # Okanogan - use a loop to plot all of them
# OKA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   OKA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Okanogan River",
#                                                                         hatchery_move_prob_array = OKA_hatchery_spillwindow_move_prob_array,
#                                                                         hatchery_covariate_experiences = OKA_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(UC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Okanogan - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("OKA_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          OKA_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }

#### Middle Columbia ####
MCH_spill_window_data <- MCH_envir$data$spill_window_data
MCW_spill_window_data <- MCW_envir$data$spill_window_data
MC_spill_window_data <- bind_rows(as.data.frame(MCH_spill_window_data), as.data.frame(MCW_spill_window_data))

# Deschutes - evaluate all fallback move probs
# Deschutes is wild only
DES_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Deschutes River", movements = mainstem_fallback_movements)

# Deschutes - get covariate experiences
DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")

# Deschutes - use a loop to plot all of them
DES_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
for (i in 1:nrow(mainstem_fallback_movements)){
  DES_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Deschutes River",
                                                                        wild_move_prob_array = DES_wild_spillwindow_move_prob_array,
                                                                        
                                                                        wild_covariate_experiences = DES_wild_covariate_experiences, 
                                                                        
                                                                        movements_evaluated = mainstem_fallback_movements,
                                                                        spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
                                                                        from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
                                                                        plot_title = paste0("Deschutes - fallback at ", mainstem_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("DES_compare_fallback_", 
                                                                                        mainstem_fallback_movements$dam[i],
                                                                                        "_spillwindow.png")),
         DES_fallback_spillwindow_plots[[i]], height = 8, width = 8)
  
}

# John Day - evaluate all fallback move probs
# JDR is wild only
JDR_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "John Day River", movements = mainstem_fallback_movements)

# John Day - get covariate experiences
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")

# John Day - use a loop to plot all of them
JDR_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
for (i in 1:nrow(mainstem_fallback_movements)){
  JDR_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "John Day River",
                                                                        wild_move_prob_array = JDR_wild_spillwindow_move_prob_array,
                                                                        
                                                                        wild_covariate_experiences = JDR_wild_covariate_experiences, 
                                                                        
                                                                        movements_evaluated = mainstem_fallback_movements,
                                                                        spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
                                                                        from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
                                                                        plot_title = paste0("John Day - fallback at ", mainstem_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("JDR_compare_fallback_", 
                                                                                        mainstem_fallback_movements$dam[i],
                                                                                        "_spillwindow.png")),
         JDR_fallback_spillwindow_plots[[i]], height = 8, width = 8)
  
}


# Fifteenmile - evaluate all fallback move probs
# Fifteenmile is wild only
FIF_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Fifteenmile Creek", movements = mainstem_fallback_movements)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")

# Fifteenmile - use a loop to plot all of them
FIF_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
for (i in 1:nrow(mainstem_fallback_movements)){
  FIF_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Fifteenmile Creek",
                                                                        wild_move_prob_array = FIF_wild_spillwindow_move_prob_array,
                                                                        
                                                                        wild_covariate_experiences = FIF_wild_covariate_experiences, 
                                                                        
                                                                        movements_evaluated = mainstem_fallback_movements,
                                                                        spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
                                                                        from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
                                                                        plot_title = paste0("Fifteenmile - fallback at ", mainstem_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("FIF_compare_fallback_", 
                                                                                        mainstem_fallback_movements$dam[i],
                                                                                        "_spillwindow.png")),
         FIF_fallback_spillwindow_plots[[i]], height = 8, width = 8)
  
}

# Umatilla - evaluate all fallback move probs
UMA_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Umatilla River", movements = mainstem_fallback_movements)
UMA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCH(origin_select= "Umatilla River", movements = mainstem_fallback_movements)

# Umatilla - get covariate experiences
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")

# Umatilla - use a loop to plot all of them
UMA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
for (i in 1:nrow(mainstem_fallback_movements)){
  UMA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Umatilla River",
                                                                        wild_move_prob_array = UMA_wild_spillwindow_move_prob_array,
                                                                        hatchery_move_prob_array = UMA_hatchery_spillwindow_move_prob_array,
                                                                        wild_covariate_experiences = UMA_wild_covariate_experiences, 
                                                                        hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                        movements_evaluated = mainstem_fallback_movements,
                                                                        spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
                                                                        from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
                                                                        plot_title = paste0("Umatilla - fallback at ", mainstem_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("UMA_compare_fallback_", 
                                                                                        mainstem_fallback_movements$dam[i],
                                                                                        "_spillwindow.png")),
         UMA_fallback_spillwindow_plots[[i]], height = 8, width = 8)
  
}

# Yakima - evaluate all fallback move probs
# Yakima is wild only
YAK_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Yakima River", movements = mainstem_fallback_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")

# Yakima - use a loop to plot all of them
YAK_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
for (i in 1:nrow(mainstem_fallback_movements)){
  YAK_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Yakima River",
                                                                        wild_move_prob_array = YAK_wild_spillwindow_move_prob_array,
                                                                        
                                                                        wild_covariate_experiences = YAK_wild_covariate_experiences, 
                                                                        
                                                                        movements_evaluated = mainstem_fallback_movements,
                                                                        spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
                                                                        from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
                                                                        plot_title = paste0("Yakima - fallback at ", mainstem_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("YAK_compare_fallback_", 
                                                                                        mainstem_fallback_movements$dam[i],
                                                                                        "_spillwindow.png")),
         YAK_fallback_spillwindow_plots[[i]], height = 8, width = 8)
  
}

# Walla Walla - evaluate all fallback move probs
WAWA_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCW(origin_select= "Walla Walla River", movements = mainstem_fallback_movements)
WAWA_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_MCH(origin_select= "Walla Walla River", movements = mainstem_fallback_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")

# Walla Walla - use a loop to plot all of them
WAWA_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
for (i in 1:nrow(mainstem_fallback_movements)){
  WAWA_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Walla Walla River",
                                                                         wild_move_prob_array = WAWA_wild_spillwindow_move_prob_array,
                                                                         hatchery_move_prob_array = WAWA_hatchery_spillwindow_move_prob_array,
                                                                         wild_covariate_experiences = WAWA_wild_covariate_experiences, 
                                                                         hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                         movements_evaluated = mainstem_fallback_movements,
                                                                         spill_predict = seq(0, max(MC_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
                                                                         plot_title = paste0("Walla Walla - fallback at ", mainstem_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("WAWA_compare_fallback_", 
                                                                                        mainstem_fallback_movements$dam[i],
                                                                                        "_spillwindow.png")),
         WAWA_fallback_spillwindow_plots[[i]], height = 8, width = 8)
  
}

# #### Snake River ####
# SRH_spill_window_data <- SRH_envir$data$spill_window_data
# SRW_spill_window_data <- SRW_envir$data$spill_window_data
# SR_spill_window_data <- bind_rows(as.data.frame(SRH_spill_window_data), as.data.frame(SRW_spill_window_data))
# 
# # Asotin Creek - evaluate all fallback move probs
# # Asotin is wild only
# ASO_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Asotin Creek", movements = mainstem_fallback_movements)
# 
# # Asotin Creek - get covariate experiences
# ASO_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Asotin Creek")
# 
# # Asotin Creek - use a loop to plot all of them
# ASO_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   ASO_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Asotin Creek",
#                                                                         wild_move_prob_array = ASO_wild_spillwindow_move_prob_array,
#                                                                         
#                                                                         wild_covariate_experiences = ASO_wild_covariate_experiences, 
#                                                                         
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Asotin Creek - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("ASO_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          ASO_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Clearwater River - evaluate all fallback move probs
# CLE_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Clearwater River", movements = mainstem_fallback_movements)
# CLE_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Clearwater River", movements = mainstem_fallback_movements)
# 
# # Clearwater River - get covariate experiences
# CLE_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Clearwater River")
# CLE_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Clearwater River")
# 
# # Clearwater River - use a loop to plot all of them
# CLE_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   CLE_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Clearwater River",
#                                                                         wild_move_prob_array = CLE_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = CLE_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = CLE_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = CLE_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Clearwater River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("CLE_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          CLE_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Grande Ronde River - evaluate all fallback move probs
# GR_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Grande Ronde River", movements = mainstem_fallback_movements)
# GR_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Grande Ronde River", movements = mainstem_fallback_movements)
# 
# # Grande Ronde River - get covariate experiences
# GR_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Grande Ronde River")
# GR_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Grande Ronde River")
# 
# # Grande Ronde River - use a loop to plot all of them
# GR_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   GR_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Grande Ronde River",
#                                                                        wild_move_prob_array = GR_wild_spillwindow_move_prob_array,
#                                                                        hatchery_move_prob_array = GR_hatchery_spillwindow_move_prob_array,
#                                                                        wild_covariate_experiences = GR_wild_covariate_experiences, 
#                                                                        hatchery_covariate_experiences = GR_hatchery_covariate_experiences,
#                                                                        movements_evaluated = mainstem_fallback_movements,
#                                                                        spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                        from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                        plot_title = paste0("Grande Ronde River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("GR_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          GR_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Imnaha River - evaluate all fallback move probs
# IMN_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Imnaha River", movements = mainstem_fallback_movements)
# IMN_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Imnaha River", movements = mainstem_fallback_movements)
# 
# # Imnaha River - get covariate experiences
# IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Imnaha River")
# IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Imnaha River")
# 
# # Imnaha River - use a loop to plot all of them
# IMN_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   IMN_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Imnaha River",
#                                                                         wild_move_prob_array = IMN_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = IMN_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = IMN_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = IMN_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Imnaha River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("IMN_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          IMN_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# # Salmon River - evaluate all fallback move probs
# SAL_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Salmon River", movements = mainstem_fallback_movements)
# SAL_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Salmon River", movements = mainstem_fallback_movements)
# 
# # Salmon River - get covariate experiences
# SAL_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Salmon River")
# SAL_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Salmon River")
# 
# # Salmon River - use a loop to plot all of them
# SAL_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   SAL_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Salmon River",
#                                                                         wild_move_prob_array = SAL_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = SAL_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = SAL_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = SAL_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Salmon River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("SAL_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          SAL_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# # Tucannon River - evaluate all fallback move probs
# TUC_wild_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRW(origin_select= "Tucannon River", movements = mainstem_fallback_movements)
# TUC_hatchery_spillwindow_move_prob_array <- estimate_spillwindow_effect_SRH(origin_select= "Tucannon River", movements = mainstem_fallback_movements)
# 
# # Tucannon River - get covariate experiences
# TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
# TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")
# 
# # Tucannon River - use a loop to plot all of them
# TUC_fallback_spillwindow_plots <- vector(mode = "list", length = nrow(mainstem_fallback_movements))
# for (i in 1:nrow(mainstem_fallback_movements)){
#   TUC_fallback_spillwindow_plots[[i]] <- plot_compare_rear_spill_effect(origin_select = "Tucannon River",
#                                                                         wild_move_prob_array = TUC_wild_spillwindow_move_prob_array,
#                                                                         hatchery_move_prob_array = TUC_hatchery_spillwindow_move_prob_array,
#                                                                         wild_covariate_experiences = TUC_wild_covariate_experiences, 
#                                                                         hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
#                                                                         movements_evaluated = mainstem_fallback_movements,
#                                                                         spill_predict = seq(0, max(SR_spill_window_data[,mainstem_fallback_movements$from[i]]),length = 100),
#                                                                         from = mainstem_fallback_movements$from[i], to = mainstem_fallback_movements$to[i], 
#                                                                         plot_title = paste0("Tucannon River - fallback at ", mainstem_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spillwindow", paste0("TUC_compare_fallback_", 
#                                                                                         mainstem_fallback_movements$dam[i],
#                                                                                         "_spillwindow.png")),
#          TUC_fallback_spillwindow_plots[[i]], height = 8, width = 8)
#   
# }
# 
# 

#### Spill days movement probability function ####
estimate_spilldays_effect_UCW <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  wild_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  wild_origin_params[subset(origin_param_map, natal_origin == origin_select)$wild] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = wild_origin_params[1]
  origin2 = wild_origin_params[2]
  origin3 = wild_origin_params[3]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- UCW_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- UCW_envir$data$movements[, "col"][UCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                                   date_numeric = as.vector(UCW_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    UCW_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> UCW_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(UCW_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(UCW_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_UCW[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_UCW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_UCW[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_UCW[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_UCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 borigin1_array_UCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_UCW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_UCW[movements$from[i],movements$to[i],iter]*origin3)/
          sum(exp(b0_array_UCW[movements$from[i],possible_movements,iter] +
                    btemp0_array_UCW[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_UCW[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_UCW[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_UCW[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_UCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    borigin1_array_UCW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_UCW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_UCW[movements$from[i],possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spilldays_effect_UCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,3)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  origin3 = hatchery_origin_params[3]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- UCH_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- UCH_envir$data$movements[, "col"][UCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    UCH_states_dates <- data.frame(state = as.vector(UCH_envir$data$y),
                                   date_numeric = as.vector(UCH_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    UCH_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> UCH_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(UCH_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(UCH_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_UCH[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_UCH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_UCH[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_UCH[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_UCH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 borigin1_array_UCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_UCH[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_UCH[movements$from[i],movements$to[i],iter]*origin3)/
          sum(exp(b0_array_UCH[movements$from[i],possible_movements,iter] +
                    btemp0_array_UCH[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_UCH[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_UCH[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_UCH[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_UCH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    borigin1_array_UCH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_UCH[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_UCH[movements$from[i],possible_movements,iter]*origin3))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}


estimate_spilldays_effect_MCW <- function(origin_select, movements){
  
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
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- MCW_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- MCW_envir$data$movements[, "col"][MCW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCW_states_dates <- data.frame(state = as.vector(MCW_envir$data$y),
                                   date_numeric = as.vector(MCW_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    MCW_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> MCW_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(MCW_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(MCW_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_MCW[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_MCW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_MCW[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_MCW[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 btemp0xorigin4_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp1xorigin4_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 btemp0xorigin5_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp1xorigin5_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 btemp0xorigin6_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 # btemp1xorigin6_array_MCW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 borigin1_array_MCW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_MCW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_MCW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_MCW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_MCW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 borigin6_array_MCW[movements$from[i],movements$to[i],iter]*origin6)/
          sum(exp(b0_array_MCW[movements$from[i],possible_movements,iter] +
                    btemp0_array_MCW[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_MCW[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_MCW[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_MCW[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    btemp0xorigin4_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp1xorigin4_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    btemp0xorigin5_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp1xorigin5_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    btemp0xorigin6_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    # btemp1xorigin6_array_MCW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    borigin1_array_MCW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_MCW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_MCW[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_MCW[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_MCW[movements$from[i],possible_movements,iter]*origin5 +
                    borigin6_array_MCW[movements$from[i],possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spilldays_effect_MCH <- function(origin_select, movements){
  
  # set up vectors to control which origin parameters are turned on
  hatchery_origin_params <- rep(0,2)
  
  # index to the right origin param to turn it on
  hatchery_origin_params[subset(origin_param_map, natal_origin == origin_select)$hatchery] <- 1
  
  # use this vector to set all of the individual origin indices
  origin1 = hatchery_origin_params[1]
  origin2 = hatchery_origin_params[2]
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- MCH_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- MCH_envir$data$movements[, "col"][MCH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    MCH_states_dates <- data.frame(state = as.vector(MCH_envir$data$y),
                                   date_numeric = as.vector(MCH_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    MCH_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> MCH_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(MCH_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(MCH_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_MCH[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_MCH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_MCH[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_MCH[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_MCH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 borigin1_array_MCH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_MCH[movements$from[i],movements$to[i],iter]*origin2)/
          sum(exp(b0_array_MCH[movements$from[i],possible_movements,iter] +
                    btemp0_array_MCH[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_MCH[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_MCH[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_MCH[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_MCH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    borigin1_array_MCH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_MCH[movements$from[i],possible_movements,iter]*origin2))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spilldays_effect_SRW <- function(origin_select, movements){
  
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
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- SRW_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- SRW_envir$data$movements[, "col"][SRW_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRW_states_dates <- data.frame(state = as.vector(SRW_envir$data$y),
                                   date_numeric = as.vector(SRW_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    SRW_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> SRW_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(SRW_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(SRW_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_SRW[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_SRW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_SRW[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_SRW[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 btemp0xorigin4_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp1xorigin4_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 btemp0xorigin5_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp1xorigin5_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 btemp0xorigin6_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 # btemp1xorigin6_array_SRW[movements$from[i],movements$to[i],iter]*med_temp*origin6 +
                                                 borigin1_array_SRW[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_SRW[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_SRW[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_SRW[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_SRW[movements$from[i],movements$to[i],iter]*origin5 +
                                                 borigin6_array_SRW[movements$from[i],movements$to[i],iter]*origin6)/
          sum(exp(b0_array_SRW[movements$from[i],possible_movements,iter] +
                    btemp0_array_SRW[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_SRW[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_SRW[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_SRW[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    btemp0xorigin4_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp1xorigin4_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    btemp0xorigin5_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp1xorigin5_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    btemp0xorigin6_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    # btemp1xorigin6_array_SRW[movements$from[i],possible_movements,iter]*med_temp*origin6 +
                    borigin1_array_SRW[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_SRW[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_SRW[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_SRW[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_SRW[movements$from[i],possible_movements,iter]*origin5 +
                    borigin6_array_SRW[movements$from[i],possible_movements,iter]*origin6))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}

estimate_spilldays_effect_SRH <- function(origin_select, movements){
  
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
  
  niter <- 4000 # this is the number of draws we have
  
  # get an array to store probabilities of movements at different spill
  # make this more memory efficient - store instead as array, with row names for from/to
  spill_move_prob_array <- array(dim = c(100, niter, nrow(movements)))
  dimnames(spill_move_prob_array)[[3]] <- paste0(movements$from, "_", movements$to)
  
  
  # spill_move_prob_array <- array(dim = c(length(model_states), length(model_states), 100, niter))
  
  
  for (i in 1:nrow(movements)){
    from <- movements$from[i]
    to <- movements$to[i]
    
    # set up spill window values to predict across
    # extract the data
    winterspill_data <- SRH_envir$data$winter_spill_days_data
    # create a sequence of values from 1 to max
    spill_days_predict <- seq(0, max(winterspill_data[,from]),length = 100)
    
    # Get the movements
    possible_movements <- SRH_envir$data$movements[, "col"][SRH_envir$data$movements[, "row"] == from]
    possible_movements <- c(possible_movements, 43)
    
    # make sure to drop all movements to upstream states (in DE years, these aren't possible)
    grep(" Upstream", model_states) -> upstream_indices
    possible_movements <- possible_movements[!(possible_movements %in% upstream_indices)]
    
    SRH_states_dates <- data.frame(state = as.vector(SRH_envir$data$y),
                                   date_numeric = as.vector(SRH_envir$data$transition_dates))
    
    # add month to states dates so that we can subset
    SRH_states_dates %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> SRH_states_dates
    
    # get median temperature for this state - only in Jan/Feb/Mar
    temp_data <- as.data.frame(SRH_envir$data$temperature)
    temp_data %>% 
      rownames_to_column("date_numeric") %>% 
      pivot_longer(cols = -c(date_numeric)) -> temp_data_long
    
    temp_data_long %>% 
      mutate(date = ymd("2005-05-31") + days(date_numeric)) %>% 
      mutate(month = month(date)) -> temp_data_long
    
    # only keep actual observed fish
    med_temp <- median(temp_data[subset(SRH_states_dates, state == from & month %in% c(1,2,3))$date_numeric,from])
    
    # do not include spill window, because code is set up to index so that it's only one or the other
    
    # Loop through all of the iterations
    for (iter in 1:niter){
      
      # Loop through a sequence of temperature values to get predicted response
      # temperature was z-scored so we can plot 2 standard deviations
      for (j in 1:length(spill_days_predict)){
        # evaluation movement 
        spill_move_prob_array[j,iter,i] <- exp(b0_array_SRH[movements$from[i],movements$to[i],iter] +
                                                 btemp0_array_SRH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # btemp1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp + # select temp0 for any winterspill estimates
                                                 # bspillwindow_array_SRH[movements$from[i],movements$to[i],iter]*med_spillwindow + # no spill window if spill days is indexed
                                                 bwinterspill_array_SRH[movements$from[i],movements$to[i],iter]*spill_days_predict[j] +
                                                 btemp0xorigin1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 # btemp1xorigin1_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin1 +
                                                 btemp0xorigin2_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin2 +
                                                 # btemp1xorigin2_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin2 + 
                                                 btemp0xorigin3_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 # btemp1xorigin3_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin3 +
                                                 btemp0xorigin4_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 # btemp1xorigin4_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin4 +
                                                 btemp0xorigin5_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 # btemp1xorigin5_array_SRH[movements$from[i],movements$to[i],iter]*med_temp*origin5 +
                                                 borigin1_array_SRH[movements$from[i],movements$to[i],iter]*origin1 +
                                                 borigin2_array_SRH[movements$from[i],movements$to[i],iter]*origin2 +
                                                 borigin3_array_SRH[movements$from[i],movements$to[i],iter]*origin3 +
                                                 borigin4_array_SRH[movements$from[i],movements$to[i],iter]*origin4 +
                                                 borigin5_array_SRH[movements$from[i],movements$to[i],iter]*origin5)/
          sum(exp(b0_array_SRH[movements$from[i],possible_movements,iter] +
                    btemp0_array_SRH[movements$from[i],possible_movements,iter]*med_temp + # select temp0 for any winterspill estimates
                    # btemp1_array_SRH[movements$from[i],possible_movements,iter]*med_temp + 
                    # bspillwindow_array_SRH[movements$from[i],possible_movements,iter]*med_spillwindow + # no spill window if spill days is indexed
                    bwinterspill_array_SRH[movements$from[i],possible_movements,iter]*spill_days_predict[j] +
                    btemp0xorigin1_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    # btemp1xorigin1_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin1 +
                    btemp0xorigin2_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin2 +
                    # btemp1xorigin2_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin2 + 
                    btemp0xorigin3_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    # btemp1xorigin3_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin3 +
                    btemp0xorigin4_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    # btemp1xorigin4_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin4 +
                    btemp0xorigin5_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    # btemp1xorigin5_array_SRH[movements$from[i],possible_movements,iter]*med_temp*origin5 +
                    borigin1_array_SRH[movements$from[i],possible_movements,iter]*origin1 +
                    borigin2_array_SRH[movements$from[i],possible_movements,iter]*origin2 +
                    borigin3_array_SRH[movements$from[i],possible_movements,iter]*origin3 +
                    borigin4_array_SRH[movements$from[i],possible_movements,iter]*origin4 +
                    borigin5_array_SRH[movements$from[i],possible_movements,iter]*origin5))
        
      }
      
      
      
    }
    
  }
  
  # return the array that contains movement probs across all movements, temps, and iter
  return(spill_move_prob_array)
  
}



#### plot spill days function ####

# create df to index to right dam temp for plotting
dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                        state = seq(2,9))

# another plot option to compare between hatchery and wild
rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")

plot_compare_rear_spilldays_effect <- function(origin_select,
                                           wild_move_prob_array = NULL, hatchery_move_prob_array = NULL,
                                           wild_covariate_experiences = NULL, hatchery_covariate_experiences = NULL,
                                           movements_evaluated, spill_predict,
                                           from, to, plot_title = NULL){
  array_index <- which(movements_evaluated$from == from & movements_evaluated$to == to)
  
  # First, determine if this origin has both a hatchery and a wild population
  
  # If hatchery is NA, run wild only
  if (is.na(subset(origin_param_map, natal_origin == origin_select)$hatchery)){
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    
    wild_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
  }
  
  # If wild is NA, run hatchery only
  else if (is.na(subset(origin_param_map, natal_origin == origin_select)$wild)){
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    hatchery_spill_move_prob_quantiles -> rear_spill_move_prob_quantiles
    
    # convert back to kcfs
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    
    hatchery_covariate_experiences %>% 
      # keep only the state that is our from state
      filter(state == from)  -> covariate_experiences
    
  } 
  # else run both
  else {
    wild_spill_move_prob <- as.data.frame(wild_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(wild_spill_move_prob) <- paste0("iter", 1:niter) 
    wild_spill_move_prob$spill <- spill_predict
    
    wild_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> wild_spill_move_prob_quantiles
    
    hatchery_spill_move_prob <- as.data.frame(hatchery_move_prob_array[,,array_index])
    
    niter <- 4000 # this is the number of draws we have
    colnames(hatchery_spill_move_prob) <- paste0("iter", 1:niter) 
    hatchery_spill_move_prob$spill <- spill_predict
    
    hatchery_spill_move_prob %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "move_prob") %>% 
      group_by(spill) %>% 
      summarise(prob = quantile(move_prob, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> hatchery_spill_move_prob_quantiles
    
    wild_spill_move_prob_quantiles %>% 
      bind_rows(., hatchery_spill_move_prob_quantiles) -> rear_spill_move_prob_quantiles
    
    # convert back to days
    rear_spill_move_prob_quantiles %>% 
      mutate(spill_actual = spill*100) -> rear_spill_move_prob_quantiles
    
    # get data organized for rug plot
    # filter to only keep jan/feb/mar
    wild_covariate_experiences %>% 
      mutate(rear = "wild") -> wild_covariate_experiences
    hatchery_covariate_experiences %>% 
      mutate(rear = "hatchery") -> hatchery_covariate_experiences
    wild_covariate_experiences %>% 
      bind_rows(., hatchery_covariate_experiences) %>% 
      # keep only the state that is our from state
      filter(state == from) -> covariate_experiences
    
    
  }
  
  # Keep only fish where they could have experienced spill conditions for rug plot
  covariate_experiences %>% 
    filter(winter_post_overshoot_vector == 1) -> covariate_experiences
  
  rear_spill_move_prob_plot <- ggplot(rear_spill_move_prob_quantiles, aes(x = spill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, 
                                                                          color = rear, fill = rear)) +
    geom_line() +
    geom_ribbon(alpha = 0.2, color = NA) +
    geom_rug(data = subset(covariate_experiences, rear == "wild"), aes(x = winter_spill*100, color = rear, y = 0), inherit.aes = FALSE,
             sides = "t", length = unit(0.3, "cm"), position = "jitter") +
    geom_rug(data = subset(covariate_experiences, rear == "hatchery"), aes(x = winter_spill*100, color = rear, y = 0), inherit.aes = FALSE,
             sides = "b", length = unit(0.3, "cm"), position = "jitter") +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(-1,max(rear_spill_move_prob_quantiles$spill_actual)+1), expand = c(0,0)) +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    xlab("Days of winter spill") +
    ylab("Movement probability") +
    ggtitle(plot_title)
  
  return(rear_spill_move_prob_plot)
}



#### Plot rear type comparison fallback spill days for post-overshot fallback movements ####
UC_postovershoot_fallback_movements <- data.frame(dam = c("RRE", "WEL"),
                                          from = c(6,7), to = c(5,6))

MC_postovershoot_fallback_movements <- data.frame(dam = c("MCN", "PRA", "ICH"),
                                                  from = c(3,4,8), to = c(2,3,3))

SR_postovershoot_fallback_movements <- data.frame(dam = c("LGR"),
                                                  from = c(9), to = c(8))

# #### Upper Columbia ####
# # only Entiat and Wenatchee
# 
# UCH_spill_days_data <- UCH_envir$data$winter_spill_days_data
# UCW_spill_days_data <- UCW_envir$data$winter_spill_days_data
# UC_spill_days_data <- bind_rows(as.data.frame(UCH_spill_days_data), as.data.frame(UCW_spill_days_data))
# 
# # Wenatchee - evaluate all fallback move probs
# WEN_wild_spilldays_move_prob_array <- estimate_spilldays_effect_UCW(origin_select= "Wenatchee River", movements = UC_postovershoot_fallback_movements)
# WEN_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_UCH(origin_select= "Wenatchee River", movements = UC_postovershoot_fallback_movements)
# 
# # Wenatchee - get covariate experiences
# WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Wenatchee River")
# WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = UCH_envir, rear = "hatchery", origin_select = "Wenatchee River")
# 
# # Wenatchee - use a loop to plot all of them
# WEN_fallback_spilldays_plots <- vector(mode = "list", length = nrow(UC_postovershoot_fallback_movements))
# for (i in 1:nrow(UC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(WEN_wild_covariate_experiences, state == UC_postovershoot_fallback_movements$from[i]) -> WEN_data_to_plot
#   
#   if(nrow(WEN_data_to_plot)>0){
#     WEN_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Wenatchee River",
#                                                                           wild_move_prob_array = WEN_wild_spilldays_move_prob_array,
#                                                                           hatchery_move_prob_array = WEN_hatchery_spilldays_move_prob_array,
#                                                                           wild_covariate_experiences = WEN_wild_covariate_experiences, 
#                                                                           hatchery_covariate_experiences = WEN_hatchery_covariate_experiences,
#                                                                           movements_evaluated = UC_postovershoot_fallback_movements,
#                                                                           spill_predict = seq(0, max(UC_spill_days_data[,UC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                           from = UC_postovershoot_fallback_movements$from[i], to = UC_postovershoot_fallback_movements$to[i], 
#                                                                           plot_title = paste0("Wenatchee - post-overshoot fallback at ", UC_postovershoot_fallback_movements$dam[i]))
#     
#     ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("WEN_compare_fallback_", 
#                                                                                         UC_postovershoot_fallback_movements$dam[i],
#                                                                                           "_spilldays.png")),
#            WEN_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
#   
# }
# 
# # Entiat - evaluate all fallback move probs
# # Only wild for entiat
# ENT_wild_spilldays_move_prob_array <- estimate_spilldays_effect_UCW(origin_select= "Entiat River", movements = UC_postovershoot_fallback_movements)
# 
# # Entiat - get covariate experiences
# ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = UCW_envir, rear = "wild", origin_select = "Entiat River")
# 
# # Entiat - use a loop to plot all of them
# ENT_fallback_spilldays_plots <- vector(mode = "list", length = nrow(UC_postovershoot_fallback_movements))
# for (i in 1:nrow(UC_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(ENT_wild_covariate_experiences, state == UC_postovershoot_fallback_movements$from[i]) -> ENT_data_to_plot
#   
#   if(nrow(ENT_data_to_plot)>0){
#     
#     ENT_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Entiat River",
#                                                                             wild_move_prob_array = ENT_wild_spilldays_move_prob_array,
#                                                                             wild_covariate_experiences = ENT_wild_covariate_experiences, 
#                                                                             movements_evaluated = UC_postovershoot_fallback_movements,
#                                                                             spill_predict = seq(0, max(UC_spill_days_data[,UC_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                             from = UC_postovershoot_fallback_movements$from[i], to = UC_postovershoot_fallback_movements$to[i], 
#                                                                             plot_title = paste0("Entiat - post-overshoot fallback at ", UC_postovershoot_fallback_movements$dam[i]))
#     
#     ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("ENT_compare_fallback_", 
#                                                                                         UC_postovershoot_fallback_movements$dam[i],
#                                                                                         "_spilldays.png")),
#            ENT_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
# 

#### Middle Columbia ####
MCH_spill_days_data <- MCH_envir$data$winter_spill_days_data
MCW_spill_days_data <- MCW_envir$data$winter_spill_days_data
MC_spill_days_data <- bind_rows(as.data.frame(MCH_spill_days_data), as.data.frame(MCW_spill_days_data))

# Deschutes - evaluate all fallback move probs
# Deschutes is wild only
DES_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Deschutes River", movements = MC_postovershoot_fallback_movements)

# Deschutes - get covariate experiences
DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")

# Deschutes - use a loop to plot all of them
DES_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(DES_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> DES_data_to_plot
  
  if(nrow(DES_data_to_plot)>0){
    DES_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Deschutes River",
                                                                            wild_move_prob_array = DES_wild_spilldays_move_prob_array,
                                                                            
                                                                            wild_covariate_experiences = DES_wild_covariate_experiences, 
                                                                            
                                                                            movements_evaluated = MC_postovershoot_fallback_movements,
                                                                            spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                            from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                            plot_title = paste0("Deschutes - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
    
    ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("DES_compare_fallback_", 
                                                                                        MC_postovershoot_fallback_movements$dam[i],
                                                                                        "_spilldays.png")),
           DES_fallback_spilldays_plots[[i]], height = 8, width = 8)
    
  } else {
    # if it is zero, do nothing
    
  }
  

  
}

# John Day - evaluate all fallback move probs
# JDR is wild only
JDR_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "John Day River", movements = MC_postovershoot_fallback_movements)

# John Day - get covariate experiences
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")

# John Day - use a loop to plot all of them
JDR_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(JDR_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> JDR_data_to_plot
  
  if(nrow(JDR_data_to_plot)>0){
  JDR_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "John Day River",
                                                                          wild_move_prob_array = JDR_wild_spilldays_move_prob_array,
                                                                          
                                                                          wild_covariate_experiences = JDR_wild_covariate_experiences, 
                                                                          
                                                                          movements_evaluated = MC_postovershoot_fallback_movements,
                                                                          spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                          from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                          plot_title = paste0("John Day - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("JDR_compare_fallback_", 
                                                                                      MC_postovershoot_fallback_movements$dam[i],
                                                                                      "_spilldays.png")),
         JDR_fallback_spilldays_plots[[i]], height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}


# Fifteenmile - evaluate all fallback move probs
# Fifteenmile is wild only
FIF_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Fifteenmile Creek", movements = MC_postovershoot_fallback_movements)

# Fifteenmile - get covariate experiences
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")

# Fifteenmile - use a loop to plot all of them
FIF_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(FIF_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> FIF_data_to_plot
  
  if(nrow(FIF_data_to_plot)>0){
  FIF_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Fifteenmile Creek",
                                                                          wild_move_prob_array = FIF_wild_spilldays_move_prob_array,
                                                                          
                                                                          wild_covariate_experiences = FIF_wild_covariate_experiences, 
                                                                          
                                                                          movements_evaluated = MC_postovershoot_fallback_movements,
                                                                          spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                          from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                          plot_title = paste0("Fifteenmile - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("FIF_compare_fallback_", 
                                                                                      MC_postovershoot_fallback_movements$dam[i],
                                                                                      "_spilldays.png")),
         FIF_fallback_spilldays_plots[[i]], height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# Umatilla - evaluate all fallback move probs
UMA_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Umatilla River", movements = MC_postovershoot_fallback_movements)
UMA_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_MCH(origin_select= "Umatilla River", movements = MC_postovershoot_fallback_movements)

# Umatilla - get covariate experiences
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")

# Umatilla - use a loop to plot all of them
UMA_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(UMA_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> UMA_wild_data_to_plot
  filter(UMA_hatchery_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> UMA_hatchery_data_to_plot
  
  if(nrow(UMA_wild_data_to_plot)>0 | nrow(UMA_hatchery_data_to_plot)>0){
  UMA_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Umatilla River",
                                                                          wild_move_prob_array = UMA_wild_spilldays_move_prob_array,
                                                                          hatchery_move_prob_array = UMA_hatchery_spilldays_move_prob_array,
                                                                          wild_covariate_experiences = UMA_wild_covariate_experiences, 
                                                                          hatchery_covariate_experiences = UMA_hatchery_covariate_experiences,
                                                                          movements_evaluated = MC_postovershoot_fallback_movements,
                                                                          spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                          from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                          plot_title = paste0("Umatilla - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("UMA_compare_fallback_", 
                                                                                      MC_postovershoot_fallback_movements$dam[i],
                                                                                      "_spilldays.png")),
         UMA_fallback_spilldays_plots[[i]], height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# Yakima - evaluate all fallback move probs
# Yakima is wild only
YAK_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Yakima River", movements = MC_postovershoot_fallback_movements)

# Yakima - get covariate experiences
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")

# Yakima - use a loop to plot all of them
YAK_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(YAK_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> YAK_data_to_plot
  
  if(nrow(YAK_data_to_plot)>0){
  YAK_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Yakima River",
                                                                          wild_move_prob_array = YAK_wild_spilldays_move_prob_array,
                                                                          
                                                                          wild_covariate_experiences = YAK_wild_covariate_experiences, 
                                                                          
                                                                          movements_evaluated = MC_postovershoot_fallback_movements,
                                                                          spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                          from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                          plot_title = paste0("Yakima - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("YAK_compare_fallback_", 
                                                                                      MC_postovershoot_fallback_movements$dam[i],
                                                                                      "_spilldays.png")),
         YAK_fallback_spilldays_plots[[i]], height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}

# Walla Walla - evaluate all fallback move probs
WAWA_wild_spilldays_move_prob_array <- estimate_spilldays_effect_MCW(origin_select= "Walla Walla River", movements = MC_postovershoot_fallback_movements)
WAWA_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_MCH(origin_select= "Walla Walla River", movements = MC_postovershoot_fallback_movements)

# Walla Walla - get covariate experiences
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")

# Walla Walla - use a loop to plot all of them
WAWA_fallback_spilldays_plots <- vector(mode = "list", length = nrow(MC_postovershoot_fallback_movements))
for (i in 1:nrow(MC_postovershoot_fallback_movements)){
  # if a movement doesn't have any data, don't plot it
  filter(WAWA_wild_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> WAWA_wild_data_to_plot
  filter(WAWA_hatchery_covariate_experiences, state == MC_postovershoot_fallback_movements$from[i]) -> WAWA_hatchery_data_to_plot
  
  if(nrow(WAWA_wild_data_to_plot)>0 | nrow(WAWA_hatchery_data_to_plot)>0){
  WAWA_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Walla Walla River",
                                                                           wild_move_prob_array = WAWA_wild_spilldays_move_prob_array,
                                                                           hatchery_move_prob_array = WAWA_hatchery_spilldays_move_prob_array,
                                                                           wild_covariate_experiences = WAWA_wild_covariate_experiences, 
                                                                           hatchery_covariate_experiences = WAWA_hatchery_covariate_experiences,
                                                                           movements_evaluated = MC_postovershoot_fallback_movements,
                                                                           spill_predict = seq(0, max(MC_spill_days_data[,MC_postovershoot_fallback_movements$from[i]]),length = 100),
                                                                           from = MC_postovershoot_fallback_movements$from[i], to = MC_postovershoot_fallback_movements$to[i], 
                                                                           plot_title = paste0("Walla Walla - post-overshoot fallback at ", MC_postovershoot_fallback_movements$dam[i]))
  
  ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("WAWA_compare_fallback_", 
                                                                                      MC_postovershoot_fallback_movements$dam[i],
                                                                                      "_spilldays.png")),
         WAWA_fallback_spilldays_plots[[i]], height = 8, width = 8)
  } else {
    # if it is zero, do nothing
    
  }
}


# #### Snake River ####
# # only one post-overshoot fallback: Tucannon
# 
# SRH_spill_days_data <- SRH_envir$data$winter_spill_days_data
# SRW_spill_days_data <- SRW_envir$data$winter_spill_days_data
# SR_spill_days_data <- bind_rows(as.data.frame(SRH_spill_days_data), as.data.frame(SRW_spill_days_data))
# 
# # Tucannon River - evaluate all fallback move probs
# TUC_wild_spilldays_move_prob_array <- estimate_spilldays_effect_SRW(origin_select= "Tucannon River", movements = SR_postovershoot_fallback_movements)
# TUC_hatchery_spilldays_move_prob_array <- estimate_spilldays_effect_SRH(origin_select= "Tucannon River", movements = SR_postovershoot_fallback_movements)
# 
# # Tucannon River - get covariate experiences
# TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = SRW_envir, rear = "wild", origin_select = "Tucannon River")
# TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = SRH_envir, rear = "hatchery", origin_select = "Tucannon River")
# 
# # Tucannon River - use a loop to plot all of them
# TUC_fallback_spilldays_plots <- vector(mode = "list", length = nrow(SR_postovershoot_fallback_movements))
# for (i in 1:nrow(SR_postovershoot_fallback_movements)){
#   # if a movement doesn't have any data, don't plot it
#   filter(TUC_wild_covariate_experiences, state == SR_postovershoot_fallback_movements$from[i]) -> TUC_wild_data_to_plot
#   filter(TUC_hatchery_covariate_experiences, state == SR_postovershoot_fallback_movements$from[i]) -> TUC_hatchery_data_to_plot
#   
#   if(nrow(TUC_wild_data_to_plot)>0 | nrow(TUC_hatchery_data_to_plot)>0){
#   TUC_fallback_spilldays_plots[[i]] <- plot_compare_rear_spilldays_effect(origin_select = "Tucannon River",
#                                                                           wild_move_prob_array = TUC_wild_spilldays_move_prob_array,
#                                                                           hatchery_move_prob_array = TUC_hatchery_spilldays_move_prob_array,
#                                                                           wild_covariate_experiences = TUC_wild_covariate_experiences, 
#                                                                           hatchery_covariate_experiences = TUC_hatchery_covariate_experiences,
#                                                                           movements_evaluated = SR_postovershoot_fallback_movements,
#                                                                           spill_predict = seq(0, max(SR_spill_days_data[,SR_postovershoot_fallback_movements$from[i]]),length = 100),
#                                                                           from = SR_postovershoot_fallback_movements$from[i], to = SR_postovershoot_fallback_movements$to[i], 
#                                                                           plot_title = paste0("Tucannon River - post-overshoot fallback at ", SR_postovershoot_fallback_movements$dam[i]))
#   
#   ggsave(here::here("stan_actual", "output_marchspilldays", "covariate_effects", "spilldays", paste0("TUC_compare_fallback_", 
#                                                                                       SR_postovershoot_fallback_movements$dam[i],
#                                                                                       "_spilldays.png")),
#          TUC_fallback_spilldays_plots[[i]], height = 8, width = 8)
#   } else {
#     # if it is zero, do nothing
#     
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
