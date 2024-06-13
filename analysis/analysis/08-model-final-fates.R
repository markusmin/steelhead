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


#### Final fates function ####

# This function takes a model fit object, a number of simulated fish,
# the identities of those fish (origin, rear), the conditions throughout the
# basin (spill, temperature, year)
# and then simulates the final fates of those fish
# By default, the conditions (spill and temperature) will be taken from the data
# itself, but note that new values of these conditions could be simulated in order
# to test different scenarios

# Function 1: Take parameters and derive the movement probability


# testing block
UCW_parameters <- UCW_fit_summary$variable
UCW_parameters[grepl("b0|btemp1_|btemp0_|bspillwindow|bwinterspill", UCW_parameters)] -> UCW_fixed_effects
UCW_fixed_effects[!(grepl("vector", UCW_fixed_effects))] -> UCW_fixed_effects

# function to take a fixed effect and store all of them in an array
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
    param_array[param_subset_indices[i, "from"], param_subset_indices[i, "to"], ] <- as.matrix(UCW_fit[,,param_subset_indices[i, "parameter"]])
  }
  
  return(param_array)
}

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

# reformat covariate values for simulation
UCW_envir$data$temperature_data
UCW_envir$data$spill_window_data
UCW_envir$data$winter_spill_days_data
UCW_envir$data$transition_dates
UCW_envir$data$y
UCW_envir$data$movements

# In order to simulate covariate values, we are going to determine the dates
# where the fish were in each state, and then sample from those dates to 
# get a representative value for temperature/spill when fish are in those states
UCW_states_dates <- data.frame(state = as.vector(UCW_envir$data$y),
                           date = as.vector(UCW_envir$data$transition_dates))



# Final states simulation
# currently we are not including a random effect of year, but we could easily
# amend the code to simulate final fates in a particular year
# origin1/origin2/origin3 are indicator variables, so set the origin you want to be 1
final_fates_simulation <- function(nsim, b0_array,btemp0_array, btemp1_array,
                                   bspillwindow_array, bwinterspill_array, btemp0xorigin1_array,
                                   btemp1xorigin1_array, btemp0xorigin2_array, btemp1xorigin2_array,
                                   btemp0xorigin3_array, btemp1xorigin3_array,
                                   borigin1_array,borigin2_array, borigin3_array, 
                                   start_state = 2, states_dates,
                                   origin1 = 0, origin2 = 0, origin3 = 0,
                                   temp_data, spillwindow_data, winterspill_data,
                                   year){
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
                                                         bwinterspill_array_UCW[i,possible_movements[j],iter]*winterspill[i] +
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
                  bwinterspill_array_UCW[i,possible_movements,iter]*winterspill[i] +
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
  
  outputs <- list(state_matrix, final_fates)
  # return(final_fates)
  return(outputs)
  
}

final_fates_wenatchee <- data.frame(state = model_states[1:(length(model_states)-1)])

  
# Now, run final fates simulation with 100 different random parameter draws, 100 fish each time
niter <- 100
for(i in 1:niter) {
  sim <- final_fates_simulation(nsim = 100, b0_array = b0_array_UCW,btemp0_array = btemp0_array_UCW, btemp1_array = btemp1_array_UCW,
                                bspillwindow_array = bspillwindow_array_UCW, bwinterspill_array = bwinterspill_array_UCW, btemp0xorigin1_array = btemp0xorigin1_array_UCW,
                                btemp1xorigin1_array = btemp1xorigin1_array_UCW, btemp0xorigin2_array = btemp0xorigin2_array_UCW, btemp1xorigin2_array = btemp1xorigin2_array_UCW,
                                btemp0xorigin3_array = btemp0xorigin3_array_UCW, btemp1xorigin3_array = btemp1xorigin3_array_UCW,
                                borigin1_array = borigin1_array_UCW,borigin2_array = borigin2_array_UCW, borigin3_array = borigin3_array_UCW, 
                                start_state = 2, states_dates = UCW_states_dates,
                                origin1 = 1, origin2 = 0, origin3 = 0,
                                temp_data = UCW_envir$data$temperature_data, spillwindow_data = UCW_envir$data$spill_window_data, winterspill_data = UCW_envir$data$winter_spill_days_data)
  final_fates_wenatchee %>% 
    bind_cols(sim[[2]]) -> final_fates_wenatchee
  
}

rownames(final_fates_wenatchee) <- NULL
column_to_rownames(final_fates_wenatchee, "state") -> final_fates_wenatchee
colnames(final_fates_wenatchee) <- paste0("iter", 1:100)

final_fates_wenatchee$total <- rowSums(final_fates_wenatchee)
rownames_to_column(final_fates_wenatchee, "state") -> final_fates_wenatchee

ggplot(final_fates_wenatchee, aes(x = state, y = total)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90))








