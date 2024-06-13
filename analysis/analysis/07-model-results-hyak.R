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

# Posterior comparisons: Snake River

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






