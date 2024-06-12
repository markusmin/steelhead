# 06-model-diagnostics

# This script takes the output from the stan model runs in 05-stan-runs and
# runs diagnostic checks

#### Load libraries, state information ####
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(tidyverse)
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
UCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization", "upper_columbia_wild","UCW_transition_counts.csv"))
UCW_movements <- paste0(UCW_transition_counts$from, "_", UCW_transition_counts$to)
UCW_movements <- UCW_movements[!(grepl("NA", UCW_movements))]

UCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization", "upper_columbia_hatchery","UCH_transition_counts.csv"))
UCH_movements <- paste0(UCH_transition_counts$from, "_", UCH_transition_counts$to)
UCH_movements <- UCH_movements[!(grepl("NA", UCH_movements))]

MCW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization", "middle_columbia_wild","MCW_transition_counts.csv"))
MCW_movements <- paste0(MCW_transition_counts$from, "_", MCW_transition_counts$to)
MCW_movements <- MCW_movements[!(grepl("NA", MCW_movements))]

MCH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization", "middle_columbia_hatchery","MCH_transition_counts.csv"))
MCH_movements <- paste0(MCH_transition_counts$from, "_", MCH_transition_counts$to)
MCH_movements <- MCH_movements[!(grepl("NA", MCH_movements))]

SRW_transition_counts <- read.csv(here::here("stan_actual", "reparameterization", "snake_river_wild","SRW_transition_counts.csv"))
SRW_movements <- paste0(SRW_transition_counts$from, "_", SRW_transition_counts$to)
SRW_movements <- SRW_movements[!(grepl("NA", SRW_movements))]

# SRH_transition_counts <- read.csv(here::here("stan_actual", "reparameterization", "snake_river_hatchery","SRH_transition_counts.csv"))
# SRH_movements <- paste0(SRH_transition_counts$from, "_", SRH_transition_counts$to)
# SRH_movements <- SRH_movements[!(grepl("NA", SRH_movements))]


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
UCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization", "upper_columbia_wild", "chain1_UCW_reparam_fit.rds"))
UCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization", "upper_columbia_wild", "chain2_UCW_reparam_fit.rds"))
UCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization", "upper_columbia_wild", "chain3_UCW_reparam_fit.rds"))
UCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization", "upper_columbia_wild", "chain4_UCW_reparam_fit.rds"))

# bind chains together
UCW_fit_raw <- bind4chains(UCW_chain1, UCW_chain2, UCW_chain3, UCW_chain4)
# thin2
thin_draws(UCW_fit_raw, thin = 2) -> UCW_fit
# summarise
UCW_fit_summary <- summarise_draws(UCW_fit)

## Upper Columbia, Hatchery
UCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization", "upper_columbia_hatchery", "chain1_UCH_reparam_fit.rds"))
UCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization", "upper_columbia_hatchery", "chain2_UCH_reparam_fit.rds"))
UCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization", "upper_columbia_hatchery", "chain3_UCH_reparam_fit.rds"))
UCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization", "upper_columbia_hatchery", "chain4_UCH_reparam_fit.rds"))

# bind chains together
UCH_fit_raw <- bind4chains(UCH_chain1, UCH_chain2, UCH_chain3, UCH_chain4)
# thin2
thin_draws(UCH_fit_raw, thin = 2) -> UCH_fit
# summarise
UCH_fit_summary <- summarise_draws(UCH_fit)

## Middle Columbia, Wild
MCW_chain1 <- readRDS(here::here("stan_actual", "reparameterization", "middle_columbia_wild", "chain1_MCW_reparam_fit.rds"))
MCW_chain2 <- readRDS(here::here("stan_actual", "reparameterization", "middle_columbia_wild", "chain2_MCW_reparam_fit.rds"))
MCW_chain3 <- readRDS(here::here("stan_actual", "reparameterization", "middle_columbia_wild", "chain3_MCW_reparam_fit.rds"))
MCW_chain4 <- readRDS(here::here("stan_actual", "reparameterization", "middle_columbia_wild", "chain4_MCW_reparam_fit.rds"))

# bind chains together
MCW_fit_raw <- bind4chains(MCW_chain1, MCW_chain2, MCW_chain3, MCW_chain4)
# thin2
thin_draws(MCW_fit_raw, thin = 2) -> MCW_fit
# summarise
MCW_fit_summary <- summarise_draws(MCW_fit)

## Middle Columbia, Hatchery
MCH_chain1 <- readRDS(here::here("stan_actual", "reparameterization", "middle_columbia_hatchery", "chain1_MCH_reparam_fit.rds"))
MCH_chain2 <- readRDS(here::here("stan_actual", "reparameterization", "middle_columbia_hatchery", "chain2_MCH_reparam_fit.rds"))
MCH_chain3 <- readRDS(here::here("stan_actual", "reparameterization", "middle_columbia_hatchery", "chain3_MCH_reparam_fit.rds"))
MCH_chain4 <- readRDS(here::here("stan_actual", "reparameterization", "middle_columbia_hatchery", "chain4_MCH_reparam_fit.rds"))

# bind chains together
MCH_fit_raw <- bind4chains(MCH_chain1, MCH_chain2, MCH_chain3, MCH_chain4)
# thin2
thin_draws(MCH_fit_raw, thin = 2) -> MCH_fit
# summarise
MCH_fit_summary <- summarise_draws(MCH_fit)

## Snake River, Wild
SRW_chain1 <- readRDS(here::here("stan_actual", "reparameterization", "snake_river_wild", "chain1_SRW_reparam_fit.rds"))
SRW_chain2 <- readRDS(here::here("stan_actual", "reparameterization", "snake_river_wild", "chain2_SRW_reparam_fit.rds"))
SRW_chain3 <- readRDS(here::here("stan_actual", "reparameterization", "snake_river_wild", "chain3_SRW_reparam_fit.rds"))
SRW_chain4 <- readRDS(here::here("stan_actual", "reparameterization", "snake_river_wild", "chain4_SRW_reparam_fit.rds"))

# bind chains together
SRW_fit_raw <- bind4chains(SRW_chain1, SRW_chain2, SRW_chain3, SRW_chain4)
# thin2
thin_draws(SRW_fit_raw, thin = 2) -> SRW_fit
# summarise
SRW_fit_summary <- summarise_draws(SRW_fit)

# ## Snake River, Hatchery
# SRH_chain1 <- readRDS(here::here("stan_actual", "reparameterization", "snake_river_hatchery", "chain1_SRH_reparam_fit.rds"))
# SRH_chain2 <- readRDS(here::here("stan_actual", "reparameterization", "snake_river_hatchery", "chain2_SRH_reparam_fit.rds"))
# SRH_chain3 <- readRDS(here::here("stan_actual", "reparameterization", "snake_river_hatchery", "chain3_SRH_reparam_fit.rds"))
# SRH_chain4 <- readRDS(here::here("stan_actual", "reparameterization", "snake_river_hatchery", "chain4_SRH_reparam_fit.rds"))
# 
# # bind chains together
# SRH_fit_raw <- bind4chains(SRH_chain1, SRH_chain2, SRH_chain3, SRH_chain4)
# # thin2
# thin_draws(SRH_fit_raw, thin = 2) -> SRH_fit
# # summarise
# SRH_fit_summary <- summarise_draws(SRH_fit)

##### Run diagnostic summaries: rhat and ess #####

# Create six-panel diagnostic plots for rhat, ess_bulk, and ess_tail

create_rhat_hist <- function(fit_summary, population){
  rhat_plot <- ggplot(fit_summary, aes(x = rhat)) +
    geom_histogram() + 
    annotate("text", label = "Range:", x = max(fit_summary$rhat, na.rm = T), y = 600, hjust = 1) +
    annotate("text", label = round(range(fit_summary$rhat, na.rm = T),3)[1], x = max(fit_summary$rhat, na.rm = T), y = 550,hjust = 1) +
    annotate("text", label = paste0(" - ", round(range(fit_summary$rhat, na.rm = T),3)[2]), x = max(fit_summary$rhat, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population)
  
  return(rhat_plot)
}

create_ess_bulk_hist <- function(fit_summary, population){
  ess_bulk_plot <- ggplot(fit_summary, aes(x = ess_bulk)) +
    geom_histogram() + 
    annotate("text", label = "Range:", x = max(fit_summary$ess_bulk, na.rm = T), y = 600, hjust = 1) +
    annotate("text", label = round(range(fit_summary$ess_bulk, na.rm = T),3)[1], x = max(fit_summary$ess_bulk, na.rm = T), y = 550,hjust = 1) +
    annotate("text", label = paste0(" - ", round(range(fit_summary$ess_bulk, na.rm = T),3)[2]), x = max(fit_summary$ess_bulk, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population)
  
  return(ess_bulk_plot)
}

create_ess_tail_hist <- function(fit_summary, population){
  ess_tail_plot <- ggplot(fit_summary, aes(x = ess_tail)) +
    geom_histogram() + 
    annotate("text", label = "Range:", x = max(fit_summary$ess_tail, na.rm = T), y = 600, hjust = 1) +
    annotate("text", label = round(range(fit_summary$ess_tail, na.rm = T),3)[1], x = max(fit_summary$ess_tail, na.rm = T), y = 550,hjust = 1) +
    annotate("text", label = paste0(" - ", round(range(fit_summary$ess_tail, na.rm = T),3)[2]), x = max(fit_summary$ess_tail, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population)
  
  return(ess_tail_plot)
}

# Run functions for all populations
UCW_rhat_plot <- create_rhat_hist(fit_summary = UCW_fit_summary, population = "UCW")
UCH_rhat_plot <- create_rhat_hist(fit_summary = UCH_fit_summary, population = "UCH")
MCW_rhat_plot <- create_rhat_hist(fit_summary = MCW_fit_summary, population = "MCW")
MCH_rhat_plot <- create_rhat_hist(fit_summary = MCH_fit_summary, population = "MCH")
SRW_rhat_plot <- create_rhat_hist(fit_summary = SRW_fit_summary, population = "SRW")
# SRH_rhat_plot <- create_rhat_hist(fit_summary = SRH_fit_summary, population = "SRH")

UCW_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = UCW_fit_summary, population = "UCW")
UCH_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = UCH_fit_summary, population = "UCH")
MCW_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = MCW_fit_summary, population = "MCW")
MCH_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = MCH_fit_summary, population = "MCH")
SRW_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = SRW_fit_summary, population = "SRW")
# SRH_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = SRH_fit_summary, population = "SRH")

UCW_ess_tail_plot <- create_ess_tail_hist(fit_summary = UCW_fit_summary, population = "UCW")
UCH_ess_tail_plot <- create_ess_tail_hist(fit_summary = UCH_fit_summary, population = "UCH")
MCW_ess_tail_plot <- create_ess_tail_hist(fit_summary = MCW_fit_summary, population = "MCW")
MCH_ess_tail_plot <- create_ess_tail_hist(fit_summary = MCH_fit_summary, population = "MCH")
SRW_ess_tail_plot <- create_ess_tail_hist(fit_summary = SRW_fit_summary, population = "SRW")
# SRH_ess_tail_plot <- create_ess_tail_hist(fit_summary = SRH_fit_summary, population = "SRH")

# arrange them
rhat_comb_plots <- ggarrange(plotlist = list(UCW_rhat_plot, UCH_rhat_plot, MCW_rhat_plot, MCH_rhat_plot, SRW_rhat_plot),
          ncol = 3, nrow = 2)

ess_tail_comb_plots <- ggarrange(plotlist = list(UCW_ess_bulk_plot, UCH_ess_bulk_plot, MCW_ess_bulk_plot, MCH_ess_bulk_plot, SRW_ess_bulk_plot),
                             ncol = 3, nrow = 2)

ess_bulk_comb_plots <- ggarrange(plotlist = list(UCW_ess_tail_plot, UCH_ess_tail_plot, MCW_ess_tail_plot, MCH_ess_tail_plot, SRW_ess_tail_plot),
                             ncol = 3, nrow = 2)

#### Run diagnostics: Parameter correlations ####

# Pairs plot for movements out of one state
one_state_param_pairs <- function(fit, fit_summary, state){
  parameters <- unique(fit_summary$variable)
  one_state_params <- parameters[grep(paste0("_", state, "_"), parameters)]
  # remove NDE parameters
  one_state_params <- one_state_params[!(grepl("NDE", one_state_params))]
  # remove year parameters
  # one_state_params <- one_state_params[!(grepl("\\[", one_state_params))]
  one_state_params <- one_state_params[!(grepl("year", one_state_params))]
  # generate the plot
  pairs_plot <- mcmc_pairs(fit, pars = one_state_params)
  return(pairs_plot)
}

# Pairs plot for only one specific movement
one_movement_param_pairs <- function(fit, fit_summary, movement){
  parameters <- unique(fit_summary$variable)
  
  if (length(grep(paste0(movement), parameters)) > 1) {
    one_movement_params <- parameters[grep(paste0(movement), parameters)]
    # remove NDE parameters
    one_movement_params <- one_movement_params[!(grepl("NDE", one_movement_params))]
    # remove year parameters
    # one_movement_params <- one_movement_params[!(grepl("\\[", one_movement_params))]
    one_movement_params <- one_movement_params[!(grepl("year", one_movement_params))]
    # generate the plot
    pairs_plot <- mcmc_pairs(fit, pars = one_movement_params,
                             off_diag_args = list(size = 0.5))
  } else {
    print(paste0("One or no parameters for movement ", movement))
  }
  
  return(pairs_plot)
}


# run pairs plot for each specific transition
UCW_movements_plotlist <- vector(mode = "list", length = length(UCW_movements))
for (i in 1:length(UCW_movements)){
  pairs_plot <- one_movement_param_pairs(fit = UCW_fit,
                                      fit_summary = UCW_fit_summary,
                                      movement = UCW_movements[i])
  
  UCW_movements_plotlist[[i]] <- pairs_plot
}

pdf(file = "stan_actual/output/UCW_onemovement_pairs_plots.pdf")
UCW_movements_plotlist
dev.off()


# run for all states, all populations

states <- 1:9
populations <- c("UCW", "UCH", "MCW", "MCH", "SRW")

UCW_plotlist <- vector(mode = "list", length = length(states))
UCH_plotlist <- vector(mode = "list", length = length(states))
MCW_plotlist <- vector(mode = "list", length = length(states))
MCH_plotlist <- vector(mode = "list", length = length(states))
SRW_plotlist <- vector(mode = "list", length = length(states))
# SRH_plotlist <- vector(mode = "list", length = length(states))

for (i in 1:length(populations)){
  for (j in 1:length(states)){
    pairs_plot <- one_state_param_pairs(fit = eval(parse(text = paste0(populations[i], "_fit"))),
                                        fit_summary = eval(parse(text = paste0(populations[i], "_fit_summary"))),
                                        state = states[j])
  
    eval(parse(text = paste0(populations[i], "_plotlist[[j]] <- pairs_plot")))
  }
}

pdf(file = "stan_actual/output/UCW_pairs_plots.pdf")
UCW_plotlist
dev.off()






