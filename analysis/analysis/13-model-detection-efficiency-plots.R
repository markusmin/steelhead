# 13-model-detection-efficiency-plots

# This script takes the output from the stan model runs in 05-stan-runs and
# plots how detection efficiency varies through time

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
# SRH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery", "chain1_SRH_reparam_v2_fit.rds"))
# temporary - chain1 didn't finish so read chain 4 twice
SRH_chain1 <- readRDS(here::here("stan_actual", "reparameterization_v2", "snake_river_hatchery", "chain4_SRH_reparam_v2_fit.rds"))
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
extract_DE_parameters <- function(fit, fit_summary){
  # extract b0 as an array
  parameters <- fit_summary$variable
  parameters[grepl("_alpha|_beta" , parameters)] -> param_subset
  
  # add in fifteenmile_beta and imnaha_beta
  DE_params <- c(param_subset[1:24], "fifteenmile_beta", param_subset[25], "imnaha_beta", param_subset[26:length(param_subset)])
  
  
  # arrange all parameter values together in a df
  DE_param_matrix <- matrix(data = 0, nrow = length(DE_params),
                                         length(as.matrix(fit[,,1])))
  
  rownames(DE_param_matrix) <- DE_params
  
  
  
  # extract the draws for each and store in a matrix
  for(i in 1:nrow(DE_param_matrix)){
    if (DE_params[i] %in% c("fifteenmile_beta", "imnaha_beta")){
      DE_param_matrix[i, ] <- 0
    } else {
      DE_param_matrix[i, ] <- as.matrix(fit[,,DE_params[i]])
    }
    
  }
  
  return(DE_param_matrix)
}

# function to reformat param matrix and calculate DE using those params + discharge
calculate_DE <- function(DE_param_matrix, tributary, tributary_design_matrices_array){
  
  # get the index for the tributary state (needs to be mouth since that's where we're correcting for DE)
  tributary_state <- intersect(grep(tributary, model_states, ignore.case = TRUE), grep("Mouth", model_states, ignore.case = TRUE))
  
  # use the state indexing to get the right design matrix
  trib_design_matrix <- tributary_design_matrices_array[,,tributary_state]
  
  # create a matrix to store DE per year, per iter
  niter <- 4000
  DE_matrix <- matrix(nrow = nrow(trib_design_matrix), ncol = niter)
  
  # for each run year, get a confidence interval around detection efficiency by using the different draws
  for (i in 1:nrow(trib_design_matrix)){
    # if there is no intercept term, that means there is no DE correction for that year - so skip and leave as NA
    if(sum(trib_design_matrix[i,1:21]) == 0){
      
      
    } else {
      for (j in 1:niter){
        eta <- sum(trib_design_matrix[i,] * DE_param_matrix[,j])
        DE_matrix[i,j] <- exp(eta)/(1 + exp(eta))
      }
      
    }
    
    
  }
  
  return(DE_matrix)
  
}

# function to plot DE plus credible intervals
plot_DE_rear <- function(DE_by_year_wild = NULL, DE_by_year_hatchery = NULL, plot_title){
  rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
  
  if(is.null(DE_by_year_hatchery)){
    niter <- 4000
    colnames(DE_by_year_wild) <- paste0("iter", 1:niter)
    
    DE_by_year_wild %>% 
      as.data.frame() %>% 
      rownames_to_column("year") %>% 
      mutate(year = as.numeric(year)) %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
      group_by(year) %>% 
      summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> DE_by_year_wild_long
    
    DE_by_year_wild_long$year_actual <- 2004:2021
    
    DE_by_year_wild_long -> DE_by_year_rear_long

  } else if (is.null(DE_by_year_wild)){
    niter <- 4000
    colnames(DE_by_year_hatchery) <- paste0("iter", 1:niter)
    
    DE_by_year_hatchery %>% 
      as.data.frame() %>% 
      rownames_to_column("year") %>% 
      mutate(year = as.numeric(year)) %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
      group_by(year) %>% 
      summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> DE_by_year_hatchery_long
    
    DE_by_year_hatchery_long$year_actual <- 2004:2021
    
    DE_by_year_hatchery_long -> DE_by_year_rear_long
    
  } else {
    niter <- 4000
    colnames(DE_by_year_wild) <- paste0("iter", 1:niter)
    colnames(DE_by_year_hatchery) <- paste0("iter", 1:niter)
    
    DE_by_year_wild %>% 
      as.data.frame() %>% 
      rownames_to_column("year") %>% 
      mutate(year = as.numeric(year)) %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
      group_by(year) %>% 
      summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "wild") -> DE_by_year_wild_long
    
    DE_by_year_wild_long$year_actual <- 2004:2021
    
    DE_by_year_hatchery %>% 
      as.data.frame() %>% 
      rownames_to_column("year") %>% 
      mutate(year = as.numeric(year)) %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
      group_by(year) %>% 
      summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> DE_by_year_hatchery_long
    
    DE_by_year_hatchery_long$year_actual <- 2004:2021
    
    DE_by_year_wild_long %>% 
      bind_rows(., DE_by_year_hatchery_long) -> DE_by_year_rear_long
    
  }
  
  
  DE_by_year_plot <- ggplot(DE_by_year_rear_long, aes(x = year_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`,
                                                      color = rear, fill = rear)) +
    geom_line() +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    geom_ribbon(alpha = 0.2, color = NA) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(2003,2022), expand = c(0,0)) +
    xlab("Year") +
    ylab("Detection probability") +
    ggtitle(plot_title)
  
  return(DE_by_year_plot)
}



### extract all params ###
UCW_DE_param_matrix <- extract_DE_parameters(fit = UCW_fit, fit_summary = UCW_fit_summary)
UCH_DE_param_matrix <- extract_DE_parameters(fit = UCH_fit, fit_summary = UCH_fit_summary)
MCW_DE_param_matrix <- extract_DE_parameters(fit = MCW_fit, fit_summary = MCW_fit_summary)
MCH_DE_param_matrix <- extract_DE_parameters(fit = MCH_fit, fit_summary = MCH_fit_summary)
SRW_DE_param_matrix <- extract_DE_parameters(fit = SRW_fit, fit_summary = SRW_fit_summary)
SRH_DE_param_matrix <- extract_DE_parameters(fit = SRH_fit, fit_summary = SRH_fit_summary)


### calculate detection efficiency by tributary ###

# first some quick sanity checks to make sure that our data is formatted properly
dim(MCW_envir$data$tributary_design_matrices_array)
# run years, parameters, states
MCW_envir$data$tributary_design_matrices_array[,,34] # this should be discharge per year for Asotin Creek


### Middle Columbia

deschutes_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
deschutes_DE_plot <- plot_DE_rear(DE_by_year_wild = deschutes_DE_by_year_wild, plot_title = "Deschutes River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "deschutes_DE_plot.png"), deschutes_DE_plot, height = 8, width = 8)

john_day_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "john", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
john_day_DE_plot <- plot_DE_rear(DE_by_year_wild = john_day_DE_by_year_wild, plot_title = "John Day River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "john_day_DE_plot.png"), john_day_DE_plot, height = 8, width = 8)

fifteenmile_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "fifteenmile", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
fifteenmile_DE_plot <- plot_DE_rear(DE_by_year_wild = fifteenmile_DE_by_year_wild, plot_title = "Fifteenmile Creek Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "fifteenmile_DE_plot.png"), fifteenmile_DE_plot, height = 8, width = 8)

umatilla_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "umatilla", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
umatilla_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "umatilla", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array)
umatilla_DE_plot <- plot_DE_rear(DE_by_year_wild = umatilla_DE_by_year_wild, 
                            DE_by_year_hatchery = umatilla_DE_by_year_hatchery, 
                            plot_title = "Umatilla River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "umatilla_DE_plot.png"), umatilla_DE_plot, height = 8, width = 8)

yakima_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "yakima", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
yakima_DE_plot <- plot_DE_rear(DE_by_year_wild = yakima_DE_by_year_wild,
                          plot_title = "Yakima River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "yakima_DE_plot.png"), yakima_DE_plot, height = 8, width = 8)

walla_walla_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "walla", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
walla_walla_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "walla", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array)
walla_walla_DE_plot <- plot_DE_rear(DE_by_year_wild = walla_walla_DE_by_year_wild,
                                    DE_by_year_hatchery = walla_walla_DE_by_year_hatchery,
                                    plot_title = "Walla Walla River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "walla_walla_DE_plot.png"), walla_walla_DE_plot, height = 8, width = 8)

### Upper Columbia
wenatchee_DE_by_year_wild <- calculate_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "wenatchee", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array)
wenatchee_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "wenatchee", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array)
wenatchee_DE_plot <- plot_DE_rear(DE_by_year_wild = wenatchee_DE_by_year_wild, 
                                  DE_by_year_hatchery = wenatchee_DE_by_year_hatchery,
                                  plot_title = "Wenatchee River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "wenatchee_DE_plot.png"), wenatchee_DE_plot, height = 8, width = 8)

entiat_DE_by_year_wild <- calculate_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "entiat", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array)
entiat_DE_plot <- plot_DE_rear(DE_by_year_wild = entiat_DE_by_year_wild, plot_title = "Entiat River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "entiat_DE_plot.png"), entiat_DE_plot, height = 8, width = 8)

okanogan_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "okanogan", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array)
okanogan_DE_plot <- plot_DE_rear(DE_by_year_hatchery = okanogan_DE_by_year_hatchery, plot_title = "Okanogan River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "okanogan_DE_plot.png"), okanogan_DE_plot, height = 8, width = 8)

methow_DE_by_year_wild <- calculate_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "methow", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array)
methow_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "methow", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array)
methow_DE_plot <- plot_DE_rear(DE_by_year_wild = methow_DE_by_year_wild,
                               DE_by_year_hatchery = methow_DE_by_year_hatchery,
                               plot_title = "Methow River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "methow_DE_plot.png"), methow_DE_plot, height = 8, width = 8)

### Snake River
asotin_DE_by_year_wild <- calculate_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "asotin", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array)
asotin_DE_plot <- plot_DE_rear(DE_by_year_wild = asotin_DE_by_year_wild, plot_title = "Asotin Creek Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "asotin_DE_plot.png"), asotin_DE_plot, height = 8, width = 8)

imnaha_DE_by_year_wild <- calculate_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "imnaha", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array)
imnaha_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = SRH_DE_param_matrix, tributary = "imnaha", tributary_design_matrices_array = SRH_envir$data$tributary_design_matrices_array)
imnaha_DE_plot <- plot_DE_rear(DE_by_year_wild = imnaha_DE_by_year_wild,
                               DE_by_year_hatchery = imnaha_DE_by_year_hatchery,
                               plot_title = "Imnaha River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "imnaha_DE_plot.png"), imnaha_DE_plot, height = 8, width = 8)

tucannon_DE_by_year_wild <- calculate_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "tucannon", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array)
tucannon_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = SRH_DE_param_matrix, tributary = "tucannon", tributary_design_matrices_array = SRH_envir$data$tributary_design_matrices_array)
tucannon_DE_plot <- plot_DE_rear(DE_by_year_wild = tucannon_DE_by_year_wild,
                                 DE_by_year_hatchery = tucannon_DE_by_year_hatchery,
                                 plot_title = "Tucannon River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "tucannon_DE_plot.png"), tucannon_DE_plot, height = 8, width = 8)





