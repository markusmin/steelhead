# 08b data simulation visualization

# the purpose of this script is to allow us to visualize the state transitions that we are simulating in script 08.
# This will help with troubleshooting especially, and make sure that our simulations are capturing
# the behavior that we are seeing in the actual dataset.

# Load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)
library(ggthemes)

# Load data
sim_list <- readRDS(here::here("simulation", "sim_600_cov_origin_hist_list.rds"))
sim_cat_data <- read.csv(here::here("simulation", "origin_rear_600.csv"))

# Load states
sim_states = c("mainstem, mouth to BON", 
               "mainstem, BON to MCN", 
               "mainstem, MCN to ICH or PRA", 
               "mainstem, PRA to RIS",
               "mainstem, ICH to LGR",
               "Deschutes River", 
               "John Day River", 
               "Tucannon River", 
               "Yakima River",
               "loss")

# Create vector of possible transitions
possible_transitions <- c("Deschutes River - loss", "Deschutes River - mainstem, BON to MCN", "John Day River - loss", "John Day River - mainstem, BON to MCN", "mainstem, BON to MCN - Deschutes River",
                          "mainstem, BON to MCN - John Day River", "mainstem, BON to MCN - loss",  "mainstem, BON to MCN - mainstem, MCN to ICH or PRA", "mainstem, BON to MCN - mainstem, mouth to BON",      
                          "mainstem, ICH to LGR - mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR - Tucannon River",  "mainstem, MCN to ICH or PRA - loss","mainstem, MCN to ICH or PRA - mainstem, BON to MCN", 
                          "mainstem, MCN to ICH or PRA - mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA - mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA - Yakima River","mainstem, mouth to BON - loss",   
                          "mainstem, mouth to BON - mainstem, BON to MCN", "mainstem, PRA to RIS - loss", "mainstem, PRA to RIS - mainstem, MCN to ICH or PRA","Tucannon River - loss",  
                          "Tucannon River - mainstem, ICH to LGR", "Yakima River - loss", "Yakima River - mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR - loss")

data_sim_summary <- function(sim_list, run_number){
  # Extract one simulation run
  sim_array <- sim_list[[run_number]]
  
  # Dimensions: rows = states, columns = occupancy in each time step, matrix slices = individual fish
  
  # Convert these arrays into counts of the number of individual transitions
  
  transitions <- data.frame(state1 = NA, state2 = NA, transition = NA, origin = NA, rear = NA, tag_code = NA)
  # transitions <- subset(transitions, !(is.na(state1)))
  
  for (z in 1:dim(sim_array)[3]){
    # Extract the det hist for an individual fish
    sim_array[,,z] -> det_hist 
    
    # Convert this into a transition history
    
    # For each occupied column in the history (this is the sum() call)
    states_vec <- vector(length = sum(det_hist))
    state_1_vec <- vector(length = sum(det_hist)-1)
    state_2_vec <- vector(length = sum(det_hist)-1)
    for (i in 1:sum(det_hist)){
      states_vec[i] <- sim_states[which(det_hist[,i] == 1)]
      if (i < sum(det_hist)) {
        state_1_vec[i] <- states_vec[i] 
      }
      if (i > 1) {
        state_2_vec[i-1] <- states_vec[i] 
      }
    }
    
    # Combine the transition info with the rear/origin info
    transition_df <- data.frame(state1 = state_1_vec, state2 = state_2_vec, transition = paste0(state_1_vec, " - ", state_2_vec), 
                                origin = sim_cat_data$natal_origin[z], rear = sim_cat_data$rear_type[z], tag_code = sim_cat_data$tag_code[z])
    
    # If it's not the first det hist, bind them together
    transitions %>% 
      bind_rows(., transition_df) -> transitions
    
    
  }
  
  # Remove the NA row
  transitions <- subset(transitions, !(is.na(state1)))
  
  # Count the number of transitions
  as.data.frame(table(transitions$transition)) %>% 
    dplyr::rename(transition = Var1, count = Freq) -> transition_counts
  
  # Figure out which aren't represented in the data
  print(paste0("Missing transitions: ", setdiff(possible_transitions, transition_counts$transition)))
  # Ah okay - this would explain why we're not getting any estimates for 5,3 and 5,8 - no fish have the "mainstem, ICH to LGR - loss" transition
  
  # Plot them
  transition_counts_barplot <- ggplot(transition_counts, aes(x = transition, y = count)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_x_discrete(limits = rev)
  
  # Look at missing states by origin & rear type
  for (i in 1:3){
    for (j in 1:2){
      transitions %>% 
        subset(., origin == i & rear == j) -> transitions_subset
      
      as.data.frame(table(transitions_subset$transition)) %>% 
        dplyr::rename(transition = Var1, count = Freq) -> transition_counts_subset
      
      
      print(paste0("origin: ", i, ", rear: ", j))
      print(paste0("missing states: ", setdiff(possible_transitions, transition_counts_subset$transition)))
      print(paste0("< 5 obs:"))
      print(paste0(subset(transition_counts_subset, count < 5)$transition, ": ", subset(transition_counts_subset, count < 5)$count))
      
    }
  }
  
  return(list(transitions, transition_counts_barplot))
  
}


##### Script to investigate counts of a single transition #####
# Start with a single from state

# create a df with names + indices of states
state_indices <- data.frame(state2 = sim_states, state2_index = seq(1,10))

from_state_counts <- function(transitions, from_state){
  state_name <- sim_states[from_state]
  transitions %>% 
    subset(., state1 == state_name) -> from_subset
  
  from_subset %>% 
    group_by(origin, rear) %>% 
    count(state2) %>% left_join(., state_indices, by = "state2") %>% 
    # mutate(indices = paste0("[",from_state,",", which(sim_states == state2),"]"))
    mutate(indices = paste0("[",from_state,",", state2_index,"]")) %>% 
    dplyr::select(-c(state2_index)) %>% 
    arrange(indices)
}


##### ORIGIN ONLY, 600 fish, run 1 #####

origin600_1 <- data_sim_summary(sim_list = sim_list, run_number = 1)
origin600_1[[1]]

from_state_counts(transitions = origin600_1[[1]], from_state = 5)
from_state_counts(transitions = origin600_1[[1]], from_state = 6) # Every rear + origin represented in each movement
from_state_counts(transitions = origin600_1[[1]], from_state = 7) # Origin 3 never is lost (??) - TUC river fish always return (but only 6 of them) # JDR
from_state_counts(transitions = origin600_1[[1]], from_state = 8) # Origins 1 and 2 are never lost # TUC
from_state_counts(transitions = origin600_1[[1]], from_state = 9) # Origin 3 is never lost # YAK







