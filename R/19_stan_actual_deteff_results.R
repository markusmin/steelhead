# 19_stan_actual_deteff_results


# This script will load and investigate the outputs from our stan models.
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(tidyverse)
library(here)
library(ggpubr)


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



# Read in 200iter runs
snake_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "snake", "200iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))
middle_columbia_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "middle_columbia", "200iter_parallel_middle_columbia_stan_actual_int_origin_stan_fit.rds"))
upper_columbia_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "upper_columbia", "200iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))

# Get each of the summary objects
snake_summary <- snake_fit$summary()
middle_columbia_summary <- middle_columbia_fit$summary()
upper_columbia_summary <- upper_columbia_fit$summary()

snake_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat))  %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable)))  %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Tucannon_River",
                         ifelse(origin == "origin2", "Asotin_Creek",
                                ifelse(origin == "origin3",  "Clearwater_River",
                                       ifelse(origin == "origin4", "Salmon_River",
                                              ifelse(origin == "origin5", "Grande_Ronde",
                                                     ifelse(origin == "origin6", "Imnaha_River", "error"))))))) -> snake_derived_probabilities 

middle_columbia_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat))  %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable)))  %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Deschutes_River",
                         ifelse(origin == "origin2", "Fifteenmile_Creek",
                                ifelse(origin == "origin3",  "John_Day_River",
                                       ifelse(origin == "origin4", "Umatilla_River",
                                              ifelse(origin == "origin5", "Yakima_River",
                                                     ifelse(origin == "origin6", "Walla_Walla_River", "error"))))))) -> middle_columbia_derived_probabilities 

upper_columbia_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat))  %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable)))  %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Wenatchee_River",
                         ifelse(origin == "origin2", "Entiat_River",
                                ifelse(origin == "origin3",  "Okanogan_River",
                                       ifelse(origin == "origin4", "Methow_River", "error"))))) -> upper_columbia_derived_probabilities 

# Split this up into DE and NDE probabilities
snake_derived_probabilities %>% 
  filter(.,grepl("_NDE", variable)) -> snake_NDE_derived_probabilities
snake_derived_probabilities %>% 
  filter(.,grepl("_DE", variable)) -> snake_DE_derived_probabilities

middle_columbia_derived_probabilities %>% 
  filter(.,grepl("_NDE", variable)) -> middle_columbia_NDE_derived_probabilities
middle_columbia_derived_probabilities %>% 
  filter(.,grepl("_DE", variable)) -> middle_columbia_DE_derived_probabilities

upper_columbia_derived_probabilities %>% 
  filter(.,grepl("_NDE", variable)) -> upper_columbia_NDE_derived_probabilities
upper_columbia_derived_probabilities %>% 
  filter(.,grepl("_DE", variable)) -> upper_columbia_DE_derived_probabilities


##### SNAKE - Export tables of parameters for report #####

# DE CORRECTED
snake_DE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> snake_DE_prob_table

# subset the probabilities shared across the DPS
snake_DE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                     "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                     "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                     "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                     "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> snake_DE_DPS_probs

snake_DE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                       "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                       "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                       "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                       "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> snake_DE_origin_probs

# export both tables
write.csv(snake_DE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "snake_DE_DPS_probs.csv"))
write.csv(snake_DE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "snake_DE_origin_probs.csv"))



# NOT DE CORRECTED
snake_NDE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> snake_NDE_prob_table

# subset the probabilities shared across the DPS
snake_NDE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                     "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                     "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                     "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                     "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> snake_NDE_DPS_probs

snake_NDE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                       "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                       "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                       "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                       "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> snake_NDE_origin_probs

# export both tables
write.csv(snake_NDE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "snake_NDE_DPS_probs.csv"))
write.csv(snake_NDE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "snake_NDE_origin_probs.csv"))

##### MIDDLE COLUMBIA - Export tables of parameters for report #####

# DE CORRECTED
middle_columbia_DE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> middle_columbia_DE_prob_table

# subset the probabilities shared across the DPS
middle_columbia_DE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                     "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                     "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                     "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                     "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> middle_columbia_DE_DPS_probs

middle_columbia_DE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                       "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                       "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                       "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                       "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> middle_columbia_DE_origin_probs

# export both tables
write.csv(middle_columbia_DE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "middle_columbia_DE_DPS_probs.csv"))
write.csv(middle_columbia_DE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "middle_columbia_DE_origin_probs.csv"))



# NOT DE CORRECTED
middle_columbia_NDE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> middle_columbia_NDE_prob_table

# subset the probabilities shared across the DPS
middle_columbia_NDE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                     "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                     "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                     "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                     "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> middle_columbia_NDE_DPS_probs

middle_columbia_NDE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                       "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                       "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                       "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                       "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> middle_columbia_NDE_origin_probs

# export both tables
write.csv(middle_columbia_NDE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "middle_columbia_NDE_DPS_probs.csv"))
write.csv(middle_columbia_NDE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "middle_columbia_NDE_origin_probs.csv"))

##### UPPER COLUMBIA - Export tables of parameters for report #####

# DE CORRECTED
upper_columbia_DE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> upper_columbia_DE_prob_table

# subset the probabilities shared across the DPS
upper_columbia_DE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                     "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                     "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                     "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                     "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> upper_columbia_DE_DPS_probs

upper_columbia_DE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                       "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                       "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                       "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                       "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> upper_columbia_DE_origin_probs

# export both tables
write.csv(upper_columbia_DE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "upper_columbia_DE_DPS_probs.csv"))
write.csv(upper_columbia_DE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "upper_columbia_DE_origin_probs.csv"))



# NOT DE CORRECTED
upper_columbia_NDE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> upper_columbia_NDE_prob_table

# subset the probabilities shared across the DPS
upper_columbia_NDE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                     "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                     "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                     "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                     "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> upper_columbia_NDE_DPS_probs

upper_columbia_NDE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                       "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River",             
                       "John Day River",  "Hood River",  "Fifteenmile Creek",  "Umatilla River", "Yakima River",
                       "Walla Walla River", "Wenatchee River", "Entiat River",     "Okanogan River",    "Methow River",                
                       "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> upper_columbia_NDE_origin_probs

# export both tables
write.csv(upper_columbia_NDE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "upper_columbia_NDE_DPS_probs.csv"))
write.csv(upper_columbia_NDE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "upper_columbia_NDE_origin_probs.csv"))

##### CONDITINONAL PROBABILITIES #####

# First - load functions
# 1) Function that reformats movement probabilities (probabilities) into a transition matrix
transition_matrix_reformat <- function(movement_probs, variable_prefix){
  movement_probs %>% 
    mutate(from = gsub(",.*$", "", variable)) %>% 
    mutate(from = gsub(paste0(variable_prefix, "\\["), "", from)) %>% 
    mutate(to = gsub(paste0(variable_prefix, "\\[.*,"), "", variable)) %>% 
    mutate(to = gsub("\\]", "", to)) -> movement_probs
  
  # distribute into a matrix
  # populate with means
  movement_probs %>% 
    dplyr::select(mean, from, to) %>% 
    pivot_wider(., names_from = to, values_from = mean) %>% 
    column_to_rownames("from") %>% 
    as.matrix() -> transition_matrix
  
  return(transition_matrix)
  
}

# 2) final fates function (if fish start in a certain state, what is their final distribution?)
non_upstream_states <- model_states[-grep("Creek Upstream|River Upstream", model_states)]

final_fates_simulation <- function(nsim, transition_matrix, start_state = 2){
  trans <- transition_matrix
  
  # Create an empty state matrix
  state_matrix <- matrix(rep(0, 29), nrow = 29, ncol = 1)
  
  transition_states <- non_upstream_states
  
  rownames(state_matrix) <- transition_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, 29), nrow = 1, ncol = 29)
  colnames(final_fate_matrix) <- transition_states
  
  # run while loop until all fish have been lost
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    state_matrix %>% 
      bind_cols(., t = rep(0, 29)) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, 29*29), ncol = 29, nrow = 29)
    rownames(movement_matrix) <- transition_states
    colnames(movement_matrix) <- transition_states
    for (j in 1:nrow(state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = state_matrix[j,i], prob = trans[j,])
    }
    # sum across rows, store in the state matrix
    state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (state_matrix[29, ncol(state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  data.frame(colSums(final_fate_matrix[,1:28])) %>% 
    dplyr::rename(count = colSums.final_fate_matrix...1.28..) %>% 
    mutate(prop = count/nsim)-> final_fates
  
  outputs <- list(state_matrix, final_fates)
  # return(final_fates)
  return(outputs)
  
}

##### SNAKE - Calculate final fates for different origins ####

subset(snake_summary, grepl("origin1_probs_DE", variable)) -> tucannon_DE_movement_probs
subset(snake_summary, grepl("origin2_probs_DE", variable)) -> asotin_DE_movement_probs
subset(snake_summary, grepl("origin3_probs_DE", variable)) -> clearwater_DE_movement_probs
subset(snake_summary, grepl("origin4_probs_DE", variable)) -> salmon_DE_movement_probs
subset(snake_summary, grepl("origin5_probs_DE", variable)) -> grande_ronde_DE_movement_probs
subset(snake_summary, grepl("origin6_probs_DE", variable)) -> imnaha_DE_movement_probs

transition_matrix_reformat(movement_probs = tucannon_DE_movement_probs, variable_prefix = "origin1_probs_DE") -> tucannon_DE_transition_matrix
transition_matrix_reformat(movement_probs = asotin_DE_movement_probs, variable_prefix = "origin2_probs_DE") -> asotin_DE_transition_matrix
transition_matrix_reformat(movement_probs = clearwater_DE_movement_probs, variable_prefix = "origin3_probs_DE") -> clearwater_DE_transition_matrix
transition_matrix_reformat(movement_probs = salmon_DE_movement_probs, variable_prefix = "origin4_probs_DE") -> salmon_DE_transition_matrix
transition_matrix_reformat(movement_probs = grande_ronde_DE_movement_probs, variable_prefix = "origin5_probs_DE") -> grande_ronde_DE_transition_matrix
transition_matrix_reformat(movement_probs = imnaha_DE_movement_probs, variable_prefix = "origin6_probs_DE") -> imnaha_DE_transition_matrix

# Take out the upstream states from each of these
upstream_indices <- grep("Creek Upstream|River Upstream", model_states)

tucannon_DE_transition_matrix <- tucannon_DE_transition_matrix[-upstream_indices, -upstream_indices]
asotin_DE_transition_matrix <- asotin_DE_transition_matrix[-upstream_indices, -upstream_indices]
clearwater_DE_transition_matrix <- clearwater_DE_transition_matrix[-upstream_indices, -upstream_indices]
salmon_DE_transition_matrix <- salmon_DE_transition_matrix[-upstream_indices, -upstream_indices]
grande_ronde_DE_transition_matrix <- grande_ronde_DE_transition_matrix[-upstream_indices, -upstream_indices]
imnaha_DE_transition_matrix <- imnaha_DE_transition_matrix[-upstream_indices, -upstream_indices]

# Manually change each so that transition probability from loss to loss is 1
tucannon_DE_transition_matrix[29,29] <- 1
asotin_DE_transition_matrix[29,29] <- 1
clearwater_DE_transition_matrix[29,29] <- 1
salmon_DE_transition_matrix[29,29] <- 1
grande_ronde_DE_transition_matrix[29,29] <- 1
imnaha_DE_transition_matrix[29,29] <- 1

# Get final fates of fish for each origin, starting at BON

tucannon_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = tucannon_DE_transition_matrix)[2]
asotin_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = asotin_DE_transition_matrix)[2]
clearwater_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = clearwater_DE_transition_matrix)[2]
salmon_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = salmon_DE_transition_matrix)[2]
grande_ronde_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = grande_ronde_DE_transition_matrix)[2]
imnaha_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = imnaha_DE_transition_matrix)[2]

# Get final fates of fish - conditional on overshoot

# For Snake River fish, we can only really look at this for Tucannon River fish
tucannon_overshoot_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = tucannon_DE_transition_matrix, start_state = 9)[2]
# How do we compare to non-overshoot fish? Start them at ICH to LGR? But then subtract out any overshoot movements manually?
tucannon_non_overshoot_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = tucannon_DE_transition_matrix, start_state = 8)[2]

##### MIDDLE COLUMBIA - Calculate final fates for different origins ####

subset(middle_columbia_summary, grepl("origin1_probs_DE", variable)) -> deschutes_DE_movement_probs
subset(middle_columbia_summary, grepl("origin2_probs_DE", variable)) -> fifteenmile_DE_movement_probs
subset(middle_columbia_summary, grepl("origin3_probs_DE", variable)) -> john_day_DE_movement_probs
subset(middle_columbia_summary, grepl("origin4_probs_DE", variable)) -> umatilla_DE_movement_probs
subset(middle_columbia_summary, grepl("origin5_probs_DE", variable)) -> yakima_DE_movement_probs
subset(middle_columbia_summary, grepl("origin6_probs_DE", variable)) -> walla_walla_DE_movement_probs

transition_matrix_reformat(movement_probs = deschutes_DE_movement_probs, variable_prefix = "origin1_probs_DE") -> deschutes_DE_transition_matrix
transition_matrix_reformat(movement_probs = fifteenmile_DE_movement_probs, variable_prefix = "origin2_probs_DE") -> fifteenmile_DE_transition_matrix
transition_matrix_reformat(movement_probs = john_day_DE_movement_probs, variable_prefix = "origin3_probs_DE") -> john_day_DE_transition_matrix
transition_matrix_reformat(movement_probs = umatilla_DE_movement_probs, variable_prefix = "origin4_probs_DE") -> umatilla_DE_transition_matrix
transition_matrix_reformat(movement_probs = yakima_DE_movement_probs, variable_prefix = "origin5_probs_DE") -> yakima_DE_transition_matrix
transition_matrix_reformat(movement_probs = walla_walla_DE_movement_probs, variable_prefix = "origin6_probs_DE") -> walla_walla_DE_transition_matrix

# Take out the upstream states from each of these
upstream_indices <- grep("Creek Upstream|River Upstream", model_states)

deschutes_DE_transition_matrix <- deschutes_DE_transition_matrix[-upstream_indices, -upstream_indices]
fifteenmile_DE_transition_matrix <- fifteenmile_DE_transition_matrix[-upstream_indices, -upstream_indices]
john_day_DE_transition_matrix <- john_day_DE_transition_matrix[-upstream_indices, -upstream_indices]
umatilla_DE_transition_matrix <- umatilla_DE_transition_matrix[-upstream_indices, -upstream_indices]
yakima_DE_transition_matrix <- yakima_DE_transition_matrix[-upstream_indices, -upstream_indices]
walla_walla_DE_transition_matrix <- walla_walla_DE_transition_matrix[-upstream_indices, -upstream_indices]

# Manually change each so that transition probability from loss to loss is 1
deschutes_DE_transition_matrix[29,29] <- 1
fifteenmile_DE_transition_matrix[29,29] <- 1
john_day_DE_transition_matrix[29,29] <- 1
umatilla_DE_transition_matrix[29,29] <- 1
yakima_DE_transition_matrix[29,29] <- 1
walla_walla_DE_transition_matrix[29,29] <- 1

# Get final fates of fish for each origin, starting at BON

deschutes_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = deschutes_DE_transition_matrix)[2]
fifteenmile_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = fifteenmile_DE_transition_matrix)[2]
john_day_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = john_day_DE_transition_matrix)[2]
umatilla_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = umatilla_DE_transition_matrix)[2]
yakima_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = yakima_DE_transition_matrix)[2]
walla_walla_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = walla_walla_DE_transition_matrix)[2]

# Get final fates of fish - conditional on overshoot

# For Snake River fish, we can only really look at this for Tucannon River fish
deschutes_overshoot_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = deschutes_DE_transition_matrix, start_state = 9)[2]
# How do we compare to non-overshoot fish? Start them at ICH to LGR? But then subtract out any overshoot movements manually?
deschutes_non_overshoot_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = deschutes_DE_transition_matrix, start_state = 8)[2]

##### UPPER COLUMBIA - Calculate final fates for different origins ####

subset(upper_columbia_summary, grepl("origin1_probs_DE", variable)) -> wenatchee_DE_movement_probs
subset(upper_columbia_summary, grepl("origin2_probs_DE", variable)) -> entiat_DE_movement_probs
subset(upper_columbia_summary, grepl("origin3_probs_DE", variable)) -> okanogan_DE_movement_probs
subset(upper_columbia_summary, grepl("origin4_probs_DE", variable)) -> methow_DE_movement_probs

transition_matrix_reformat(movement_probs = wenatchee_DE_movement_probs, variable_prefix = "origin1_probs_DE") -> wenatchee_DE_transition_matrix
transition_matrix_reformat(movement_probs = entiat_DE_movement_probs, variable_prefix = "origin2_probs_DE") -> entiat_DE_transition_matrix
transition_matrix_reformat(movement_probs = okanogan_DE_movement_probs, variable_prefix = "origin3_probs_DE") -> okanogan_DE_transition_matrix
transition_matrix_reformat(movement_probs = methow_DE_movement_probs, variable_prefix = "origin4_probs_DE") -> methow_DE_transition_matrix

# Take out the upstream states from each of these
upstream_indices <- grep("Creek Upstream|River Upstream", model_states)

wenatchee_DE_transition_matrix <- wenatchee_DE_transition_matrix[-upstream_indices, -upstream_indices]
entiat_DE_transition_matrix <- entiat_DE_transition_matrix[-upstream_indices, -upstream_indices]
okanogan_DE_transition_matrix <- okanogan_DE_transition_matrix[-upstream_indices, -upstream_indices]
methow_DE_transition_matrix <- methow_DE_transition_matrix[-upstream_indices, -upstream_indices]

# Manually change each so that transition probability from loss to loss is 1
wenatchee_DE_transition_matrix[29,29] <- 1
entiat_DE_transition_matrix[29,29] <- 1
okanogan_DE_transition_matrix[29,29] <- 1
methow_DE_transition_matrix[29,29] <- 1

# Get final fates of fish for each origin, starting at BON

wenatchee_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = wenatchee_DE_transition_matrix)[2]
entiat_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = entiat_DE_transition_matrix)[2]
okanogan_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = okanogan_DE_transition_matrix)[2]
methow_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = methow_DE_transition_matrix)[2]

# Get final fates of fish - conditional on overshoot

# For Snake River fish, we can only really look at this for Tucannon River fish
wenatchee_overshoot_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = wenatchee_DE_transition_matrix, start_state = 9)[2]
# How do we compare to non-overshoot fish? Start them at ICH to LGR? But then subtract out any overshoot movements manually?
wenatchee_non_overshoot_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = wenatchee_DE_transition_matrix, start_state = 8)[2]

##### Plot detection efficiency #####

# Step 1: Extract correct parameters
det_eff_params <- c("asotin_alpha1", "asotin_alpha2", "deschutes_alpha1", "entiat_alpha1"           ,
                    "fifteenmile_alpha1", "hood_alpha1", "imnaha_alpha1", "john_day_alpha1", "methow_alpha1", "methow_alpha2"           ,
                    "okanogan_alpha1", "tucannon_alpha1", "tucannon_alpha2", "umatilla_alpha1", "umatilla_alpha2", "walla_walla_alpha1",
                    "walla_walla_alpha2", "walla_walla_alpha3", "wenatchee_alpha1", "yakima_alpha1", "asotin_beta", "deschutes_beta"         ,
                    "entiat_beta", "hood_beta", "john_day_beta", "methow_beta"  ,
                    "okanogan_beta", "tucannon_beta", "umatilla_beta", "walla_walla_beta", "wenatchee_beta", "yakima_beta")

subset(snake_summary, variable %in% det_eff_params) -> snake_DE_params
subset(middle_columbia_summary, variable %in% det_eff_params) -> middle_columbia_DE_params
subset(upper_columbia_summary, variable %in% det_eff_params) -> upper_columbia_DE_params


# Step 2: Make a function to predict DE from parameters + discharge
predict_DE <- function(alpha_term, beta_term, params_df){
  subset(params_df, variable %in% c(alpha_term, beta_term)) -> trib_DE
  
  seq(-2,2, 0.01) -> zscore_discharge_values
  
  predicted_DE <- data.frame(discharge = zscore_discharge_values, 
                                 mean_DE = plogis(subset(params_df, variable == alpha_term)$mean + 
                                                    subset(params_df, variable == beta_term)$mean * zscore_discharge_values),
                                 upper_DE = plogis(subset(params_df, variable == alpha_term)$q95 + 
                                                     subset(params_df, variable == beta_term)$mean * zscore_discharge_values),
                                 lower_DE = plogis(subset(params_df, variable == alpha_term)$q5 + 
                                                     subset(params_df, variable == beta_term)$mean * zscore_discharge_values))
  
  predicted_DE %>% 
    pivot_longer(., cols = c("mean_DE", "upper_DE", "lower_DE")) -> predicted_DE_long
  
  return(predicted_DE_long)
  
}

# Step 3: Make a function to plot them
DE_plot_function <- function(df, trib_name){
  plot <- ggplot(df, aes(x = discharge, y = value, color = name, linetype = name)) +
    geom_line(size = 1) +
    xlab("Discharge (Z-scored)") +
    ylab("Detection Efficiency") +
    scale_color_manual(values = c("grey50", "black", "grey50")) +
    scale_linetype_manual(values = c(2,1,2)) +
    ggtitle(trib_name) +
    theme(legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 14)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  return(plot)
}




# Step 4: Plot them
# Middle Columbia Tributaries
predict_DE(alpha_term = "umatilla_alpha1", beta_term = "umatilla_beta", params_df = middle_columbia_DE_params) -> umatilla_predicted_DE1
umatilla_DE_plot1 <- DE_plot_function(df = umatilla_predicted_DE1, trib_name = "Umatilla River (07/08-13/14)")

predict_DE(alpha_term = "umatilla_alpha2", beta_term = "umatilla_beta", params_df = middle_columbia_DE_params) -> umatilla_predicted_DE2
umatilla_DE_plot2 <- DE_plot_function(df = umatilla_predicted_DE2, trib_name = "Umatilla River (14/15-21/22)")

predict_DE(alpha_term = "walla_walla_alpha1", beta_term = "walla_walla_beta", params_df = middle_columbia_DE_params) -> walla_walla_predicted_DE1
walla_walla_DE_plot1 <- DE_plot_function(df = walla_walla_predicted_DE1, trib_name = "Walla Walla River (05/06-11/12)")

predict_DE(alpha_term = "walla_walla_alpha2", beta_term = "walla_walla_beta", params_df = middle_columbia_DE_params) -> walla_walla_predicted_DE2
walla_walla_DE_plot2 <- DE_plot_function(df = walla_walla_predicted_DE2, trib_name = "Walla Walla River (12/13-18/19)")

predict_DE(alpha_term = "walla_walla_alpha3", beta_term = "walla_walla_beta", params_df = middle_columbia_DE_params) -> walla_walla_predicted_DE3
walla_walla_DE_plot3 <- DE_plot_function(df = walla_walla_predicted_DE3, trib_name = "Walla Walla River (19/20-21/22)")

predict_DE(alpha_term = "yakima_alpha1", beta_term = "yakima_beta", params_df = middle_columbia_DE_params) -> yakima_predicted_DE1
yakima_DE_plot1 <- DE_plot_function(df = yakima_predicted_DE1, trib_name = "Yakima River (04/05-21/22)")

predict_DE(alpha_term = "deschutes_alpha1", beta_term = "deschutes_beta", params_df = middle_columbia_DE_params) -> deschutes_predicted_DE1
deschutes_DE_plot1 <- DE_plot_function(df = deschutes_predicted_DE1, trib_name = "Deschutes River (12/13-19/20)")

predict_DE(alpha_term = "john_day_alpha1", beta_term = "john_day_beta", params_df = middle_columbia_DE_params) -> john_day_predicted_DE1
john_day_DE_plot1 <- DE_plot_function(df = john_day_predicted_DE1, trib_name = "John Day River (12/13-21/22)")

# Hood River is technically lower Columbia
predict_DE(alpha_term = "hood_alpha1", beta_term = "hood_beta", params_df = middle_columbia_DE_params) -> hood_predicted_DE1
hood_DE_plot1 <- DE_plot_function(df = hood_predicted_DE1, trib_name = "Hood River (15/16-21/22)")


# Upper Columbia Tributaries
predict_DE(alpha_term = "okanogan_alpha1", beta_term = "okanogan_beta", params_df = upper_columbia_DE_params) -> okanogan_predicted_DE
okanogan_DE_plot <- DE_plot_function(df = okanogan_predicted_DE, trib_name = "Okanogan (12/13-21/22)")

predict_DE(alpha_term = "wenatchee_alpha1", beta_term = "wenatchee_beta", params_df = upper_columbia_DE_params) -> wenatchee_predicted_DE1
wenatchee_DE_plot1 <- DE_plot_function(df = wenatchee_predicted_DE1, trib_name = "Wenatchee River (2010-2020)")

predict_DE(alpha_term = "entiat_alpha1", beta_term = "entiat_beta", params_df = upper_columbia_DE_params) -> entiat_predicted_DE1
entiat_DE_plot1 <- DE_plot_function(df = entiat_predicted_DE1, trib_name = "Entiat River (07/08-21/22)")

predict_DE(alpha_term = "methow_alpha1", beta_term = "methow_beta", params_df = upper_columbia_DE_params) -> methow_predicted_DE1
methow_DE_plot1 <- DE_plot_function(df = methow_predicted_DE1, trib_name = "Methow River (09/10-16/17)")

predict_DE(alpha_term = "methow_alpha2", beta_term = "methow_beta", params_df = upper_columbia_DE_params) -> methow_predicted_DE2
methow_DE_plot2 <- DE_plot_function(df = methow_predicted_DE2, trib_name = "Methow River (17/18-21/22)")


# Snake River Tributaries
predict_DE(alpha_term = "asotin_alpha1", beta_term = "asotin_beta", params_df = snake_DE_params) -> asotin_predicted_DE1
asotin_DE_plot1 <- DE_plot_function(df = asotin_predicted_DE1, trib_name = "Asotin Creek (11/12-17/18)")

predict_DE(alpha_term = "asotin_alpha2", beta_term = "asotin_beta", params_df = snake_DE_params) -> asotin_predicted_DE2
asotin_DE_plot2 <- DE_plot_function(df = asotin_predicted_DE2, trib_name = "Asotin Creek (18/19-21/22)")

predict_DE(alpha_term = "tucannon_alpha1", beta_term = "tucannon_beta", params_df = snake_DE_params) -> tucannon_predicted_DE1
tucannon_DE_plot1 <- DE_plot_function(df = tucannon_predicted_DE1, trib_name = "Tucannon River (10/11-19/20)")

predict_DE(alpha_term = "tucannon_alpha2", beta_term = "tucannon_beta", params_df = snake_DE_params) -> tucannon_predicted_DE2
tucannon_DE_plot2 <- DE_plot_function(df = tucannon_predicted_DE2, trib_name = "Tucannon River (20/21-21/22)")

# predict_DE(alpha_term = "imnaha_alpha1", beta_term = "imnaha_beta", params_df = snake_DE_params) -> imnaha_predicted_DE1
# imnaha_DE_plot1 <- DE_plot_function(df = imnaha_predicted_DE1, trib_name = "Imnaha River (12/13-21/22)")
# 
# predict_DE(alpha_term = "fifteenmile_alpha1", beta_term = "fifteenmile_beta", params_df = middle_columbia_DE_params) -> fifteenmile_predicted_DE1
# fifteenmile_DE_plot1 <- DE_plot_function(df = fifteenmile_predicted_DE1, trib_name = "Fifteenmile Creek (12/13-18/19)")

# Fifteenmile and Imnaha have to be treated differently because they don't have discharge components
# Fifteenmile
subset(middle_columbia_DE_params, variable %in% c("fifteenmile_alpha1")) -> fifteenmile_DE
  
  seq(-2,2, 0.01) -> zscore_discharge_values
  
  fifteenmile_predicted_DE <- data.frame(discharge = zscore_discharge_values, 
                             mean_DE = plogis(subset(middle_columbia_DE_params, variable == "fifteenmile_alpha1")$mean),
                             upper_DE = plogis(subset(middle_columbia_DE_params, variable == "fifteenmile_alpha1")$q95),
                             lower_DE = plogis(subset(middle_columbia_DE_params, variable == "fifteenmile_alpha1")$q5))
  
  fifteenmile_predicted_DE %>% 
    pivot_longer(., cols = c("mean_DE", "upper_DE", "lower_DE")) -> fifteenmile_predicted_DE1
  
fifteenmile_DE_plot1 <- DE_plot_function(df = fifteenmile_predicted_DE1, trib_name = "Fifteenmile Creek (12/13-18/19)")

# Imnaha
subset(snake_DE_params, variable %in% c("imnaha_alpha1")) -> imnaha_DE

seq(-2,2, 0.01) -> zscore_discharge_values

imnaha_predicted_DE <- data.frame(discharge = zscore_discharge_values, 
                                       mean_DE = plogis(subset(middle_columbia_DE_params, variable == "imnaha_alpha1")$mean),
                                       upper_DE = plogis(subset(middle_columbia_DE_params, variable == "imnaha_alpha1")$q95),
                                       lower_DE = plogis(subset(middle_columbia_DE_params, variable == "imnaha_alpha1")$q5))

imnaha_predicted_DE %>% 
  pivot_longer(., cols = c("mean_DE", "upper_DE", "lower_DE")) -> imnaha_predicted_DE1

imnaha_DE_plot1 <- DE_plot_function(df = imnaha_predicted_DE1, trib_name = "Imnaha River (12/13-21/22)")


# ARRANGE ALL DE PLOTS

ggarrange(umatilla_DE_plot1, 
          umatilla_DE_plot2,
                       walla_walla_DE_plot1,
                       walla_walla_DE_plot2,
                       walla_walla_DE_plot3,
                       yakima_DE_plot1,
                       deschutes_DE_plot1,
                       john_day_DE_plot1,
                       hood_DE_plot1,
                       okanogan_DE_plot,
                       wenatchee_DE_plot1,
                       entiat_DE_plot1,
                       methow_DE_plot1,
                       methow_DE_plot2,
                       asotin_DE_plot1,
                       asotin_DE_plot2,
                       tucannon_DE_plot1,
                       tucannon_DE_plot2,
          fifteenmile_DE_plot1,
                       imnaha_DE_plot1, ncol = 5, nrow = 4) -> DE_plots

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "DE_plots.pdf"), DE_plots, height = 15, width = 20)




