### 18 - final fates v2

library(tidyverse)
library(lubridate)

##### STEP 1: Create the transition matrix #####

# So we will be able to just input the derived movement probabilities from the stan model, which should make it easier

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

nstates <- length(model_states)

transition_matrix <- matrix(0, nrow = nstates, ncol = nstates)
rownames(transition_matrix) <- model_states
colnames(transition_matrix) <- model_states

# Populate every possible option with a 1
# 1: mainstem, mouth to BON
transition_matrix["mainstem, mouth to BON", "loss"] <- 1
transition_matrix["mainstem, mouth to BON", "mainstem, BON to MCN"] <- 1


# 2: mainstem, BON to MCN
transition_matrix["mainstem, BON to MCN", "loss"] <- 1
transition_matrix["mainstem, BON to MCN", "mainstem, mouth to BON"] <- 1
transition_matrix["mainstem, BON to MCN", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- 1
transition_matrix["mainstem, BON to MCN", "Deschutes River Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "Deschutes River Upstream"] <- 1
# transition_matrix["mainstem, BON to MCN", "John Day River"] <- 1
transition_matrix["mainstem, BON to MCN", "John Day River Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "John Day River Upstream"] <- 1
# transition_matrix["mainstem, BON to MCN", "Hood River"] <- 1
transition_matrix["mainstem, BON to MCN", "Hood River Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "Hood River Upstream"] <- 1
# transition_matrix["mainstem, BON to MCN", "Fifteenmile Creek"] <- 1
transition_matrix["mainstem, BON to MCN", "Fifteenmile Creek Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "Fifteenmile Creek Upstream"] <- 1
# transition_matrix["mainstem, BON to MCN", "Umatilla River"] <- 1
transition_matrix["mainstem, BON to MCN", "Umatilla River Mouth"] <- 1
transition_matrix["mainstem, BON to MCN", "Umatilla River Upstream"] <- 1
transition_matrix["mainstem, BON to MCN", "BON to MCN other tributaries"] <- 1


# 3: mainstem, MCN to ICH or PRA
transition_matrix["mainstem, MCN to ICH or PRA", "loss"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, BON to MCN"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR"] <- 1
# transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River Mouth"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River Upstream"] <- 1
# transition_matrix["mainstem, MCN to ICH or PRA", "Walla Walla River"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Walla Walla River Mouth"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Walla Walla River Upstream"] <- 1


# 4: mainstem, PRA to RIS
transition_matrix["mainstem, PRA to RIS", "loss"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, RIS to RRE"] <- 1


# 5: mainstem, RIS to RRE
transition_matrix["mainstem, RIS to RRE", "loss"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, RRE to WEL"] <- 1
# transition_matrix["mainstem, RIS to RRE", "Wenatchee River"] <- 1
transition_matrix["mainstem, RIS to RRE", "Wenatchee River Mouth"] <- 1
transition_matrix["mainstem, RIS to RRE", "Wenatchee River Upstream"] <- 1


# 6: mainstem, RRE to WEL
transition_matrix["mainstem, RRE to WEL", "loss"] <- 1
transition_matrix["mainstem, RRE to WEL", "mainstem, RIS to RRE"] <- 1
transition_matrix["mainstem, RRE to WEL", "mainstem, upstream of WEL"] <- 1
# transition_matrix["mainstem, RRE to WEL", "Entiat River"] <- 1
transition_matrix["mainstem, RRE to WEL", "Entiat River Mouth"] <- 1
transition_matrix["mainstem, RRE to WEL", "Entiat River Upstream"] <- 1


# 7: mainstem, upstream of WEL
transition_matrix["mainstem, upstream of WEL", "loss"] <- 1
transition_matrix["mainstem, upstream of WEL", "mainstem, RRE to WEL"] <- 1
# transition_matrix["mainstem, upstream of WEL", "Okanogan River"] <- 1
transition_matrix["mainstem, upstream of WEL", "Okanogan River Mouth"] <- 1
transition_matrix["mainstem, upstream of WEL", "Okanogan River Upstream"] <- 1
# transition_matrix["mainstem, upstream of WEL", "Methow River"] <- 1
transition_matrix["mainstem, upstream of WEL", "Methow River Mouth"] <- 1
transition_matrix["mainstem, upstream of WEL", "Methow River Upstream"] <- 1
transition_matrix["mainstem, upstream of WEL", "Upstream WEL other tributaries"] <- 1


# 8: mainstem, ICH to LGR
transition_matrix["mainstem, ICH to LGR", "loss"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, upstream of LGR"] <- 1
# transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- 1
transition_matrix["mainstem, ICH to LGR", "Tucannon River Mouth"] <- 1
transition_matrix["mainstem, ICH to LGR", "Tucannon River Upstream"] <- 1


# 9: mainstem, upstream of LGR
transition_matrix["mainstem, upstream of LGR", "loss"] <- 1
transition_matrix["mainstem, upstream of LGR", "mainstem, ICH to LGR"] <- 1
# transition_matrix["mainstem, upstream of LGR", "Asotin Creek"] <- 1
transition_matrix["mainstem, upstream of LGR", "Asotin Creek Mouth"] <- 1
transition_matrix["mainstem, upstream of LGR", "Asotin Creek Upstream"] <- 1
transition_matrix["mainstem, upstream of LGR", "Clearwater River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Salmon River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Grande Ronde River"] <- 1
# transition_matrix["mainstem, upstream of LGR", "Imnaha River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Imnaha River Mouth"] <- 1
transition_matrix["mainstem, upstream of LGR", "Imnaha River Upstream"] <- 1


# Deschutes River
# 10: Deschutes River Mouth
# 11: Deschutes River Upstream
transition_matrix["Deschutes River Mouth", "loss"] <- 1
transition_matrix["Deschutes River Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["Deschutes River Mouth", "Deschutes River Upstream"] <- 1
transition_matrix["Deschutes River Upstream", "loss"] <- 1
transition_matrix["Deschutes River Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["Deschutes River Upstream", "Deschutes River Mouth"] <- 1


# John Day River
# 12: John Day River Mouth
# 13: John Day River Upstream
transition_matrix["John Day River Mouth", "loss"] <- 1
transition_matrix["John Day River Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["John Day River Mouth", "John Day River Upstream"] <- 1
transition_matrix["John Day River Upstream", "loss"] <- 1
transition_matrix["John Day River Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["John Day River Upstream", "John Day River Mouth"] <- 1


# Hood River
# 14: Hood River Mouth
# 15: Hood River Upstream
transition_matrix["Hood River Mouth", "loss"] <- 1
transition_matrix["Hood River Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["Hood River Mouth", "Hood River Upstream"] <- 1
transition_matrix["Hood River Upstream", "loss"] <- 1
transition_matrix["Hood River Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["Hood River Upstream", "Hood River Mouth"] <- 1


# Fifteenmile Creek
# 16: Fifteenmile Creek Mouth
# 17: Fifteenmile Creek Upstream
transition_matrix["Fifteenmile Creek Mouth", "loss"] <- 1
transition_matrix["Fifteenmile Creek Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["Fifteenmile Creek Mouth", "Fifteenmile Creek Upstream"] <- 1
transition_matrix["Fifteenmile Creek Upstream", "loss"] <- 1
transition_matrix["Fifteenmile Creek Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["Fifteenmile Creek Upstream", "Fifteenmile Creek Mouth"] <- 1


# Umatilla River
# 18: Umatilla River Mouth
# 19: Umatilla River Upstream
transition_matrix["Umatilla River Mouth", "loss"] <- 1
transition_matrix["Umatilla River Mouth", "mainstem, BON to MCN"] <- 1
# transition_matrix["Umatilla River Mouth", "Umatilla River Upstream"] <- 1
transition_matrix["Umatilla River Upstream", "loss"] <- 1
transition_matrix["Umatilla River Upstream", "mainstem, BON to MCN"] <- 1
# transition_matrix["Umatilla River Upstream", "Umatilla River Mouth"] <- 1


# Yakima River
# 20: Yakima River Mouth
# 21: Yakima River Upstream
transition_matrix["Yakima River Mouth", "loss"] <- 1
transition_matrix["Yakima River Mouth", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["Yakima River Mouth", "Yakima River Upstream"] <- 1
transition_matrix["Yakima River Upstream", "loss"] <- 1
transition_matrix["Yakima River Upstream", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["Yakima River Upstream", "Yakima River Mouth"] <- 1


# Walla Walla River
# 22: Walla Walla River Mouth
# 23: Walla Walla River Upstream
transition_matrix["Walla Walla River Mouth", "loss"] <- 1
transition_matrix["Walla Walla River Mouth", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["Walla Walla River Mouth", "Walla Walla River Upstream"] <- 1
transition_matrix["Walla Walla River Upstream", "loss"] <- 1
transition_matrix["Walla Walla River Upstream", "mainstem, MCN to ICH or PRA"] <- 1
# transition_matrix["Walla Walla River Upstream", "Walla Walla River Mouth"] <- 1


# Wenatchee River
# 24: Wenatchee River Mouth
# 25: Wenatchee River Upstream
transition_matrix["Wenatchee River Mouth", "loss"] <- 1
transition_matrix["Wenatchee River Mouth", "mainstem, RIS to RRE"] <- 1
# transition_matrix["Wenatchee River Mouth", "Wenatchee River Upstream"] <- 1
transition_matrix["Wenatchee River Upstream", "loss"] <- 1
transition_matrix["Wenatchee River Upstream", "mainstem, RIS to RRE"] <- 1
# transition_matrix["Wenatchee River Upstream", "Wenatchee River Mouth"] <- 1


# Entiat River
# 26: Entiat River Mouth
# 27: Entiat River Upstream
transition_matrix["Entiat River Mouth", "loss"] <- 1
transition_matrix["Entiat River Mouth", "mainstem, RRE to WEL"] <- 1
# transition_matrix["Entiat River Mouth", "Entiat River Upstream"] <- 1
transition_matrix["Entiat River Upstream", "loss"] <- 1
transition_matrix["Entiat River Upstream", "mainstem, RRE to WEL"] <- 1
# transition_matrix["Entiat River Upstream", "Entiat River Mouth"] <- 1


# Okanogan River
# 28: Okanogan River Mouth
# 29: Okanogan River Upstream
transition_matrix["Okanogan River Mouth", "loss"] <- 1
transition_matrix["Okanogan River Mouth", "mainstem, upstream of WEL"] <- 1
# transition_matrix["Okanogan River Mouth", "Okanogan River Upstream"] <- 1
transition_matrix["Okanogan River Upstream", "loss"] <- 1
transition_matrix["Okanogan River Upstream", "mainstem, upstream of WEL"] <- 1
# transition_matrix["Okanogan River Upstream", "Okanogan River Mouth"] <- 1


# Methow River
# 30: Methow River Mouth
# 31: Methow River Upstream
transition_matrix["Methow River Mouth", "loss"] <- 1
transition_matrix["Methow River Mouth", "mainstem, upstream of WEL"] <- 1
# transition_matrix["Methow River Mouth", "Methow River Upstream"] <- 1
transition_matrix["Methow River Upstream", "loss"] <- 1
transition_matrix["Methow River Upstream", "mainstem, upstream of WEL"] <- 1
# transition_matrix["Methow River Upstream", "Methow River Mouth"] <- 1


# Tucannon River
# 32: Tucannon River Mouth
# 33: Tucannon River Upstream
transition_matrix["Tucannon River Mouth", "loss"] <- 1
transition_matrix["Tucannon River Mouth", "mainstem, ICH to LGR"] <- 1
# transition_matrix["Tucannon River Mouth", "Tucannon River Upstream"] <- 1
transition_matrix["Tucannon River Upstream", "loss"] <- 1
transition_matrix["Tucannon River Upstream", "mainstem, ICH to LGR"] <- 1
# transition_matrix["Tucannon River Upstream", "Tucannon River Mouth"] <- 1


# Asotin Creek
# 34: Asotin Creek Mouth
# 35: Asotin Creek Upstream
transition_matrix["Asotin Creek Mouth", "loss"] <- 1
transition_matrix["Asotin Creek Mouth", "mainstem, upstream of LGR"] <- 1
# transition_matrix["Asotin Creek Mouth", "Asotin Creek Upstream"] <- 1
transition_matrix["Asotin Creek Upstream", "loss"] <- 1
transition_matrix["Asotin Creek Upstream", "mainstem, upstream of LGR"] <- 1
# transition_matrix["Asotin Creek Upstream", "Asotin Creek Mouth"] <- 1


# 36: Clearwater River
transition_matrix["Clearwater River", "loss"] <- 1
transition_matrix["Clearwater River", "mainstem, upstream of LGR"] <- 1


# 37: Salmon River
transition_matrix["Salmon River", "loss"] <- 1
transition_matrix["Salmon River", "mainstem, upstream of LGR"] <- 1


# 38: Grande Ronde River
transition_matrix["Grande Ronde River", "loss"] <- 1
transition_matrix["Grande Ronde River", "mainstem, upstream of LGR"] <- 1


# Imnaha River
# 39: Imnaha River Mouth
# 40: Imnaha River Upstream
transition_matrix["Imnaha River Mouth", "loss"] <- 1
transition_matrix["Imnaha River Mouth", "mainstem, upstream of LGR"] <- 1
# transition_matrix["Imnaha River Mouth", "Imnaha River Upstream"] <- 1
transition_matrix["Imnaha River Upstream", "loss"] <- 1
transition_matrix["Imnaha River Upstream", "mainstem, upstream of LGR"] <- 1
# transition_matrix["Imnaha River Upstream", "Imnaha River Mouth"] <- 1

# 41: BON to MCN other tributaries
transition_matrix["BON to MCN other tributaries", "loss"] <- 1
transition_matrix["BON to MCN other tributaries", "mainstem, BON to MCN"] <- 1

# 42: Upstream WEL other tributaries
transition_matrix["Upstream WEL other tributaries", "loss"] <- 1
transition_matrix["Upstream WEL other tributaries", "mainstem, upstream of WEL"] <- 1


# Now - remove all upstream states
transition_matrix2 <- transition_matrix[seq(1,43,1)[-grep("Creek Upstream|River Upstream", rownames(transition_matrix))], seq(1,43,1)[-grep("Creek Upstream|River Upstream", colnames(transition_matrix))]]

upstream_indices <- grep("Creek Upstream|River Upstream", rownames(transition_matrix))

##### now do final fates (?) #####

# We previously did it via simulation. That's a ton of simulations, would take a long time.

# Can we use an analytic solution? I think so - pick a number of transitions (say 50, or 100), then multiply out all parameter values

# First - they all have to start in mainstem, BON to MCN, since we first see them at BON

# We need to create a df, which contains the total probability of ending up in any state

# But don't include any upstream states, since we don't care about those (they don't get a separate probability from the RM states)
final_fate_probabilities <- data.frame(final_state = rownames(transition_matrix2), total_prob = 0)

# Load the movement estimates from the stan model output
# Temporarily, use estimates from this intermediate stan model:
snake_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "diagnostics", "snake20", "100iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))

# Get the summary object (this takes a while)
snake_summary <- snake_fit$summary()

# Subset all of the different origins (since they each have different probabilities)
# For the Snake:
# Origin 1 = Tucannon River
# Origin 2 = Asotin Creek
# Origin 3 = Clearwater River
# Origin 4 = Salmon River
# Origin 5 = Grande Ronde River
# Origin 6 = Imnaha River

subset(snake_summary, grepl("origin1_probs_DE", variable)) -> origin1_DE_movement_probs
subset(snake_summary, grepl("origin2_probs_DE", variable)) -> origin2_DE_movement_probs
subset(snake_summary, grepl("origin3_probs_DE", variable)) -> origin3_DE_movement_probs
subset(snake_summary, grepl("origin4_probs_DE", variable)) -> origin4_DE_movement_probs
subset(snake_summary, grepl("origin5_probs_DE", variable)) -> origin5_DE_movement_probs
subset(snake_summary, grepl("origin6_probs_DE", variable)) -> origin6_DE_movement_probs

# Create a function to reformat these into a transition matrix

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

transition_matrix_reformat(movement_probs = origin1_DE_movement_probs, variable_prefix = "origin1_probs_DE") -> origin1_DE_transition_matrix

# Take out the upstream states
origin1_DE_transition_matrix <- origin1_DE_transition_matrix[-upstream_indices, -upstream_indices] 
# Manually change it so that transition probability from loss to loss is 1
origin1_DE_transition_matrix[29,29] <- 1
# rownames(origin1_DE_transition_matrix) <- seq(1,29)
# colnames(origin1_DE_transition_matrix) <- seq(1,29)
# check that the rows still sum to 1 - yep!
rowSums(origin1_DE_transition_matrix)

##### Final fates - simulation-based approach function #####

# Start state argument is default to 2 (BON to MCN, ie fish once we sae them at BON), 
# but can be others (so for example, can start fish in an overshoot state)
final_fates_simulation <- function(nsim, transition_matrix, start_state = 2){
  trans <- transition_matrix
  
  # Create an empty state matrix
  state_matrix <- matrix(rep(0, 29), nrow = 29, ncol = 1)
  
  transition_states <- rownames(transition_matrix2)
  
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

##### Calculate final fates #####
print(paste0("final fate start time: ", Sys.time()))
tuc_river_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = origin1_DE_transition_matrix)[2]
print(paste0("final fate end time: ", Sys.time()))
# for one million fish, this only takes four seconds
