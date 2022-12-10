# UPPER COLUMBIA - 03_stan_actual_int_origin - for mox!

# in this version of the model, we now have tributary detection efficiencies

# this model will fit to the following populations/origins:
# "Wenatchee_River", "Entiat_River", "Okanogan_River","Methow_River"


# This script will fit an intercept + origin model using stan to the actual dataset

# FOR TESTING: setwd
# setwd("/Users/markusmin/Documents/CBR/steelhead/stan_actual/deteff_ESU_models/upper_columbia/")

# library("rstan")
library(cmdstanr)
library(posterior)
library(tidyverse)
library(lubridate)
# library(bayesplot)

##### Get all of the data in order #####

# Create a transition matrix of 1s and 0s for movements from (rows) to (columns)

# Create a starting matrix of sites, with the first occasion populated (BON to MCN)
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
# 43 states

# Create a rows = from, columns = to matrix for movement probabilities


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


##### for this ESU - cut out state transitions for which we have no data #####
# Remove the states that no fish ever enter
# SEE END OF SCRIPT WHERE WE CHECK
# Check to see if every transition in our model is represented

# When we remove some of the upstream transitions, we're left with a lot of upstream states
# that no fish have ever been in.
# There are also never any upstream to RM or vice versa transitions, because we removed those. I set
# those all to zeros in the above chunk.

# Now that we have separate river mouth and upstream states, no fish have ever been in:
# 1) "Deschutes River Upstream"     
# 2) "John Day River Mouth"         
# 3) "John Day River Upstream"      
# 4) "Hood River Upstream"          
# 5) "Fifteenmile Creek Mouth"     
# 6) "Fifteenmile Creek Upstream"   
# 7) "Umatilla River Mouth"         
# 8) "Umatilla River Upstream"      
# 9) "Yakima River Mouth"           
# 10) "Yakima River Upstream"       
# 11) "Walla Walla River Mouth"      
# 12) "Walla Walla River Upstream"   
# 13) "Tucannon River Mouth"         
# 14) "Tucannon River Upstream"      
# 15) "Asotin Creek Mouth"          
# 16)  "Asotin Creek Upstream"        
# 17) "Salmon River"                 
# 18) "Grande Ronde River"           
# 19) "Imnaha River Mouth"           
# 20) "Imnaha River Upstream"       
# 21) "BON to MCN other tributaries" 
# set these all to zero
transition_matrix["Deschutes River Upstream", ] <- 0
transition_matrix["John Day River Mouth", ] <- 0
transition_matrix["John Day River Upstream", ] <- 0
transition_matrix["Hood River Upstream", ] <- 0
transition_matrix["Fifteenmile Creek Mouth", ] <- 0
transition_matrix["Fifteenmile Creek Upstream", ] <- 0
transition_matrix["Umatilla River Mouth", ] <- 0
transition_matrix["Umatilla River Upstream", ] <- 0
transition_matrix["Yakima River Mouth", ] <- 0
transition_matrix["Yakima River Upstream", ] <- 0
transition_matrix["Walla Walla River Mouth" , ] <- 0
transition_matrix["Walla Walla River Upstream" , ] <- 0
transition_matrix["Tucannon River Mouth", ] <- 0
transition_matrix["Tucannon River Upstream", ] <- 0
transition_matrix["Asotin Creek Mouth", ] <- 0
transition_matrix["Asotin Creek Upstream", ] <- 0
transition_matrix["Salmon River", ] <- 0
transition_matrix["Grande Ronde River", ] <- 0
transition_matrix["Imnaha River Mouth", ] <- 0
transition_matrix["Imnaha River Upstream", ] <- 0
transition_matrix["BON to MCN other tributaries", ] <- 0

##### Data continued #####
# Use the transition matrix to calculate the possible movements
possible_movements <- rowSums(transition_matrix)

# Get the indices that are 1s (except loss, since that's 1 minus the others)
movements <- which(transition_matrix[,1:(nstates-1)] == 1, arr.ind = TRUE)
nmovements <- dim(movements)[1]

# Now get all of the movements which are fixed to zero
not_movements <- which(transition_matrix[,1:(nstates-1)] == 0, arr.ind = TRUE)
n_notmovements <- dim(not_movements)[1]

# Paste movements to name parameters
movements %>% 
  as.data.frame() %>% 
  mutate(b0_matrix_name = paste0("b0_matrix_", row, "_", col)) -> b0_matrix_names

movements %>% 
  as.data.frame() -> movements_df

# Now, note which movements are movements from mainstem into tributaries (these will need two versions of parameters)
print(paste0(seq(1,43,1), " - ", model_states)) # this just makes it easier to tell state indexing
trib_det_possible <- c(seq(10,35,1), c(39,40))
mainstem_indices <- seq(1,9,1)



# We also need to note what the indices are of river mouth vs. upstream tributary sites. In years where
# we can calculate detection efficiency at the river mouth, we need to ignore movements into to upstream sites
# In years where we can't, we don't ignore them, but the same parameters govern movements into both upstream
# and river mouth sites.

# We also don't care about any movements from upstream to river mouth states and vice versa
upstream_indices <- c(11,13,15,17,19,21,23,25,27,29,31,33,35,40)
river_mouth_indices <- upstream_indices - 1

b0_matrix_names %>% 
  mutate(mainstem_upstream_movement = ifelse(row %in% mainstem_indices & col %in% upstream_indices, 1, 0)) -> b0_matrix_names

b0_matrix_names %>% 
  mutate(mainstem_river_mouth_movement = ifelse(row %in% mainstem_indices & col %in% river_mouth_indices, 1, 0)) -> b0_matrix_names

b0_matrix_names %>% 
  mutate(within_trib_movement = ifelse(row %in% upstream_indices & col %in% river_mouth_indices |
                                       col %in% upstream_indices & row %in% river_mouth_indices, 1, 0)) -> b0_matrix_names

b0_matrix_names %>% 
  mutate(upstream_mainstem_movement = ifelse(row %in% upstream_indices & col %in% mainstem_indices, 1, 0)) -> b0_matrix_names

b0_matrix_names %>% 
  mutate(river_mouth_mainstem_movement = ifelse(row %in% river_mouth_indices & col %in% mainstem_indices, 1, 0)) -> b0_matrix_names


movements_df %>% 
  mutate(mainstem_upstream_movement = ifelse(row %in% mainstem_indices & col %in% upstream_indices, 1, 0)) -> movements_df

movements_df %>% 
  mutate(mainstem_river_mouth_movement = ifelse(row %in% mainstem_indices & col %in% river_mouth_indices, 1, 0)) -> movements_df

movements_df %>% 
  mutate(within_trib_movement = ifelse(row %in% upstream_indices & col %in% river_mouth_indices |
                                         col %in% upstream_indices & row %in% river_mouth_indices, 1, 0)) -> movements_df

movements_df %>% 
  mutate(upstream_mainstem_movement = ifelse(row %in% upstream_indices & col %in% mainstem_indices, 1, 0)) -> movements_df

movements_df %>% 
  mutate(river_mouth_mainstem_movement = ifelse(row %in% river_mouth_indices & col %in% mainstem_indices, 1, 0)) -> movements_df

# upper_columbia_movements <- subset(movements_df, row %in% c(8,9,21,22,23,24,25,26) |
#                                      col %in% c(8,9,21,22,23,24,25,26))
upper_columbia_movements <- subset(movements_df, row %in% c(4,5,6,7,seq(24,31,1),42) |
                            col %in% c(4,5,6,7,seq(24,31,1),42))

non_upper_columbia_movements  <- subset(movements_df, !(row %in% c(4,5,6,7,seq(24,31,1),42)) &
                                 !(col %in% c(4,5,6,7,seq(24,31,1),42)))

# I think these both need to be order by from (rows), to (columns)
upper_columbia_movements %>% 
  arrange(row, col) -> upper_columbia_movements

b0_matrix_names %>% 
  arrange(row,col) -> b0_matrix_names




# print declaration of b0_matrix parameters - distinguish which need two versions based on mainstem_trib_movement column
# Note that this block doesn't change between models for different ESUs
for (i in 1:(nrow(b0_matrix_names))){
  # Paste the numerator
  
  # If it is a movement from the mainstem into the river mouth state, then it gets two versions of the parameter.
  # If it's not, then it just gets one.
  # If its a movement from the upstream to the river mouth state or vice versa, it doesn't get a parameter at all - 
  # we don't care about these movements.
  
  # Movements from mainstem to river mouth: print two versions
  if (b0_matrix_names$mainstem_river_mouth_movement[i] == 1){
    cat("real ", b0_matrix_names$b0_matrix_name[i], "_DE", ";", "\n", sep = "")
    cat("real ", b0_matrix_names$b0_matrix_name[i], "_NDE", ";", "\n", sep = "")
  }
  
  # If it's a within tributary movement - we don't want it. These have been
  # removed from the data (we don't care about them!)
  else if (b0_matrix_names$within_trib_movement[i] == 1){
    # do nothing!
  }
  # If it's a movement from mainstem to upstream - we don't want a parameter for it,
  # because it'll share the same parameter (see transformed parameters block) 
  # as the movement from mainstem to river mouth
  else if (b0_matrix_names$mainstem_upstream_movement[i] == 1){
    # do nothing!
  }
  
  # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
  # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
  # these transitions will get the same parameters as the river mouth to mainstem
  else if (b0_matrix_names$upstream_mainstem_movement[i] == 1){
    # do nothing!
  } else {
  # Finally - if it's any other movement, just treat it like normal!
  # If it's not, print just the one version
    cat("real ", b0_matrix_names$b0_matrix_name[i],";", "\n", sep = "")
  }
  
}

# In this ESU, there are four origins, so we will have three origin parameters
# We will only allow an origin effect within the ESU (so after PRA). Before, they will only have an intercept term
# print declaration of borigin parameters
for (i in 1:3){
  for (j in 1:nrow(upper_columbia_movements)){
    
    # Movements from mainstem to river mouth: print two versions
    if (upper_columbia_movements$mainstem_river_mouth_movement[j] == 1){
      cat("real ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_DE", ";", "\n", sep = "")
      cat("real ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_NDE", ";", "\n", sep = "")
    }
    
    # If it's a within tributary movement - we don't want it
    else if (upper_columbia_movements$within_trib_movement[j] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream - we also don't want it
    else if (upper_columbia_movements$mainstem_upstream_movement[j] == 1){
      # do nothing!
    } 
    # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
    # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
    # these transitions will get the same parameters as the river mouth to mainstem
    else if (upper_columbia_movements$upstream_mainstem_movement[j] == 1){
      # do nothing!
    } else {
      # Finally - if it's any other movement, just treat it like normal!
      # If it's not, print just the one version
      cat("real ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], ";", "\n", sep = "")
    }
    
  }
  cat ("\n")
}


# write it out for the transformed parameters section

# Write out two separate matrices - an NDE (no detection efficiency correction) and a DE (detection efficiency) matrix
for (i in 1:(nrow(b0_matrix_names))){
  # Movements from mainstem to river mouth: store the two different versions
  if (b0_matrix_names$mainstem_river_mouth_movement[i] == 1){
      cat("b0_matrix_DE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i], "_DE",";", "\n", sep = "")
      cat("b0_matrix_NDE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i], "_NDE",";", "\n", sep = "")
    
  }
  
  # If it's a within tributary movement - we don't want it
  else if (b0_matrix_names$within_trib_movement[i] == 1){
    # do nothing!
  }
  # If it's a movement from mainstem to upstream:
  # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
  # For DE - we need to give these all -100000 values, so that in the logit they come out as zeros (since they're not transitions that are allowed)
  else if (b0_matrix_names$mainstem_upstream_movement[i] == 1){
    # cat("b0_matrix_DE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i-1], "_DE",";", "\n", sep = "")
    cat("b0_matrix_DE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", "-100000;", "\n", sep = "")
    cat("b0_matrix_NDE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i-1], "_NDE",";", "\n", sep = "")
  }
  # If it's a movement from the upstream state back to the mainstem, give it the same parameter as movement from the river mouth to the upstream
  # And note - that these are not DE or NDE parameters
  # Because of order, we need -2 here instead of -1, to pick the right from state
  # NO WE DON'T, just use -1
  else if (b0_matrix_names$upstream_mainstem_movement[i] == 1){
    cat("b0_matrix_DE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i-1],";", "\n", sep = "")
    cat("b0_matrix_NDE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i-1],";", "\n", sep = "")
  # If it's not, store the same parameter in both matrices
  } else {
    cat("b0_matrix_DE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i], ";", "\n", sep = "")
    cat("b0_matrix_NDE[", b0_matrix_names$row[i], ",", b0_matrix_names$col[i], "]", " = ", b0_matrix_names$b0_matrix_name[i], ";", "\n", sep = "")
  }
    
} 

# There are four origins, so we need three parameters (origins - 1)
for (i in 1:3){
  for (j in 1:nrow(upper_columbia_movements)){
    # Movements from mainstem to river mouth: store the two different versions
    if (upper_columbia_movements$mainstem_river_mouth_movement[j] == 1){
      cat("borigin", i, "_matrix_DE[", upper_columbia_movements$row[j], ",", upper_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_DE", ";", "\n", sep = "")
      cat("borigin", i, "_matrix_NDE[", upper_columbia_movements$row[j], ",", upper_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_NDE", ";", "\n", sep = "")
      
    }
    
    # If it's a within tributary movement - we don't want it
    else if (upper_columbia_movements$within_trib_movement[j] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream:
    # For NDE - give it the same parameter as the one for the mainstem to the river mouth site (which by index, is the site before)
    # For DE - we don't care about this, because we're removing the upstream states
    else if (upper_columbia_movements$mainstem_upstream_movement[j] == 1){
      # cat("borigin", i, "_matrix_DE[", upper_columbia_movements$row[j], ",", upper_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j-1], "_", upper_columbia_movements$col[j-1], "_DE", ";", "\n", sep = "")
      cat("borigin", i, "_matrix_NDE[", upper_columbia_movements$row[j], ",", upper_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j-1], "_", upper_columbia_movements$col[j-1], "_NDE", ";", "\n", sep = "")
      
    }
    # If it's a movement from the upstream state back to the mainstem, give it the same parameter as movement from the river mouth to the upstream
    # And note - that these are not DE or NDE parameters
    else if (upper_columbia_movements$upstream_mainstem_movement[j] == 1){
      cat("borigin", i, "_matrix_DE[", upper_columbia_movements$row[j], ",", upper_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j-1], "_", upper_columbia_movements$col[j],  ";", "\n", sep = "")
      cat("borigin", i, "_matrix_NDE[", upper_columbia_movements$row[j], ",", upper_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j-1], "_", upper_columbia_movements$col[j], ";", "\n", sep = "")
      # If it's not, store the same parameter in both matrices
    } else {
      cat("borigin", i, "_matrix_DE[", upper_columbia_movements$row[j], ",", upper_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], ";", "\n", sep = "")
      cat("borigin", i, "_matrix_NDE[", upper_columbia_movements$row[j], ",", upper_columbia_movements$col[j], "]", " = ", "borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j],  ";", "\n", sep = "")
    }
    
  }
  
  cat ("\n")
}
  


# write out the priors
# print b0 matrix priors
for (i in 1:(nrow(b0_matrix_names))){
  
  # Paste the numerator
  
  # If it is a movement from the mainstem into the river mouth state, then it gets two versions of the parameter.
  # If it's not, then it just gets one.
  # If its a movement from the upstream to the river mouth state or vice versa, it doesn't get a parameter at all - 
  # we don't care about these movements.
  
  # Movements from mainstem to river mouth: print two versions
  if (b0_matrix_names$mainstem_river_mouth_movement[i] == 1){
    cat(b0_matrix_names$b0_matrix_name[i], "_DE", " ~ normal(0,10);", "\n", sep = "")
    cat(b0_matrix_names$b0_matrix_name[i], "_NDE", " ~ normal(0,10);", "\n", sep = "")
  }
  
  # If it's a within tributary movement - we don't want it (no prior)
  else if (b0_matrix_names$within_trib_movement[i] == 1){
    # do nothing!
  }
  # If it's a movement from mainstem to upstream - we also don't want it (no prior)
  else if (b0_matrix_names$mainstem_upstream_movement[i] == 1){
    # do nothing!
  }
  
  # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
  # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter (so no prior).
  # these transitions will get the same parameters as the river mouth to mainstem
  else if (b0_matrix_names$upstream_mainstem_movement[i] == 1){
    # do nothing!
  } else {
    # Finally - if it's any other movement, just treat it like normal!
    # If it's not, print just the one version
    cat(b0_matrix_names$b0_matrix_name[i], " ~ normal(0,10);", "\n", sep = "")
  }
  

}

# There are four origins in this model, so we will have three origin parameters.
for (i in 1:3){
  for (j in 1:nrow(upper_columbia_movements)){
    
    # Movements from mainstem to river mouth: print two versions
    if (upper_columbia_movements$mainstem_river_mouth_movement[j] == 1){
      cat("borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_DE", " ~ normal(0,10);", "\n", sep = "")
      cat("borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], "_NDE", " ~ normal(0,10);", "\n", sep = "")
    }
    
    # If it's a within tributary movement - we don't want it
    else if (upper_columbia_movements$within_trib_movement[j] == 1){
      # do nothing!
    }
    # If it's a movement from mainstem to upstream - we also don't want it
    else if (upper_columbia_movements$mainstem_upstream_movement[j] == 1){
      # do nothing!
    } 
    # If we move from the river mouth state to the mainstem, we want a parameter - this is captured in the below
    # BUT - if it's a movement from the upstream state to the mainstem, we don't want a parameter.
    # these transitions will get the same parameters as the river mouth to mainstem
    else if (upper_columbia_movements$upstream_mainstem_movement[j] == 1){
      # do nothing!
    } else {
      # Finally - if it's any other movement, just treat it like normal!
      # If it's not, print just the one version
      cat("borigin", i, "_matrix_", upper_columbia_movements$row[j], "_", upper_columbia_movements$col[j], " ~ normal(0,10);", "\n", sep = "")
    }
    
  }
  cat ("\n")
}

##### LOAD and REFORMAT STATES DATA #####

# Load states complete
# states_complete <- read.csv(here::here("from_hyak_transfer", "2022-07-21-complete_det_hist", "states_complete.csv"), row.names = 1)
states_complete <- read.csv("upper_columbia_adults_states_complete.csv", row.names = 1)

##### FIRST STEP - EDIT STATES DATA - REMOVE UPSTREAM DETECTIONS IN DE YEARS #####

# first create the run year df
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2022-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2023-05-31 23:59:59"), by = "years")
run_year_numeric = seq(4, 22, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)

# First, let's join the run year df with the states complete 
# check how long this takes:
# takes about 4-5 minutes for this dataset
print(Sys.time())
states_complete %>% 
  rowwise() %>% # this is apparently crucial to getting this to work
  dplyr::mutate(date_time = ymd_hms(date_time)) %>% 
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_time & run_year_end >= date_time)$run_year) -> states_complete
print(Sys.time())

# write this to csv, to speed up code later
# write.csv(states_complete, "states_complete_run_year.csv")
# states_complete <- read.csv("states_complete_run_year.csv", row.names = 1)


# SAVE THIS AS A MODEL INPUT - to index which matrix to use
run_year_indexing <- data.frame(run_year = run_year[1:18],
                                run_year_index = seq(1,18,1))

states_complete %>% 
  left_join(., run_year_indexing, by = "run_year") -> states_complete

# Remove three fish that have detections in the 22/23 run year
run_year_2223_fish <- unique(subset(states_complete, is.na(run_year_index))$tag_code)

states_complete %>% 
  subset(., !(tag_code %in% run_year_2223_fish)) -> states_complete


# Create a data frame for when each tributary has detection efficiency capability

# Make them first each individually, then join them
deschutes_river_trib_det_eff_capability <- data.frame(state = "Deschutes River Upstream",
                                                      run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                      DE = c(rep(0,9), rep(1,6), rep(0,3)))

john_day_river_trib_det_eff_capability <- data.frame(state = "John Day River Upstream",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,8), rep(1,10)))

hood_river_trib_det_eff_capability <- data.frame(state = "Hood River Upstream",
                                                 run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                 DE = c(rep(0,8), rep(1,10)))

fifteenmile_creek_trib_det_eff_capability <- data.frame(state = "Fifteenmile Creek Upstream",
                                                        run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                        DE = c(rep(0,7), rep(1,11)))

umatilla_river_trib_det_eff_capability <- data.frame(state = "Umatilla River Upstream",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,10), rep(1,8)))

yakima_river_trib_det_eff_capability <- data.frame(state = "Yakima River Upstream",     
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   # DE = c(rep(1,18)))
                                                   # For now - remove 04/05 run year
                                                   DE = c(0, rep(1,17)))

walla_walla_river_trib_det_eff_capability <- data.frame(state = "Walla Walla River Upstream",
                                                        run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                        DE = c(rep(0,1), rep(1,17)))

wenatchee_river_trib_det_eff_capability <- data.frame(state = "Wenatchee River Upstream",
                                                      run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                      DE = c(rep(0,6), rep(1,12)))

entiat_river_trib_det_eff_capability <- data.frame(state = "Entiat River Upstream",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,3), rep(1,15)))

okanogan_river_trib_det_eff_capability <- data.frame(state = "Okanogan River Upstream",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,9), rep(1,9)))

methow_river_trib_det_eff_capability <- data.frame(state = "Methow River Upstream",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,5), rep(1,13)))

tucannon_river_trib_det_eff_capability <- data.frame(state = "Tucannon River Upstream",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,6), rep(1,12)))

asotin_creek_trib_det_eff_capability <- data.frame(state = "Asotin Creek Upstream",      
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,7), rep(1,11)))

imnaha_river_trib_det_eff_capability <- data.frame(state = "Imnaha River Upstream",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,6), rep(1,12)))


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

# join by state and run year
states_complete %>% 
  left_join(.,trib_det_eff_capability, by = c("state", "run_year")) -> states_complete

# table(subset(states_complete, DE == 1)$state)
# table(subset(states_complete, DE == 0)$state)

# check on our Snake River tribs with detection efficiency
# table(subset(states_complete, state == "Tucannon River Upstream")$run_year_index) # no individuals at upstream sites in run years without DE
# table(subset(states_complete, state == "Asotin Creek Upstream")$run_year_index) # definitely individuals in years without DE
# table(subset(states_complete, state == "Asotin Creek Upstream")$DE) # so then why are they all coming out as DE = 1?
# table(subset(states_complete, state == "Imnaha River Upstream")$run_year_index)

# remove transitions within tributaries (we don't care about these!)
# states_complete %>% 
#   # group_by(tag_code_2) %>% 
#   subset(state == "Asotin Creek Upstream" & lead(state) == "Asotin Creek Mouth" |
#            lag(state) == "Asotin Creek Upstream" & state == "Asotin Creek Mouth") -> AC_test

# Here's the root of our problem: this fish goes mouth - upstream - mouth - upstream.
# If we were to only remove upstream states, you would be left with mouth-mouth transitions, which are not allowed
# subset(states_complete, tag_code == "384.3B23B18C05")


# Solution: 1) Remove any transitions that occur between mouth and upstream or upstream and mouth; 2) remove all only upstream transitions
# that occur in years where DE = 0

# KEY STEP - in all run years in upstream states where DE = 1, remove the upstream state, and remove all transitions between mouth-upstream and upstream-mouth
states_complete_orig <- states_complete # this for testing
# states_complete <- states_complete_orig

# get all of the upstream states in a vector
upstream_states <- model_states[upstream_indices]

# get all of the river mouth states in a vector
river_mouth_states <- model_states[river_mouth_indices]

# get all mainstem states in a vector
mainstem_states <- model_states[mainstem_indices]


states_complete %>% 
  # step 1): remove states where the tag code and the next tag code are the same (therefore same fish)
  subset(., !(tag_code == lag(tag_code) & 
                # and where the next state is an upstream state and the last state is a river mouth state
                state %in% upstream_states & lag(state) %in% river_mouth_states |
                # OR where the next state is a river mouth state and the last state is an upstream state
                tag_code == lag(tag_code) & 
                state %in% river_mouth_states & lag(state) %in% upstream_states)) %>% 
  # step 1.5): for fish that go mainstem - upstream - mainstem - remove the second two states
  subset(., !(tag_code == lag(tag_code) & 
                tag_code == lead(tag_code) &
                state %in% upstream_states & 
                lag(state) %in% mainstem_states & 
                lead(state) %in% mainstem_states |
           lag(tag_code) == tag_code & 
           lag(tag_code, n = 2) == tag_code & 
           state %in% mainstem_states & 
           lag(state) %in% upstream_states &
           lag(state, 2) %in% mainstem_states)) %>% 
  # step 2): remove only upstream states in years where DE == 1
  # keep NAs and zeros, which will remove 1s
  subset(., is.na(DE) | DE == 0 ) -> states_complete

# this introduces an issue - if a fish goes straight to an upstream state and then comes back to the mainstem in a DE = 1 year,
# we end up with consecutive mainstem states because the intervening state got removed.
# Solution: See step 1.5 above. If a fish goes mainstem - upstream - mainstem, remove the upstream state and the next mainstem state

# Now - we have to check to make sure everything makes sense
states_complete %>% 
  ungroup() %>% 
  mutate(current_tag_code = tag_code_2) %>% 
  mutate(next_tag_code = lead(tag_code_2)) %>% 
  mutate(current_state = state) %>% 
  mutate(next_state = lead(state)) %>% 
  mutate(transition = paste0(current_state, " - ", next_state)) %>% 
  subset(current_tag_code == next_tag_code) -> transitions_df
  # dplyr::select(tag_code_2, state, transition, run_year) -> transitions_df

table(transitions_df$transition)

# table(subset(transitions_df, state %in% mainstem_states)$transition)


# here's a problem: mainstem, BON to MCN - mainstem, BON to MCN. Also other consecutive mainstem transitions
problem_tag_codes <-  subset(transitions_df, transition %in% c("mainstem, BON to MCN - mainstem, BON to MCN",
                                                               "mainstem, ICH to LGR - mainstem, ICH to LGR",
                                                               "mainstem, MCN to ICH or PRA - mainstem, MCN to ICH or PRA"))$tag_code_2

# subset(states_complete, tag_code_2 %in% problem_tag_codes) -> investigate
# subset(states_complete_orig, tag_code_2 %in% problem_tag_codes) -> investigate_orig
# subset(transitions_df, tag_code %in% problem_tag_codes)


##### Get values to pass as data for stan model #####

unique_tag_code_2s <- unique(states_complete$tag_code_2)

# first, get the maximum number of visits by any fish
states_complete %>% 
  # group_by(tag_code) %>% 
  # We need to use tag_code_2 because this splits apart repeat spawners
  group_by(tag_code_2) %>% 
  dplyr::count() %>% 
  as.data.frame() -> site_visits_by_tag_code

# + 1 here because we don't yet have a loss state
max_visits <- max(site_visits_by_tag_code$n) + 1

unique_tag_codes <- unique(states_complete$tag_code_2)
nfish <- length(unique_tag_codes)

# second, create an array with the dimensions number of fish, maximum number of site visits, number of states
state_data <- array(data = 0, dim = c(nstates, max_visits, nfish))

# third, populate the state data

# Add a column for the observation
states_complete %>% 
  # group_by(tag_code) %>% 
  group_by(tag_code_2) %>% 
  dplyr::mutate(order = row_number()) -> states_complete

# add a column to index by fish
states_complete %>% 
  # group_by(tag_code) %>% 
  group_by(tag_code_2) %>% 
  dplyr::mutate(tag_code_number = cur_group_id()) -> states_complete


##### Prepare final states data for model #####

# Loop through to convert df to array
# This takes about 30 seconds
print(Sys.time())
for (i in 1:nrow(states_complete)){
# for (i in 1:100){
  # We need to skip all upstream detections in DE years
  # index by 1) which state, 2) which visit, 3) which fish
  state_data[which(model_states == states_complete[i, "state", drop = TRUE]),states_complete[i, "order", drop = TRUE], states_complete[i, "tag_code_number", drop = TRUE]] <- 1 
  
}
print(Sys.time())

# Now, add the loss state
# First, get the number of non-loss states
n_not_loss_states <- vector(length = nfish)
# Now, get the indices of the fish that already have a loss state
states_complete %>% 
  subset(state == "loss") -> trapped_fish

trapped_fish_tags <- trapped_fish$tag_code_2
# Get the indices of these
which(unique_tag_codes %in% trapped_fish_tags) -> trapped_fish_indices

for (i in 1:nfish){
  # 43 is the loss state
  n_not_loss_states[i] <- sum(state_data[,,i])
  state_data[43, n_not_loss_states[i]+1,i] <- 1
}

# remove the extra loss state for the trapped fish
for (i in 1:length(trapped_fish_indices)){
  j <- trapped_fish_indices[i]
  state_data[43, n_not_loss_states[j]+1,j] <-0
}

saveRDS(state_data, "upper_columbia_state_data.csv")

##### Load and reformat tributary data for detection efficiency #####

# Load tributary discharge data
# tributary_discharge_data <- read.csv("tributary_discharge_data.csv")
tributary_discharge_data <- read.csv("tributary_discharge_data_zscore.csv")

# Get tributary categorical data (equipment eras)
# this by run year
# we will have multiple eras, then NAs for run years without an ability to calculate detection efficiency

# Create a design matrix for all tributary data
# first columns = intercepts for eras
# last columns = discharge values
# rows = run years

# This creates the master file; we then subset it to make a design matrix for each tributary, where we replace every value except
# the relevant columns with zeros. We will then store this in one big array; each slice of the array will be one tributary,
# then rows will be run years, and columns will be the design matrix for the parameters.
# Then, in the actual model code, we will use the tributary a fish was detected in to select a slice of the array;
# then that slice (which is the design matrix for that tributary) will then be indexed by the run year to get a row vector.
# The row vector will be multiplied by the full parameter vector (all alphas and betas) to get the eta (linear predictor) for the detection efficiency GLM





# initialize the matrix, with dimensions rows = run years, columns = N eras (for alphas) + N tributaries (for betas)
# NOTE: for the Imnaha and for Fifteenmile Creek, we will not have a relationship with discharge. So we will only have an intercept, and put
# in fake data (all zeros, so that the slope term doesn't matter) for the discharge data in order to allow it to run.

# Tributaries that we have no capability to estimate detection capability: Clearwater River, Grande Ronde River, Salmon River

# Note that we are choosing the start run year depending on when in the year the site came online (e.g. 12/2011 would be
# 11/12 start year). If it's active in the spring of that run year then we accept it

# Here's how many eras we have (do it alphabetically)
# Asotin Creek 1: 11/12 - 17/18
# Asotin Creek 2: 18/19 - 21/22
# Deschutes River 1: 13/14 - 18/19
# Entiat River 1: 07/08 - 21/22
# Fifteenmile Creek: 11/12 - 18/19; river mouth site still active after this point but no upstream sites
# Hood River 1: 12/13 - 21/22
# Imnaha River 1: 10/11 - 21/22
# John Day River 1: 12/13 - 21/22
# Methow River 1: 09/10 - 16/17
# Methow River 2: 17/18 - 21/22
# Okanogan River 1: 13/14 - 21/22
# Tucannon River 1: 10/11 - 19/20
# Tucannon River 2: 20/21 - 21/22
# Umatilla River 1: 06/07 - 13/14
# Umatilla River 2: 14/15 - 21/22 (see email from Stacy Remple for info on this - not in tributary detection efficiency RMD)
# Walla Walla River 1 (ORB): 05/06 - 14/15
# Walla Walla River 2 (PRV): 12/13 - 18/19
# Walla Walla River 3 (WWB): 19/20 - 21/22
# We're going to have to deal with this overlap in the Walla Walla; since PRV is 1 RKM downstream (closer to mouth) of 
# ORB, we can call ORB in years of overlap as an upstream site
# Wenatchee River 1: 10/11 (starts January 2011) - 21/22
# Yakima River 1: 04/05 - 21/22

# Total eras: 20
# Total tribs: 14 (17 total - Clearwater, Grande Ronde, and Salmon)

# So the first 20 columns are the era/intercept terms, then the next 14 are the slope terms for discharge

# We then have 18 run years (04/05 through 21/22)

# Our parameter vector will also be a column vector of length 34

tributary_design_matrix <- matrix(0, nrow = 18, ncol = 34)

# This is to make it easier to tell indexing
rownames(tributary_design_matrix) <- paste0(run_year_df$run_year[1:(nrow(run_year_df)-1)], "-", seq(1,18,1))

# Now let's populate it with 1s/discharge values

# (1) Asotin Creek 1: 11/12 - 17/18
tributary_design_matrix[8:14,1] <- 1

# (2) Asotin Creek 2: 18/19 - 21/22
tributary_design_matrix[15:18,2] <- 1

# (3) Deschutes River 1: 13/14 - 18/19
tributary_design_matrix[10:15,3] <- 1

# (4) Entiat River 1: 07/08 - 21/22
tributary_design_matrix[4:18,4] <- 1

# (5) Fifteenmile Creek: 11/12 - 18/19; river mouth site still active after this point but no upstream sites
# tributary_design_matrix[8:15,5] <- 1
tributary_design_matrix[8:18,5] <- 1

# (6) Hood River 1: 12/13 - 21/22
tributary_design_matrix[9:18,6] <- 1

# (7) Imnaha River 1: 10/11 - 21/22
tributary_design_matrix[7:18,7] <- 1

# (8) John Day River 1: 12/13 - 21/22
tributary_design_matrix[9:18,8] <- 1

# (9) Methow River 1: 09/10 - 16/17
tributary_design_matrix[6:13,9] <- 1

# (10) Methow River 2: 17/18 - 21/22
tributary_design_matrix[14:18,10] <- 1

# (11) Okanogan River 1: 13/14 - 21/22
tributary_design_matrix[10:18,11] <- 1

# (12) Tucannon River 1: 10/11 - 19/20
tributary_design_matrix[7:16,12] <- 1

# (13) Tucannon River 2: 20/21 - 21/22
tributary_design_matrix[17:18,13] <- 1

# (14) Umatilla River 1: 06/07 - 13/14
tributary_design_matrix[3:10,14] <- 1

# (15) Umatilla River 2: 14/15 - 21/22 (see email from Stacy Remple for info on this - not in tributary detection efficiency RMD or on PTAGIS, but there was an antenna installation in 2014)
tributary_design_matrix[11:18,15] <- 1

# (16) Walla Walla River 1 (ORB): 05/06 - 11/12 (this site continues to 14/15, but in 12/13 the PRV site comes online which is closer to the river mouth and also seems to have better detection efficiency)
tributary_design_matrix[2:8,16] <- 1

# (17) Walla Walla River 2 (PRV): 12/13 - 18/19
tributary_design_matrix[9:15,17] <- 1

# (18) Walla Walla River 3 (WWB): 19/20 - 21/22
tributary_design_matrix[16:18,18] <- 1

# (19) Wenatchee River 1: 10/11 (starts January 2011) - 21/22
tributary_design_matrix[7:18,19] <- 1

# (20) Yakima River 1: 04/05 - 21/22
# tributary_design_matrix[1:18,20] <- 1
# for now - remove 04/05 run year
tributary_design_matrix[2:18,20] <- 1


# Now, add discharge values in the remaining columns
tributary_design_matrix[2:18,21] <-subset(tributary_discharge_data, tributary == "Asotin Creek")$mean_discharge_zscore

tributary_design_matrix[2:18,22] <-subset(tributary_discharge_data, tributary == "Deschutes River")$mean_discharge_zscore

tributary_design_matrix[2:18,23] <-subset(tributary_discharge_data, tributary == "Entiat River")$mean_discharge_zscore

# tributary_design_matrix[2:18,24] <-subset(tributary_discharge_data, tributary == "Fifteenmile Creek")$mean_discharge_zscore
# No discharge data for Fifteenmile Creek - keep as zeros, since we don't have a complete time series

tributary_design_matrix[2:18,25] <-subset(tributary_discharge_data, tributary == "Hood River")$mean_discharge_zscore

# tributary_design_matrix[2:10,26] <-subset(tributary_discharge_data, tributary == "Imnaha River")$mean_discharge_zscore
# Indexing is different for the Imnaha because we only have discharge data through 12/13 (13/14 data is not from the complete run year)
# Keep Imnaha as all zeros - we only have one year (12/13) with overlap between discharge data and detection capability

tributary_design_matrix[2:18,27] <-subset(tributary_discharge_data, tributary == "John Day River")$mean_discharge_zscore

tributary_design_matrix[2:18,28] <-subset(tributary_discharge_data, tributary == "Methow River")$mean_discharge_zscore

tributary_design_matrix[2:18,29] <-subset(tributary_discharge_data, tributary == "Okanogan River")$mean_discharge_zscore

tributary_design_matrix[2:18,30] <-subset(tributary_discharge_data, tributary == "Tucannon River")$mean_discharge_zscore

tributary_design_matrix[2:18,31] <-subset(tributary_discharge_data, tributary == "Umatilla River")$mean_discharge_zscore

tributary_design_matrix[2:18,32] <-subset(tributary_discharge_data, tributary == "Walla Walla River")$mean_discharge_zscore

tributary_design_matrix[2:18,33] <-subset(tributary_discharge_data, tributary == "Wenatchee River")$mean_discharge_zscore

tributary_design_matrix[2:18,34] <-subset(tributary_discharge_data, tributary == "Yakima River")$mean_discharge_zscore




# change row names just to run years
rownames(tributary_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]

# Remove the 04/05 run year
# tributary_design_matrix <- tributary_design_matrix[2:18,]
# Don't do this - just leave it in for indexing purposes

# Take the master tributary design matrix, and create design matrices for each of the tributaries individually
asotin_creek_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(asotin_creek_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
asotin_creek_design_matrix[,1] <- tributary_design_matrix[,1]
asotin_creek_design_matrix[,2] <- tributary_design_matrix[,2]
asotin_creek_design_matrix[,21] <- tributary_design_matrix[,21]

deschutes_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(deschutes_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
deschutes_river_design_matrix[,3] <- tributary_design_matrix[,3]
deschutes_river_design_matrix[,22] <- tributary_design_matrix[,22]

entiat_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(entiat_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
entiat_river_design_matrix[,4] <- tributary_design_matrix[,4]
entiat_river_design_matrix[,23] <- tributary_design_matrix[,23]

fifteenmile_creek_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(fifteenmile_creek_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
fifteenmile_creek_design_matrix[,5] <- tributary_design_matrix[,5]
fifteenmile_creek_design_matrix[,24] <- tributary_design_matrix[,24]

hood_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(hood_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
hood_river_design_matrix[,6] <- tributary_design_matrix[,6]
hood_river_design_matrix[,25] <- tributary_design_matrix[,25]

imnaha_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(imnaha_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
imnaha_river_design_matrix[,7] <- tributary_design_matrix[,7]
imnaha_river_design_matrix[,26] <- tributary_design_matrix[,26]

john_day_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(john_day_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
john_day_river_design_matrix[,8] <- tributary_design_matrix[,8]
john_day_river_design_matrix[,27] <- tributary_design_matrix[,27]

methow_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(methow_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
methow_river_design_matrix[,9] <- tributary_design_matrix[,9]
methow_river_design_matrix[,10] <- tributary_design_matrix[,10]
methow_river_design_matrix[,28] <- tributary_design_matrix[,28]

okanogan_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(okanogan_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
okanogan_river_design_matrix[,11] <- tributary_design_matrix[,11]
okanogan_river_design_matrix[,29] <- tributary_design_matrix[,29]

tucannon_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(tucannon_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
tucannon_river_design_matrix[,12] <- tributary_design_matrix[,12]
tucannon_river_design_matrix[,13] <- tributary_design_matrix[,13]
tucannon_river_design_matrix[,30] <- tributary_design_matrix[,30]

umatilla_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(umatilla_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
umatilla_river_design_matrix[,14] <- tributary_design_matrix[,14]
umatilla_river_design_matrix[,15] <- tributary_design_matrix[,15]
umatilla_river_design_matrix[,31] <- tributary_design_matrix[,31]

walla_walla_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(walla_walla_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
walla_walla_river_design_matrix[,16] <- tributary_design_matrix[,16]
walla_walla_river_design_matrix[,17] <- tributary_design_matrix[,17]
walla_walla_river_design_matrix[,18] <- tributary_design_matrix[,18]
walla_walla_river_design_matrix[,32] <- tributary_design_matrix[,32]

wenatchee_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(wenatchee_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
wenatchee_river_design_matrix[,19] <- tributary_design_matrix[,19]
wenatchee_river_design_matrix[,33] <- tributary_design_matrix[,33]

yakima_river_design_matrix <- matrix(0, nrow = 18, ncol = 34)
rownames(yakima_river_design_matrix) <- run_year_df$run_year[1:(nrow(run_year_df)-1)]
yakima_river_design_matrix[,20] <- tributary_design_matrix[,20]
yakima_river_design_matrix[,34] <- tributary_design_matrix[,34]

# Now take these matrices and store them in one big array
# indices: rows = run years, columns = parameters/design matrix, slices = tributaries
# Note: this will have the same number of slices as possible state transitions (42 - which is all transitions minus loss).
# But only 14 of these slices (corresponding to the 14 tributaries that have detection efficiency capability)
# will have any non-zero values.
tributary_design_matrices_array <- array(0, dim = c(18,34,42))

# Note that these are indices of only the mouth states (since these are the states that we are calculating detection efficiency at)
tributary_design_matrices_array[,,34] <- asotin_creek_design_matrix
tributary_design_matrices_array[,,10] <- deschutes_river_design_matrix
tributary_design_matrices_array[,,26] <- entiat_river_design_matrix
tributary_design_matrices_array[,,16] <- fifteenmile_creek_design_matrix
tributary_design_matrices_array[,,14] <- hood_river_design_matrix
tributary_design_matrices_array[,,39] <- imnaha_river_design_matrix
tributary_design_matrices_array[,,12] <- john_day_river_design_matrix
tributary_design_matrices_array[,,30] <- methow_river_design_matrix
tributary_design_matrices_array[,,28] <- okanogan_river_design_matrix
tributary_design_matrices_array[,,32] <- tucannon_river_design_matrix
tributary_design_matrices_array[,,18] <- umatilla_river_design_matrix
tributary_design_matrices_array[,,22] <- walla_walla_river_design_matrix
tributary_design_matrices_array[,,24] <- wenatchee_river_design_matrix
tributary_design_matrices_array[,,20] <- yakima_river_design_matrix

# Create a vector for the order of these tributaries
trib_det_eff_order <- c("Asotin_Creek", 
"Deschutes_River", 
"Entiat_River", 
"Fifteenmile_Creek", 
"Hood_River",
"Imnaha_River",
"John_Day_River", 
"Methow_River", 
"Okanogan_River", 
"Tucannon_River", 
"Umatilla_River",
"Walla_Walla_River",
"Wenatchee_River", 
"Yakima_River")


##### Convert the dates into a numeric index #####
# for now, comment this whole block out

# 
# # First, populate the dates for implicit site visits
# # Create a column to indicate if the date_time was interpolated
# # Note the arrival date
# states_complete %>% 
#   mutate(date_source = ifelse(!is.na(date_time), "known_arrival", "interpolated")) %>% 
#   mutate(date = NA) -> states_complete
# 
# # I think we have to loop this
# 
# # First, figure out indices of missing date_time
# missing_date_time <- is.na(states_complete$date_time)
# date_times <- states_complete$date_time
# 
# # For efficiency, only loop through the missing dates
# missing_date_time_indices <- which(missing_date_time == TRUE)
# 
# # Extract date for all known date times
# states_complete %>% 
#   mutate(date =  format(as.Date(date(states_complete[i,"date_time", drop = TRUE]), origin = "1970-01-01"))) -> states_complete
# 
# # Loop to get arrival dates in state
# # We currently don't need this, because we interpolated the time earlier
# # takes about 90 minutes to run
# # for (j in 1:length(missing_date_time_indices)){
# # # for (i in 1:100){
# #   i <- missing_date_time_indices[j]
# # 
# #     # Figure out how far back the last known time was
# #     # Truncate the missing date times to only ones prior to current state
# #     missing_date_time_subset <- missing_date_time[1:(i-1)]
# #     prev_time_index <- max(which(missing_date_time_subset %in% FALSE))
# #     
# #     # Truncate the missing date times to only ones after current state
# #     missing_date_time_subset <- missing_date_time[(i+1):nrow(states_complete)]
# #     next_time_index <- min(which(missing_date_time_subset %in% FALSE)) + i
# #     
# #     # Now, interpolate the missing time
# #     # First, figure out how long between the two known times
# #     prev_time <- ymd_hms(date_times[prev_time_index])
# #     next_time <- ymd_hms(date_times[next_time_index])
# #     time_diff <- next_time - prev_time
# #     
# #     # Get the missing time - add the time difference divided by the number 
# #     # of missing steps plus 1, multiply by which number step it is
# #     missing_time <- prev_time + (time_diff/(next_time_index - prev_time_index) * (i - prev_time_index))
# #     
# #     # Extract just the date
# #     missing_date <- date(missing_time)
# #     
# #     # populate the missing date_time
# #     states_complete[i, "date"] <- format(as.Date(missing_date, origin = "1970-01-01"))
# # }
#   
# 
# # create an empty matrix
# transition_date_matrix <- matrix(data = NA, nrow = nfish, ncol = max_visits)
# # populate it with the dates
# for (i in 1:nrow(states_complete)){
#   transition_date_matrix[states_complete[i,"tag_code_number", drop = TRUE], states_complete[i, "order", drop = TRUE]] <- states_complete[i, "date_time", drop = TRUE]
# }
# 
# 
# # We'll use 2005-05-31 as day 0, so that day 1 is 2005-06-01, which is the first day in run year 05/06 (first run year in our dataset)
# date_numeric <- function(x, na.rm = FALSE) as.numeric(ymd(date(x)) - ymd("2005-05-31"))
# transition_date_matrix %>% 
#   as_tibble() %>% 
#   mutate_all(date_numeric) -> transition_date_numeric


##### Reformat origin and rear information #####
tag_code_metadata <- read.csv("tag_code_metadata.csv")
# keep only the fish that are in the dataset
tag_code_metadata <- subset(tag_code_metadata, tag_code %in% states_complete$tag_code)

# Convert origins into numbers
# All relative to Yakima, since it's the last one alphabetically
origin_numeric <- data.frame(natal_origin = c("Asotin_Creek", 
                                        "Clearwater_River",
                                        "Deschutes_River", 
                                        "Entiat_River", 
                                        "Fifteenmile_Creek", 
                                        "Grande_Ronde_River", 
                                        "Hood_River",
                                        "Imnaha_River",
                                        "John_Day_River", 
                                        "Methow_River", 
                                        "Okanogan_River", 
                                        "Salmon_River", 
                                        "Tucannon_River", 
                                        "Umatilla_River",
                                        "Walla_Walla_River",
                                        "Wenatchee_River", 
                                        "Yakima_River"),
                             natal_origin_numeric = seq(1,17,1))

natal_origin_table <- read.csv("natal_origin_table.csv")
tag_code_metadata %>% 
  left_join(., natal_origin_table, by = "release_site_name") %>% 
  left_join(., origin_numeric, by = "natal_origin") %>% 
  mutate(rear_type_numeric = ifelse(rear_type_code %in% c("H", "U"), 2, 1))-> tag_code_metadata

# At this point, we need to recreate tag_code_metadata but with the tag_code_2 field
states_complete %>% 
  distinct(tag_code_2, .keep_all = TRUE) %>% 
  dplyr::select(tag_code, tag_code_2) -> tag_codes_2

# reformat this into origin_rear info
tag_codes_2 %>%
  left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin_numeric, rear_type_numeric), by = "tag_code") %>% 
  dplyr::rename(natal_origin = natal_origin_numeric, rear_type = rear_type_numeric) %>% 
  dplyr::select(-tag_code) -> origin_rear_actual

write.csv(origin_rear_actual,"upper_columbia_origin_rear_actual.csv")
  

fish_sim_cat_data_actual <- origin_rear_actual
  
  
# Store quantities for loop
# Store the total number of individuals
n.ind <- dim(state_data)[3]
  
  # Store the number of observations per individual
  # -1 because the last state is loss, which isn't actually an observation
  n.obs <- vector(length = n.ind)
  for (i in 1:n.ind){
    n.obs[i] <- sum(state_data[,,i]) - 1
  }
  
  
  
  # Get the state that each fish was in at each n.obs
  
  # First initialize an empty list
  states_list <- list()
  
  for (i in 1:n.ind){
    vec <- vector(length = (n.obs[i]-1))
    
    # Store in list
    states_list[[i]] <- vec
  }
  
  # Store with the states
  for (i in 1:n.ind){
    for (j in 1:(n.obs[i])){
      # states_list[[i]][j] <- rownames(as.data.frame(which(state_data[[i]][,j] == 1)))
      states_list[[i]][j] <- which(state_data[,j,i] == 1) # Get the index of the site instead of the name
    }
  }
  
  # Turn into matrix for stan
  states_mat <- matrix(nrow = n.ind, ncol = max(n.obs))
  for (i in 1:n.ind){
    states_mat[i,1:(n.obs[i])] <- states_list[[i]]
  }
  
  
  # Create the design matrix for categorical variables
  # new for detection efficiency: add a run year field to allow for glm calculation
  # of detection efficiency, as well as to note when certain movements are not possible
  cat_X_mat_actual <- matrix(0, nrow = n.ind, ncol = 5)
  # Start it so that they're all 0s
  # The first column everyone gets a 1 (this is beta 0/grand mean mu)
  cat_X_mat_actual[,1] <- 1
  
  # This is for origin + rear + run year
  for (i in 1:n.ind){
    # Rear type
    if (fish_sim_cat_data_actual$rear_type[i] == 1){
      cat_X_mat_actual[i,2] <- 1
    }
    else {
      cat_X_mat_actual[i,2] <- -1
    }
    # Natal origin
    if (fish_sim_cat_data_actual$natal_origin[i] == 16){ # Wenatchee River
      cat_X_mat_actual[i,3] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == 4){ # Entiat River
      cat_X_mat_actual[i,4] <- 1
    }
    else if (fish_sim_cat_data_actual$natal_origin[i] == 11){ # Okanogan River
      cat_X_mat_actual[i,5] <- 1
    }
    else { # for Methow River
      cat_X_mat_actual[i,3:5] <- -1
    }
  }
  # "Tucannon_River", "Asotin_Creek", "Clearwater_River", "Salmon_River", "Grande_Ronde_River", "Imnaha_River"
  # One tweak: We have to replace all NAs in our input data for stan to accept it
  states_mat[is.na(states_mat)] <- -999
  
  
  # We are going to transform our detection histories to instead be vectors containing
  # the number of each site that was visited, instead of a matrix of sites with 0s for
  # not visited and 1s for visited
  
  # First, create an empty matrix to store the newly formatted detection histories
  state_data_2 <- matrix(0, nrow = dim(state_data)[3], ncol = dim(state_data)[2])
  
  # Now fill it in
  for (i in 1:dim(state_data)[3]){
    det_hist <- state_data[,,i]
    
    # Count the number of site visits
    nsite_visits <- sum(det_hist)
    
    for (j in 1:nsite_visits){
      state_data_2[i,j] <- which(det_hist[,j] == 1, arr.ind = TRUE)
    }
    
  }
    
  
  # Get vector of run years of fish
  # need the -1 correction for indices because we're ignoring the 04/05 run year
  fish_run_years <- vector(length = n.ind)
  
  states_complete %>% 
    distinct(tag_code_2, .keep_all = TRUE) %>% 
    mutate(run_year = subset(run_year_df, run_year_start <= date_time & run_year_end >= date_time)$run_year) %>% 
    dplyr::select(tag_code_2, run_year) -> tag_code_2_run_years
  
  for (i in 1:n.ind){
    fish_run_years[i] <- which(run_year_df$run_year == tag_code_2_run_years$run_year[i])-1
  }
  
  # Get the state numbers that are mainstem states that connect to tributaries for which we can estimate detection efficiency
  mainstem_trib_states <- c(2,3,5,6,7,8,9)
  
  # Now create a list of vectors of which tributaries connect to each of these
  # To make 
  mainstem_trib_connections <- list(c(2,4,5,7,11),
                                    c(12,14,0,0,0),
                                    c(13,0,0,0,0),
                                    c(3,0,0,0,0),
                                    c(9,8,0,0,0),
                                    c(10,0,0,0,0),
                                    c(1,6,0,0,0))
  
  # Now, figure out how many detection efficiencies you need to calculate for each mainstem stem
  # make this the length of number of states for ease of indexing
  n_detection_efficiencies <- rep(0, 43)
  n_detection_efficiencies[2] <- 5
  n_detection_efficiencies[3] <- 2
  n_detection_efficiencies[5] <- 1
  n_detection_efficiencies[6] <- 1
  n_detection_efficiencies[7] <- 2
  n_detection_efficiencies[8] <- 1
  n_detection_efficiencies[9] <- 2
  
  ##### Create the run year/detection efficiency array #####
  # we need an array that we can index - current (from, rows) x k (to, columns) x run years (slices) - and the values
  # indicate which of the two matrices to use
  # 0 means no detection efficiency correction; 1 means use detection efficiency matrix to pull parameter value from
  # 18 run years: 04/05 through 21/22 (this is the same ad the detection efficiency GLM indexing)
  
  # Also note - we only need to see transitions from mainstem to river mouth
  run_year_DE_array <- array(0, dim = c(nstates, nstates, 18))
  
  # figure out what the indices are of the various tributaries with detection efficiency estimation capability
  print(paste0(seq(1,43,1), " - ", model_states))
  
  # Asotin Creek
  # From Upstream LGR (9) to mouth (34) or upstream (35)
  run_year_DE_array[9,34,8:18] <- 1
  
  # Deschutes River
  # From mainstem BON to MCN (2) to mouth (10) or upstream (11)
  run_year_DE_array[2,10,10:15] <- 1
  
  # Entiat River
  # From RRE to WEL (6) to mouth (26) or upstream (27)
  run_year_DE_array[6,26,4:18] <- 1
  
  # Fifteenmile Creek
  # From BON to MCN (2) to mouth (16) or upstream (17)
  # Note here that in years 16:18, we no longer have upstream sites, but still have a river mouth site
  run_year_DE_array[2,16,8:18] <- 1
  
  # Hood River
  # From BON to MCN (2) to mouth (14) or upstream (15)
  run_year_DE_array[2,14,9:18] <- 1
  
  # Imnaha River
  # From Upstream LGR (9) to mouth (39) or upstream (40)
  run_year_DE_array[9,39,7:18] <- 1
  
  # John Day River
  # From BON to MCN (2) to mouth (12) or upstream (13)
  run_year_DE_array[2,12,9:18] <- 1
  
  # Methow River
  # From Upstream WEL (7) to mouth (30) or upstream (31)
  run_year_DE_array[7,30,6:18] <- 1
  
  # Okanogan River
  # From Upstream WEL (7) to mouth (28) or upstream (29)
  run_year_DE_array[7,28,10:18] <- 1
  
  # Tucannon River
  # From ICH to LGR (8) to mouth (32) or upstream (33)
  run_year_DE_array[8,32,7:18] <- 1
  
  # Umatilla River
  # From BON to MCN (2) to mouth (18) or upstream (19)
  run_year_DE_array[2,18,3:18] <- 1
  
  # Walla Walla River
  # From MCN to ICH or PRA (3) to mouth (22) or upstream (23)
  run_year_DE_array[3,22,2:18] <- 1
  
  # Wenatchee River
  # From RIS to RRE (5) to mouth (24) or upstream (25)
  run_year_DE_array[5,24,7:18] <- 1
  
  # Yakima River
  # From MCN to ICH or PRA (3) to mouth (20) or upstream (21)
  # run_year_DE_array[3,20,1:18] <- 1
  # so for now - remove run year 04/05, otherwise it'll all get confused because there's no discharge data for that year
  run_year_DE_array[3,20,2:18] <- 1
  
  
  
  
  # Load the data from the detection efficiency model to use as priors in this model
  det_eff_param_posteriors <- as.matrix(read.csv("det_eff_param_posteriors.csv", row.names = 1))
  
  
  # get a vector of run years that transitions occurred within - need to do this at the end, once you've made all changes to states_complete
  transition_run_years <- states_complete$run_year_index
  
  ntransitions = length(transition_run_years)
  
  
  ##### Run stan model #####
  
  # step 0: data in a list #
  data <- list(y = state_data_2, n_ind = n.ind, n_obs = n.obs, possible_movements = possible_movements,
               states_mat = states_mat, max_visits = dim(state_data_2)[2],
               movements = movements, not_movements = not_movements,
               nmovements = nmovements, # dates = dates,
               n_notmovements = n_notmovements, possible_states = transition_matrix, cat_X_mat = cat_X_mat_actual,
               grainsize = 1, N = dim(state_data_2)[1],
               # New data for detection efficiency
               tributary_design_matrices_array = tributary_design_matrices_array,
               ntransitions = ntransitions,
               transition_run_years = transition_run_years,
               run_year_DE_array = run_year_DE_array,
               det_eff_param_posteriors = det_eff_param_posteriors)
  
  
  print(Sys.time())
  
  # Fit stan model
  
  # fit <- stan(file = '01_stan_sim_int_only.stan', data = data)
  
  # Fit stan model using cmdstan
  # Step 1: load the model
  # mod <- cmdstan_model("01_stan_sim_int_only.stan", compile = FALSE)
  mod <- cmdstan_model("parallel_upper_columbia_03_stan_actual_int_origin_deteff.stan", compile = FALSE)
  
  print("Compiling model")
  # Step 2: Compile the model, set up to run in parallel
  mod$compile(cpp_options = list(stan_threads = TRUE))
  
  print("model compiled at")
  print(Sys.time())
  
  # Step 3: Run MCMC (HMC)
  fit <- mod$sample(
    data = data, 
    seed = 101, 
    # seed = 456,
    # chains = 3, 
    chains = 1,
    parallel_chains = 1,
    # parallel_chains = 3,
    refresh = 10, # print update every 10 iter
    # iter_sampling = 1000,
    # iter_warmup = 1000,
    iter_warmup = 500,
    iter_sampling = 500,
    # adapt_delta = 0.95,
    init = 0,
    threads_per_chain = 28
  )
  
# saveRDS(fit, "100iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds")
fit$save_object(file = "seed101_500iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds")

# Troubleshoot our data
# Check to see if every transition in our model is represented
model_states %>% 
  as.data.frame() %>% 
  dplyr::rename(state = ".") %>% 
  mutate(index = row_number()) -> model_states_df

states_complete %>% 
  left_join(., model_states_df, by = "state") %>% 
  mutate(transition = ifelse(tag_code_2 == lag(tag_code_2), paste0(lag(state), " - ", state), NA)) %>% 
  mutate(transition_numeric = ifelse(tag_code_2 == lag(tag_code_2), paste0(lag(index), " - ", index), NA)) -> state_transitions

table(state_transitions$transition)
table(state_transitions$transition_numeric)
setdiff(model_states, unique(state_transitions$state))
setdiff(unique(state_transitions$state), model_states)
