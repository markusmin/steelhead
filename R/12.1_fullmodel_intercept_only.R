# 12.1 FULL MODEL

# For testing:


# In this initial fit to the full dataset, we will be fitting an intercept-only model.

# This script will write the JAGS code and run the model.

# Load libraries
library(tidyverse)
library(lubridate)
library(janitor)
library(rjags)
library(jagsUI)
library(R2jags)
library(coda)


##### Get model data inputs setup #####


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
  
  # Tributary sites (17)
  "Deschutes River",
  "John Day River",
  "Hood River",
  "Fifteenmile Creek",
  "Umatilla River",
  "Yakima River",
  "Walla Walla River",
  "Wenatchee River",
  "Entiat River",
  "Okanogan River",
  "Methow River",
  "Tucannon River",
  "Asotin Creek",
  "Clearwater River",
  "Salmon River",
  "Grande Ronde River",
  "Imnaha River",
  "BON to MCN other trib",
  "Upstream WEL other trib",
  
  # Loss
  "loss"
)

# 28 states, plus loss
nstates <- length(model_states)


##### Create the transition matrix #####
# Need to update this with the additional states

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
transition_matrix["mainstem, BON to MCN", "Deschutes River"] <- 1
transition_matrix["mainstem, BON to MCN", "John Day River"] <- 1
transition_matrix["mainstem, BON to MCN", "Hood River"] <- 1
transition_matrix["mainstem, BON to MCN", "Fifteenmile Creek"] <- 1
transition_matrix["mainstem, BON to MCN", "Umatilla River"] <- 1


# 3: mainstem, MCN to ICH or PRA
transition_matrix["mainstem, MCN to ICH or PRA", "loss"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, BON to MCN"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Yakima River"] <- 1
transition_matrix["mainstem, MCN to ICH or PRA", "Walla Walla River"] <- 1


# 4: mainstem, PRA to RIS
transition_matrix["mainstem, PRA to RIS", "loss"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, PRA to RIS", "mainstem, RIS to RRE"] <- 1


# 5: mainstem, RIS to RRE
transition_matrix["mainstem, RIS to RRE", "loss"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, PRA to RIS"] <- 1
transition_matrix["mainstem, RIS to RRE", "mainstem, RRE to WEL"] <- 1
transition_matrix["mainstem, RIS to RRE", "Wenatchee River"] <- 1


# 6: mainstem, RRE to WEL
transition_matrix["mainstem, RRE to WEL", "loss"] <- 1
transition_matrix["mainstem, RRE to WEL", "mainstem, RIS to RRE"] <- 1
transition_matrix["mainstem, RRE to WEL", "mainstem, upstream of WEL"] <- 1
transition_matrix["mainstem, RRE to WEL", "Entiat River"] <- 1


# 7: mainstem, upstream of WEL
transition_matrix["mainstem, upstream of WEL", "loss"] <- 1
transition_matrix["mainstem, upstream of WEL", "mainstem, RRE to WEL"] <- 1
transition_matrix["mainstem, upstream of WEL", "Okanogan River"] <- 1
transition_matrix["mainstem, upstream of WEL", "Methow River"] <- 1


# 8: mainstem, ICH to LGR
transition_matrix["mainstem, ICH to LGR", "loss"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, MCN to ICH or PRA"] <- 1
transition_matrix["mainstem, ICH to LGR", "mainstem, upstream of LGR"] <- 1
transition_matrix["mainstem, ICH to LGR", "Tucannon River"] <- 1


# 9: mainstem, upstream of LGR
transition_matrix["mainstem, upstream of LGR", "loss"] <- 1
transition_matrix["mainstem, upstream of LGR", "mainstem, ICH to LGR"] <- 1
transition_matrix["mainstem, upstream of LGR", "Asotin Creek"] <- 1
transition_matrix["mainstem, upstream of LGR", "Clearwater River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Salmon River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Grande Ronde River"] <- 1
transition_matrix["mainstem, upstream of LGR", "Imnaha River"] <- 1


# 10: Deschutes River
transition_matrix["Deschutes River", "loss"] <- 1
transition_matrix["Deschutes River", "mainstem, BON to MCN"] <- 1


# 11: John Day River
transition_matrix["John Day River", "loss"] <- 1
transition_matrix["John Day River", "mainstem, BON to MCN"] <- 1


# 12: Hood River
transition_matrix["Hood River", "loss"] <- 1
transition_matrix["Hood River", "mainstem, BON to MCN"] <- 1


# 13: Fifteenmile Creek
transition_matrix["Fifteenmile Creek", "loss"] <- 1
transition_matrix["Fifteenmile Creek", "mainstem, BON to MCN"] <- 1


# 14: Umatilla River
transition_matrix["Umatilla River", "loss"] <- 1
transition_matrix["Umatilla River", "mainstem, BON to MCN"] <- 1


# 15: Yakima River
transition_matrix["Yakima River", "loss"] <- 1
transition_matrix["Yakima River", "mainstem, MCN to ICH or PRA"] <- 1


# 16: Walla Walla River
transition_matrix["Walla Walla River", "loss"] <- 1
transition_matrix["Walla Walla River", "mainstem, MCN to ICH or PRA"] <- 1


# 17: Wenatchee River
transition_matrix["Wenatchee River", "loss"] <- 1
transition_matrix["Wenatchee River", "mainstem, RIS to RRE"] <- 1


# 18: Entiat River
transition_matrix["Entiat River", "loss"] <- 1
transition_matrix["Entiat River", "mainstem, RRE to WEL"] <- 1


# 19: Okanogan River
transition_matrix["Okanogan River", "loss"] <- 1
transition_matrix["Okanogan River", "mainstem, upstream of WEL"] <- 1


# 20: Methow River
transition_matrix["Methow River", "loss"] <- 1
transition_matrix["Methow River", "mainstem, upstream of WEL"] <- 1


# 21: Tucannon River
transition_matrix["Tucannon River", "loss"] <- 1
transition_matrix["Tucannon River", "mainstem, ICH to LGR"] <- 1


# 22: Asotin Creek
transition_matrix["Asotin Creek", "loss"] <- 1
transition_matrix["Asotin Creek", "mainstem, upstream of LGR"] <- 1


# 23: Clearwater River
transition_matrix["Clearwater River", "loss"] <- 1
transition_matrix["Clearwater River", "mainstem, upstream of LGR"] <- 1


# 24: Salmon River
transition_matrix["Salmon River", "loss"] <- 1
transition_matrix["Salmon River", "mainstem, upstream of LGR"] <- 1


# 25: Grande Ronde River
transition_matrix["Grande Ronde River", "loss"] <- 1
transition_matrix["Grande Ronde River", "mainstem, upstream of LGR"] <- 1


# 26: Imnaha River
transition_matrix["Imnaha River", "loss"] <- 1
transition_matrix["Imnaha River", "mainstem, upstream of LGR"] <- 1

# 27: BON to MCN other trib
transition_matrix["BON to MCN other trib", "loss"] <- 1
transition_matrix["BON to MCN other trib", "mainstem, BON to MCN"] <- 1

# 28: BON to MCN other trib
transition_matrix["Upstream WEL other trib", "loss"] <- 1
transition_matrix["Upstream WEL other trib", "mainstem, upstream of WEL"] <- 1

# Export the full model transition matrix for import into other scripts

write.csv(transition_matrix, here::here("figures", "full_model_transition_matrix.csv"))


# Get the indices that are 1s (except loss, since that's 1 minus the others)
movements <- which(transition_matrix[,1:(nstates-1)] == 1, arr.ind = TRUE)
nmovements <- dim(movements)[1]
# Now get all of the movements which are fixed to zero
not_movements <- which(transition_matrix[,1:(nstates-1)] == 0, arr.ind = TRUE)
n_notmovements <- dim(not_movements)[1]




##### Convert complete detection history into an array #####
# This array will have dimensions nstates x noccasions (max) x nfish

# Read in states
read.csv(here::here("from_hyak_transfer", "2022-05-25-complete_det_hist", "states_complete.csv")) %>% 
  dplyr::select(-X) -> states_complete

# Keep only individuals from the natal origins we're interested in

# Read in natal origin data
natal_origin_table <- read.csv(here::here("covariate_data", "natal_origin_table.csv"))
# Read in tag code metadata
tag_code_metadata <- read.csv(here::here("covariate_data", "tag_code_metadata.csv"))

tag_code_metadata %>% 
  dplyr::select(tag_code, release_site_name) %>% 
  left_join(., natal_origin_table, by = "release_site_name") -> tag_code_origins

# Remove fish from Klickitat and Wind River
KLIC_WIND_tag_codes <- subset(tag_code_origins, natal_origin %in% c("Klickitat_River", "Wind_River"))$tag_code
# Subset out those
states_complete %>% 
  subset(., !(tag_code %in% KLIC_WIND_tag_codes)) -> states_complete

# Keep fish from only years that arrays were active
tag_code_metadata %>% 
  left_join(., natal_origin_table, by = "release_site_name") %>% 
  dplyr::select(tag_code, run_year, natal_origin, rear_type_code) %>% 
  subset(., tag_code %in% BON_arrival$tag_code) -> origin_metadata

# Remove individuals from run years where arrays were not active
origin_metadata %>% 
  mutate(year_numeric = as.numeric(sub("\\/.*", "", run_year))) %>% 
  mutate(array_active = ifelse(natal_origin == "Methow_River" & year_numeric < 9, "inactive",
                               ifelse(natal_origin == "Klickitat_River" & year_numeric < 12, "inactive",
                                      ifelse(natal_origin == "Okanogan_River" & year_numeric < 14, "inactive",
                                             ifelse(natal_origin == "Wind_River" & year_numeric < 13, "inactive",
                                                    ifelse(natal_origin == "Asotin_Creek" & year_numeric < 12, "inactive", "active")))))) -> origin_metadata

# Fish from inactive array years
inactive_array_fish_tag_codes <- subset(origin_metadata, array_active == "inactive")$tag_code
# Subset out those
states_complete %>% 
  subset(., !(tag_code %in% inactive_array_fish_tag_codes)) -> states_complete



# complete_det_hist <- read.csv(here::here("from_hyak_transfer", "2022-05-23-complete_det_hist", "complete_det_hist.csv"), row.names = 1)

# Convert states_complete into an array
n.ind <- length(unique(states_complete$tag_code))
tag_codes <- unique(states_complete$tag_code)
# Figure out the maximum number of state visits
states_complete %>% 
  count(tag_code) -> states_by_tagcode_count

n.obs <- states_by_tagcode_count$n

max_visits <- max(states_by_tagcode_count$n)

# Make nfish smaller for testing
# n.ind <- 100
# Create an empty array
movement_array <- array(0, c(nstates,max_visits,n.ind))
# Holy cow this is huge

# Start a row counter for indexing
counter <- 1

for (i in 1:n.ind){ # loop through each fish
  # Loop through the observations for each fish
  for (j in 1:n.obs[i]){
    movement_array[which(model_states == states_complete[counter, 'state']),j,i] <- 1
    counter <- counter + 1
  }
  # Once the detection history is done, then move it to loss state
  movement_array[nstates, n.obs[i] + 1,i] <- 1
}



##### Get other data for JAGS model #####
# Get the state that each fish was in at each n.obs

# Turn into matrix for JAGS
states_mat <- matrix(nrow = n.ind, ncol = max(n.obs))
for (i in 1:n.ind){
  # states_mat[i,1:(n.obs[i])] <- states_list[[i]]
  for (j in 1:(n.obs[i])){
    states_mat[i,j] <- which(movement_array[,j,i] == 1)
  }
}


###### Parameters monitored #####
parameters <- c(
  "b0_matrix"
)

##### Write function to write out JAGS model #####

# To make this more efficient, we will loop this to write out the JAGS model
# Nstates is all of physical states plus the loss state - so therefore -1 to remove loss state

sink(file = here::here("full_JAGS_models", "full_JAGS_nocov.txt"), type = "output")
cat("model {", "\n", 
    "for(i in 1:n.ind){ # Loop through the detection matrices for each individual", "\n", 
    "for(j in 1:(n.obs[i])){ # Loop through each of the observations", "\n",
    "y[,j+1,i] ~ dmulti(c(", "\n", "\n" )

for (p in 1:(nstates-1)){
  # Paste the numerator
  cat("# State", p, "\n")
  cat(paste0("possible_states[states_mat[i,j], ", p, "] * exp(b0_matrix[states_mat[i,j], ", p, "])", "/", "\n", "\n")) # this is the numerator
  # Paste the denominator
  cat("(1 + ", "\n")
  for (q in 1:(nstates-2)){ # Need to stop 1 before last state, since we need to close parentheses and not add a "+" on the last line
    cat(paste0("possible_states[states_mat[i,j], ", q, "] * exp(b0_matrix[states_mat[i,j], ", q, "])", " +", "\n"))

  }
  # Add the last line of the denominator
  cat(paste0("possible_states[states_mat[i,j], ", nstates-1, "] * exp(b0_matrix[states_mat[i,j], ", nstates-1, "])", "),", "\n"))
  
  cat("\n")
}

# Now write out the loss (1 - everything else)
cat("# Loss", "\n")
cat("1 - sum(")
for (p in 1:(nstates-2)){ # the last state is loss, so that's -1, and the second to last state is
  # different because we need to close the parentheses and not add a comma
  # Paste the numerator
  cat("# State", p, "\n")
  cat(paste0("possible_states[states_mat[i,j], ", p, "] * exp(b0_matrix[states_mat[i,j], ", p, "])", "/", "\n", "\n")) # this is the numerator
  # Paste the denominator
  cat("(1 + ", "\n")
  for (q in 1:(nstates-2)){ # Need to stop 1 before last state, since we need to close parentheses and not add a "+" on the last line
    cat(paste0("possible_states[states_mat[i,j], ", q, "] * exp(b0_matrix[states_mat[i,j], ", q, "])", " +", "\n"))
    
  }
  # Add the last line of the denominator
  cat(paste0("possible_states[states_mat[i,j], ", nstates-1, "] * exp(b0_matrix[states_mat[i,j], ", nstates-1, "]", ")),", "\n"))
  
  cat("\n")
}

# The last line of the loss
cat("# State", nstates-1, "\n")
cat(paste0("possible_states[states_mat[i,j], ", nstates-1, "] * exp(b0_matrix[states_mat[i,j], ", nstates-1, "])", "/", "\n", "\n")) # this is the numerator
# Paste the denominator
cat("(1 + ", "\n")
for (q in 1:(nstates-2)){ # Need to stop 1 before last state, since we need to close parentheses and not add a "+" on the last line
  cat(paste0("possible_states[states_mat[i,j], ", q, "] * exp(b0_matrix[states_mat[i,j], ", q, "])", " +", "\n"))
  
}
# Add the last line of the denominator
cat(paste0("possible_states[states_mat[i,j], ", nstates-1, "] * exp(b0_matrix[states_mat[i,j], ", nstates-1, "]", ")))", "\n"))

cat("\n")




# Close the sum() call for loss
cat (")")

# Close the dmulti call
cat(", 1)", "\n")

# Close the `for(i in 1:n.ind){` loop
cat("}", "\n")

# Close the `for(j in 1:(n.obs[i])){` loop
cat("}", "\n", "\n")

# Now write out the priors

cat("#### PRIORS ####",  "\n",  "\n")
cat("# Set priors for all possible movements", "\n")
cat("for (i in 1:nmovements){", "\n")
cat("b0_matrix[movements[i,1], movements[i,2]] ~ dnorm(0,0.001)", "\n")
cat("}", "\n")

cat("# Set all not possible movements to -9999 - these aren't used", "\n")
cat("for (i in 1:n_notmovements){", "\n")
cat("b0_matrix[not_movements[i,1], not_movements[i,2]] <- -9999", "\n")
cat("}", "\n", "\n")

# Now close the original bracket for model
cat("}")

sink()


##### Data #####


data <- list(y = movement_array, n.ind = n.ind, n.obs = n.obs,
             states_mat = states_mat, nstates = length(model_states),
             movements = movements, not_movements = not_movements,
             nmovements = nmovements, 
             n_notmovements = n_notmovements, possible_states = transition_matrix)


###### Initial values #####
inits <- function(){
  b0_matrix <- matrix(NA, nrow = nstates, ncol = (nstates-1))
  
  for (j in 1:dim(movements)[1]){
    b0_matrix[movements[j,1], movements[j,2]] <- runif(1,-1,1)
  }
  
  return(list(
    b0_matrix = b0_matrix
  ))
}

start_time <- Sys.time()
print(paste0("Start time: ", start_time))


##### Run JAGS model #####
out.jags <- 
  jags.parallel(
    data = data,
    inits = inits,
    model.file="full_JAGS_nocov.txt",
    # model.file=here::here("full_JAGS_models", "full_JAGS_nocov.txt"),
    parameters.to.save = parameters,
    n.chains = 3, n.iter = 10000, n.burnin = 5000,
    n.thin = 10,
    jags.seed = 123
  )

# Save JAGS model output
saveRDS(out.jags, "fullmodel_intercept_only_JAGS.rds")



end_time <- Sys.time()
print(paste0("End time: ", end_time))
print(paste0("Total run time: ", end_time - start_time))


