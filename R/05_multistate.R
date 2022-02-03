### 05 - Fit multistate

# This R script fits the multistate, bidirectional model

# Note: You must run scripts 01 to 03 to generate the detections used in this model


### Load libraries
library(here)
library(plyr)
library(tidyverse)
library(lubridate)
library(janitor)


# Load data
JDR_states <- read.csv(here("model_files", "JDR_states.csv"), row.names = 1)

# Break up state history into individual detections
JDR_tag_codes <- unique(JDR_states$tag_code)
# Loop through each fish
# for (i in 1:length(unique(JDR_tag_codes))){
# Make a df to store values
JDR_stepwise_detections <- data.frame(tag_code = character(),
                                      state_1 = character(),
                                      state_2 = character(),
                                      date_time_1 = as.POSIXct(character()),
                                      date_time_2 = as.POSIXct(character()),
                                      pathway = character())

for (i in 1:(nrow(JDR_states)-1)){
  
  # If it's the same fish:
  if (JDR_states[i,'tag_code'] == JDR_states[i+1,'tag_code']){
    # Store the tag code
    JDR_stepwise_detections[i,'tag_code'] <- JDR_states[i,'tag_code']
    # Store the current site, and the next site
    JDR_stepwise_detections[i,'state_1'] <- JDR_states[i,'state']
    JDR_stepwise_detections[i,'state_2'] <- JDR_states[i+1,'state']
    # Store the time leaving this state, and the time entering the next state
    JDR_stepwise_detections[i,'date_time_1'] <- JDR_states[i,'date_time']
    JDR_stepwise_detections[i,'date_time_2'] <- JDR_states[i+1,'date_time']
    # Store the pathway between states
    JDR_stepwise_detections[i,'pathway'] <- JDR_states[i+1, 'pathway']
  }
  
  # If it's a different fish, end the entry (note that it's lost from the detection history)
  # This may be because a fish spawned, or it could be undetermined loss at the end
  else {
    # End the previous entry
    # Store the tag code
    JDR_stepwise_detections[i,'tag_code'] <- JDR_states[i,'tag_code']
    # Store the current state, and the next state
    JDR_stepwise_detections[i,'state_1'] <- JDR_states[i,'state']
    JDR_stepwise_detections[i,'state_2'] <- "lost"
    # Store the time leaving this state, and the time entering the next state
    JDR_stepwise_detections[i,'date_time_1'] <- JDR_states[i,'date_time']
    JDR_stepwise_detections[i,'date_time_2'] <- NA
  }
  
}


# Remove all of the NA rows
JDR_stepwise_detections <- JDR_stepwise_detections[!is.na(JDR_stepwise_detections$tag_code),]

##### Data summaries #####

# We can fairly quickly assess rates of fallback using this data (via a tally-based approach)

# Over McNary Dam: state_1 == "mainstem, MCN to ICH or PRA", state_2 == "mainstem, BON to MCN"
JDR_stepwise_detections %>% 
  subset(., state_1 == "mainstem, MCN to ICH or PRA" & state_2 == "mainstem, BON to MCN") -> JDR_MCN_fallback

JDR_MCN_fallback_individuals <- unique(JDR_MCN_fallback$tag_code)

# Count number of overshooting individuals at McNary
JDR_stepwise_detections %>% 
  subset(., state_1 == "mainstem, BON to MCN"  & state_2 == "mainstem, MCN to ICH or PRA") -> JDR_MCN_overshoot

JDR_MCN_overshoot_individuals <- unique(JDR_MCN_overshoot$tag_code)

length(JDR_MCN_fallback_individuals)/length(JDR_MCN_overshoot_individuals)

# This would be a good individual to illustrate how my code works
subset(JDR_det_hist_original, tag_code == "3D9.1C2C89DE1F")
subset(JDR_states, tag_code == "3D9.1C2C89DE1F")









##### Multistate model #####

# After detection in the adult fishways at each dam, there are four or five probabilities
# that must sum to 1, for all of the options for a fish.
# Let's use an individual that was just detected at Bonneville as an example:
# 1) Overshooting MCN
# 2) Falling back over BON (can either be detected at a fallback array, or seen at the same dam)
# 3) Straying to a tributary between BON and MCN
# 4) Undetermined loss (not seen again)
# 5) Homing to the natal tributary (between BON and MCN) - only for individuals between BON and MCN

# For those downstream of Bonneville, we can only estimate overshooting and loss probabilities


# We can define these as o, f, s, l, and h, respectively. Each will have a subscript for each 
# reach/dam.

# Let's turn the state history for each tag code into a detection history
JDR_states %>% 
  # First make a column to order the observations by tag_code (individual)
  group_by(tag_code) %>% 
  mutate(order = paste0("state_", row_number())) %>% 
  pivot_wider(., id_cols = tag_code, names_from = order, values_from = state) -> test
  
subset(test, !is.na(state_23)) -> what_is_this_fish

# This would also be a good individual to illustrate how my code works
longest_det_hist <- subset(JDR_det_hist, tag_code == "3D9.1C2C31C103")
# This guy is wild! Looks like two spawning events, then return to sea


# You could then input these into a multinomial likelihood, but boy would
# it struggle to calculate it

# Turn those stepwise detections into a parameter for each

# Make a dataframe that contains this information
JDR_state_model_probabilities <- 
  
  # For fish that are downstream of Bonneville
  data.frame(state_1 = "mainstem, mouth to BON", state_2 = "lost", probability = "l_bon") %>% # loss
  bind_rows(., data.frame(state_1 = "mainstem, mouth to BON", state_2 = "mainstem, BON to MCN", probability = "o_bon")) %>% # overshoot
  
  # Fish that are BON -> MCN
  bind_rows(., data.frame(state_1 = "mainstem, BON to MCN", state_2 = "mainstem, MCN to ICH or PRA", probability = "o_mcn")) %>%  #overshoot
  bind_rows(., data.frame(state_1 = "mainstem, BON to MCN", state_2 = "mainstem, mouth to BON", probability = "f_bon")) %>% # fallback
  bind_rows(., data.frame(state_1 = "mainstem, BON to MCN", state_2 = "BON to MCN tributaries", probability = "s_bon_mcn")) %>% # stray
  bind_rows(., data.frame(state_1 = "mainstem, BON to MCN", state_2 = "natal tributaries", probability = "h_bon_mcn")) %>% # home to natal
  bind_rows(., data.frame(state_1 = "mainstem, BON to MCN", state_2 = "lost", probability = "l_bon_mcn")) %>% # loss
  
  # Fish that are MCN to ICH or PRA
  # Extra overshoot probability, because we're branching here
  bind_rows(., data.frame(state_1 = "mainstem, MCN to ICH or PRA", state_2 = "mainstem, PRA to RIS", probability = "o_pra")) %>%  #overshoot PRA
  bind_rows(., data.frame(state_1 = "mainstem, MCN to ICH or PRA", state_2 = "mainstem, ICH to LGR", probability = "o_ich")) %>%  #overshoot ICH
  bind_rows(., data.frame(state_1 = "mainstem, MCN to ICH or PRA", state_2 = "mainstem, BON to MCN", probability = "f_mcn")) %>% # fallback
  bind_rows(., data.frame(state_1 = "mainstem, MCN to ICH or PRA", state_2 = "MCN to PRA or ICH tributaries", probability = "s_mcn_ich_pra")) %>% # stray
  bind_rows(., data.frame(state_1 = "mainstem, MCN to ICH or PRA", state_2 = "lost", probability = "l_mcn_ich_pra")) %>% # loss
  
  # UPPER COLUMBIA
  
  # Fish that are PRA to RIS
  bind_rows(., data.frame(state_1 = "mainstem, PRA to RIS", state_2 = "mainstem, RIS to RRE", probability = "o_ris")) %>%  #overshoot RIS
  bind_rows(., data.frame(state_1 = "mainstem, PRA to RIS", state_2 = "mainstem, MCN to ICH or PRA", probability = "f_pra")) %>% # fallback 
  # bind_rows(., data.frame(state_1 = "mainstem, PRA to RIS", state_2 = "PRA to RIS tributaries", probability = "s_pra_ris")) %>% # stray - no individuals strayed to these
  bind_rows(., data.frame(state_1 = "mainstem, PRA to RIS", state_2 = "lost", probability = "l_pra_ris")) %>% # loss
 
  # Fish that are RIS to RRE
  bind_rows(., data.frame(state_1 = "mainstem, RIS to RRE", state_2 = "mainstem, RRE to WEL", probability = "o_rre")) %>%  #overshoot RRE
  bind_rows(., data.frame(state_1 = "mainstem, RIS to RRE", state_2 = "mainstem, PRA to RIS", probability = "f_ris")) %>% # fallback 
  # bind_rows(., data.frame(state_1 = "mainstem, RIS to RRE", state_2 = "RIS to RRE tributaries", probability = "s_ris_rre")) %>% # stray - no individuals strayed to these
  bind_rows(., data.frame(state_1 = "mainstem, RIS to RRE", state_2 = "lost", probability = "l_ris_rre")) %>% # loss
  
  # Fish that are RRE to WEL
  bind_rows(., data.frame(state_1 = "mainstem, RRE to WEL", state_2 = "mainstem, upstream of WEL", probability = "o_wel")) %>%  #overshoot WEL
  bind_rows(., data.frame(state_1 = "mainstem, RRE to WEL", state_2 = "mainstem, RIS to RRE", probability = "f_rre")) %>% # fallback 
  # bind_rows(., data.frame(state_1 = "mainstem, RRE to WEL", state_2 = "RRE to WEL tributaries", probability = "s_rre_wel")) %>% # stray - no individuals strayed to these
  bind_rows(., data.frame(state_1 = "mainstem, RRE to WEL", state_2 = "lost", probability = "l_rre_wel")) %>% # loss
  
  # Fish that are upstream of WEL
  # nothing to overshoot
  bind_rows(., data.frame(state_1 = "mainstem, upstream of WEL", state_2 = "mainstem, RIS to RRE", probability = "f_rre")) %>% # fallback 
  # bind_rows(., data.frame(state_1 = "mainstem, RRE to WEL", state_2 = "RRE to WEL tributaries", probability = "s_rre_wel")) %>% # stray - no individuals strayed to these
  bind_rows(., data.frame(state_1 = "mainstem, upstream of WEL", state_2 = "lost", probability = "l_wel")) %>% # loss
  
  # SNAKE RIVER
  
  # Fish that are ICH to LGR
  bind_rows(., data.frame(state_1 = "mainstem, ICH to LGR", state_2 = "mainstem, upstream of LGR", probability = "o_lgr")) %>%  #overshoot RRE
  bind_rows(., data.frame(state_1 = "mainstem, ICH to LGR", state_2 = "mainstem, MCN to ICH or PRA", probability = "f_ich")) %>% # fallback 
  bind_rows(., data.frame(state_1 = "mainstem, ICH to LGR", state_2 = "ICH to LGR tributaries", probability = "s_ich_lgr")) %>% # stray - no individuals strayed to these
  bind_rows(., data.frame(state_1 = "mainstem, ICH to LGR", state_2 = "lost", probability = "l_ich_lgr")) %>% # loss
  
  # Fish that are upstream of LGR
  # nothing to overshoot
  bind_rows(., data.frame(state_1 = "mainstem, upstream of LGR", state_2 = "mainstem, ICH to LGR", probability = "f_lgr")) %>% # fallback 
  bind_rows(., data.frame(state_1 = "mainstem, upstream of LGR", state_2 = "Upstream LGR tributaries", probability = "s_lgr")) %>% # stray - no individuals strayed to these
  bind_rows(., data.frame(state_1 = "mainstem, upstream of LGR", state_2 = "lost", probability = "l_lgr")) %>%  # loss
  
  # Fish that are in natal tributaries
  bind_rows(., data.frame(state_1 = "natal tributaries", state_2 = "mainstem, BON to MCN", probability = "r_nat_trib")) %>% # fallback 
  bind_rows(., data.frame(state_1 = "natal tributaries", state_2 = "lost", probability = "l_nat_trib")) %>% # lost (likely spawned)
  
  # Fish in BON to MCN tributaries (stray)
  bind_rows(., data.frame(state_1 = "BON to MCN tributaries", state_2 = "mainstem, BON to MCN", probability = "r_bon_mcn_trib")) %>% # fallback 
  bind_rows(., data.frame(state_1 = "BON to MCN tributaries", state_2 = "lost", probability = "l_bon_mcn_trib")) %>% # lost (likely spawned)
  
  # Fish in MCN to PRA or ICH tributaries (stray)
  bind_rows(., data.frame(state_1 = "MCN to PRA or ICH tributaries", state_2 = "mainstem, MCN to ICH or PRA", probability = "r_mcn_pra_ich_trib")) %>% # fallback 
  bind_rows(., data.frame(state_1 = "MCN to PRA or ICH tributaries", state_2 = "lost", probability = "l_mcn_pra_ich_trib")) %>% # lost (likely spawned)
  
  # Fish in ICH to LGR tributaries (stray)
  bind_rows(., data.frame(state_1 = "ICH to LGR tributaries", state_2 = "mainstem, ICH to LGR", probability = "r_lgr_trib")) %>% # fallback 
  bind_rows(., data.frame(state_1 = "ICH to LGR tributaries", state_2 = "lost", probability = "l_ich_lgr_trib")) %>% # lost (likely spawned)
  
  # Fish in tributaries upstream of LGR (stray)
  bind_rows(., data.frame(state_1 = "Upstream LGR tributaries", state_2 = "mainstem, upstream of LGR", probability = "r_lgr_trib")) %>% # fallback 
  bind_rows(., data.frame(state_1 = "Upstream LGR tributaries", state_2 = "lost", probability = "l_lgr_trib")) # lost (likely spawned)



  # Add model probabilities to stepwise detections
JDR_stepwise_detections %>% 
  left_join(., JDR_state_model_probabilities, by = c("state_1", "state_2")) -> JDR_stepwise_probabilities

# Export this
write.csv(JDR_stepwise_probabilities, here("model_files", "JDR_stepwise_probabilities"))


# check for any missing probabilities
subset(JDR_stepwise_probabilities, is.na(probability)) -> missing_prob


##### Convert into probabilities for each detection history

# create a new df to store these
JDR_complete_probs <- data.frame(tag_code = unique(JDR_stepwise_probabilities$tag_code), prob = rep(NA, length(unique(JDR_stepwise_probabilities$tag_code))))

# Loop through tag codes
for (i in 1:dim(JDR_complete_probs)[1]){
  # Get the unique tag code
  ind_tag_code <- unique(JDR_complete_probs$tag_code)[i]
  
  # Subset the stepwise probabilities for that tag code
  ind_det_hist <- subset(JDR_stepwise_probabilities, tag_code == ind_tag_code)
  
  # Concatenate the stepwise probabilities
  ind_prob <- paste(ind_det_hist$probability, collapse = " * ")
  
  JDR_complete_probs$prob[i] <- ind_prob
  
}

# Tally the unique combinations of probabilities
as.data.frame(table(JDR_complete_probs$prob)) %>% 
  dplyr::rename(prob = Var1, count = Freq) -> JDR_unique_probs

#### Write out NLL function ####

# Write parameters to text

JDR_params <- data.frame(param = sort(unique(JDR_stepwise_probabilities$probability)),
                            params = rep(NA, length(unique(JDR_stepwise_probabilities$probability))))

# Remove all of the loss parameters, since these are 1 - the other params
JDR_params <- JDR_params[!(grepl("l_", JDR_params$param)),]

JDR_params %>% 
  mutate(params = paste0("params[",row_number(JDR_params$param),"]")) %>% 
  mutate(to_text = paste0(param, " <- ", params)) -> JDR_params

# Write all loss parameters out
JDR_loss_params <- data.frame(param = sort(unique(JDR_stepwise_probabilities$probability)),
                              eqn = rep(NA, length(unique(JDR_stepwise_probabilities$probability))))
JDR_loss_params <- JDR_loss_params[grepl("l_", JDR_loss_params$param),]

JDR_loss_params$eqn <- c("1 - o_bon", "1 - (o_mcn + s_bon_mcn + h_bon_mcn + f_bon)", 
                         "1 - (r_bon_mcn_trib)", "1 - (o_lgr + f_ich + s_ich_lgr)", 
                         "1  - r_ich_lgr_trib" , "1 - (f_lgr + s_lgr)", "1 - (r_lgr_trib)", 
                         "1 - (f_mcn + o_pra + o_ich + s_mcn_ich_pra)", 
                         "1 - (r_mcn_pra_ich_trib)", "1 - r_nat_trib", 
                         "1 - (f_pra + o_ris)", "1 - (o_rre + f_ris)", "1 - f_wel")

JDR_loss_params %>% 
  mutate(to_text = paste0(param, " <- ", eqn)) -> JDR_loss_params


# Write out probabilities for each detection history

JDR_unique_probs %>% 
  mutate(to_text = paste0(paste0("# ", row_number(JDR_unique_probs$prob)), # Comment the probability
                          "\n",  paste0("n", row_number(JDR_unique_probs$prob), " <- ", prob), "\n")) -> JDR_unique_probs

# Write out vector of probabilities for multinomial likelihood
# write out n values as a vector
n_vec_length <- length(unique(JDR_complete_probs$prob))
n_vec <- paste0("n", 1:n_vec_length)
# Add breaks in to make it legible
# Turn vector into matrix
# First, augment n_vec to make it a multiple of 6
augmented_values <- rep("", round_any(n_vec_length, 6, f = ceiling)-n_vec_length)
n_vec_aug <- c(n_vec, augmented_values)
n_vec_matrix <- matrix(data = n_vec_aug, nrow = 6)

# Now, create a vector of line breaks
line_breaks <- rep("\n", ncol(n_vec_matrix))

# Combine line breaks and n_vec
rbind(n_vec_matrix, line_breaks) -> n_vec_matrix_breaks
n_vec_2 <- c(as.matrix(n_vec_matrix_breaks))
# remove the augmented values
n_vec_final <- n_vec_final[1:(length(n_vec_2) - (length(augmented_values)+1))]
# collapse vector into string of probabilities
n_vec_string <- paste(n_vec_final, collapse = ", ")
# Remove commas after line breaks
n_vec_string <- str_replace_all(n_vec_string, "\n,", "\n")
prob_vec <- paste0("p <- c(",n_vec_string, ")")

# Multinomial likelihood function - Write to file: start of function call, params (from optim), 
# detection history probabilities, prob vector, optim call
writeLines(c("multistate_model <- function(data) {", "\n", "data <- data", "\n",
             "negLL = function(params, data){", "\n", "# Get parameters from param vector",
             JDR_params$to_text, "\n", "# Solve for all of the loss parameters as 1 - sum of the other parameters", "\n",
              JDR_loss_params$to_text, "\n",  "# Get probabilities for each movement history",
             JDR_unique_probs$to_text, "\n", "# Get vector of probabilities for multinomial likelihood",
             prob_vec, "\n", "# dmultinom call", "negLL <- -1* dmultinom(x = data$count, prob = p, log = TRUE)", "}", "\n",
             "optim_results <- optim(par = optim_inits, data = data, fn = negLL, method = 'L-BFGS-B', hessian = FALSE, lower = 0.0001, upper = 0.9999)",
             "\n", "return(optim_results)", "\n", "}"),
           file(here::here("model_files", "JDR_multinomial_model_fit_function.txt")), sep = "\n")

# Paste the multinomial likelihood MLE fit function you just wrote

##### Multinomial MLE function v2 - use convenience function, 1 - other params, to constrain optimization #####
multistate_model <- function(data, optim_inits) {
  
  
  data <- data
  
  
  negLL = function(params, data){
    
    
    # Get parameters from param vector
    f_bon <- params[1]
    f_ich <- params[2]
    f_lgr <- params[3]
    f_mcn <- params[4]
    f_pra <- params[5]
    f_ris <- params[6]
    h_bon_mcn <- params[7]
    # Solve for all of the loss parameters as 1 - sum of the other parameters
    # l_bon <- params[8]
    # l_bon_mcn <- params[9]
    # l_bon_mcn_trib <- params[10]
    # l_ich_lgr <- params[11]
    # l_ich_lgr_trib <- params[12]
    # l_lgr <- params[13]
    # l_lgr_trib <- params[14]
    # l_mcn_ich_pra <- params[15]
    # l_mcn_pra_ich_trib <- params[16]
    # l_nat_trib <- params[17]
    # l_pra_ris <- params[18]
    # l_ris_rre <- params[19]
    # l_wel <- params[20]
    
    # We need to reorder these now
    o_bon <- params[8]
    o_ich <- params[9]
    o_lgr <- params[10]
    o_mcn <- params[11]
    o_pra <- params[12]
    o_ris <- params[13]
    o_rre <- params[14]
    o_wel <- params[15]
    r_bon_mcn_trib <- params[16]
    r_lgr_trib <- params[17]
    r_mcn_pra_ich_trib <- params[18]
    r_nat_trib <- params[19]
    s_bon_mcn <- params[20]
    s_ich_lgr <- params[21]
    s_lgr <- params[22]
    s_mcn_ich_pra <- params[23]
    
    # Solve for loss probabilities
    l_bon <- 1 - o_bon
    l_bon_mcn <- 1 - (o_mcn + s_bon_mcn + h_bon_mcn + f_bon)
    l_bon_mcn_trib <- 1 - (r_bon_mcn_trib)
    l_ich_lgr <- 1 - (o_lgr + f_ich + s_ich_lgr)
    l_ich_lgr_trib <- 1 # Note: This would be 1 - r_ich_lgr_trib, but no fish returned after straying into these tributaries
    l_lgr <- 1 - (f_lgr + s_lgr)
    l_lgr_trib <- 1 - (r_lgr_trib)
    l_mcn_ich_pra <- 1 - (f_mcn + o_pra + o_ich + s_mcn_ich_pra)
    l_mcn_pra_ich_trib <- 1 - (r_mcn_pra_ich_trib)
    l_nat_trib <- 1 - r_nat_trib
    l_pra_ris <- 1 - (f_pra + o_ris)
    l_ris_rre <- 1 - (o_rre + f_ris)
    l_wel <- 1 # Note: This would be 1 - f_wel, but no fish fell back over Well's dam
    
    
    # Get probabilities for each movement history
    # 1
    n1 <- f_bon * l_bon
    
    # 2
    n2 <- f_bon * o_bon * f_bon * l_bon
    
    # 3
    n3 <- f_bon * o_bon * f_bon * o_bon * f_bon * o_bon * f_bon * o_bon * l_bon_mcn
    
    # 4
    n4 <- f_bon * o_bon * f_bon * o_bon * f_bon * o_bon * f_bon * o_bon * o_mcn * l_mcn_ich_pra
    
    # 5
    n5 <- f_bon * o_bon * f_bon * o_bon * f_bon * o_bon * h_bon_mcn * l_nat_trib
    
    # 6
    n6 <- f_bon * o_bon * f_bon * o_bon * f_bon * o_bon * l_bon_mcn
    
    # 7
    n7 <- f_bon * o_bon * f_bon * o_bon * h_bon_mcn * l_nat_trib
    
    # 8
    n8 <- f_bon * o_bon * f_bon * o_bon * l_bon_mcn
    
    # 9
    n9 <- f_bon * o_bon * f_bon * o_bon * o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 10
    n10 <- f_bon * o_bon * f_bon * o_bon * o_mcn * f_mcn * o_mcn * o_ich * l_ich_lgr
    
    # 11
    n11 <- f_bon * o_bon * f_bon * o_bon * o_mcn * l_mcn_ich_pra
    
    # 12
    n12 <- f_bon * o_bon * h_bon_mcn * l_nat_trib
    
    # 13
    n13 <- f_bon * o_bon * h_bon_mcn * r_nat_trib * f_bon * l_bon
    
    # 14
    n14 <- f_bon * o_bon * h_bon_mcn * r_nat_trib * h_bon_mcn * l_nat_trib
    
    # 15
    n15 <- f_bon * o_bon * l_bon_mcn
    
    # 16
    n16 <- f_bon * o_bon * o_mcn * f_mcn * f_bon * l_bon
    
    # 17
    n17 <- f_bon * o_bon * o_mcn * f_mcn * f_bon * o_bon * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 18
    n18 <- f_bon * o_bon * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 19
    n19 <- f_bon * o_bon * o_mcn * f_mcn * h_bon_mcn * r_nat_trib * l_bon_mcn
    
    # 20
    n20 <- f_bon * o_bon * o_mcn * f_mcn * l_bon_mcn
    
    # 21
    n21 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * f_mcn * f_bon * l_bon
    
    # 22
    n22 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 23
    n23 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * f_mcn * l_bon_mcn
    
    # 24
    n24 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 25
    n25 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * r_nat_trib * f_bon * l_bon
    
    # 26
    n26 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 27
    n27 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 28
    n28 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 29
    n29 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 30
    n30 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * o_ich * f_ich * o_ich * l_ich_lgr
    
    # 31
    n31 <- f_bon * o_bon * o_mcn * f_mcn * o_mcn * o_ich * l_ich_lgr
    
    # 32
    n32 <- f_bon * o_bon * o_mcn * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 33
    n33 <- f_bon * o_bon * o_mcn * f_mcn * s_bon_mcn * r_bon_mcn_trib * h_bon_mcn * l_nat_trib
    
    # 34
    n34 <- f_bon * o_bon * o_mcn * l_mcn_ich_pra
    
    # 35
    n35 <- f_bon * o_bon * o_mcn * o_ich * f_ich * f_mcn * f_bon * l_bon
    
    # 36
    n36 <- f_bon * o_bon * o_mcn * o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 37
    n37 <- f_bon * o_bon * o_mcn * o_ich * f_ich * f_mcn * o_mcn * o_ich * o_lgr * l_lgr
    
    # 38
    n38 <- f_bon * o_bon * o_mcn * o_ich * f_ich * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 39
    n39 <- f_bon * o_bon * o_mcn * o_ich * f_ich * l_mcn_ich_pra
    
    # 40
    n40 <- f_bon * o_bon * o_mcn * o_ich * f_ich * o_ich * f_ich * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 41
    n41 <- f_bon * o_bon * o_mcn * o_ich * l_ich_lgr
    
    # 42
    n42 <- f_bon * o_bon * o_mcn * o_ich * o_lgr * f_lgr * f_ich * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 43
    n43 <- f_bon * o_bon * o_mcn * o_ich * o_lgr * f_lgr * l_ich_lgr
    
    # 44
    n44 <- f_bon * o_bon * o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * f_ich * f_mcn * f_bon * l_bon
    
    # 45
    n45 <- f_bon * o_bon * o_mcn * o_ich * o_lgr * f_lgr * o_lgr * l_lgr
    
    # 46
    n46 <- f_bon * o_bon * o_mcn * o_ich * o_lgr * f_lgr * s_ich_lgr * l_ich_lgr_trib
    
    # 47
    n47 <- f_bon * o_bon * o_mcn * o_ich * o_lgr * l_lgr
    
    # 48
    n48 <- f_bon * o_bon * o_mcn * o_ich * o_lgr * s_ich_lgr * l_lgr_trib
    
    # 49
    n49 <- f_bon * o_bon * o_mcn * o_ich * o_lgr * s_ich_lgr * r_lgr_trib * f_lgr * l_ich_lgr
    
    # 50
    n50 <- f_bon * o_bon * o_mcn * o_pra * f_pra * f_mcn * h_bon_mcn * l_nat_trib
    
    # 51
    n51 <- f_bon * o_bon * o_mcn * o_pra * f_pra * o_pra * o_ris * f_ris * f_pra * o_ich * o_lgr * f_lgr * o_lgr * l_lgr
    
    # 52
    n52 <- f_bon * o_bon * o_mcn * o_pra * o_ris * l_ris_rre
    
    # 53
    n53 <- f_bon * o_bon * o_mcn * s_mcn_ich_pra * l_mcn_pra_ich_trib
    
    # 54
    n54 <- f_bon * o_bon * s_bon_mcn * l_bon_mcn_trib
    
    # 55
    n55 <- f_bon * o_bon * s_bon_mcn * r_bon_mcn_trib * o_mcn * l_mcn_ich_pra
    
    # 56
    n56 <- h_bon_mcn * l_nat_trib
    
    # 57
    n57 <- h_bon_mcn * r_nat_trib * f_bon * l_bon
    
    # 58
    n58 <- h_bon_mcn * r_nat_trib * l_bon_mcn
    
    # 59
    n59 <- h_bon_mcn * r_nat_trib * o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 60
    n60 <- h_bon_mcn * r_nat_trib * o_mcn * o_ich * o_lgr * l_lgr
    
    # 61
    n61 <- l_bon_mcn
    
    # 62
    n62 <- o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 63
    n63 <- o_mcn * f_mcn * f_bon * l_bon
    
    # 64
    n64 <- o_mcn * f_mcn * f_bon * o_bon * h_bon_mcn * l_nat_trib
    
    # 65
    n65 <- o_mcn * f_mcn * f_bon * o_bon * h_bon_mcn * r_nat_trib * f_bon * l_bon
    
    # 66
    n66 <- o_mcn * f_mcn * f_bon * o_bon * l_bon_mcn
    
    # 67
    n67 <- o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 68
    n68 <- o_mcn * f_mcn * h_bon_mcn * r_nat_trib * f_bon * l_bon
    
    # 69
    n69 <- o_mcn * f_mcn * h_bon_mcn * r_nat_trib * h_bon_mcn * r_nat_trib * f_bon * o_bon * l_bon_mcn
    
    # 70
    n70 <- o_mcn * f_mcn * h_bon_mcn * r_nat_trib * l_bon_mcn
    
    # 71
    n71 <- o_mcn * f_mcn * l_bon_mcn
    
    # 72
    n72 <- o_mcn * f_mcn * o_mcn * f_mcn * f_bon * l_bon
    
    # 73
    n73 <- o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn
    
    # 74
    n74 <- o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 75
    n75 <- o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * r_nat_trib * h_bon_mcn * l_nat_trib
    
    # 76
    n76 <- o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * r_nat_trib * l_bon_mcn
    
    # 77
    n77 <- o_mcn * f_mcn * o_mcn * f_mcn * l_bon_mcn
    
    # 78
    n78 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 79
    n79 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 80
    n80 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 81
    n81 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 82
    n82 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 83
    n83 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 84
    n84 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * o_ich * f_ich * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 85
    n85 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 86
    n86 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 87
    n87 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * o_ich * l_ich_lgr
    
    # 88
    n88 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * o_ich * o_lgr * f_lgr * f_ich * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 89
    n89 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * o_ich * o_lgr * f_lgr * s_ich_lgr * l_ich_lgr_trib
    
    # 90
    n90 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * o_ich * o_lgr * l_lgr
    
    # 91
    n91 <- o_mcn * f_mcn * o_mcn * f_mcn * o_mcn * s_mcn_ich_pra * l_mcn_pra_ich_trib
    
    # 92
    n92 <- o_mcn * f_mcn * o_mcn * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 93
    n93 <- o_mcn * f_mcn * o_mcn * f_mcn * s_bon_mcn * r_bon_mcn_trib * f_bon * l_bon
    
    # 94
    n94 <- o_mcn * f_mcn * o_mcn * f_mcn * s_bon_mcn * r_bon_mcn_trib * h_bon_mcn * l_nat_trib
    
    # 95
    n95 <- o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 96
    n96 <- o_mcn * f_mcn * o_mcn * o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 97
    n97 <- o_mcn * f_mcn * o_mcn * o_ich * f_ich * f_mcn * l_bon_mcn
    
    # 98
    n98 <- o_mcn * f_mcn * o_mcn * o_ich * f_ich * f_mcn * o_mcn * o_ich * o_lgr * l_lgr
    
    # 99
    n99 <- o_mcn * f_mcn * o_mcn * o_ich * f_ich * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 100
    n100 <- o_mcn * f_mcn * o_mcn * o_ich * f_ich * o_ich * l_ich_lgr
    
    # 101
    n101 <- o_mcn * f_mcn * o_mcn * o_ich * l_ich_lgr
    
    # 102
    n102 <- o_mcn * f_mcn * o_mcn * o_ich * o_lgr * f_lgr * f_ich * f_mcn * l_bon_mcn
    
    # 103
    n103 <- o_mcn * f_mcn * o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * f_ich * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 104
    n104 <- o_mcn * f_mcn * o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * l_ich_lgr
    
    # 105
    n105 <- o_mcn * f_mcn * o_mcn * o_ich * o_lgr * f_lgr * o_lgr * l_lgr
    
    # 106
    n106 <- o_mcn * f_mcn * o_mcn * o_ich * o_lgr * f_lgr * o_lgr * s_ich_lgr * l_lgr_trib
    
    # 107
    n107 <- o_mcn * f_mcn * o_mcn * o_ich * o_lgr * l_lgr
    
    # 108
    n108 <- o_mcn * f_mcn * o_mcn * o_ich * o_lgr * s_ich_lgr * l_lgr_trib
    
    # 109
    n109 <- o_mcn * f_mcn * o_mcn * o_ich * o_lgr * s_ich_lgr * r_lgr_trib * f_lgr * l_ich_lgr
    
    # 110
    n110 <- o_mcn * f_mcn * o_mcn * o_ich * s_ich_lgr * r_lgr_trib * o_lgr * l_lgr
    
    # 111
    n111 <- o_mcn * f_mcn * o_mcn * o_pra * f_pra * f_mcn * h_bon_mcn * l_nat_trib
    
    # 112
    n112 <- o_mcn * f_mcn * o_mcn * o_pra * o_ris * l_ris_rre
    
    # 113
    n113 <- o_mcn * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 114
    n114 <- o_mcn * f_mcn * s_bon_mcn * r_bon_mcn_trib * h_bon_mcn * l_nat_trib
    
    # 115
    n115 <- o_mcn * f_mcn * s_bon_mcn * r_bon_mcn_trib * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 116
    n116 <- o_mcn * f_mcn * s_bon_mcn * r_bon_mcn_trib * o_mcn * l_mcn_ich_pra
    
    # 117
    n117 <- o_mcn * l_mcn_ich_pra
    
    # 118
    n118 <- o_mcn * o_ich * f_ich * f_mcn * f_bon * o_bon * o_mcn * l_mcn_ich_pra
    
    # 119
    n119 <- o_mcn * o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 120
    n120 <- o_mcn * o_ich * f_ich * f_mcn * h_bon_mcn * r_nat_trib * f_bon * l_bon
    
    # 121
    n121 <- o_mcn * o_ich * f_ich * f_mcn * l_bon_mcn
    
    # 122
    n122 <- o_mcn * o_ich * f_ich * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 123
    n123 <- o_mcn * o_ich * f_ich * f_mcn * o_mcn * o_ich * o_lgr * l_lgr
    
    # 124
    n124 <- o_mcn * o_ich * f_ich * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 125
    n125 <- o_mcn * o_ich * f_ich * o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 126
    n126 <- o_mcn * o_ich * f_ich * o_ich * f_ich * f_mcn * l_bon_mcn
    
    # 127
    n127 <- o_mcn * o_ich * f_ich * o_ich * f_ich * f_mcn * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    # 128
    n128 <- o_mcn * o_ich * f_ich * o_ich * f_ich * f_mcn * o_mcn * f_mcn * h_bon_mcn * r_nat_trib * f_bon * l_bon
    
    # 129
    n129 <- o_mcn * o_ich * f_ich * o_ich * f_ich * f_mcn * o_mcn * f_mcn * o_mcn * f_mcn * s_bon_mcn * r_bon_mcn_trib * f_bon * l_bon
    
    # 130
    n130 <- o_mcn * o_ich * f_ich * o_ich * f_ich * o_ich * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 131
    n131 <- o_mcn * o_ich * f_ich * o_ich * f_ich * o_ich * l_ich_lgr
    
    # 132
    n132 <- o_mcn * o_ich * f_ich * o_ich * f_ich * o_ich * o_lgr * l_lgr
    
    # 133
    n133 <- o_mcn * o_ich * f_ich * o_ich * l_ich_lgr
    
    # 134
    n134 <- o_mcn * o_ich * f_ich * o_ich * o_lgr * f_lgr * l_ich_lgr
    
    # 135
    n135 <- o_mcn * o_ich * f_ich * o_ich * o_lgr * f_lgr * o_lgr * l_lgr
    
    # 136
    n136 <- o_mcn * o_ich * f_ich * o_ich * o_lgr * l_lgr
    
    # 137
    n137 <- o_mcn * o_ich * l_ich_lgr
    
    # 138
    n138 <- o_mcn * o_ich * o_lgr * f_lgr * f_ich * f_mcn * h_bon_mcn * l_nat_trib
    
    # 139
    n139 <- o_mcn * o_ich * o_lgr * f_lgr * f_ich * f_mcn * l_bon_mcn
    
    # 140
    n140 <- o_mcn * o_ich * o_lgr * f_lgr * f_ich * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 141
    n141 <- o_mcn * o_ich * o_lgr * f_lgr * l_ich_lgr
    
    # 142
    n142 <- o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * f_ich * o_ich * o_lgr * l_lgr
    
    # 143
    n143 <- o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * f_ich * s_mcn_ich_pra * l_mcn_pra_ich_trib
    
    # 144
    n144 <- o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * l_ich_lgr
    
    # 145
    n145 <- o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * o_lgr * f_lgr * o_lgr * s_ich_lgr * l_lgr_trib
    
    # 146
    n146 <- o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * o_lgr * f_lgr * s_ich_lgr * r_lgr_trib * l_ich_lgr
    
    # 147
    n147 <- o_mcn * o_ich * o_lgr * f_lgr * o_lgr * f_lgr * o_lgr * l_lgr
    
    # 148
    n148 <- o_mcn * o_ich * o_lgr * f_lgr * o_lgr * l_lgr
    
    # 149
    n149 <- o_mcn * o_ich * o_lgr * f_lgr * s_ich_lgr * l_ich_lgr_trib
    
    # 150
    n150 <- o_mcn * o_ich * o_lgr * l_lgr
    
    # 151
    n151 <- o_mcn * o_ich * o_lgr * s_ich_lgr * l_lgr_trib
    
    # 152
    n152 <- o_mcn * o_ich * o_lgr * s_ich_lgr * r_lgr_trib * f_lgr * f_ich * f_mcn * f_bon * o_bon * f_bon * o_bon * o_mcn * o_ich * o_lgr * s_ich_lgr * r_lgr_trib * f_lgr * f_ich * o_ich * f_ich * f_mcn * f_bon * l_bon
    
    # 153
    n153 <- o_mcn * o_ich * o_lgr * s_ich_lgr * r_lgr_trib * f_lgr * f_ich * f_mcn * f_bon * o_bon * o_mcn * o_ich * o_lgr * s_ich_lgr * l_lgr_trib
    
    # 154
    n154 <- o_mcn * o_ich * o_lgr * s_ich_lgr * r_lgr_trib * s_ich_lgr * l_lgr_trib
    
    # 155
    n155 <- o_mcn * o_ich * s_ich_lgr * l_ich_lgr_trib
    
    # 156
    n156 <- o_mcn * o_pra * f_pra * f_mcn * h_bon_mcn * l_nat_trib
    
    # 157
    n157 <- o_mcn * o_pra * f_pra * o_ich * l_ich_lgr
    
    # 158
    n158 <- o_mcn * o_pra * f_pra * o_ich * o_lgr * l_lgr
    
    # 159
    n159 <- o_mcn * o_pra * f_pra * o_pra * f_pra * f_mcn * h_bon_mcn * l_nat_trib
    
    # 160
    n160 <- o_mcn * o_pra * f_pra * o_pra * f_pra * o_pra * l_pra_ris
    
    # 161
    n161 <- o_mcn * o_pra * l_pra_ris
    
    # 162
    n162 <- o_mcn * o_pra * o_ris * f_ris * o_ris * f_ris * f_pra * f_mcn * o_mcn * f_mcn * o_mcn * l_mcn_ich_pra
    
    # 163
    n163 <- o_mcn * o_pra * o_ris * o_rre * o_wel * l_wel
    
    # 164
    n164 <- o_mcn * s_mcn_ich_pra * l_mcn_pra_ich_trib
    
    # 165
    n165 <- o_mcn * s_mcn_ich_pra * r_mcn_pra_ich_trib * f_mcn * h_bon_mcn * l_nat_trib
    
    # 166
    n166 <- o_mcn * s_mcn_ich_pra * r_mcn_pra_ich_trib * f_mcn * s_bon_mcn * l_bon_mcn_trib
    
    # 167
    n167 <- s_bon_mcn * l_bon_mcn_trib
    
    # 168
    n168 <- s_bon_mcn * r_bon_mcn_trib * h_bon_mcn * l_nat_trib
    
    # 169
    n169 <- s_bon_mcn * r_bon_mcn_trib * o_mcn * f_mcn * h_bon_mcn * l_nat_trib
    
    
    
    # Get vector of probabilities for multinomial likelihood
    p <- c(n1, n2, n3, n4, n5, n6, 
           n7, n8, n9, n10, n11, n12, 
           n13, n14, n15, n16, n17, n18, 
           n19, n20, n21, n22, n23, n24, 
           n25, n26, n27, n28, n29, n30, 
           n31, n32, n33, n34, n35, n36, 
           n37, n38, n39, n40, n41, n42, 
           n43, n44, n45, n46, n47, n48, 
           n49, n50, n51, n52, n53, n54, 
           n55, n56, n57, n58, n59, n60, 
           n61, n62, n63, n64, n65, n66, 
           n67, n68, n69, n70, n71, n72, 
           n73, n74, n75, n76, n77, n78, 
           n79, n80, n81, n82, n83, n84, 
           n85, n86, n87, n88, n89, n90, 
           n91, n92, n93, n94, n95, n96, 
           n97, n98, n99, n100, n101, n102, 
           n103, n104, n105, n106, n107, n108, 
           n109, n110, n111, n112, n113, n114, 
           n115, n116, n117, n118, n119, n120, 
           n121, n122, n123, n124, n125, n126, 
           n127, n128, n129, n130, n131, n132, 
           n133, n134, n135, n136, n137, n138, 
           n139, n140, n141, n142, n143, n144, 
           n145, n146, n147, n148, n149, n150, 
           n151, n152, n153, n154, n155, n156, 
           n157, n158, n159, n160, n161, n162, 
           n163, n164, n165, n166, n167, n168, 
           n169)
    
    
    # dmultinom call
    
    # Add an if statement: If any parameters are negative, give negLL a crazy high value rather than call dmultinom
    if(all(p >= 0)){
      negLL <- -1* dmultinom(x = data$count, prob = p, log = TRUE)
    }
    else{
      negLL <- 9999999
    }
    
    return(negLL)
  }
  
  
  optim_results <- optim(par = optim_inits, data = data, fn = negLL, method = 'L-BFGS-B', hessian = FALSE, lower = 0.0001, upper = 0.9999)
  
  return(optim_results)
}






##### Fit multinomial MLE function #####

# Get inits for optim
set.seed(123)
# optim_inits <- runif(length(unique(JDR_params$param)), min = 0, max = 1)
optim_inits <- runif(length(unique(JDR_params$param)), min = 0, max = 0.2)

# Run model

# Run optim multiple times
optim_repeat <- function(niter){
  optim_outputs <- vector(mode = "list", length = niter)
  for (i in 1:niter){
    # Generate new starting values
    # Start in a place that you know won't be rejected
    optim_inits <- runif(length(unique(JDR_params$param)), min = 0, max = 0.2)
    optim_outputs[[i]] <- multistate_model(data = JDR_unique_probs, optim_inits = optim_inits)
  }
  
  return(optim_outputs)
}

optim_outputs <- optim_repeat(niter = 10)

# Examine all runs
# initialize an empty DF to store
optim_results <- data.frame(run = seq(1, 10, 1),
                            counts = rep(NA, 10),
                            NLL = rep(NA, 10))

# Store info on number of runs, NLL
for (i in 1:nrow(optim_results)){
  optim_output <- optim_outputs[[i]]
  # Store run number
  optim_results$run[i] <- i
  # Store number of times optim ran
  optim_results$counts[i] <- optim_output$counts[1]
  # Store NLL
  optim_results$NLL[i] <- optim_output$value
  # Store parameter values
  # for (j in 1:length(optim_output$par)){
  #   optim_results <- cbind(optim_results, optim_output$par[j])
  # }
}

# Store parameter estimates

# Initialize a df to store parameter estimates
optim_par_est <- data.frame(row.names = seq(1, 10, 1))

# For each parameter
for (i in 1:23){
  # For each optim run
  
  # Create an empty vector
  par_est <- rep(NA, 10)
  for (j in 1:nrow(optim_results)){
    optim_output <- optim_outputs[[j]]
    par_est[j] <- optim_output$par[i]
  }
  
  # Bind it to the df
  optim_par_est %>% 
    bind_cols(., par_est) -> optim_par_est
  
}

# Rename columns
colnames(optim_par_est) <- JDR_params$param

# Get the loss parameters
optim_par_est %>% 
  mutate(l_bon = 1 - o_bon) %>% 
  mutate(l_bon_mcn = 1 - (o_mcn + s_bon_mcn + h_bon_mcn + f_bon)) %>% 
  mutate(l_bon_mcn_trib = 1 - (r_bon_mcn_trib)) %>% 
  mutate(l_ich_lgr = 1 - (o_lgr + f_ich + s_ich_lgr)) %>% 
  mutate(l_ich_lgr_trib = 1) %>% 
  mutate(l_lgr = 1 - (f_lgr + s_lgr)) %>% 
  mutate(l_lgr_trib = 1 - (r_lgr_trib)) %>% 
  mutate(l_mcn_ich_pra = 1 - (f_mcn + o_pra + o_ich + s_mcn_ich_pra)) %>% 
  mutate(l_mcn_pra_ich_trib = 1 - (r_mcn_pra_ich_trib)) %>% 
  mutate(l_nat_trib = 1 - r_nat_trib) %>% 
  mutate(l_pra_ris = 1 - (f_pra + o_ris)) %>% 
  mutate(l_ris_rre = 1 - (o_rre + f_ris)) %>% 
  mutate(l_wel = 1) -> optim_par_est

# Check that parameter estimates make sense
sum(optim_par_est$f_bon, optim_par_est$o_mcn, optim_par_est$h_bon_mcn, optim_par_est$s_bon_mcn, optim_par_est$l_bon_mcn)
# Looks good so far

# Get mean values for each
as.data.frame(colMeans(optim_par_est)) %>% 
  rownames_to_column("Parameter") %>% 
  dplyr::rename(mean = `colMeans(optim_par_est)`)-> mean_optim_par_est

# Inspect PRA probs - estimate seems too high
pra_overshoot_ind <- JDR_complete_probs[grep("o_pra", JDR_complete_probs$prob),]$tag_code

# PRA overshoot inspect
subset(JDR_states, tag_code %in% pra_overshoot_ind) -> pra_overshoot_det_hist
# Nothing looks off here

# Look at the stepwise states
subset(JDR_stepwise_probabilities, state_1 == "mainstem, MCN to ICH or PRA") -> MCN_ICH_PRA_ind
sum(as.data.frame(table(MCN_ICH_PRA_ind$state_2))$Freq)
table(MCN_ICH_PRA_ind$state_2)
table(MCN_ICH_PRA_ind$probability)
table(MCN_ICH_PRA_ind$probability)/sum(as.data.frame(table(MCN_ICH_PRA_ind$state_2))$Freq)

24/1616

sum(optim_par_est$f_mcn, optim_par_est$o_pra, optim_par_est$o_ich, optim_par_est$s_mcn_ich_pra, optim_par_est$l_mcn_ich_pra)
# The model estimates do sum to 1
# However - the overshoot ICH value appears too high, and the loss parameter appears much too low (because o_pra and o_ich are too high)
# Well, I'm not sure that the o_ich parameter is much too high. It's pretty similar to Shelby's values

# Okay, 24 of 1616 overshot PRA. That percentage lines up more with what Shelby estimated (around 1.4 - 1.6 percent)
# So what explains the discrepancy with our numbers?
# Check model code
# I don't see anything obvious... could it be something to do with the branching?

