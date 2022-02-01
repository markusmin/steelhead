### 05 - Fit multistate

# This R script fits the multistate, bidirectional model

# Note: You must run scripts 01 to 03 to generate the detections used in this model


### Load libraries
library(here)
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
  bind_rows(., data.frame(state_1 = "mainstem, upstream of LGR", state_2 = "Upstream LGR tributaries", probability = "s_ich_lgr")) %>% # stray - no individuals strayed to these
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


