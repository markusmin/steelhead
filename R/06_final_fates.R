### 06 - final fates

# Description: This R script uses a transition matrix to convert movement probabilities into final fates.

library(here)
library(tidyverse)

##### Step 1: Populate the transition matrix #####

param_table <- read.csv(here::here("model_files", "JDR_param_estimates.csv"), row.names = 1)
# Assign median estimates to parameter names
for (i in 1:nrow(param_table)){
  assign(param_table$Param[i], param_table$median[i])
}

# populate with zeros
JDR_transition_matrix <- matrix(rep(0, 15*15), nrow = 15, ncol = 15)
  
JDR_transition_states <- c("mouth_bon", "bon_mcn", "bon_mcn_trib", 
                           "ich_lgr", "ich_lgr_trib", "ups_lgr", 
                           "ups_lgr_trib", "mcn_ich_pra", "mcn_ich_pra_trib",
                           "nat_trib", "pra_ris", "ris_rre",
                           "rre_wel", "ups_wel", "loss")

colnames(JDR_transition_matrix) <- JDR_transition_states
rownames(JDR_transition_matrix) <- JDR_transition_states

# fill in the non-zero probabilities

### 1: Mouth to BON
JDR_transition_matrix['mouth_bon','bon_mcn'] <- o_bon
JDR_transition_matrix['mouth_bon','loss'] <- l_bon


### 2: BON to MCN
JDR_transition_matrix['bon_mcn','mcn_ich_pra'] <- o_mcn
JDR_transition_matrix['bon_mcn','bon_mcn_trib'] <- s_bon_mcn
JDR_transition_matrix['bon_mcn','nat_trib'] <- h_bon_mcn
JDR_transition_matrix['bon_mcn','mouth_bon'] <- f_bon
JDR_transition_matrix['bon_mcn','loss'] <- l_bon_mcn


### 3: BON to MCN tributaries
JDR_transition_matrix['bon_mcn_trib','bon_mcn'] <- r_bon_mcn_trib
JDR_transition_matrix['bon_mcn_trib','loss'] <- l_bon_mcn_trib


### 4: ICH to LGR
JDR_transition_matrix['ich_lgr','ups_lgr'] <- o_lgr
JDR_transition_matrix['ich_lgr','mcn_ich_pra'] <- f_ich
JDR_transition_matrix['ich_lgr','ich_lgr_trib'] <- s_ich_lgr
JDR_transition_matrix['ich_lgr','loss'] <- l_ich_lgr


### 5: ICH to LGR tributaries
# l_ich_lgr_trib <- 1 # Note: This would be 1 - r_ich_lgr_trib, but no fish returned after straying into these tributaries
JDR_transition_matrix['ich_lgr_trib','loss'] <- l_ich_lgr_trib


### 6: Upstream of LGR
JDR_transition_matrix['ups_lgr','ich_lgr'] <- f_lgr
JDR_transition_matrix['ups_lgr','ups_lgr_trib'] <- s_lgr
JDR_transition_matrix['ups_lgr','loss'] <- l_ich_lgr_trib


### 7: Upstream of LGR tributaries
JDR_transition_matrix['ups_lgr_trib','ups_lgr'] <- r_lgr_trib
JDR_transition_matrix['ups_lgr_trib','loss'] <- l_lgr_trib


### 8: MCN to ICH or PRA
JDR_transition_matrix['mcn_ich_pra','bon_mcn'] <- f_mcn
JDR_transition_matrix['mcn_ich_pra','pra_ris'] <- o_pra
JDR_transition_matrix['mcn_ich_pra','ich_lgr'] <- o_ich
JDR_transition_matrix['mcn_ich_pra','mcn_ich_pra_trib'] <- s_mcn_ich_pra
JDR_transition_matrix['mcn_ich_pra','loss'] <- l_mcn_ich_pra


### 9: MCN to ICH or PRA tributaries
JDR_transition_matrix['mcn_ich_pra_trib','mcn_ich_pra'] <- r_mcn_pra_ich_trib
JDR_transition_matrix['mcn_ich_pra_trib','loss'] <- l_mcn_pra_ich_trib


### 10: Natal tributaries
JDR_transition_matrix['nat_trib','bon_mcn'] <- r_nat_trib
JDR_transition_matrix['nat_trib','loss'] <- l_nat_trib


### 11: PRA to RIS
JDR_transition_matrix['pra_ris','mcn_ich_pra'] <- f_pra
JDR_transition_matrix['pra_ris','ris_rre'] <- o_ris
JDR_transition_matrix['pra_ris','loss'] <- l_pra_ris


### 12: RIS to RRE
JDR_transition_matrix['ris_rre','rre_wel'] <- o_rre
JDR_transition_matrix['ris_rre','pra_ris'] <- f_ris
JDR_transition_matrix['ris_rre','loss'] <- l_ris_rre


### 13: RRE to WEL
JDR_transition_matrix['rre_wel','ups_wel'] <- o_wel
JDR_transition_matrix['rre_wel','loss'] <- l_rre_wel
# l_rre_wel <- 1 - (o_wel + f_rre)
# No individuals ever fell back over RRE. So we will tweak this to remove it


### 14: Upstream of WEL
# l_wel <- 1 # Note: This would be 1 - f_wel, but no fish fell back over Well's dam
JDR_transition_matrix['ups_wel','loss'] <- l_wel

# Loss to loss = 1
JDR_transition_matrix['loss', 'loss'] <- 1


##### Step 2: Simulation testing #####

# Start 1 million individuals in the mouth to BON state. Where do they end up?

JDR_transition_simulation <- function(nsim, transition_matrix){
  trans <- transition_matrix
  
  # Create an empty state matrix
  JDR_state_matrix <- matrix(rep(0, 15), nrow = 15, ncol = 1)
  
  JDR_transition_states <- c("mouth_bon", "bon_mcn", "bon_mcn_trib", 
                             "ich_lgr", "ich_lgr_trib", "ups_lgr", 
                             "ups_lgr_trib", "mcn_ich_pra", "mcn_ich_pra_trib",
                             "nat_trib", "pra_ris", "ris_rre",
                             "rre_wel", "ups_wel", "loss")
  
  rownames(JDR_state_matrix) <- JDR_transition_states
  
  JDR_state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> JDR_state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  JDR_state_matrix[2,2] <- nsim
    
  # Loop through 100 iterations, they should all reach a final state  
  for (i in 2:100){
    # First, append a row for the next states
    JDR_state_matrix %>% 
      bind_cols(., t = rep(0, 15)) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> JDR_state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, 15*15), ncol = 15, nrow = 15)
    rownames(movement_matrix) <- JDR_transition_states
    colnames(movement_matrix) <- JDR_transition_states
    for (j in 1:nrow(JDR_state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = JDR_state_matrix[j,i], prob = trans[j,])
    }
    # sum across rows
    JDR_state_matrix[,i+1] <- rowSums(movement_matrix)
    
  }
  # or: use a while loop, keep looping until condition of all being in the loss state is reached
  
}


##### V2 ######

JDR_transition_simulation <- function(nsim, transition_matrix){
  trans <- transition_matrix
  
  # Create an empty state matrix
  JDR_state_matrix <- matrix(rep(0, 15), nrow = 15, ncol = 1)
  
  JDR_transition_states <- c("mouth_bon", "bon_mcn", "bon_mcn_trib", 
                             "ich_lgr", "ich_lgr_trib", "ups_lgr", 
                             "ups_lgr_trib", "mcn_ich_pra", "mcn_ich_pra_trib",
                             "nat_trib", "pra_ris", "ris_rre",
                             "rre_wel", "ups_wel", "loss")
  
  rownames(JDR_state_matrix) <- JDR_transition_states
  
  JDR_state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> JDR_state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  JDR_state_matrix[2,2] <- nsim
  
  # create a final fate matrix
  JDR_final_fate_matrix <- matrix(rep(0, 15), nrow = 1, ncol = 15)
  colnames(JDR_final_fate_matrix) <- JDR_transition_states
  
  # run while loop until all fish have been lost
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    JDR_state_matrix %>% 
      bind_cols(., t = rep(0, 15)) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> JDR_state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, 15*15), ncol = 15, nrow = 15)
    rownames(movement_matrix) <- JDR_transition_states
    colnames(movement_matrix) <- JDR_transition_states
    for (j in 1:nrow(JDR_state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = JDR_state_matrix[j,i], prob = trans[j,])
    }
    # sum across rows, store in the state matrix
    JDR_state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    JDR_final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> JDR_final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (JDR_state_matrix[15, ncol(JDR_state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  data.frame(colSums(JDR_final_fate_matrix[,1:14])) %>% 
    dplyr::rename(count = colSums.JDR_final_fate_matrix...1.14..) %>% 
    mutate(prop = count/nsim)-> final_fates
  
  return(final_fates)
  
}


# run simulation with 1 million fish
set.seed(123)
JDR_transition_simulation(nsim = 1000000, transition_matrix = JDR_transition_matrix)
