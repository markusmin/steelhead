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
  
  outputs <- list(JDR_state_matrix, final_fates)
  # return(final_fates)
  return(outputs)
  
}


# run simulation with 1 million fish
set.seed(123)
msm_est <- JDR_transition_simulation(nsim = 10000000, transition_matrix = JDR_transition_matrix)

##### Eigenvectors - analytical solution #####
# transpose it
JDR_transition_matrix <- t(JDR_transition_matrix)

final_matrix <- JDR_transition_matrix^100
init_vector <- c(0, 1000000, rep(0, 13))
final_matrix * init_vector

eigenvalues <- eigen(JDR_transition_matrix)$values
eigenvectors <- eigen(JDR_transition_matrix)$vectors
abs(eigenvalues)
max(abs(eigenvalues))
max(Re(eigenvalues))
Re(eigenvectors[,1])

# Try again with popbio library
library(popbio)
JDR_mat1 <- eigen.analysis(JDR_transition_matrix)
JDR_mat1$stable.stage
# This just gives that all individuals end up in the loss - not useful!
# How do you do eigen analysis with an absorbing state?



##### Compare results to tallies from detection history #####

JDR_det_hist <- read.csv(here::here("model_files", "JDR_det_hist.csv"), row.names = 1)

# Get the metadata
JDR_det_hist %>% 
  dplyr::select(-c(tag_code, start_time, end_time)) %>% 
  distinct(event_site_name, .keep_all = TRUE) -> JDR_event_site_metadata

columbia_dams <- c("Bonneville Adult Fishways (combined)", "McNary Adult Fishways (combined)", 
                   "PRA - Priest Rapids Adult","RIA - Rock Island Adult", 
                   "RRF - Rocky Reach Fishway", "WEA - Wells Dam, DCPUD Adult Ladders")
snake_dams <- c("ICH - Ice Harbor Dam (Combined)",  "Lower Granite Dam Adult Fishways (combined)")

##### Tributary detection sites

# Natal tributary detection sites
BON_MCN_natal_sites <- JDR_event_site_metadata$event_site_name[grep("John Day", JDR_event_site_metadata$event_site_basin_name)]

# Straying sites between BON and MCN:
# Basins: Descuhtes River
# Subbasins: Klickitat, Umatilla, Fifteenmile Creek
# grep: Hood in event_site_name, since "Middle Columbia-Hood" contains mainstem sites
# Note: When using grep, you need to take out the river mouth arrays
BON_MCN_stray_sites <- c(JDR_event_site_metadata$event_site_name[grep("Hood|Fifteenmile|Deschutes",
                                                                      JDR_event_site_metadata$event_site_basin_name)],
                         JDR_event_site_metadata$event_site_name[grep("Klickitat|Umatilla",
                                                                      JDR_event_site_metadata$event_site_subbasin_name)],
                         JDR_event_site_metadata$event_site_name[grep("Hood",
                                                                      JDR_event_site_metadata$event_site_name)]
)
# Take out river mouth arrays
BON_MCN_stray_sites <- BON_MCN_stray_sites[!grepl("mouth|Mouth", BON_MCN_stray_sites)]


# Straying sites between MCN and branches (PRA or ICH):
# Basins: Yakima
# Subbasins: Walla Walla
MCN_PRA_ICH_stray_sites <- c(JDR_event_site_metadata$event_site_name[grep("Yakima",
                                                                          JDR_event_site_metadata$event_site_basin_name)],
                             JDR_event_site_metadata$event_site_name[grep("Walla Walla",
                                                                          JDR_event_site_metadata$event_site_subbasin_name)])
# Take out river mouth arrays
MCN_PRA_ICH_stray_sites <- MCN_PRA_ICH_stray_sites[!grepl("mouth|Mouth", MCN_PRA_ICH_stray_sites)]

# Straying sites between ICH and LGR:
# Basins: NA
# Subbasins: NA
# Here, you have to grep the site name, because Little Goose and Lower Monumental
# are considered part of the "Lower Snake-Tucannon" subbasin
ICH_LGR_stray_sites <- c(JDR_event_site_metadata$event_site_name[grep("Tucannon",
                                                                      JDR_event_site_metadata$event_site_name)])

# Take out river mouth arrays
# Here I'm not sure if the Lower Tucannon River Array is a similar distance as a
# river mouth array
ICH_LGR_stray_sites <- ICH_LGR_stray_sites[!grepl("mouth|Mouth", ICH_LGR_stray_sites)]

# Straying sites upstream of LGR:
# Basins: Salmon
# Subbasins: Upper Grande Ronde, Clearwater, Lower Grande Ronde
LGR_upstream_stray_sites <- c(JDR_event_site_metadata$event_site_name[grep("Salmon",
                                                                           JDR_event_site_metadata$event_site_basin_name)],
                              JDR_event_site_metadata$event_site_name[grep("Upper Grande Ronde|Clearwater|Lower Grande Ronde",
                                                                           JDR_event_site_metadata$event_site_subbasin_name)])

# Take out river mouth arrays
# Don't do this for these sites - because there are creek mouths in here,
# i.e., "LAP - Lapwai Creek, near its mouth" and "SWT - Sweetwater Cr. near its mouth"  
# LGR_upstream_stray_sites <- LGR_upstream_stray_sites[!grepl("mouth|Mouth", LGR_upstream_stray_sites)]

##### In-river detection sites #####
# Figure out sites that are "in river" (in this I include river mouths), that 
# can indicate fallback depending on the order

pre_BON_inriver <- c("ESANIS - East Sand Island, Columbia River", "TWX - Estuary Towed Array (Exp.)")
# This includes some dams that we are ignoring
BON_MCN_inriver <- c("COLR4 - Columbia River - Bonneville Dam to John Day Dam (km 234-347)",
                     "The Dalles Adult Fishways (combined)", "DRM - Deschutes River mouth",
                     "HRM - Hood River Mouth", "JDJ - John Day Dam Juvenile")
MCN_ICH_PRA_inriver <- c("BADGEI - Badger Island, Columbia River")
# This includes dams (LMO and LGO) that we are ignoring
# These are some fallback routes - but because we aren't looking at these
# dams specifically, we can consider them in river detection sites
ICH_LGR_inriver <- c("GRJ - Lower Granite Dam Juvenile", "GOJ - Little Goose Dam Juvenile",
                     "LMJ - Lower Monumental Dam Juvenile",
                     "LGRTAL - LGR - Release into the Tailrace within 0.5 km downstream of Dam",
                     "LMA - Lower Monumental Adult Ladders", "GOA - Little Goose Fish Ladder")

##### Sites that are confirmed fallback (only certain dams) ####
BON_fallback_arrays <- c("BCC - BON PH2 Corner Collector", "B2J - Bonneville PH2 Juvenile")
MCN_fallback_arrays <- c("MCJ - McNary Dam Juvenile")
LGR_fallback_arrays <- c("GRJ - Lower Granite Dam Juvenile")


# Confirm that all sites have been categorized
setdiff(JDR_event_site_metadata$event_site_name, 
        c(columbia_dams, snake_dams, # Adult fishways at dams
          LGR_upstream_stray_sites, ICH_LGR_stray_sites,MCN_PRA_ICH_stray_sites, 
          BON_MCN_stray_sites,BON_MCN_natal_sites, # Tributary sites
          pre_BON_inriver, BON_MCN_inriver, MCN_ICH_PRA_inriver, ICH_LGR_inriver, #in river arrays
          BON_fallback_arrays, MCN_fallback_arrays, LGR_fallback_arrays)) # fallback arrays

# Turn this into a dataframe (this was a silly way to do this)
JDR_site_classification <- data.frame(event_site_name = c(columbia_dams, snake_dams, # Adult fishways at dams
                                                          LGR_upstream_stray_sites, ICH_LGR_stray_sites,MCN_PRA_ICH_stray_sites, 
                                                          BON_MCN_stray_sites,BON_MCN_natal_sites, # Tributary sites
                                                          pre_BON_inriver, BON_MCN_inriver, MCN_ICH_PRA_inriver, ICH_LGR_inriver, #in river arrays
                                                          BON_fallback_arrays, MCN_fallback_arrays, LGR_fallback_arrays, "lost"),
                                      site_class = c(c("BON (adult)", "MCN (adult)", "PRA (adult)", "RIS (adult)",
                                                       "RRE (adult)", "WEL (adult)"),
                                                     c("ICH (adult)", "LGR (adult)"), 
                                                     rep("LGR_upstream_stray_sites", length(LGR_upstream_stray_sites)),
                                                     rep("ICH_LGR_stray_sites", length(ICH_LGR_stray_sites)),
                                                     rep("MCN_PRA_ICH_stray_sites", length(MCN_PRA_ICH_stray_sites)),
                                                     rep("BON_MCN_stray_sites", length(BON_MCN_stray_sites)),
                                                     rep("BON_MCN_natal_sites", length(BON_MCN_natal_sites)),
                                                     rep("pre_BON_inriver", length(pre_BON_inriver)), 
                                                     rep("BON_MCN_inriver", length(BON_MCN_inriver)),
                                                     rep("MCN_ICH_PRA_inriver", length(MCN_ICH_PRA_inriver)),
                                                     rep("ICH_LGR_inriver", length(ICH_LGR_inriver)), 
                                                     rep("BON_fallback_arrays", length(BON_fallback_arrays)),
                                                     rep("MCN_fallback_arrays", length(MCN_fallback_arrays)),
                                                     rep("LGR_fallback_arrays", length(LGR_fallback_arrays)), "lost"))

# Create a new dataframe for state (location in the system)
JDR_site_classification %>%
  mutate(
    state = ifelse(
      site_class == "BON (adult)", "mainstem, BON to MCN",
      ifelse(
        site_class == "MCN (adult)", "mainstem, MCN to ICH or PRA",
        ifelse(
          site_class == "PRA (adult)", "mainstem, PRA to RIS",
          ifelse(
            site_class == "RIS (adult)", "mainstem, RIS to RRE",
            ifelse(
              site_class == "RRE (adult)", "mainstem, RRE to WEL",
              ifelse(
                site_class == "WEL (adult)", "mainstem, upstream of WEL",
                ifelse(
                  site_class == "ICH (adult)", "mainstem, ICH to LGR",
                  ifelse(
                    site_class == "LGR (adult)", "mainstem, upstream of LGR",
                    ifelse(
                      site_class == "LGR_upstream_stray_sites", "Upstream LGR tributaries",
                      ifelse(
                        site_class == "ICH_LGR_stray_sites", "ICH to LGR tributaries",
                        ifelse(
                          site_class == "MCN_PRA_ICH_stray_sites", "MCN to PRA or ICH tributaries",
                          ifelse(
                            site_class == "BON_MCN_stray_sites", "BON to MCN tributaries",
                            ifelse(
                              site_class == "BON_MCN_natal_sites", "natal tributaries",
                              ifelse(
                                site_class == "pre_BON_inriver", "mainstem, mouth to BON",
                                ifelse(
                                  site_class == "BON_MCN_inriver", "mainstem, BON to MCN",
                                  ifelse(
                                    site_class == "MCN_ICH_PRA_inriver", "mainstem, MCN to ICH or PRA",
                                    ifelse(
                                      site_class == "ICH_LGR_inriver", "mainstem, ICH to LGR",
                                      ifelse(
                                        site_class == "BON_fallback_arrays", "mainstem, mouth to BON",
                                        ifelse(
                                          site_class == "MCN_fallback_arrays", "mainstem, BON to MCN",
                                          ifelse(site_class == "LGR_fallback_arrays", "mainstem, ICH to LGR", 
                                                 ifelse(site_class == "lost", "lost", NA)
                                          ))))))))))))))))))))) -> JDR_site_classification

# Add info to detection history
JDR_det_hist %>% 
  left_join(., JDR_site_classification, by = "event_site_name") -> JDR_det_hist

JDR_det_hist %>% 
  group_by(tag_code) %>% 
  filter(row_number() == n()) -> JDR_final_state

as.data.frame(table(JDR_final_state$state)) %>% 
  dplyr::rename(state = Var1, count = Freq) %>% 
  mutate(prop = count/nrow(JDR_final_state)) -> JDR_final_state_counts


##### Compare tallies and multistate estimates

msm_est %>% 
  rownames_to_column("state_abb") %>% 
  mutate(state = ifelse(state_abb == "mouth_bon", "mainstem, mouth to BON",
                        ifelse(state_abb == "bon_mcn", "mainstem, BON to MCN",
                               ifelse(state_abb == "bon_mcn_trib", "BON to MCN tributaries",
                                      ifelse(state_abb == "ich_lgr", "mainstem, ICH to LGR",
                                             ifelse(state_abb == "ich_lgr_trib", "ICH to LGR tributaries",
                                                    ifelse(state_abb == "ups_lgr", "mainstem, upstream of LGR",
                                                           ifelse(state_abb == "ups_lgr_trib", "Upstream LGR tributaries",
                                                                  ifelse(state_abb == "mcn_ich_pra", "mainstem, MCN to ICH or PRA",
                                                                         ifelse(state_abb == "mcn_ich_pra_trib", "MCN to PRA or ICH tributaries",
                                                                                ifelse(state_abb == "nat_trib", "natal tributaries",
                                                                                       ifelse(state_abb == "pra_ris", "mainstem, PRA to RIS",
                                                                                              ifelse(state_abb == "ris_rre", "mainstem, RIS to RRE",
                                                                                                     ifelse(state_abb == "rre_wel", "mainstem, RRE to WEL",
                                                                                                            ifelse(state_abb == "ups_wel", "mainstem, upstream of WEL", NA))))))))))))))) %>% 
  dplyr::rename(model_prop = prop) %>% 
  dplyr::select(-count) %>% 
  left_join(., JDR_final_state_counts, by = "state") %>% 
  dplyr::select(-c(state_abb, count)) %>% 
  relocate(state) %>% 
  dplyr::rename( tally_prop = prop) %>% 
  # fill in the zeros from the tally counts 
  mutate(tally_prop = ifelse(is.na(tally_prop), 0, tally_prop)) %>% 
  mutate(model_v_tally = model_prop - tally_prop) -> state_comp
state_comp

# See which ones are more than 1% off
subset(state_comp, abs(model_v_tally) > 0.01)


##### Troubleshoot differences #####

# Load state data
JDR_stepwise_probabilities <- read.csv(here::here("model_files", "JDR_stepwise_probabilities.csv"), row.names = 1)

# Group by state 1 and see where individuals moved
# JDR_stepwise_probabilities %>% 
#   group_by(state_1) 




