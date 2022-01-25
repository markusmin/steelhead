### 03 - Fit multistate model

# NOTE: To fit this model, you will first need to run 02_detection_histories.R
# for each of your river systems, with each output named "<river_system>_det_hist"


### Load libraries
library(here)
library(tidyverse)
library(lubridate)
library(janitor)

# Look at counts by site
JDR_site_count <- as.data.frame(table(JDR_det_hist$event_site_name))

##### Assign detection sites to parts of river ####

##### Dams in order - separate Columbia and Snake #####
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

# Inspect detection histories after tributary detections
# Do we want to remove this? Probably indicates post-spawning behavior
JDR_det_hist %>% 
  group_by(tag_code) %>% 
  subset(state %in% tributary_mainstem$tributary) %>% 
  filter(end_time == max(end_time)) %>% 
  dplyr::select(tag_code, end_time) %>% 
  dplyr::rename(spawning = end_time) -> spawning_times

JDR_det_hist %>% 
  left_join(., spawning_times, by = "tag_code") %>%  
  group_by(tag_code) %>% 
  filter(end_time > spawning) -> post_spawn_history

# A lot of these do seem like iteroparous individuals migrating back downstream. But not all of them.
# Would bias our fallback estimates to include them! But we need to be flexible enough to 
# include this behavior, since some of it is tributary dip-ins


# Turn this into state history, not detection history.
# This needs to include that is implicit in the detection history; i.e., if
# a fish was seen in the BON adult ladders in consecutive detection events,
# it must have fallen back over and thus entered the pre BON mainstem state

# So this code currently works to detect some implicit fallback, i.e., detection
# at the same dam twice in a row. However, it isn't yet able to capture implicit
# fallback when it is detected in a tributary. For example, if a fish is seen
# at McNary, and then in the John Day, it doesn't register that it was in the
# mainstem between BON and MCN between those times. But perhaps we can
# change that after the fact?
# Between two dams, it should always be BON -> mainstem -> (tributary) -> mainstem -> MCN
# Don't have to go to tributary, but if they do they must go to mainstem before or after

# Need a function that can insert rows into the DF for the implicit states

# Order of sites (no tributaries), Columbia River
site_order_notrib_columbia <- c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                                "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                "mainstem, upstream of WEL")
# Order of sites (including tributaries), Columbia River, straying tribs
site_order_tribs_columbia_stray <- c("mainstem, mouth to BON", "mainstem, BON to MCN",
                                     "BON to MCN tributaries", "mainstem, BON to MCN",
                                     "mainstem, MCN to ICH or PRA", "MCN to PRA or ICH tributaries",
                               "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                               "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                               "mainstem, upstream of WEL")

site_order_tribs_columbia_natal <- c("mainstem, mouth to BON", "mainstem, BON to MCN",
                                     "natal tributaries", "mainstem, BON to MCN",
                                     "mainstem, MCN to ICH or PRA", "MCN to PRA or ICH tributaries",
                                     "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                     "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                     "mainstem, upstream of WEL")

# Order of sites (no tributaries), Snake
# site_order_notrib_snake <- c("mainstem, ICH to LGR", "mainstem, upstream of LGR")
# Contains all non-tributary sites from upstream of LGR to downstream of BON
site_order_notrib_snake <- c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                            "mainstem, MCN to ICH or PRA",
                            "mainstem, ICH to LGR",
                            "mainstem, upstream of LGR")
# Order of sites (includes tributaries), Snake
site_order_tribs_snake <- c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                            "mainstem, MCN to ICH or PRA", "mainstem, ICH to LGR", 
                             "ICH to LGR tributaries", "mainstem, ICH to LGR",
                             "mainstem, upstream of LGR", "Upstream LGR tributaries",
                             "mainstem, upstream of LGR")


# Site order, no tributaries, Snake directly to upper Columbia
site_order_notrib_columbia_snake <- c("mainstem, upstream of WEL",
                                      "mainstem, RRE to WEL", 
                                      "mainstem, RIS to RRE",
                                      "mainstem, PRA to RIS",
                             "mainstem, MCN to ICH or PRA",
                             "mainstem, ICH to LGR",
                             "mainstem, upstream of LGR")


# Upper Columbia sites
upper_columbia_sites <- c("mainstem, PRA to RIS",
                                     "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                     "mainstem, upstream of WEL")

# Snake River sites
# snake_sites <- c("mainstem, ICH to LGR", 
#                             "ICH to LGR tributaries", "mainstem, ICH to LGR",
#                             "mainstem, upstream of LGR", "Upstream LGR tributaries",
#                             "mainstem, upstream of LGR")
snake_sites <- c("mainstem, ICH to LGR", 
                 "mainstem, upstream of LGR")


# Shared sites (lower Columbia, before confluence of Snake and Columbia)
# shared_sites <- c("mainstem, mouth to BON", "mainstem, BON to MCN",
#                                      "natal tributaries",
#                                      "BON to MCN tributaries",
#                   "MCN to PRA or ICH tributaries",
#                                      "mainstem, MCN to ICH or PRA")
shared_sites <- c("mainstem, mouth to BON", "mainstem, BON to MCN",
                                     "mainstem, MCN to ICH or PRA")

# Tributary df - tributaries and what part of the mainstem they're in
# We will use this to determine where they need to move to next, and then
# use the no tributary site order to fill in the rest of the history
tributary_mainstem <- data.frame(tributary = c("natal tributaries", 
                                               "BON to MCN tributaries", 
                                               "MCN to PRA or ICH tributaries",
                                               "ICH to LGR tributaries",
                                               "Upstream LGR tributaries"),
                                 mainstem = c("mainstem, BON to MCN",
                                              "mainstem, BON to MCN",
                                              "mainstem, MCN to ICH or PRA",
                                              "mainstem, ICH to LGR",
                                              "mainstem, upstream of LGR"))


# Make a function to insert rows for implicit detections
insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  return(existingDF)
}


# Create empty df to store states
JDR_stepwise_states <- data.frame(tag_code = character(),
                                      state = character(),
                                      date_time = as.POSIXct(character()))


# Before running this, save the original detection history
# JDR_det_hist_original <- JDR_det_hist
# Reset detection history
# JDR_det_hist <- JDR_det_hist_original


# Add a counter to keep track of the number of added rows.
# This isn't actually needed for the code after all, but is still interesting
added_rows <- 0

# Our for loop is too long, and will give us an error. However, we need to
# make sure that we are accounting for the rows that we are inserting into
# our dataframe during this for loop, so to be safe we are doubling the 
# length of the for loop. It will run fine, it will just give an error once
# the dataframe ends, but that doesn't matter for our purposes.
# NOTE: If you don't get an error, that means your for loop isn't long enough
# and you need to extend it
# The error should read: 
# "Error in if (JDR_det_hist[i, "tag_code"] == JDR_det_hist[i - 1, "tag_code"]) { : 
# missing value where TRUE/FALSE needed

for (i in 1:(2 * nrow(JDR_det_hist) - 1)) {
  # Update i to keep up with number of added rows
  # NOT NECESSARY
  # i <- i + added_rows
  
  # If it's the first entry, store it
  if (i == 1) {
    # Store the tag code
    JDR_stepwise_states[i, 'tag_code'] <- JDR_det_hist[i, 'tag_code']
    # Store the current state
    JDR_stepwise_states[i, 'state'] <- JDR_det_hist[i, 'state']
    # Store the time entering this state
    JDR_stepwise_states[i, 'date_time'] <-
      JDR_det_hist[i, 'end_time']
  }
  
  # For all other entries, look at the previous entry.
  # If it's the same fish as the previous entry:
  else if (JDR_det_hist[i, 'tag_code'] == JDR_det_hist[i - 1, 'tag_code']) {
    # If it's in a different state than it was previously, record it
    if (JDR_det_hist[i, 'state'] != JDR_det_hist[i - 1, 'state']) {
      ### NEW CONDITIONS: COMPLETE DETECTION HISTORY (NOT REPEAT DAM DET)
      # We have to account for what happens if we see a fish at an
      # site, and then next at a different site while
      # skipping the sites in between (could be upstream or downstream).
      
      
      # if it's in the right part of the mainstem before, then we skip it.
      # else (anything else) means it's not. In that case, we have to determine
      # what part of the mainstem it needs to have been in to get to that tributary,
      # then use that information to fill in the missing sites from the order of
      # sites DFs
      
      ### TRIBUTARIES
      # If it's in a tributary and wasn't in the right part of the mainstem prev:
      if (JDR_det_hist[i, 'state'] %in% tributary_mainstem$tributary) {
        if (JDR_det_hist[i, 'state'] %in% tributary_mainstem$tributary &
            (JDR_det_hist[i - 1, 'state'] !=  subset(tributary_mainstem, tributary ==
                                                    JDR_det_hist[i, 'state'])$mainstem)) {
          mainstem_site <-
            subset(tributary_mainstem, tributary == JDR_det_hist[i, 'state'])$mainstem
          
          # INSERT ROW FOR THE MAINSTEM SITE
          # Insert a new row into stepwise states, with the implicit detection site
          # Tag code, state, and time (which is NA)
          implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                              mainstem_site,
                              NA)
          
          JDR_stepwise_states <-
            insertRow(JDR_stepwise_states, implicit_state, i)
          
          # Insert a new row into original detection history, with
          # implicit detection site info
          implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                                  NA,
                                  NA,
                                  NA,
                                  NA,
                                  NA,
                                  NA,
                                  NA,
                                  "implicit site",
                                  mainstem_site)
          
          # Also, change the value in the original dataframe
          JDR_det_hist <-
            insertRow(JDR_det_hist, implicit_detection, i)
          
          # Add 1 to the number of added rows
          added_rows <- added_rows + 1
          
          # After this, the code should recognize if there are any missing sites in between

        }          # If it was in the right part of the mainstem previously, store it like normal
        else {
          # Store the tag code
          JDR_stepwise_states[i, 'tag_code'] <-
            JDR_det_hist[i, 'tag_code']
          # Store the current state
          JDR_stepwise_states[i, 'state'] <-
            JDR_det_hist[i, 'state']
          # Store the time entering this state
          JDR_stepwise_states[i, 'date_time'] <-
            JDR_det_hist[i, 'end_time']
        }
        
        
      }
      
      ### MAINSTEM SITES
      else{
        
        # If it was previously in a tributary and this observation isn't
        # in the corresponding mainstem segment, add a line to the detection
        # history
        if (JDR_det_hist[i-1, 'state'] %in% tributary_mainstem$tributary){
          if (JDR_det_hist[i, 'state'] !=  subset(tributary_mainstem, tributary ==
                                                  JDR_det_hist[i-1, 'state'])$mainstem){
            
            mainstem_site <-
              subset(tributary_mainstem, tributary == JDR_det_hist[i-1, 'state'])$mainstem
            
            # INSERT ROW FOR THE MAINSTEM SITE
            # Insert a new row into stepwise states, with the implicit detection site
            # Tag code, state, and time (which is NA)
            implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                                mainstem_site,
                                NA)
            
            JDR_stepwise_states <-
              insertRow(JDR_stepwise_states, implicit_state, i)
            
            # Insert a new row into original detection history, with
            # implicit detection site info
            implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    "implicit site",
                                    mainstem_site)
            
            # Also, change the value in the original dataframe
            JDR_det_hist <-
              insertRow(JDR_det_hist, implicit_detection, i)
            
            # Add 1 to the number of added rows
            added_rows <- added_rows + 1
            
          }
        } 
        
        # Zeroth: Jumping between Snake and Columbia, without visiting 
        # shared sites in between. We see this with some fish, if they are
        # seen at PRA and then ICH or vice versa, without anything in between
        else if (JDR_det_hist[i-1, 'state'] %in% upper_columbia_sites &
                 JDR_det_hist[i, 'state'] %in% snake_sites |
                 JDR_det_hist[i-1, 'state'] %in% snake_sites &
                 JDR_det_hist[i, 'state'] %in% upper_columbia_sites){
          
          # These inherently have to skip, since they're not seen
          # at a shared site (MCN to ICH/PRA is the only route, and by
          # definition they're not seen there)
          # Get the index of the current site
          current_index <-
            which(site_order_notrib_columbia_snake %in% JDR_det_hist[i, 'state'])
          # Get the index of the previous site
          previous_index <-
            which(site_order_notrib_columbia_snake %in% JDR_det_hist[i - 1, 'state'])
          
          # sequence from current to previous index
          if (current_index < previous_index) {
            index_order <- seq(current_index, previous_index, by = 1)
            
          } else {
            index_order <- seq(previous_index, current_index, by = 1)
          }
          
          # Count the number of sites you need to add and loop through
          for (j in 1:(length(index_order) - 2)) {
            # Add a row for each of them
            
            # Insert a new row into stepwise states, with the implicit detection site
            # Tag code, state, and time (which is NA)
            missing_site <-
              site_order_notrib_columbia_snake[index_order[1 + j]]
            implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                                missing_site,
                                NA)
            
            JDR_stepwise_states <-
              insertRow(JDR_stepwise_states, implicit_state, i)
            
            # Insert a new row into original detection history, with
            # implicit detection site info
            implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    NA,
                                    "implicit site",
                                    missing_site)
            
            # Also, change the value in the original dataframe
            JDR_det_hist <-
              insertRow(JDR_det_hist, implicit_detection, i)
            
            # Add 1 to the number of added rows
            added_rows <- added_rows + 1
          }
          
          
        }
        
        
        # First: Shared sites - need two options, depending on where it previously was
        # Shared sites can use either Columbia or Snake sites
        # Condition: If the indices are off by more than one (plus or minus one)
        else if (JDR_det_hist[i, 'state'] %in% shared_sites) {
          if (JDR_det_hist[i-1, 'state'] %in% c(shared_sites, upper_columbia_sites)){
            if (which(site_order_notrib_columbia %in% JDR_det_hist[i, 'state']) >
                which(site_order_notrib_columbia %in% JDR_det_hist[i - 1, 'state']) + 1 |
                which(site_order_notrib_columbia %in% JDR_det_hist[i, 'state']) <
                which(site_order_notrib_columbia %in% JDR_det_hist[i - 1, 'state']) - 1){
              # If we see it at two non-consecutive mainstem locations
              
              # If the next site skips sites, insert the missing sites
              # Get the index of the current site
              current_index <-
                which(site_order_notrib_columbia %in% JDR_det_hist[i, 'state'])
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_columbia %in% JDR_det_hist[i - 1, 'state'])
              
              # sequence from current to previous index
              if (current_index < previous_index) {
                index_order <- seq(current_index, previous_index, by = 1)
                
              } else {
                index_order <- seq(previous_index, current_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order) - 2)) {
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_columbia[index_order[1 + j]]
                implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                                    missing_site,
                                    NA)
                
                JDR_stepwise_states <-
                  insertRow(JDR_stepwise_states, implicit_state, i)
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit site",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                JDR_det_hist <-
                  insertRow(JDR_det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            # Finally, if we're not missing any sites in between, just store it
            else {
              # Store the tag code
              JDR_stepwise_states[i, 'tag_code'] <-
                JDR_det_hist[i, 'tag_code']
              # Store the current state
              JDR_stepwise_states[i, 'state'] <- JDR_det_hist[i, 'state']
              # Store the time entering this state
              JDR_stepwise_states[i, 'date_time'] <-
                JDR_det_hist[i, 'end_time']
            }
            
          }
          # Snake Sites
          else if (JDR_det_hist[i-1, 'state'] %in% snake_sites){
            if (which(site_order_notrib_snake %in% JDR_det_hist[i, 'state']) >
                which(site_order_notrib_snake %in% JDR_det_hist[i - 1, 'state']) + 1 |
                which(site_order_notrib_snake %in% JDR_det_hist[i, 'state']) <
                which(site_order_notrib_snake %in% JDR_det_hist[i - 1, 'state']) - 1){
              # If we see it at two non-consecutive mainstem locations
              
              # If the next site skips sites, insert the missing sites
              # Get the index of the current site
              current_index <-
                which(site_order_notrib_snake %in% JDR_det_hist[i, 'state'])
              # Get the index of the previous site
              previous_index <-
                which(site_order_notrib_snake %in% JDR_det_hist[i - 1, 'state'])
              
              # sequence from current to previous index
              if (current_index < previous_index) {
                index_order <- seq(current_index, previous_index, by = 1)
                
              } else {
                index_order <- seq(previous_index, current_index, by = 1)
              }
              
              # Count the number of sites you need to add and loop through
              for (j in 1:(length(index_order) - 2)) {
                # Add a row for each of them
                
                # Insert a new row into stepwise states, with the implicit detection site
                # Tag code, state, and time (which is NA)
                missing_site <-
                  site_order_notrib_snake[index_order[1 + j]]
                implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                                    missing_site,
                                    NA)
                
                JDR_stepwise_states <-
                  insertRow(JDR_stepwise_states, implicit_state, i)
                
                # Insert a new row into original detection history, with
                # implicit detection site info
                implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        NA,
                                        "implicit site",
                                        missing_site)
                
                # Also, change the value in the original dataframe
                JDR_det_hist <-
                  insertRow(JDR_det_hist, implicit_detection, i)
                
                # Add 1 to the number of added rows
                added_rows <- added_rows + 1
              }
            }
            
            # Finally, if we're not missing any sites in between, just store it
            else {
              # Store the tag code
              JDR_stepwise_states[i, 'tag_code'] <-
                JDR_det_hist[i, 'tag_code']
              # Store the current state
              JDR_stepwise_states[i, 'state'] <- JDR_det_hist[i, 'state']
              # Store the time entering this state
              JDR_stepwise_states[i, 'date_time'] <-
                JDR_det_hist[i, 'end_time']
            }
            
          }
          
          
          
             
          
        }
        
        # Second: Upper Columbia sites
        else if (JDR_det_hist[i, 'state'] %in% upper_columbia_sites){
          if (which(site_order_notrib_columbia %in% JDR_det_hist[i, 'state']) >
              which(site_order_notrib_columbia %in% JDR_det_hist[i -
                                                                 1, 'state']) + 1 |
              which(site_order_notrib_columbia %in% JDR_det_hist[i, 'state']) <
              which(site_order_notrib_columbia %in% JDR_det_hist[i -
                                                                 1, 'state']) - 1){
            # If we see it at two non-consecutive mainstem locations
            
            # If the next site skips sites, insert the missing sites
            # Get the index of the current site
            current_index <-
              which(site_order_notrib_columbia %in% JDR_det_hist[i, 'state'])
            # Get the index of the previous site
            previous_index <-
              which(site_order_notrib_columbia %in% JDR_det_hist[i - 1, 'state'])
            
            # sequence from current to previous index
            if (current_index < previous_index) {
              index_order <- seq(current_index, previous_index, by = 1)
              
            } else {
              index_order <- seq(previous_index, current_index, by = 1)
            }
            
            # Count the number of sites you need to add and loop through
            for (j in 1:(length(index_order) - 2)) {
              # Add a row for each of them
              
              # Insert a new row into stepwise states, with the implicit detection site
              # Tag code, state, and time (which is NA)
              missing_site <-
                site_order_notrib_columbia[index_order[1 + j]]
              implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                                  missing_site,
                                  NA)
              
              JDR_stepwise_states <-
                insertRow(JDR_stepwise_states, implicit_state, i)
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      "implicit site",
                                      missing_site)
              
              # Also, change the value in the original dataframe
              JDR_det_hist <-
                insertRow(JDR_det_hist, implicit_detection, i)
              
              # Add 1 to the number of added rows
              added_rows <- added_rows + 1
            }
          } 
          
          # Finally, if we're not missing any sites in between, just store it
          else {
            # Store the tag code
            JDR_stepwise_states[i, 'tag_code'] <-
              JDR_det_hist[i, 'tag_code']
            # Store the current state
            JDR_stepwise_states[i, 'state'] <- JDR_det_hist[i, 'state']
            # Store the time entering this state
            JDR_stepwise_states[i, 'date_time'] <-
              JDR_det_hist[i, 'end_time']
          }
                   
          
        }
        
        # Third: Snake River sites
        else if (JDR_det_hist[i, 'state'] %in% snake_sites){
          if (which(site_order_notrib_snake %in% JDR_det_hist[i, 'state']) >
              which(site_order_notrib_snake %in% JDR_det_hist[i -
                                                                 1, 'state']) + 1 |
              JDR_det_hist[i, 'state'] %in% snake_sites &
              which(site_order_notrib_snake %in% JDR_det_hist[i, 'state']) <
              which(site_order_notrib_snake %in% JDR_det_hist[i -
                                                                 1, 'state']) - 1){
            # If we see it at two non-consecutive mainstem locations
            
            # If the next site skips sites, insert the missing sites
            # Get the index of the current site
            current_index <-
              which(site_order_notrib_snake %in% JDR_det_hist[i, 'state'])
            # Get the index of the previous site
            previous_index <-
              which(site_order_notrib_snake %in% JDR_det_hist[i - 1, 'state'])
            
            # sequence from current to previous index
            if (current_index < previous_index) {
              index_order <- seq(current_index, previous_index, by = 1)
              
            } else {
              index_order <- seq(previous_index, current_index, by = 1)
            }
            
            # Count the number of sites you need to add and loop through
            for (j in 1:(length(index_order) - 2)) {
              # Add a row for each of them
              
              # Insert a new row into stepwise states, with the implicit detection site
              # Tag code, state, and time (which is NA)
              missing_site <-
                site_order_notrib_snake[index_order[1 + j]]
              implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                                  missing_site,
                                  NA)
              
              JDR_stepwise_states <-
                insertRow(JDR_stepwise_states, implicit_state, i)
              
              # Insert a new row into original detection history, with
              # implicit detection site info
              implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      NA,
                                      "implicit site",
                                      missing_site)
              
              # Also, change the value in the original dataframe
              JDR_det_hist <-
                insertRow(JDR_det_hist, implicit_detection, i)
              
              # Add 1 to the number of added rows
              added_rows <- added_rows + 1
            }
          }
          
          # Finally, if we're not missing any sites in between, just store it
          else {
            # Store the tag code
            JDR_stepwise_states[i, 'tag_code'] <-
              JDR_det_hist[i, 'tag_code']
            # Store the current state
            JDR_stepwise_states[i, 'state'] <- JDR_det_hist[i, 'state']
            # Store the time entering this state
            JDR_stepwise_states[i, 'date_time'] <-
              JDR_det_hist[i, 'end_time']
          }
                 
          
        }
        
        # Finally, if we're not missing any sites in between, just store it
        else {
          # Store the tag code
          JDR_stepwise_states[i, 'tag_code'] <-
            JDR_det_hist[i, 'tag_code']
          # Store the current state
          JDR_stepwise_states[i, 'state'] <- JDR_det_hist[i, 'state']
          # Store the time entering this state
          JDR_stepwise_states[i, 'date_time'] <-
            JDR_det_hist[i, 'end_time']
        }
      }
      
    }
    
    
    
    
    
    
    
    
    ### REPEAT DETECTION AT THE SAME DAM
    # If it's next seen at the same dam, add a line to indicate it fell back
    # I think we need a statement for each dam, because each is unique
    else if (JDR_det_hist[i, 'site_class'] == JDR_det_hist[i - 1, 'site_class'] &
             JDR_det_hist[i, 'site_class'] == "BON (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                          "mainstem, mouth to BON",
                          NA)
      
      JDR_stepwise_states <-
        insertRow(JDR_stepwise_states, implicit_state, i)
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit site",
                              "mainstem, mouth to BON")
      
      # Also, change the value in the original dataframe
      JDR_det_hist <- insertRow(JDR_det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (JDR_det_hist[i, 'site_class'] == JDR_det_hist[i - 1, 'site_class'] &
             JDR_det_hist[i, 'site_class'] == "MCN (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                          "mainstem, BON to MCN",
                          NA)
      
      JDR_stepwise_states <-
        insertRow(JDR_stepwise_states, implicit_state, i)
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit site",
                              "mainstem, BON to MCN")
      
      # Also, change the value in the original dataframe
      JDR_det_hist <- insertRow(JDR_det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (JDR_det_hist[i, 'site_class'] == JDR_det_hist[i - 1, 'site_class'] &
             JDR_det_hist[i, 'site_class'] == "PRA (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                          "mainstem, MCN to ICH or PRA",
                          NA)
      
      JDR_stepwise_states <-
        insertRow(JDR_stepwise_states, implicit_state, i)
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(
        JDR_det_hist[i, 'tag_code'],
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        "implicit site",
        "mainstem, MCN to ICH or PRA"
      )
      
      # Also, change the value in the original dataframe
      JDR_det_hist <- insertRow(JDR_det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (JDR_det_hist[i, 'site_class'] == JDR_det_hist[i - 1, 'site_class'] &
             JDR_det_hist[i, 'site_class'] == "RIS (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                          "mainstem, PRA to RIS",
                          NA)
      
      JDR_stepwise_states <-
        insertRow(JDR_stepwise_states, implicit_state, i)
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit site",
                              "mainstem, PRA to RIS")
      
      # Also, change the value in the original dataframe
      JDR_det_hist <- insertRow(JDR_det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (JDR_det_hist[i, 'site_class'] == JDR_det_hist[i - 1, 'site_class'] &
             JDR_det_hist[i, 'site_class'] == "RRE (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                          "mainstem, RIS to RRE",
                          NA)
      
      JDR_stepwise_states <-
        insertRow(JDR_stepwise_states, implicit_state, i)
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit site",
                              "mainstem, RIS to RRE")
      
      # Also, change the value in the original dataframe
      JDR_det_hist <- insertRow(JDR_det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (JDR_det_hist[i, 'site_class'] == JDR_det_hist[i - 1, 'site_class'] &
             JDR_det_hist[i, 'site_class'] == "WEL (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                          "mainstem, RRE to WEL",
                          NA)
      
      JDR_stepwise_states <-
        insertRow(JDR_stepwise_states, implicit_state, i)
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit site",
                              "mainstem, RRE to WEL")
      
      # Also, change the value in the original dataframe
      JDR_det_hist <- insertRow(JDR_det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (JDR_det_hist[i, 'site_class'] == JDR_det_hist[i - 1, 'site_class'] &
             JDR_det_hist[i, 'site_class'] == "ICH (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                          "mainstem, MCN to ICH or PRA",
                          NA)
      
      JDR_stepwise_states <-
        insertRow(JDR_stepwise_states, implicit_state, i)
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(
        JDR_det_hist[i, 'tag_code'],
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        NA,
        "implicit site",
        "mainstem, MCN to ICH or PRA"
      )
      
      # Also, change the value in the original dataframe
      JDR_det_hist <- insertRow(JDR_det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    else if (JDR_det_hist[i, 'site_class'] == JDR_det_hist[i - 1, 'site_class'] &
             JDR_det_hist[i, 'site_class'] == "LGR (adult)") {
      # Insert a new row into stepwise states, with the implicit detection site
      # Tag code, state, and time (which is NA)
      implicit_state <- c(JDR_det_hist[i, 'tag_code'],
                          "mainstem, ICH to LGR",
                          NA)
      
      JDR_stepwise_states <-
        insertRow(JDR_stepwise_states, implicit_state, i)
      
      # Insert a new row into original detection history, with
      # implicit detection site info
      implicit_detection <- c(JDR_det_hist[i, 'tag_code'],
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              NA,
                              "implicit site",
                              "mainstem, ICH to LGR")
      
      # Also, change the value in the original dataframe
      JDR_det_hist <- insertRow(JDR_det_hist, implicit_detection, i)
      
      # Add 1 to the number of added rows
      added_rows <- added_rows + 1
    }
    
    # If it's in the same state, don't record it
    else {
      # Nothing!
    }
  }
  # If it's a different fish:
  else {
    # Start a new entry
    # Store the tag code
    JDR_stepwise_states[i, 'tag_code'] <- JDR_det_hist[i, 'tag_code']
    # Store the current state
    JDR_stepwise_states[i, 'state'] <- JDR_det_hist[i, 'state']
    # Store the time entering this state
    JDR_stepwise_states[i, 'date_time'] <-
      JDR_det_hist[i, 'end_time']
  }
}















# # Remove all NAs
JDR_stepwise_states_noNA <- JDR_stepwise_states[!is.na(JDR_stepwise_states$tag),]

# Compare them to original detection histories

subset(JDR_stepwise_states_noNA, tag_code == "3D9.1BF1989388")
subset(JDR_det_hist_original, tag_code == "3D9.1BF1989388")

subset(JDR_stepwise_states_noNA, tag_code == "3D9.1BF192A09F")
subset(JDR_det_hist_original, tag_code == "3D9.1BF192A09F")

# These look okay!


# 
# ##### Sort individuals into run years #####
# run_year <- c("05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15")
# run_year_start <- seq(ymd("2005-06-01"), ymsd("2014-06-01"), by = "years")
# run_year_end <- seq(ymd("2006-05-31"), ymd("2015-05-31"), by = "years")
# 
# run_year_df <- data.frame(run_year, run_year_start, run_year_end)
# 
# JDR_det_hist %>% 
#   group_by(tag_code) %>% 
#   filter(row_number() == 1) %>% 
#   mutate(run_year = subset(run_year_df, run_year_start < start_time & run_year_end > start_time)$run_year) %>% 
#   dplyr::select(tag_code, run_year) -> tag_codes_run_year
# 
# JDR_det_hist %>% 
#   left_join(., tag_codes_run_year, by = "tag_code") -> JDR_det_hist
# 
# ##### Break up detection history into individual transitions #####
# JDR_tag_codes <- unique(JDR_det_hist$tag_code)
# # Loop through each fish
# # for (i in 1:length(unique(JDR_tag_codes))){
# # Make a df to store values
# JDR_stepwise_detections <- data.frame(tag_code = character(),
#                                       site_1 = character(),
#                                       site_2 = character(),
#                                       date_time_1 = as.POSIXct(character()),
#                                       date_time_2 = as.POSIXct(character()))
# 
# for (i in 1:(nrow(JDR_det_hist)-1)){
#   
#   # If it's the same fish:
#   if (JDR_det_hist[i,'tag_code'] == JDR_det_hist[i+1,'tag_code']){
#     # Store the tag code
#     JDR_stepwise_detections[i,'tag_code'] <- JDR_det_hist[i,'tag_code']
#     # Store the current site, and the next site
#     JDR_stepwise_detections[i,'site_1'] <- JDR_det_hist[i,'event_site_name']
#     JDR_stepwise_detections[i,'site_2'] <- JDR_det_hist[i+1,'event_site_name']
#     # Store the time leaving this site, and the time entering the next site
#     JDR_stepwise_detections[i,'date_time_1'] <- JDR_det_hist[i,'end_time']
#     JDR_stepwise_detections[i,'date_time_2'] <- JDR_det_hist[i+1,'start_time']
#   }
# 
#   # If it's a different fish, end the entry (note that it's lost from the detection history)
#   # This may be because a fish spawned, or it could be undetermined loss at the end
#   else {
#     # End the previous entry
#     # Store the tag code
#     JDR_stepwise_detections[i,'tag_code'] <- JDR_det_hist[i,'tag_code']
#     # Store the current site, and the next site
#     JDR_stepwise_detections[i,'site_1'] <- JDR_det_hist[i,'event_site_name']
#     JDR_stepwise_detections[i,'site_2'] <- "lost"
#     # Store the time leaving this site, and the time entering the next site
#     JDR_stepwise_detections[i,'date_time_1'] <- JDR_det_hist[i,'end_time']
#     JDR_stepwise_detections[i,'date_time_2'] <- NA
#   }
# 
# }
# 
# 
# # Remove all of the NA rows
# # JDR_stepwise_detections <- JDR_stepwise_detections[!is.na(JDR_stepwise_detections$tag_code),]
# 
# 
# ##### Multistate model #####
# 
# # After detection in the adult fishways at each dam, there are four or five probabilities
# # that must sum to 1, for all of the options for a fish.
# # Let's use an individual detected at Bonneville as an example:
# # 1) Overshooting MCN
# # 2) Falling back over BON (can either be detected at a fallback array, or seen at the same dam)
# # 3) Straying to a tributary between BON and MCN
# # 4) Homing to the natal tributary (between BON and MCN)
# # 5) Undetermined loss (not seen again)
# 
# # We can use a tally-based approach to count these
# JDR_stepwise_detections %>% 
#   subset(site_1 == "Bonneville Adult Fishways (combined)") %>% 
#   left_join(dplyr::rename(JDR_site_classification, site_1 = event_site_name, 
#                           site_1_class = site_class), by = "site_1") %>% 
#   left_join(dplyr::rename(JDR_site_classification, site_2 = event_site_name, 
#                         site_2_class = site_class), by = "site_2") -> JDR_BON_site_1
# 
# table(JDR_BON_site_1$site_2_class)
# # Dealing with implicit site use - for example, if a fish is seen again 
# # at the same dam without any other sites in between, it means it 
# # fell back. OR if a fish is seen at the next dam but not at any of the in
# # river sites, it means it must have been in the river
# 
# # What about those that are next detected at an upstream site?
# JDR_stepwise_detections %>% 
#   left_join(dplyr::rename(JDR_site_classification, site_1 = event_site_name, 
#                           site_1_class = site_class), by = "site_1") %>% 
#   left_join(dplyr::rename(JDR_site_classification, site_2 = event_site_name, 
#                           site_2_class = site_class), by = "site_2") %>% 
#   subset(site_1_class == "BON_MCN_inriver") -> JDR_BON_MCN_inriver_site_1
# 
# table(JDR_BON_MCN_inriver_site_1$site_2_class)
# # Define parameters
# 
# ##### Update stepwise detection histories to indicate state transitions #####
# JDR_stepwise_detections %>% 
#   left_join(dplyr::rename(JDR_site_classification, site_1 = event_site_name, 
#                           site_1_class = site_class), by = "site_1") %>% 
#   left_join(dplyr::rename(JDR_site_classification, site_2 = event_site_name, 
#                           site_2_class = site_class), by = "site_2") -> JDR_stepwise_detections
# 
# 
# 
