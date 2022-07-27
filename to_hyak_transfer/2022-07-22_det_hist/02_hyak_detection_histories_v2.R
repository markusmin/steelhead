### 02 - Detection history generation

### Load libraries
library(tidyverse)
library(lubridate)
library(janitor)

# For testing - setwd
setwd("/Users/markusmin/Documents/CBR/steelhead/to_hyak_transfer/2022-07-22_det_hist/")

##### All tributaries #####


### Load complete detection history files
CTH_1 <- clean_names(read.csv("CTH_tag_codes_1.csv"))
CTH_2 <- clean_names(read.csv("CTH_tag_codes_2.csv"))
CTH_3 <- clean_names(read.csv("CTH_tag_codes_3.csv"))
CTH_4 <- clean_names(read.csv("CTH_tag_codes_4.csv"))
CTH_5 <- clean_names(read.csv("CTH_tag_codes_5.csv"))
CTH_6 <- clean_names(read.csv("CTH_tag_codes_6.csv"))
CTH_7 <- clean_names(read.csv("CTH_tag_codes_7.csv"))
CTH_8 <- clean_names(read.csv("CTH_tag_codes_8.csv"))
CTH_9 <- clean_names(read.csv("CTH_tag_codes_9.csv"))
CTH_10 <- clean_names(read.csv("CTH_tag_codes_10.csv"))
CTH_11 <- clean_names(read.csv("CTH_tag_codes_11.csv"))
CTH_12 <- clean_names(read.csv("CTH_tag_codes_12.csv"))
CTH_13 <- clean_names(read.csv("CTH_tag_codes_13.csv"))
CTH_14 <- clean_names(read.csv("CTH_tag_codes_14.csv"))


### Combine files, fix some column data types, rename some columns with long names
CTH_1 %>% 
  # dplyr::select(-event_site_subbasin_code) %>% 
  bind_rows(., dplyr::select(CTH_2, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_3, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_4, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_5, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_6, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_7, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_8, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_9, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_10, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_11, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_12, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_13, -event_site_subbasin_code)) %>% 
  bind_rows(., dplyr::select(CTH_14, -event_site_subbasin_code)) %>% 
  mutate(event_date_time_value = mdy_hms(event_date_time_value)) %>% 
  dplyr::rename(ant_config = antenna_group_configuration_value) -> CTH_complete

# FOR TESTING - read in just the first file
CTH_1 <- clean_names(read.csv(here::here("PTAGIS_queries", "complete_tag_histories", "antenna_info", "2022-07-25-CTH1.csv")))
CTH_1 %>% 
  mutate(event_date_time_value = mdy_hms(event_date_time_value)) %>% 
  dplyr::rename(ant_config = antenna_group_configuration_value) -> CTH_complete


# Remove the Klickitat and Wind Rivers
origin_table <- read.csv("natal_origin_table.csv")
# Read in the metadata
tag_code_metadata <- read.csv("tag_code_metadata.csv")
tag_code_metadata %>% 
  left_join(., origin_table, by = "release_site_name") -> tag_code_origin_metadata
KLIC_WIND_tag_codes <- subset(tag_code_origin_metadata, natal_origin %in% c("Klickitat_River", "Wind_River"))$tag_code
# Subset out those
CTH_complete %>% 
  subset(., !(tag_code %in% KLIC_WIND_tag_codes)) -> CTH_complete


# Get the date of arrival at BON - take min in case it fell back over Bonneville
CTH_complete %>% 
  group_by(tag_code) %>% 
  subset(event_site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                                "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF", 
                                "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) %>% 
  filter(event_date_time_value == min(event_date_time_value)) %>% 
  dplyr::select(tag_code, event_date_time_value) %>% 
  dplyr::rename(BON_arrival = event_date_time_value) -> BON_arrival

# Add run year info
# Add run year info
run_year <- c("05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21", "21/22")
run_year_start <- seq(ymd("2005-06-01"), ymd("2021-06-01"), by = "years")
run_year_end <- seq(ymd("2006-05-31"), ymd("2022-05-31"), by = "years")

run_year_df <- data.frame(run_year, run_year_start, run_year_end)

BON_arrival %>% 
  mutate(dummy =TRUE) %>% 
  left_join(run_year_df %>% mutate(dummy=TRUE), by = "dummy") %>% 
  filter(run_year_start <= BON_arrival, run_year_end >= BON_arrival) %>% 
  select(-c(dummy, run_year_start, run_year_end)) -> BON_arrival

write.csv(BON_arrival, "complete_BON_arrival.csv", row.names = FALSE)


# Subset only the adult migration history
CTH_complete %>% 
  left_join(., BON_arrival, by = "tag_code") %>% 
  # Keep only observations, no mark or recapture events (there is some weird stuff going on with those)
  subset(., event_type_name == "Observation") %>% 
  group_by(tag_code) %>% 
  filter(event_date_time_value >= BON_arrival) %>% 
  dplyr::select(-BON_arrival) %>% 
  # Make sure that they are in chronological order
  arrange(tag_code, event_date_time_value) -> CTH_adult

CTH_adult %>% 
  group_by(tag_code) %>% 
  filter(any(event_site_name %in% c("BHL - Adult Fishway at BONH", "BON - Bonneville Dam Complex"))) -> BON_non_ladders

# So it would appear that there wasn't a single detection at "BON - Bonneville Dam Complex" in our detection histories - unsure of where this came from, then

# EDIT CTH_adult
# We have to edit this so that the code recognizes that different detectors
# are at the same site

# HOWEVER: We want to keep the original event site names as well so that we don't lose that detail in later files
CTH_adult %>% 
  mutate(event_site_name_original = event_site_name) -> CTH_adult

# Let's check where the BO4 fish were previously
CTH_adult %>% 
  mutate(last_event_site_name = lag(event_site_name)) %>% 
  subset(event_site_name == "BO4 - Bonneville WA Ladder Slots") -> BO4_fish

BO4_fish %>% 
  subset(last_event_site_name != "BO4 - Bonneville WA Ladder Slots") -> BO4_fish

# table(BO4_fish$last_event_site_name)

# Look at fish that weren't previously at BO2 or BO3
subset(BO4_fish, !(last_event_site_name %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF")))$tag_code -> nondirect_BO4s
subset(BO4_fish, last_event_site_name == "BO1 - Bonneville Bradford Is. Ladder")$tag_code -> BO1_BO4s

subset(BO4_fish, last_event_site_name == "DRM - Deschutes River mouth")$tag_code -> DRM_BO4s

# Let's look at fish that were released into the adult fish ladder
subset(CTH_adult, event_site_name == "LGRLDR - LGR - Release into the Adult Fish Ladder")$tag_code -> LGR_ladder_held_fish

# Okay, some edits:
# BHL is the Bonneville Hatchery ladder, which was only active 2012-2014, and is actually before the dam itself (and is not part of the adult fishways)
# Combine the arrays for the dams that don't separate states in our model (because they weren't operational for the full study period)
CTH_adult %>%
#   mutate(event_site_name = ifelse(event_site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder",
#                                                          "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
#                                                          "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility",
#                                                          "BON - Bonneville Dam Complex"),
#                                   "Bonneville Adult Fishways (combined)", event_site_name)) %>%
#   mutate(event_site_name = ifelse(event_site_name %in% c("MC1 - McNary Oregon Shore Ladder", "MC2 - McNary Washington Shore Ladder"),
#                                   "McNary Adult Fishways (combined)", event_site_name)) %>% 
  mutate(event_site_name = ifelse(event_site_name %in% c("TD1 - The Dalles East Fish Ladder", "TD2 - The Dalles North Fish Ladder"),
                                  "The Dalles Adult Fishways (combined)", event_site_name)) %>%
  mutate(event_site_name = ifelse(event_site_name %in% c("JO1 - John Day South Fish Ladder", "JO2 - John Day North Fish Ladder"),
                                  "John Day Dam Adult Fishways (combined)", event_site_name)) -> CTH_adult
#   mutate(event_site_name = ifelse(event_site_name %in% c("ICH - Ice Harbor Dam (Combined)", "IHR - Ice Harbor Dam"),
#                                   "Ice Harbor Adult Fishways (combined)", event_site_name)) %>% 
#   mutate(event_site_name = ifelse(event_site_name %in% c("LGRLDR - LGR - Release into the Adult Fish Ladder", "GRA - Lower Granite Dam Adult",
#                                                          "LGR - Lower Granite Dam"),
#                                   "Lower Granite Dam Adult Fishways (combined)", event_site_name)) %>% 
#   mutate(event_site_name = ifelse(event_site_name %in% c("RIA - Rock Island Adult", "RIS - Rock Island Dam"),
#                                   "Rock Island Adult Fishways (combined)", event_site_name)) %>% 
#   mutate(event_site_name = ifelse(event_site_name %in% c("WEA - Wells Dam, DCPUD Adult Ladders", 
#                                                          "WELLD2 - WEL - Release into the West Adult Fish Ladder",
#                                                          "WELLD1 - WEL - Release into the East Adult Fish Ladder",
#                                                          "WEL - Wells Dam"), "Wells Dam Adult Fishways (combined)", event_site_name)) %>% 
#   mutate(event_site_name = ifelse(event_site_name %in% c("PRA - Priest Rapids Adult", "PRD - Priest Rapids Dam",
#                                                          "PRDLD1 - PRD - Release into the Left Bank (facing downstream) Adult Fish Ladder"), 
#                                   "Priest Rapids Adult Fishways (combined)", event_site_name)) -> CTH_adult


# BON array names
BON_arrays <- c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility",
                "BON - Bonneville Dam Complex")

# McNary array names
MCN_arrays <- c("MC1 - McNary Oregon Shore Ladder", "MC2 - McNary Washington Shore Ladder")

# ICH array names
ICH_arrays <- c("ICH - Ice Harbor Dam (Combined)", "IHR - Ice Harbor Dam")

# LGR array names
LGR_arrays <- c("LGRLDR - LGR - Release into the Adult Fish Ladder", "GRA - Lower Granite Dam Adult", "LGR - Lower Granite Dam")

# RIA array names
RIS_arrays <- c("RIA - Rock Island Adult", "RIS - Rock Island Dam")

# WEL array names
WEL_arrays <- c("WEA - Wells Dam, DCPUD Adult Ladders", "WELLD2 - WEL - Release into the West Adult Fish Ladder",
                "WELLD1 - WEL - Release into the East Adult Fish Ladder", "WEL - Wells Dam")

# PRA array names
PRA_arrays <- c("PRA - Priest Rapids Adult", "PRD - Priest Rapids Dam", 
                "PRDLD1 - PRD - Release into the Left Bank (facing downstream) Adult Fish Ladder")


##### Store the antenna configurations #####






##### Convert complete tag histories into history of detection events #####


# Get the first and last time of an event, store those


# Get a list of the unique tag IDs
unique_tag_IDs <- unique(CTH_adult$tag_code) 

# Create a subset of data to test code - every 20th individual
# unique_tag_IDs <- unique_tag_IDs[seq(1, length(unique_tag_IDs), 10)]

# Create a new dataframe that will store our new detection history
# det_hist <- data.frame(tag_code = character(), event_site_name = character(), 
#                            start_time = as.POSIXct(character()), end_time = as.POSIXct(character()), 
#                            event_site_basin_name = character(), event_site_subbasin_name = character(),
#                            event_site_latitude = numeric(), event_site_longitude = numeric())



# Make a function to store all of the fields we want
# 
# store_det_info <- function(counter, i, j){
#   
#   # store the tag code
#   ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
#   
#   # store the event type name
#   ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
#   
#   # store the antenna info
#   # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
#   # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
#   ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
#   
#   # store the location fields
#   ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
#   # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
#   # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
#   # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
#   # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
#   
#   # store the end time
#   ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
#   # Store the end antenna
#   ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
#   
# }
# 
# # second version
# store_det_info <- function(counter, i, j){
#   
#   # store all of these values in a vector
#   value_vec <- vector(length = 9)
#   
#   # store the tag code
#   # ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
#   value_vec[1] <- unique_tag_IDs[i]
#   
#   # store the event type name
#   # ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
#   value_vec[2] <- as.character(tag_hist[j,'event_type_name'])
#   
#   # store the event site name
#   # ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
#   value_vec[3] <- as.character(tag_hist[j,'event_site_name'])
#   
#   # no start time in this function
#   
#   # store the end time
#   # ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
#   # value_vec[5] <- ymd_hms(as.character(as_datetime(as.numeric(tag_hist[j,'event_date_time_value']))))
#   end_time <- ymd_hms(as.character(as_datetime(as.numeric(tag_hist[j,'event_date_time_value']))))
# 
#   # no start antenna group
#   
#   # Store the end antenna group
#   # ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
#   value_vec[7] <- as.character(tag_hist[j,'antenna_group_name'])
#   
#   
#   # store the antenna configuration
#   value_vec[8] <- as.character(tag_hist[j,'ant_config'])
#   
#   # return the values in a vector
#   return(list(value_vec, end_time))
#   
# }


det_hist <- data.frame(tag_code = character(), event_type_name = character(), event_site_name = character(), 
                       start_time = as.POSIXct(character()), end_time = as.POSIXct(character()), 
                       # antenna_group_name = character(), antenna_id = character(), ant_config = numeric())
                       start_antenna = character(), end_antenna = character(), ant_config = numeric(), aborted = character())


# Loop through the unique tags
# for (i in 1:length(unique_tag_IDs)){
for (i in 1:50){
  # Get the start time
  if (i == 1){
    start_time <- Sys.time() 
  }
  
  # Make a new dataframe to store the history for each fish
  ind_det_hist <- data.frame(tag_code = character(), event_type_name = character(), event_site_name = character(), 
                             start_time = as.POSIXct(character()), end_time = as.POSIXct(character()), 
                             # antenna_group_name = character(), antenna_id = character(), ant_config = numeric())
                             start_ant_group = character(), end_ant_group = character(),
                             ant_config = numeric(), aborted = character())
  
  # Make a new dataframe - we don't know how big it needs to be, but we can make it 
  # # large and then cut it down at the end
  # Didn't really help with speed
  # ind_det_hist <- data.frame(tag_code = rep(NA, 40), event_site_name = rep(NA, 40), 
  #                            start_time = as.POSIXct(rep(NA, 40)), end_time = as.POSIXct(rep(NA, 40)), 
  #                            event_site_basin_name = rep(NA, 40), event_site_subbasin_name = rep(NA, 40),
  #                            event_site_latitude = as.numeric(rep(NA, 40)), event_site_longitude = as.numeric(rep(NA, 40)))
  
  # subset the complete dataset to only this fish
  tag_hist <- subset(CTH_adult, tag_code == unique_tag_IDs[i])
  
  
# Loop through the rows of the tag history
for (j in 1:nrow(tag_hist)){
    
    
    ### For the first entry, just store these values plus the time as the start time
    if (j == 1){
      # store the tag code
      ind_det_hist[1,'tag_code'] <- unique_tag_IDs[i]
        
      # store the event type name
      ind_det_hist[1,'event_type_name'] <- tag_hist[j,'event_type_name']
    
      
      # store the antenna info
      # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
      # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
      ind_det_hist[1,'ant_config'] <- tag_hist[j,'ant_config']
      
      # If it's BO2-BO3-BO4, store the site name that way
      if (tag_hist[j,'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                               "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")){
        ind_det_hist[1,'event_site_name'] <- 'BO2-BO3-BO4' 
        
      }
      # otherwise, store the site name like normal
      else {
        ind_det_hist[1,'event_site_name'] <- tag_hist[j,'event_site_name']
      }  
      
      
      
      # Store the start time
      ind_det_hist[1, 'start_time'] <- tag_hist[j,'event_date_time_value']
      # Store the start antenna group
      ind_det_hist[1,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
      
      # START THE COUNTER
      counter <- 1
      
      # SPECIAL CASE: If there is only a single detection in the entire history, then end it
      if(nrow(tag_hist) == 1){
        ind_det_hist[1,'end_time'] <- tag_hist[j,'event_date_time_value']
        ind_det_hist[1,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
      }
      
      # SPECIAL CASE: If there is only one detection at the first site, store the 
      # end time as well

      else if (tag_hist[j+1, 'event_site_name'] != tag_hist[j, 'event_site_name'] &
               # Need to make sure they're not in the BO2-BO3-BO4 complex
               !(tag_hist[j, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                        "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility"))){
        
        # Store the end time
        ind_det_hist[counter, 'end_time'] <- tag_hist[j,'event_date_time_value']
        # Store the end antenna group
        ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
        
        # UPDATE THE COUNTER
        # every time we store an end time, we update the counter. This allows
        # us to move through the detection history df
        counter <- counter + 1
      }
      

    }
    
    ### If it's the last entry, store those values
    
    else if (j == nrow(tag_hist)){
      # Store the end time
      ind_det_hist[counter, 'end_time'] <- tag_hist[j,'event_date_time_value']
      # don't need to update the counter here because it's the last observation for the fish
      
      # SPECIAL CASE: If it's the last entry AND the only detection at a site
      if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
          tag_hist[j, 'event_date_time_value'] -
          tag_hist[j-1, 'event_date_time_value'] >= hours(x = 6)){
        
        # store the tag code
        ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
        
        # store the event type name
        ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
        
        # store the antenna info
        # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
        # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
        ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
        
        # store the location field
        ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
        
        # Store the start time
        ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
        # Store the start antenna
        ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
        
      }
    }
    
    ### If there was only one detection at a site, then store the start and end times at the same time
    
    # See if this event site name is different than the one ahead and behind
    else if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] &
             tag_hist[j+1, 'event_site_name'] != tag_hist[j, 'event_site_name'] &
             # Need to make sure they're not in the BO2-BO3-BO4 complex
             !(tag_hist[j, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                        "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility"))){
      
      # store the tag code
      ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
      
      # store the event type name
      ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
      
      # store the antenna info
      # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
      # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
      ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
      
      # store the location fields
      ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
      
      # Store the start time
      ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
      # Store the start antenna
      ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
      
      # Store the end time
      ind_det_hist[counter, 'end_time'] <- tag_hist[j,'event_date_time_value']
      # Store the end antenna
      ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
      # UPDATE THE COUNTER
      # every time we store an end time, we update the counter. This allows
      # us to move through the detection history df
      counter <- counter + 1
      
    }
    
    
    
    
    # Look at the adult ladders specifically
    
    ##### BONNEVILLE #####
    
    else if (tag_hist[j, 'event_site_name'] %in% BON_arrays){
      
      # BON ROUTE 1 - BO1
      # This one only has a single array
      if(tag_hist[j, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder"){
        
        # If it's the first detection in this route, store the start time and antenna
        if (tag_hist[j-1, 'event_site_name'] != "BO1 - Bonneville Bradford Is. Ladder") {
          # Store the event site name
          ind_det_hist[counter, 'event_site_name'] <- tag_hist[j,'event_site_name']
          # Store the start time
          ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
          # Store the start antenna
          ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
        }
        
        
        # Different if statements for varying antenna configurations
        
        # BO1 100
        else if(tag_hist[j, 'ant_config'] == 100){
          
          # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
          if (tag_hist[j+1, 'event_site_name'] != "BO1 - Bonneville Bradford Is. Ladder" | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's within 48 hours
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            # store the antenna info
            # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']

            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
          
          
          }
          
          # if it's within 48 hours at BO1, then don't record it
          else {
            
          }
        
        }
        
        # BO1 110 & 120
        # These are the same for our purposes, since the only difference is the installation of lamprey monitors for BO1 120
        if(tag_hist[j, 'ant_config'] %in% c(110, 120)){
          
          
          # If the next detection is at a different site, or is more than 48 hours later, store the current time as the end time
          if (tag_hist[j+1, 'event_site_name'] != "BO1 - Bonneville Bradford Is. Ladder" | # If it's at a different site
              tag_hist[j+1, 'event_date_time_value'] - 
              tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's within 48 hours
            
            # store the info
            # store the tag code
            ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
            
            # store the event type name
            ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
            
            # store the antenna info
            # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
            # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
            ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
            
            # store the location fields
            ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
            # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
            # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
            # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
            # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
            
            # store the end time
            ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
            # Store the end antenna
            ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
            
            # UPDATE THE COUNTER
            # every time we store an end time, we update the counter. This allows
            # us to move through the detection history df
            counter <- counter + 1
            
            # However, if the last antenna seen is not an exit coil, note that it was an aborted attempt
            if (tag_hist[j, 'antenna_id'] %in% BO1_110_lower){
              
              ind_det_hist[counter,'aborted'] <- "aborted"
            
            }
            
          }
          
          # if it's within 48 hours at BO1, then don't record it
          
          else {
            
          }
          
        }
        
      } 
      # ROUTE 2 - BO2-BO3-BO4 complex
        # For this route, we need to make sure that the fish came out through BO4, otherwise it's an aborted attempt.
        # We will use the same code as for BO1, but use BO2 and BO3 as the "lower" portions of the ladder
      else if(tag_hist[j, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                    "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")){
        
        # if it's the first detection in this complex, store the start time 
        if (!(tag_hist[j-1, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                       "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility"))) {
          # Store the event site name
          ind_det_hist[counter, 'event_site_name'] <- 'BO2-BO3-BO4'
          # Store the start time
          ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
          # Store the start antenna
          ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
        }
        
        # If the next detection is at a site not in BO2, BO3, or BO4,
        # or is more than 48 hours later, store the current time as the end time
        else if (!(tag_hist[j+1, 'event_site_name']  %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                       "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) | # If it's at a different site
            tag_hist[j+1, 'event_date_time_value'] - 
            tag_hist[j, 'event_date_time_value'] >= hours(x = 48)){ # or it's over 48 hours later
        
          # store the info
          # store the tag code
          ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
          
          # store the event type name
          ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
          
          # store the antenna info
          # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
          # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
          ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
          
          # store the location fields
          # ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
          # ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
          # ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
          # ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
          # ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
          
          # store the end time
          ind_det_hist[counter,'end_time'] <- tag_hist[j,'event_date_time_value']  
          # Store the end antenna
          ind_det_hist[counter,'end_ant_group'] <- tag_hist[j,'antenna_group_name']
          
          # Save event site name
          ind_det_hist[counter,'event_site_name'] <- 'BO2-BO3-BO4'
          


          # However, if the last antenna seen is not in BO4, then note it was an aborted attempt
          if (tag_hist[j, 'event_site_name'] %in% c("BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF",
                                                    "BONAFF - BON - Adult Fish Facility")){
            
            ind_det_hist[counter,'aborted'] <- "aborted"
            
          }
          
          # UPDATE THE COUNTER
          # every time we store an end time, we update the counter. This allows
          # us to move through the detection history df
          counter <- counter + 1
        
        
        }
        
        
      
      }
      
    }
    
    
    ##### McNary #####
    
    
    
    
    
    
    # For every other entry that is not at a dam, look at the previous entry to see if it
    # was the same site < 6 hours ago
    else {
      
      # If the next entry isn't the same site OR is >=6 hours ahead, store
      # the current time as the end time
      if (tag_hist[j+1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
          tag_hist[j+1, 'event_date_time_value'] -
          tag_hist[j, 'event_date_time_value'] >= hours(x = 6)){
        
        # Store the end time
        ind_det_hist[counter, 'end_time'] <- tag_hist[j,'event_date_time_value']
        
        # SPECIAL CASE: If there is only one detection at a site, store the 
        # start time as well
        # Here we will look at both the previous and next entry
        # If the previous detection is at a different site OR is more than 6 hours ago,
        # store the start time
        if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
            tag_hist[j, 'event_date_time_value'] -
            tag_hist[j-1, 'event_date_time_value'] >= hours(x = 6)){
          
          # store the tag code
          ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
          
          # store the event type name
          ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
          
          # store the antenna info
          # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
          # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
          ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
          
          # store the location fields
          ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
          
          # Store the start time
          ind_det_hist[counter, 'start_time'] <- tag_hist[j,'event_date_time_value']
          # Store the start antenna
          ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
          
        }
        
        # UPDATE THE COUNTER
        # every time we store an end time, we update the counter. This allows
        # us to move through the detection history df
        counter <- counter + 1
        
      }
      
      # If the previous site was a different site OR was >6 hours ago, 
      # start a new entry and start time
      else if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
               tag_hist[j, 'event_date_time_value'] -
               tag_hist[j-1, 'event_date_time_value'] >= hours(x = 6)){
      
        
        # store the tag code
        ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
        
        # store the event type name
        ind_det_hist[counter,'event_type_name'] <- tag_hist[j,'event_type_name']
        
        # store the antenna info
        # ind_det_hist[counter,'antenna_id'] <- tag_hist[j,'antenna_id']
        # ind_det_hist[counter,'antenna_group_name'] <- tag_hist[j,'antenna_group_name']
        ind_det_hist[counter,'ant_config'] <- tag_hist[j,'ant_config']
        
        # store the location fields
        ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
        
        # store the start time
        ind_det_hist[counter,'start_time'] <- tag_hist[j,'event_date_time_value']
        # Store the start antenna
        ind_det_hist[counter,'start_ant_group'] <- tag_hist[j,'antenna_group_name']
        
      }
      
      # If the previous entry was the same site <=6 hours ago AND
      # The next entry is the same site <=6 hours ago, skip it
      else {
        
      }
      
    }
  }
  
  # Cut out extra rows
  ind_det_hist <- ind_det_hist[!(is.na(ind_det_hist$tag_code)),]
  
  # Append the individual tag history to the complete tag history
  det_hist %>% 
    bind_rows(., ind_det_hist) -> det_hist
  
  # Print the run time
  if (i == length(unique_tag_IDs)){
    end_time <- Sys.time()
    
    print(paste0("Total tags: ", length(unique_tag_IDs)))
    print(paste0("Total time: ", end_time - start_time))
  }
  
  # Make sure that it's running
  print(paste0("Tag ", i))
}

# Export this detection history
# write.csv(det_hist, "complete_det_hist.csv")
write.csv(det_hist, "test_hist.csv")


##### 

# det_hist <- read.csv(here("model_files", "complete_det_hist.csv"), row.names = 1)


# Inspect sites - for mapping

det_hist %>% 
  # Make a correction for JDA - shouldn't be necessary later once the script is re-run
  mutate(event_site_name = ifelse(event_site_name %in% c("JO1 - John Day South Fish Ladder", "JO2 - John Day North Fish Ladder"),
                                  "John Day Dam Adult Fishways (combined)", event_site_name)) %>% 
  # mutate(event_site_latitude = ifelse(event_site_name == "John Day Dam Adult Fishways (combined)", 45.71866, event_site_latitude)) %>% 
  # mutate(event_site_longitude = ifelse(event_site_name == "John Day Dam Adult Fishways (combined)", -120.6978, event_site_longitude)) %>% 
  group_by(event_site_name) %>% 
  summarise(n()) %>% 
  dplyr::rename(count = `n()`)-> event_det_counts

# Get the metadata
det_hist %>% 
  dplyr::select(-c(tag_code, start_time, end_time)) %>% 
  distinct(event_site_name, .keep_all = TRUE) -> event_site_metadata

write.csv(event_site_metadata, "complete_event_site_metadata.csv")

event_site_metadata %>% 
  # Make a correction for JDA - shouldn't be necessary later once the script is re-run
  mutate(event_site_name = ifelse(event_site_name %in% c("JO1 - John Day South Fish Ladder", "JO2 - John Day North Fish Ladder"),
                                  "John Day Dam Adult Fishways (combined)", event_site_name)) -> event_site_metadata

# # Get the JDR event det counts - this somehow has some sites that aren't in the complete data (because the complete data isn't truly complete yet)
# JDR_event_det_counts <- read.csv(here::here("model_files", "JDR_event_det_counts.csv"), row.names = 1)
# 
# JDR_event_det_counts %>% 
#   rownames_to_column("event_site_name") %>% 
#   dplyr::select(event_site_name, count) -> JDR_event_det_counts
# 
# # Add these together
# event_det_counts %>% 
#   bind_rows(., JDR_event_det_counts) %>% 
#   subset(., !duplicated(event_site_name)) -> event_det_counts



event_det_counts %>% 
  left_join(., event_site_metadata, by = "event_site_name") -> event_det_counts

# Get another df for dams
# Note: this is outdated, but also not important to the output
event_det_counts %>% 
  mutate(dam = ifelse(event_site_name %in% c("Bonneville Adult Fishways (combined)",
                                             "McNary Adult Fishways (combined)",
                                             "Priest Rapids Adult Fishways (combined)",
                                             "Rock Island Adult Fishways (combined)", 
                                             "RRF - Rocky Reach Fishway", 
                                             "Wells Dam Adult Fishways (combined)",
                                             "Ice Harbor Adult Fishways (combined)",  
                                             "Lower Granite Dam Adult Fishways (combined)",
                                             "John Day Dam Adult Fishways (combined)",
                                             # Dams without consistent PIT tag detectors
                                             # Missing John Day for this dataset (installed 2017)
                                             "The Dalles Adult Fishways (combined)",
                                             "LMA - Lower Monumental Adult Ladders",
                                             "GOA - Little Goose Fish Ladder"), "dam",
                      "in stream")) %>% 
  # Get a field for dam abbreviations
  mutate(dam_abbr = ifelse(event_site_name == "Bonneville Adult Fishways (combined)", "BON",
                           ifelse(event_site_name == "McNary Adult Fishways (combined)", "MCN",
                                  ifelse(event_site_name == "Priest Rapids Adult Fishways (combined)", "PRA",
                                         ifelse(event_site_name == "Rock Island Adult Fishways (combined)", "RIS",
                                                ifelse(event_site_name == "RRF - Rocky Reach Fishway", "RRE",
                                                       ifelse(event_site_name == "Wells Dam Adult Fishways (combined)", "WEL",
                                                              ifelse(event_site_name == "Ice Harbor Adult Fishways (combined)", "ICH",
                                                                     ifelse(event_site_name == "Lower Granite Dam Adult Fishways (combined)", "LGR",
                                                                            # Dams without consistent PIT tag detectors
                                                                            ifelse(event_site_name == "John Day Dam Adult Fishways (combined)", "JDA",
                                                                            ifelse(event_site_name == "The Dalles Adult Fishways (combined)", "TDA",
                                                                                   ifelse(event_site_name == "LMA - Lower Monumental Adult Ladders", "LMO",
                                                                                          ifelse(event_site_name == "GOA - Little Goose Fish Ladder", "LGO", NA))))))))))))) -> event_det_counts

write.csv(event_det_counts,  "complete_event_det_counts.csv", row.names = FALSE)






