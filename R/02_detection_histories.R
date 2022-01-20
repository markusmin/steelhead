### 02 - Detection history generation

### Load libraries
library(here)
library(tidyverse)
library(lubridate)
library(janitor)

##### John Day River - wild #####

### Load complete detection history files
JDR_CTH_1 <- clean_names(read.csv(here::here("PTAGIS_queries", "complete_tag_histories", "john_day_river", "2022-01-14-john_day_river_CTH_2005_2015_1.csv")))
JDR_CTH_2 <- clean_names(read.csv(here::here("PTAGIS_queries", "complete_tag_histories", "john_day_river", "2022-01-14-john_day_river_CTH_2005_2015_2.csv")))

### Combine files, fix some column data types
JDR_CTH_1 %>% 
  bind_rows(., JDR_CTH_2) %>% 
  mutate(event_date_time_value = mdy_hms(event_date_time_value))-> JDR_CTH

# Get the date of arrival at BON - take min in case it fell back over Bonneville
JDR_CTH %>% 
  group_by(tag_code) %>% 
  subset(event_site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                                              "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF", 
                                              "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) %>% 
  filter(event_date_time_value == min(event_date_time_value)) %>% 
  dplyr::select(tag_code, event_date_time_value) %>% 
  dplyr::rename(BON_arrival = event_date_time_value) -> JDR_BON_arrival

# Subset only the adult migration history
JDR_CTH %>% 
  left_join(., JDR_BON_arrival, by = "tag_code") %>% 
  group_by(tag_code) %>% 
  filter(event_date_time_value >= BON_arrival) %>% 
  dplyr::select(-BON_arrival) %>% 
  # Make sure that they are in chronological order
  arrange(tag_code, event_date_time_value) -> JDR_CTH_adult

# EDIT JDR_CTH_adult
# We have to edit this so that the code recognizes that different detectors
# are at the same site
JDR_CTH_adult %>% 
  mutate(event_site_name = ifelse(event_site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                                                        "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF", 
                                                        "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility"),
                                  "Bonneville Adult Fishways (combined)", event_site_name)) %>% 
  mutate(event_site_name = ifelse(event_site_name %in% c("MC1 - McNary Oregon Shore Ladder", "MC2 - McNary Washington Shore Ladder"),
                                  "McNary Adult Fishways (combined)", event_site_name))-> JDR_CTH_adult



### Convert complete tag histories into history of detection events

# Use cutoff of six hours between events to describe different events

# Get the first and last time of an event, store those

# Get a list of the unique tag IDs
unique_tag_IDs <- unique(JDR_CTH_adult$tag_code) 

# Create a new dataframe that will store our new detection history
JDR_det_hist <- data.frame(tag_code = character(), event_site_name = character(), 
                           start_time = as.POSIXct(character()), end_time = as.POSIXct(character()), 
                           event_site_basin_name = character(), event_site_subbasin_name = character(),
                           event_site_latitude = numeric(), event_site_longitude = numeric())

# Loop through the unique tags
for (i in 1:length(unique_tag_IDs)){
# for (i in 1:10){
  # Get the start time
  if (i == 1){
    start_time <- Sys.time() 
  }
  
  # Make a new dataframe to store the history for each fish
  ind_det_hist <- data.frame(tag_code = character(), event_site_name = character(), 
                             start_time = as.POSIXct(character()), end_time = as.POSIXct(character()), 
                             event_site_basin_name = character(), event_site_subbasin_name = character(),
                             event_site_latitude = numeric(), event_site_longitude = numeric())
  
  # subset the complete dataset to only this fish
  tag_hist <- subset(JDR_CTH_adult, tag_code == unique_tag_IDs[i])
  
  
  # Loop through the rows of the tag history
  for (j in 1:nrow(tag_hist)){
    
    
    # For the first entry, just store these values
    if (j == 1){
      # store the tag code
      ind_det_hist[1,'tag_code'] <- unique_tag_IDs[i]
      
      # store the location fields
      ind_det_hist[1,'event_site_name'] <- tag_hist[j,'event_site_name']
      ind_det_hist[1,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
      ind_det_hist[1,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
      ind_det_hist[1,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
      ind_det_hist[1,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
      
      # store the start time
      ind_det_hist[1,'start_time'] <- tag_hist[[j,'event_date_time_value']]
      
      # START THE COUNTER
      counter <- 1
    }
    
    # If it's the last entry, store those values
    
    else if (j == nrow(tag_hist)){
      # Store the end time
      ind_det_hist[counter, 'end_time'] <- tag_hist[[j,'event_date_time_value']]
      
      # SPECIAL CASE: If it's the last entry AND the only detection at a site
      if (tag_hist[j-1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
          tag_hist[j, 'event_date_time_value'] -
          tag_hist[j-1, 'event_date_time_value'] >= hours(x = 6)){
        
        # store the tag code
        ind_det_hist[counter,'tag_code'] <- unique_tag_IDs[i]
        
        # store the location fields
        ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
        ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
        ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
        ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
        ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
        
        # Store the start time
        ind_det_hist[counter, 'start_time'] <- tag_hist[[j,'event_date_time_value']]
        
      }
    }
    
    # For every other entry, look at the previous entry to see if it
    # was the same site < 6 hours ago
    else {
      
      # If the next entry isn't the same site OR is >=6 hours ahead, store it as
      # the end time
      if (tag_hist[j+1, 'event_site_name'] != tag_hist[j, 'event_site_name'] |
               tag_hist[j+1, 'event_date_time_value'] -
               tag_hist[j, 'event_date_time_value'] >= hours(x = 6)){
        
        # Store the end time
        ind_det_hist[counter, 'end_time'] <- tag_hist[[j,'event_date_time_value']]
        
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
          
          # store the location fields
          ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
          ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
          ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
          ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
          ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
          
          # Store the start time
          ind_det_hist[counter, 'start_time'] <- tag_hist[[j,'event_date_time_value']]
          
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
        
        # store the location fields
        ind_det_hist[counter,'event_site_name'] <- tag_hist[j,'event_site_name']
        ind_det_hist[counter,'event_site_basin_name'] <- tag_hist[j,'event_site_basin_name']
        ind_det_hist[counter,'event_site_subbasin_name'] <- tag_hist[j,'event_site_subbasin_name']
        ind_det_hist[counter,'event_site_latitude'] <- tag_hist[j,'event_site_latitude_value']
        ind_det_hist[counter,'event_site_longitude'] <- tag_hist[j,'event_site_longitude_value']
        
        # store the start time
        ind_det_hist[counter,'start_time'] <- tag_hist[[j,'event_date_time_value']]
        
      }
      
      # If the previous entry was the same site <=6 hours ago AND
      # The next entry is the same site <=6 hours ago, skip it
      else {
        
      }
      
    }
  }
  
  # Append the individual tag history to the complete tag history
  JDR_det_hist %>% 
    bind_rows(., ind_det_hist) -> JDR_det_hist
  
  # Print the run time
  if (i == length(unique_tag_IDs)){
    end_time <- Sys.time()
    
    print(paste0("Total tags: ", length(unique_tag_IDs)))
    print(paste0("Total time: ", end_time - start_time))
  }
  
  
}
 
