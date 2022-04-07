### 07 - covariates

# This script will reformat covariate data for inclusion in the multistate model


# Load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)

##### Data exploration #####

# Load data from McNary

# Tailrace temperature
MCN_t_tailrace <- clean_names(read.csv(here::here("covariate_data","MCN","basin_tempc_tailrace_MCN_allyears.csv")))

# fix column names
colnames(MCN_t_tailrace) <- gsub("x", "", colnames(MCN_t_tailrace))

# Pivot longer
MCN_t_tailrace %>% 
  pivot_longer(., cols = colnames(MCN_t_tailrace)[2:ncol(MCN_t_tailrace)]) %>% 
  dplyr::rename(year = name, temp = value) %>% 
  mutate(day = as.numeric(day), year = as.factor(year)) %>% 
  arrange(year, day) -> MCN_t_tailrace_long

# Plot it

ggplot(MCN_t_tailrace_long, aes(x = day, y = temp, color = year)) +
  geom_line()

# Add a data filtering step to remove what appears to be instrument error

# Data filtering function
data_filter <- function(data, limit, var){
  # Flag values that are more than the limit different from the previous
  data %>% 
    mutate(drop = ifelse(is.na(eval(parse(text = var))), FALSE,
                         ifelse(abs(eval(parse(text = var)) - lag(eval(parse(text = var)))) > limit, TRUE, FALSE))) -> data
  
  clean_data <- list(clean = subset(data, drop == FALSE), dropped = subset(data, drop == TRUE))
  
  return(clean_data)
}

subset(data, drop == TRUE)

# this isn't quite working right and we don't need it, so we will just manually subset


# Store clean data
MCN_t_tailrace_long_clean <- subset(MCN_t_tailrace_long, temp < 28)

ggplot(MCN_t_tailrace_long_clean, aes(x = day, y = temp, color = year)) +
  geom_line()

#### Reformat into mean temperature for certain portions of year

MCN_t_tailrace_long_clean %>% 
  mutate(window = ifelse(day <= 90 & day >= 0, "winter",
                         ifelse(day <= 250 & day >= 160, "summer", "none"))) -> MCN_t_tailrace_long_clean

MCN_t_tailrace_long_clean %>% 
  group_by(year, window) %>% 
  summarise(mean(temp)) %>% 
  dplyr::rename(mean_temp = `mean(temp)`)-> MCN_t_tailrace_window_means
  
##### Look at distribution of run timing

adult_returns <- clean_names(read.csv(here::here("PTAGIS_queries", "intermediate_files", "all_adult_returns.csv")))

adult_returns %>% 
  mutate(event_date = mdy(first_obs_date_max))-> adult_returns

# Get the date of arrival at BON - take min in case it fell back over Bonneville
adult_returns %>% 
  group_by(tag_code) %>% 
  subset(site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                                "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF", 
                                "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) %>% 
  filter(event_date == min(event_date)) %>% 
  dplyr::select(tag_code, event_date) %>% 
  dplyr::rename(BON_arrival = event_date) -> adult_BON_arrival

# Add run year info

# Sort into run years
run_year <- c("05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21", "21/22")
run_year_start <- seq(ymd("2005-06-01"), ymd("2021-06-01"), by = "years")
run_year_end <- seq(ymd("2006-05-31"), ymd("2022-05-31"), by = "years")

run_year_df <- data.frame(run_year, run_year_start, run_year_end)

adult_BON_arrival %>% 
  mutate(dummy =TRUE) %>% 
  left_join(run_year_df %>% mutate(dummy=TRUE), by = "dummy") %>% 
  filter(run_year_start <= BON_arrival, run_year_end >= BON_arrival) %>% 
  select(-c(dummy, run_year_start, run_year_end)) -> adult_BON_arrival

# Get date timing
adult_BON_arrival %>% 
  mutate(DayMonth = format(as.Date(BON_arrival), "%m-%d")) %>% 
  mutate(DayMonth = yday(BON_arrival)) -> adult_BON_arrival


# Plot date
ggplot(data = adult_BON_arrival, aes(x = BON_arrival)) +
  geom_bar(stat = "count")

ggplot(data = adult_BON_arrival, aes(x = DayMonth)) +
  geom_bar(stat = "count") +
  xlab("Day of year") +
  ggtitle("Arrival at Bonneville, all adult Steelhead, 2005-2022")


# Plot John Day River steelhead, arrival at natal tributary
JDR_CTH_1 <- clean_names(read.csv(here::here("PTAGIS_queries", "complete_tag_histories", "john_day_river", "2022-01-14-john_day_river_CTH_2005_2015_1.csv")))
JDR_CTH_2 <- clean_names(read.csv(here::here("PTAGIS_queries", "complete_tag_histories", "john_day_river", "2022-01-14-john_day_river_CTH_2005_2015_2.csv")))

### Combine files, fix some column data types
JDR_CTH_1 %>% 
  bind_rows(., JDR_CTH_2) %>% 
  mutate(event_date_time_value = mdy_hms(event_date_time_value))-> JDR_CTH

# Get names of natal tributary sites
JDR_det_hist <- read.csv(file = here::here("model_files", "JDR_det_hist.csv"))
JDR_det_hist %>% 
  group_by(event_site_name) %>% 
  summarise(n()) %>% 
  dplyr::rename(count = `n()`)-> JDR_event_det_counts

# Get the metadata
JDR_det_hist %>% 
  dplyr::select(-c(tag_code, start_time, end_time)) %>% 
  distinct(event_site_name, .keep_all = TRUE) -> JDR_event_site_metadata

BON_MCN_natal_sites <- JDR_event_site_metadata$event_site_name[grep("John Day", JDR_event_site_metadata$event_site_basin_name)]

# Get the date of arrival at BON - take min in case it fell back over Bonneville
JDR_CTH %>% 
  group_by(tag_code) %>% 
  subset(event_site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                                "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF", 
                                "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) %>% 
  filter(event_date_time_value == min(event_date_time_value)) %>% 
  dplyr::select(tag_code, event_date_time_value) %>% 
  dplyr::rename(BON_arrival = event_date_time_value) -> JDR_BON_arrival

JDR_CTH %>% 
  group_by(tag_code) %>% 
  subset(event_site_name %in% BON_MCN_natal_sites) %>% 
  filter(event_date_time_value == max(event_date_time_value)) %>% # take max, since min would be release
  dplyr::select(tag_code, event_date_time_value) %>% 
  dplyr::rename(nat_trib_arrival = event_date_time_value) -> JDR_nat_trib_arrival

# Join the two, keep only those where the nat_trib_arrival is greater than the BON arrival

JDR_nat_trib_arrival %>% 
  left_join(., JDR_BON_arrival, by = "tag_code") %>% 
  filter(nat_trib_arrival > BON_arrival) %>% 
  # Get day and month
  mutate(DayMonth = yday(nat_trib_arrival)) -> JDR_nat_trib_final

ggplot(data = JDR_nat_trib_final, aes(x = DayMonth)) +
  geom_bar(stat = "count") +
  xlab("Day of year") +
  ggtitle("Arrival at natal tributary, JDR steelhead 05/06-14/15")


##### Window of time for in-stream covariates #####

# Step 1: Find median travel time between previous dam and next dam
# Import the shortened detection history
det_hist <- read.csv(here::here("model_files", "complete_det_hist.csv"), row.names = 1)
BON_arrival <- read.csv(here::here("covariate_data", "complete_BON_arrival.csv"))
origin_table <- read.csv(here::here("covariate_data", "natal_origin_table.csv"))
# Read in the metadata
tag_code_metadata <- read.csv(here::here("covariate_data", "tag_code_metadata.csv"))
tag_code_metadata %>% 
  left_join(., origin_table, by = "release_site_name") %>% 
  dplyr::select(tag_code, run_year, natal_origin, rear_type_code, release_site_name) %>% 
  subset(., tag_code %in% BON_arrival$tag_code) -> origin_metadata

# Load the states from the 03 script
read.csv(here::here("model_files", "states.csv")) %>% 
  dplyr::select(-X) -> states

# Check to see distribution of origins in this subset of the data
tag_code_subset <- unique(states$tag_code)
subset_origin_metadata <- subset(origin_metadata, tag_code %in% tag_code_subset)
table(subset_origin_metadata$natal_origin)
# Remove the Klickitat and Wind Rivers - shouldn't have to do this in the future, since I made sure to subset earlier
KLIC_WIND_tag_codes <- subset(origin_metadata, natal_origin %in% c("Klickitat_River", "Wind_River"))$tag_code

det_hist %>% 
  left_join(., origin_metadata, by = "tag_code") %>% 
  mutate(start_time = ymd_hms(start_time),
         end_time = ymd_hms(end_time)) %>% 
  # add a month field
  mutate(start_month = month(start_time))-> det_hist


# Create a function for this
median_travel_time_calc <- function(prev_dam, next_dam, data){
  # Keep only the parts of the df where the detection ns are prev_dam, then next dam
  data %>% 
    mutate(prev_site = lag(event_site_name),
           prev_time = lag(end_time),
           prev_tag_code = lag(tag_code)) %>% 
    mutate(prev_time = ymd_hms(prev_time),
           end_time = ymd_hms(end_time)) %>% 
    filter(., event_site_name == next_dam & prev_site == prev_dam) -> two_dam_data
  
  # New calculate the travel times and get the median
  two_dam_data %>% 
    mutate(travel_time = ifelse(tag_code == prev_tag_code, time_length(as.period(end_time - prev_time), unit = "days"),NA)) %>% 
    subset(travel_time >= 0) -> two_dam_movement_data
  
  median_travel_time <- median(two_dam_movement_data$travel_time)
  
  return(median_travel_time)
    
}

# Make an alternative function to calculate travel times but not median, so you can subset
travel_time_calc <- function(prev_dam, next_dam, data){
  # Keep only the parts of the df where the detection ns are prev_dam, then next dam
  data %>% 
    mutate(prev_site = lag(event_site_name),
           prev_time = lag(end_time),
           prev_tag_code = lag(tag_code)) %>% 
    mutate(prev_time = ymd_hms(prev_time),
           end_time = ymd_hms(end_time)) %>% 
    filter(., event_site_name == next_dam & prev_site == prev_dam) -> two_dam_data
  
  # New calculate the travel times and get the median
  two_dam_data %>% 
    mutate(travel_time = ifelse(tag_code == prev_tag_code, time_length(as.period(end_time - prev_time), unit = "days"),NA)) %>% 
    subset(travel_time >= 0) -> two_dam_movement_data
  
  return(two_dam_movement_data)
  
}

ggplot(test2, aes(x = travel_time)) +
  geom_histogram()

# Loop through this function for each set of sequential dams, each origin, and each month
dams <-
  c(
    "Bonneville Adult Fishways (combined)",
    "McNary Adult Fishways (combined)",
    "ICH - Ice Harbor Dam (Combined)",
    "LMA - Lower Monumental Adult Ladders",
    "GOA - Little Goose Fish Ladder",
    "Lower Granite Dam Adult Fishways (combined)",
    "PRA - Priest Rapids Adult",
    "RIA - Rock Island Adult",
    "RRF - Rocky Reach Fishway",
    "WEA - Wells Dam, DCPUD Adult Ladders"
  )

# make a table with dam, upstream dam, downstream dam
dam_relationships <- data.frame(dam = dams, 
                                downstream_dam = c(NA, 
                                                   "Bonneville Adult Fishways (combined)",
                                                   "McNary Adult Fishways (combined)",
                                                   "ICH - Ice Harbor Dam (Combined)",
                                                   "LMA - Lower Monumental Adult Ladders",
                                                   "GOA - Little Goose Fish Ladder",
                                                   "McNary Adult Fishways (combined)",
                                                   "PRA - Priest Rapids Adult",
                                                   "RIA - Rock Island Adult",
                                                   "RRF - Rocky Reach Fishway"),
                                upstream_dam = c("McNary Adult Fishways (combined)",
                                                 "ICH - Ice Harbor Dam (Combined)",
                                                 "LMA - Lower Monumental Adult Ladders",
                                                 "GOA - Little Goose Fish Ladder",
                                                 "Lower Granite Dam Adult Fishways (combined)",
                                                 NA,
                                                 "RIA - Rock Island Adult",
                                                 "RRF - Rocky Reach Fishway",
                                                 "WEA - Wells Dam, DCPUD Adult Ladders",
                                                 NA))

# Get natal origins
unique_origins <- unique(det_hist$natal_origin)[!is.na(unique(det_hist$natal_origin))]
# get months
unique_months <- unique(det_hist$start_month)

# Create an empty data frame to store all of the median travel times
median_travel_time_df <- crossing(unique_origins, unique_months, dam_relationships$dam, dam_relationships$downstream_dam, dam_relationships$upstream_dam)
median_travel_time_df <- crossing(unique_origins, unique_months, dam_relationships$dam)
colnames(median_travel_time_df) <- c("natal_origin","month","dam")
median_travel_time_df %>% 
  left_join(., dam_relationships, by = "dam") -> median_travel_time_df

# Add variables to store upstream and downstream travel times
median_travel_time_df %>% 
  mutate(median_upstream_travel_time = NA,
         median_downstream_travel_time = NA) -> median_travel_time_df


for (i in 1:length(unique_origins)){
  # continually update unique months for each natal origin
  unique_months <- unique(subset(det_hist, natal_origin == unique_origins[i])$start_month)
  
  for (j in 1:length(unique_months)){
    origin <- unique_origins[i]
    current_month <- unique_months[j]
    origin_month_data <- subset(det_hist, start_month == current_month & natal_origin == origin)
    
    # get the unique dam combinations for this dataset
    unique_sites <- unique(origin_month_data$event_site_name)
    unique_dams <- unique_sites[unique_sites %in% dams]
    
    # Create an if else statement to only loop through if there are any dams in the detection history
    for (k in 1:length(unique_dams)){
      current_dam <- unique_dams[k]
      upstream_dam <- subset(dam_relationships, dam == current_dam)$upstream_dam
      downstream_dam <- subset(dam_relationships, dam == current_dam)$downstream_dam
      
      print(paste0("Current Dam: ", current_dam,
                   "; Upstream Dam: ", upstream_dam,
                   "; Downstream Dam: ", downstream_dam,
                   "; Origin: ", origin,
                   "; Current Month: ", current_month))
      
      median_upstream_travel_time <- median_travel_time_calc(prev_dam = current_dam,
                              next_dam = upstream_dam,
                              data = origin_month_data)
      median_downstream_travel_time <- median_travel_time_calc(prev_dam = current_dam,
                                                               next_dam = downstream_dam,
                                                               data = origin_month_data)
      
      median_travel_time_df$median_upstream_travel_time[median_travel_time_df$natal_origin == origin &
                                                          median_travel_time_df$month == current_month & 
                                                          median_travel_time_df$dam == current_dam]<- median_upstream_travel_time
      
      
      median_travel_time_df$median_downstream_travel_time[median_travel_time_df$natal_origin == origin &
                                                          median_travel_time_df$month == current_month & 
                                                          median_travel_time_df$dam == current_dam]<- median_downstream_travel_time
      
      }
  }
}


## v2

for (i in 1:length(unique_origins)){
  origin <- unique_origins[i]
  origin_data <- subset(det_hist, natal_origin == origin)
  
  # get the unique dam combinations for this dataset
  unique_sites <- unique(origin_data$event_site_name)
  unique_dams <- unique_sites[unique_sites %in% dams]

  for (j in 1:length(unique_dams)){
      current_dam <- unique_dams[j]
      upstream_dam <- subset(dam_relationships, dam == current_dam)$upstream_dam
      downstream_dam <- subset(dam_relationships, dam == current_dam)$downstream_dam
      

      # get a df with individual travel times
      upstream_travel_time <- travel_time_calc(prev_dam = current_dam,
                                                             next_dam = upstream_dam,
                                                             data = origin_data)
      downstream_travel_time <- travel_time_calc(prev_dam = current_dam,
                                                               next_dam = downstream_dam,
                                                               data = origin_data)
        
        
        # get unique months
        unique_months <- unique(origin_data$start_month)
        
        # Loop through the unique months, subset to get median travel time for these
        for (k in 1:length(unique_months)){
          current_month <- unique_months[k]
          
          print(paste0("Current Dam: ", current_dam,
                       "; Upstream Dam: ", upstream_dam,
                       "; Downstream Dam: ", downstream_dam,
                       "; Origin: ", origin,
                       "; Current Month: ", current_month))
          
          # get the upstream travel time
          origin_month_upstream_travel_time <- subset(upstream_travel_time, start_month == current_month)
          median_upstream_travel_time <- median(origin_month_upstream_travel_time$travel_time)
          
          # get the downstream travel time
          origin_month_downstream_travel_time <- subset(downstream_travel_time, start_month == current_month)
          median_downstream_travel_time <- median(origin_month_downstream_travel_time$travel_time)
      
        median_travel_time_df$median_upstream_travel_time[median_travel_time_df$natal_origin == origin &
                                                          median_travel_time_df$month == current_month & 
                                                          median_travel_time_df$dam == current_dam]<- median_upstream_travel_time
      
      
        median_travel_time_df$median_downstream_travel_time[median_travel_time_df$natal_origin == origin &
                                                            median_travel_time_df$month == current_month & 
                                                            median_travel_time_df$dam == current_dam]<- median_downstream_travel_time
      
    }
  }
}



# Step 2: For each individual fish, create a weeklong window center around the median travel time, i.e. +/- 3 days

# rename the median travel time df to prepare for join
median_travel_time_df %>% 
  dplyr::rename(start_month = month, event_site_name = dam) -> median_travel_time_df_2

det_hist %>% 
  # add the travel time to each detection
  left_join(., median_travel_time_df_2, by = c("natal_origin", "start_month", "event_site_name")) -> det_hist_2

# ah, so this won't work. We need to add this to the complete detection history, to match by state

##### Interpolate times at states #####

# The upstream vs. downstream travel time approach won't work, because we can't assume in what direction
# they were traveling. Instead, we will use the time halfway between the time arriving in the current state
# and the time arriving in the next state as the time for covariates
# Now, when we have implicit site usage, we don't know when they arrived. So in this
# case, we will interpolate these times by taking the halfway point between two known times. If
# there are multiple (N) implict site visits between known times, then we will take N 
# breakpoints; i.e., if we saw an individual at state 1 at t = 1, and at state 4 at t = 10,
# we will say that they were at state 2 at t = 4 and state 3 at t = 7

# Load the states from the 03 script
read.csv(here::here("model_files", "states.csv")) %>% 
  dplyr::select(-X) -> states

# Create a column to indicate if the date_time was interpolated
# Note the arrival date
states %>% 
  mutate(date_source = ifelse(!is.na(date_time), "known_arrival", "interpolated")) %>% 
  # get an empty date column
  mutate(arrival_date = NA) -> states

# I think we have to loop this

# First, figure out indices of missing date_time
missing_date_time <- is.na(states$date_time)
date_times <- states$date_time

# Loop to get arrival dates in state
for (i in 1:nrow(states)){
  # if it's NA
  print(paste0("Row #", i))
  if (is.na(states[i,"date_time"])){
    # Figure out how far back the last known time was
    # Truncate the missing date times to only ones prior to current state
    missing_date_time_subset <- missing_date_time[1:(i-1)]
    prev_time_index <- max(which(missing_date_time_subset %in% FALSE))
    
    # Truncate the missing date times to only ones after current state
    missing_date_time_subset <- missing_date_time[(i+1):nrow(states)]
    next_time_index <- min(which(missing_date_time_subset %in% FALSE)) + i
    
    # Now, interpolate the missing time
    # First, figure out how long between the two known times
    prev_time <- ymd_hms(date_times[prev_time_index])
    next_time <- ymd_hms(date_times[next_time_index])
    time_diff <- next_time - prev_time
    
    # Get the missing time - add the time difference divided by the number 
    # of missing steps plus 1, multiply by which number step it is
    missing_time <- prev_time + (time_diff/(next_time_index - prev_time_index) * (i - prev_time_index))
    
    # Extract just the date
    missing_date <- date(missing_time)
    
    # populate the missing date_time
    states[i, "arrival_date"] <- format(as.Date(missing_date, origin = "1970-01-01"))
    
  }
  # If it's not, leave it alone and just extract date
  else {
    states[i, "arrival_date"] <- format(as.Date(date(states[i,"date_time"]), origin = "1970-01-01"))
  }
}

# Loop again to get the "experience date" - halfway between arrival 
# date in state and arrival date in next state

states %>% 
  mutate(date = NA) %>% 
  mutate(state_time = NA) %>% 
  mutate(arrival_date = ymd(arrival_date))-> states

for (i in 1:(nrow(states) - 1)) {
  # If the next fish is the same fish
  if (i == 1 | (states[i, "tag_code"] == states[i+1, "tag_code"])){
    movement_time <- states[i+1, "arrival_date"] - states[i, "arrival_date"]
    
    states[i, "date"] <- format(as.Date(date(states[i, "arrival_date"] + movement_time/2), origin = "1970-01-01"))
    states[i, "state_time"] <- movement_time
    
  }
  
  # If it's the last entry for a fish, just store the arrival time
  else {
    states[i, "date"] <- format(as.Date(date(states[i, "arrival_date"]), origin = "1970-01-01"))
    states[i, "state_time"] <- NA # populate time in state with NA, since we don't know
  }
}

# Get the next state - useful to determine time in state (possible covariate)
states %>% 
  mutate(next_state = NA) -> states

# Get next state, prev state
for (i in 1:(nrow(states) - 1)) {
  # If the next fish is the same fish
  if (i == 1 | (states[i, "tag_code"] == states[i+1, "tag_code"])){
    states[i, "next_state"] <- states[i+1, "state"]
    
  }
  
  # If it's the last entry for a fish, just store the arrival time
  else {
    states[i, "next_state"] <- NA
  }
}

# Export this as states_time

write.csv(states, here::here("model_files", "states_times.csv"))




##### Import dam covariate data #####
# Pivot longer and reformat function
cov_reformat <- function(input_data, variable_name){
  input_data %>% 
    pivot_longer(., cols = colnames(.)[2:ncol(.)]) %>% 
    dplyr::rename(year = name,
                  !!quo_name(variable_name) := value) %>% 
    mutate(year = gsub("X","",year)) %>% 
    mutate(year = as.numeric(year), day = as.numeric(day)) %>% 
    # Add a new column for date (combine year and day)
    mutate(date = format(as.Date(day, origin = paste0(year, "-01-01")))) %>% 
    arrange(year, day) -> output_data
  
  return(output_data)
}

# BON
BON_flow <- read.csv(here::here("covariate_data", "BON", "basin_outflow_BON_allyears.csv"))[1:366,]
BON_flow_long <- cov_reformat(input_data = BON_flow, variable_name = "flow")
BON_spill <- read.csv(here::here("covariate_data", "BON", "basin_spill_BON_allyears.csv"))[1:366,]
BON_spill_long <- cov_reformat(input_data = BON_spill, variable_name = "spill")
BON_temp <- read.csv(here::here("covariate_data", "BON", "basin_tempc_tailrace_BON_allyears.csv"))[1:366,]
BON_temp_long <- cov_reformat(input_data = BON_temp, variable_name = "temp")


# ICH
ICH_flow <- read.csv(here::here("covariate_data", "ICH", "basin_outflow_IHR_allyears.csv"))[1:366,]
ICH_flow_long <- cov_reformat(input_data = ICH_flow, variable_name = "flow")
ICH_spill <- read.csv(here::here("covariate_data", "ICH", "basin_spill_IHR_allyears.csv"))[1:366,]
ICH_spill_long <- cov_reformat(input_data = ICH_spill, variable_name = "spill")
ICH_temp <- read.csv(here::here("covariate_data", "ICH", "basin_tempc_tailrace_IHR_allyears.csv"))[1:366,]
ICH_temp_long <- cov_reformat(input_data = ICH_temp, variable_name = "temp")


# JDA
JDA_flow <- read.csv(here::here("covariate_data", "JDA", "basin_outflow_JDA_allyears.csv"))[1:366,]
JDA_flow_long <- cov_reformat(input_data = JDA_flow, variable_name = "flow")
JDA_spill <- read.csv(here::here("covariate_data", "JDA", "basin_spill_JDA_allyears.csv"))[1:366,]
JDA_spill_long <- cov_reformat(input_data = JDA_spill, variable_name = "spill")
JDA_temp <- read.csv(here::here("covariate_data", "JDA", "basin_tempc_tailrace_JDA_allyears.csv"))[1:366,]
JDA_temp_long <- cov_reformat(input_data = JDA_temp, variable_name = "temp")


# LGR
LGR_flow <- read.csv(here::here("covariate_data", "LGR", "basin_outflow_LWG_allyears.csv"))[1:366,]
LGR_flow_long <- cov_reformat(input_data = LGR_flow, variable_name = "flow")
LGR_spill <- read.csv(here::here("covariate_data", "LGR", "basin_spill_LWG_allyears.csv"))[1:366,]
LGR_spill_long <- cov_reformat(input_data = LGR_spill, variable_name = "spill")
LGR_temp <- read.csv(here::here("covariate_data", "LGR", "basin_tempc_tailrace_LWG_allyears.csv"))[1:366,]
LGR_temp_long <- cov_reformat(input_data = LGR_temp, variable_name = "temp")


# MCN
MCN_flow <- read.csv(here::here("covariate_data", "MCN", "basin_outflow_MCN_allyears.csv"))[1:366,]
MCN_flow_long <- cov_reformat(input_data = MCN_flow, variable_name = "flow")
MCN_spill <- read.csv(here::here("covariate_data", "MCN", "basin_spill_MCN_allyears.csv"))[1:366,]
MCN_spill_long <- cov_reformat(input_data = MCN_spill, variable_name = "spill")
MCN_temp <- read.csv(here::here("covariate_data", "MCN", "basin_tempc_tailrace_MCN_allyears.csv"))[1:366,]
MCN_temp_long <- cov_reformat(input_data = MCN_temp, variable_name = "temp")


# PRA
PRA_flow <- read.csv(here::here("covariate_data", "PRA", "basin_outflow_PRD_allyears.csv"))[1:366,]
PRA_flow_long <- cov_reformat(input_data = PRA_flow, variable_name = "flow")
PRA_spill <- read.csv(here::here("covariate_data", "PRA", "basin_spill_PRD_allyears.csv"))[1:366,]
PRA_spill_long <- cov_reformat(input_data = PRA_spill, variable_name = "spill")
PRA_temp <- read.csv(here::here("covariate_data", "PRA", "basin_tempc_tailrace_PRD_allyears.csv"))[1:366,]
PRA_temp_long <- cov_reformat(input_data = PRA_temp, variable_name = "temp")


# RIS
RIS_flow <- read.csv(here::here("covariate_data", "RIS", "basin_outflow_RIS_allyears.csv"))[1:366,]
RIS_flow_long <- cov_reformat(input_data = RIS_flow, variable_name = "flow")
RIS_spill <- read.csv(here::here("covariate_data", "RIS", "basin_spill_RIS_allyears.csv"))[1:366,]
RIS_spill_long <- cov_reformat(input_data = RIS_spill, variable_name = "spill")
RIS_temp <- read.csv(here::here("covariate_data", "RIS", "basin_tempc_tailrace_RIS_allyears.csv"))[1:366,]
RIS_temp_long <- cov_reformat(input_data = RIS_temp, variable_name = "temp")


# RRE
RRE_flow <- read.csv(here::here("covariate_data", "RRE", "basin_outflow_RRH_allyears.csv"))[1:366,]
RRE_flow_long <- cov_reformat(input_data = RRE_flow, variable_name = "flow")
RRE_spill <- read.csv(here::here("covariate_data", "RRE", "basin_spill_RRH_allyears.csv"))[1:366,]
RRE_spill_long <- cov_reformat(input_data = RRE_spill, variable_name = "spill")
RRE_temp <- read.csv(here::here("covariate_data", "RRE", "basin_tempc_tailrace_RRH_allyears.csv"))[1:366,]
RRE_temp_long <- cov_reformat(input_data = RRE_temp, variable_name = "temp")


# WEL
WEL_flow <- read.csv(here::here("covariate_data", "WEL", "basin_outflow_WEL_allyears.csv"))[1:366,]
WEL_flow_long <- cov_reformat(input_data = WEL_flow, variable_name = "flow")
WEL_spill <- read.csv(here::here("covariate_data", "WEL", "basin_spill_WEL_allyears.csv"))[1:366,]
WEL_spill_long <- cov_reformat(input_data = WEL_spill, variable_name = "spill")
WEL_temp <- read.csv(here::here("covariate_data", "WEL", "basin_tempc_tailrace_WEL_allyears.csv"))[1:366,]
WEL_temp_long <- cov_reformat(input_data = WEL_temp, variable_name = "temp")

# Reformat this into one big DF - date, plus value at each state

# Three: one for flow, one for spill, one for temp

##### Temperature #####

dplyr::rename(BON_temp_long, "mainstem, mouth to BON" = temp) %>% 
  full_join(., dplyr::rename(ICH_temp_long, "mainstem, MCN to ICH or PRA (ICH)" = temp), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(LGR_temp_long, "mainstem, ICH to LGR" = temp), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(MCN_temp_long, "mainstem, BON to MCN" = temp), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(PRA_temp_long, "mainstem, MCN to ICH or PRA (PRA)" = temp), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(RIS_temp_long, "mainstem, PRA to RIS" = temp), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(RRE_temp_long, "mainstem, RIS to RRE" = temp), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(WEL_temp_long, "mainstem, RRE to WEL" = temp), by = c("day", "year", "date")) -> temp_cov_df

temp_cov_df

# Shorten column names
colnames(temp_cov_df) <- gsub("mainstem, ", "", colnames(temp_cov_df))
# Reorder
temp_cov_df %>% 
  dplyr::select(-c(day, year)) %>% 
  dplyr::relocate(., date) -> temp_cov_df

# Save this
write.csv(temp_cov_df, here::here("covariate_data", "temperature_by_state.csv"))



# How do we deal with MCN to ICH or PRA state?
# I figure that maybe they could be different covariates, for overshoot probability at each dam?


##### Spill #####


dplyr::rename(BON_spill_long, "mainstem, mouth to BON" = spill) %>% 
  full_join(., dplyr::rename(ICH_spill_long, "mainstem, MCN to ICH or PRA (ICH)" = spill), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(LGR_spill_long, "mainstem, ICH to LGR" = spill), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(MCN_spill_long, "mainstem, BON to MCN" = spill), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(PRA_spill_long, "mainstem, MCN to ICH or PRA (PRA)" = spill), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(RIS_spill_long, "mainstem, PRA to RIS" = spill), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(RRE_spill_long, "mainstem, RIS to RRE" = spill), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(WEL_spill_long, "mainstem, RRE to WEL" = spill), by = c("day", "year", "date")) -> spill_cov_df


##### Flow #####

dplyr::rename(BON_flow_long, "mainstem, mouth to BON" = flow) %>% 
  full_join(., dplyr::rename(ICH_flow_long, "mainstem, MCN to ICH or PRA (ICH)" = flow), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(LGR_flow_long, "mainstem, ICH to LGR" = flow), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(MCN_flow_long, "mainstem, BON to MCN" = flow), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(PRA_flow_long, "mainstem, MCN to ICH or PRA (PRA)" = flow), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(RIS_flow_long, "mainstem, PRA to RIS" = flow), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(RRE_flow_long, "mainstem, RIS to RRE" = flow), by = c("day", "year", "date")) %>% 
  full_join(., dplyr::rename(WEL_flow_long, "mainstem, RRE to WEL" = flow), by = c("day", "year", "date")) -> flow_cov_df


##### Index to match date at state to covariates




