# stepwise_states_cleaning

library(here)
library(tidyverse)
library(lubridate)
library(janitor)

# Dealing with juveniles, kelts, and repeat spawners

# Input file: complete detection history from 03 script



# Here we will remove any juvenile movements and identify kelt and repeat spawner movements

# Combine the tag codes so that we can figure out what file to look in for each file
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_1.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 1) -> tag_codes_1
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_2.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 2) -> tag_codes_2
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_3.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 3) -> tag_codes_3
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_4.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 4) -> tag_codes_4
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_5.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 5) -> tag_codes_5
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_6.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 6) -> tag_codes_6
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_7.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 7) -> tag_codes_7
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_8.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 8) -> tag_codes_8
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_9.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 9) -> tag_codes_9
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_10.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 10) -> tag_codes_10
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_11.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 11) -> tag_codes_11
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_12.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 12) -> tag_codes_12
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_13.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 13) -> tag_codes_13
read.table(here::here("PTAGIS_queries", "intermediate_files", "tag_codes", "tag_codes_14.txt")) %>% 
  dplyr::rename(tag_code = V1) %>% 
  mutate(file = 14) -> tag_codes_14

tag_codes_1 %>% 
  bind_rows(., tag_codes_2) %>% 
  bind_rows(., tag_codes_3) %>% 
  bind_rows(., tag_codes_4) %>% 
  bind_rows(., tag_codes_5) %>% 
  bind_rows(., tag_codes_6) %>% 
  bind_rows(., tag_codes_7) %>% 
  bind_rows(., tag_codes_8) %>% 
  bind_rows(., tag_codes_9) %>% 
  bind_rows(., tag_codes_10) %>% 
  bind_rows(., tag_codes_11) %>% 
  bind_rows(., tag_codes_12) %>% 
  bind_rows(., tag_codes_13) %>% 
  bind_rows(., tag_codes_14) -> tag_code_file_finder

# Read in states complete
# read.csv(here::here("from_hyak_transfer", "2022-05-25-complete_det_hist", "states_complete.csv")) %>%
# read.csv(here::here("from_hyak_transfer", "2022-07-16-complete_det_hist", "states_complete.csv")) %>%
# read.csv(here::here("from_hyak_transfer", "2022-07-18-complete_det_hist", "states_complete.csv")) %>%
# read.csv(here::here("from_hyak_transfer", "2022-07-21-complete_det_hist", "states_complete.csv")) %>%
#   dplyr::select(-X) -> states_complete

# states_complete_part1 <- read.csv(here::here("from_hyak_transfer", "2022-08-08-complete_det_hist_ckpt", "states_complete_part1.csv"), row.names = 1)
# states_complete_part2 <- read.csv(here::here("from_hyak_transfer", "2022-08-08-complete_det_hist_ckpt", "states_complete_part2.csv"), row.names = 1)
# states_complete_part3 <- read.csv(here::here("from_hyak_transfer", "2022-08-08-complete_det_hist_ckpt", "states_complete_part3.csv"), row.names = 1)
# states_complete_part4 <- read.csv(here::here("from_hyak_transfer", "2022-08-08-complete_det_hist_ckpt", "states_complete_part4.csv"), row.names = 1)

states_complete_part1 <- read.csv(here::here("from_hyak_transfer", "2022-10-27-complete_det_hist_ckpt", "states_complete_part1.csv"), row.names = 1)
states_complete_part2 <- read.csv(here::here("from_hyak_transfer", "2022-10-27-complete_det_hist_ckpt", "states_complete_part2.csv"), row.names = 1)
states_complete_part3 <- read.csv(here::here("from_hyak_transfer", "2022-10-27-complete_det_hist_ckpt", "states_complete_part3.csv"), row.names = 1)
states_complete_part4 <- read.csv(here::here("from_hyak_transfer", "2022-10-27-complete_det_hist_ckpt", "states_complete_part4.csv"), row.names = 1)

states_complete_part1 %>% 
  bind_rows(., states_complete_part2) %>% 
  bind_rows(., states_complete_part3) %>% 
  bind_rows(., states_complete_part4) -> states_complete

# Get rid of dummy fish
states_complete %>% 
  subset(tag_code != "dummy_fish") -> states_complete

# Read in tag code metadata
read.csv(here::here("covariate_data", "tag_code_metadata.csv")) -> tag_code_metadata

# Fix the release dates
tag_code_metadata %>% 
  mutate(release_date = mdy(release_date_mmddyyyy)) -> tag_code_metadata

# Join release date with the states complete
states_complete %>% 
  left_join(dplyr::select(tag_code_metadata, tag_code, release_date), by = "tag_code") %>% 
  mutate(date_time = ymd_hms(date_time)) -> states_complete

# First, we're going to have to interpolate the times for implicit visits for the code to work
# states_complete %>% 
#   # mutate(date_time = ifelse(pathway == "implicit", format(as_datetime(lag(date_time) + (lead(date_time) - lag(date_time))/2, origin = "1970-01-01 00:00:00")),
#   #                           format(as_datetime(date_time, origin = "1970-01-01 00:00:00")))) -> states_complete
#   mutate(date_time = ifelse(pathway == "implicit", format(as_datetime(lag(date_time) + (lead(date_time) - lag(date_time))/2)),
#                             format(as_datetime(date_time)))) %>% 
#   mutate(date_time = ymd_hms(date_time)) -> states_complete

# First, figure out indices of missing date_time
missing_date_time <- is.na(states_complete$date_time)
date_times <- states_complete$date_time

notime_tags <- c("3D9.1BF15C4665", "3D9.1BF1608BA7", "3D9.1BF16B2E82", "3D9.1BF17A3856", "3D9.1BF1937E5E", "3D9.1BF1CAB788",
                 "3D9.1BF1CBCB02", "3D9.1BF20BDB85", "3D9.1BF2803D80", "3D9.1C2C4C0362",
                 "3D9.1C2C594B70", "3D9.1C2C7FE403", "3D9.1C2CE3EEAD", "3D9.1C2CFFD8AE", "3D9.1C2D681A4E", "3D9.1C2D7C28C1",
                 "3D9.1C2DAC38B6", "3D9.1C2DCAC7C3", "3D9.1C2DEC13C7", "3D9.257C5CDF73",
                 "3D9.257C5D22AA", "3DD.003BC73B1A", "3DD.00778B221C", "3DD.0077BD17FF")

# notime_tags_times <- subset(complete_det_hist, tag_code %in% notime_tags & event_site_name == "Bonneville Adult Fishways (combined)" & is.na(end_time))

# for (i in 1:nrow(states_complete)){
#   # this is the temporary fix
#   if (states_complete[i,"tag_code"] %in% notime_tags & states_complete[i,"pathway"] == "BON (adult)"){
#     states_complete[i, "date_time"] <- subset(notime_tags_times, tag_code == states_complete[i,"tag_code"])$start_time
#   }
#   # end temporary fix
#   
#   else if (is.na(states_complete[i,"date_time"])){
#     # Figure out how far back the last known time was
#     # Truncate the missing date times to only ones prior to current state
#     missing_date_time_subset <- missing_date_time[1:(i-1)]
#     prev_time_index <- max(which(missing_date_time_subset %in% FALSE))
#     
#     # Truncate the missing date times to only ones after current state
#     missing_date_time_subset <- missing_date_time[(i+1):nrow(states_complete)]
#     next_time_index <- min(which(missing_date_time_subset %in% FALSE)) + i
#     
#     # Now, interpolate the missing time
#     # First, figure out how long between the two known times
#     prev_time <- ymd_hms(date_times[prev_time_index])
#     next_time <- ymd_hms(date_times[next_time_index])
#     time_diff <- next_time - prev_time
#     
#     # Get the missing time - add the time difference divided by the number 
#     # of missing steps plus 1, multiply by which number step it is
#     missing_time <- prev_time + (time_diff/(next_time_index - prev_time_index) * (i - prev_time_index))
#     
#     # populate the missing date_time
#     states_complete[i, "date_time"] <- missing_time
#     
#   }
#   # If it's not, leave it alone (do nothing)
#   else {
#     
#   }
# }


# V2: don't get middle time, instead just pick next actual time
for (i in 1:nrow(states_complete)){
  # this is the temporary fix
  # if (states_complete[i,"tag_code"] %in% notime_tags & states_complete[i,"pathway"] == "BON (adult)"){
  #   states_complete[i, "date_time"] <- subset(notime_tags_times, tag_code == states_complete[i,"tag_code"])$start_time
  # }
  # end temporary fix

  # else if (is.na(states_complete[i,"date_time"])){
  if (is.na(states_complete[i,"date_time"])){
    # Find the next known time
    # Truncate the missing date times to only ones after current state
    missing_date_time_subset <- missing_date_time[(i+1):nrow(states_complete)]
    next_time_index <- min(which(missing_date_time_subset %in% FALSE)) + i

    # Now, interpolate the missing time
    # next_time <- ymd_hms(date_times[next_time_index])
    next_time <- format(as.POSIXct(date_times[next_time_index], tz = "UTC"), "%Y-%m-%d %H:%M:%S")

    # Get the missing time - add the time difference divided by the number
    # of missing steps plus 1, multiply by which number step it is
    missing_time <- ymd_hms(next_time)

    # populate the missing date_time
    # Edit - force PST
    # apparently it's not actually necessary
    states_complete[i, "date_time"] <- missing_time
    # states_complete[i, "date_time"] <- force_tz(missing_time, "PST")

  }
  # If it's not, leave it alone (do nothing)
  else {

  }
}



# Put this checkpoint in so we don't have to re-run the for loop
write.csv(states_complete, here::here("from_hyak_transfer", "2022-10-27-complete_det_hist_ckpt", "states_complete_times_interpolated.csv"))

# read.csv(here::here("from_hyak_transfer", "2022-07-16-complete_det_hist", "states_complete_times_interpolated.csv")) %>%
# read.csv(here::here("from_hyak_transfer", "2022-07-18-complete_det_hist", "states_complete_times_interpolated.csv")) %>%
read.csv(here::here("from_hyak_transfer", "2022-10-27-complete_det_hist_ckpt", "states_complete_times_interpolated.csv")) %>%
  dplyr::select(-X) -> states_complete

# Okay, I'm not exactly sure why these fish are missing date/times, but manually fix them
# states_complete[174991,"date_time"] <- states_complete[174992,"date_time"]
# states_complete[270821,"date_time"] <- states_complete[270822,"date_time"]


states_complete %>% 
  subset(tag_code != "dummy_fish") -> states_complete

# Fix date time
states_complete %>% 
  mutate(date_time = ymd_hms(date_time)) -> states_complete

# Figure out how long detections were after release time
states_complete %>% 
  mutate(time_after_release = date_time - as.POSIXct(release_date)) -> states_complete

# Take apart the date time into y, m, d for subsetting
states_complete %>% 
  mutate(release_year = year(release_date)) %>% 
  mutate(release_month = month(release_date)) %>% 
  mutate(release_day = day(release_date)) %>% 
  mutate(event_year = year(date_time)) %>% 
  mutate(event_month = month(date_time)) %>% 
  mutate(event_day = day(date_time)) -> states_complete

##### Identify juvenile movements #####

# Identify juvenile movements as:
# 1) within 30 days of release; or
# 2) On or before June 15 of the release year

states_complete %>% 
  mutate(life_stage = "Adult") %>% #first, classify all movement as adult (we will then change this to kelt, repeat spawner, and juvenile later if certain conditions are met)
  mutate(life_stage = ifelse(time_after_release <= days(x = 90) | # fish that were seen right after release
  # mutate(life_stage = ifelse(
                             # or fish that were seen the same year as release and before June 15
                             release_year == event_year & event_month <= 5 | 
                               release_year == event_year & event_month == 6 & event_day <= 15 |
                            # or fish that were released the year before (on or after July 1) and seen the subsequent year on or before June 15
                               release_year == event_year-1 & release_month >= 7 & event_month <= 5 | 
                               release_year == event_year-1 & release_month >= 7 & event_month == 6 & event_day <= 15, "Juvenile",  
  life_stage)) -> states_complete

# Remove all detection histories that are only juvenile (keep only those detection histories with at least one adult movement)
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "Adult")) -> states_complete

# Okay, so we need to fix our arrival times to not pick the juvenile arrival
# This tag code (3D9.1BF2803D80) somehow has three of the same observation - BON adult on 2014-09-10.
# I think this was a problem tag from a previous script too, so it's coming from before

# subset(states_complete, tag_code == "3D9.1BF2803D80") -> bad_boy
# Looks good now!

# Let's figure out the first ADULT arrival at BON

# EDIT 08-11-2022: We need to allow the adult arrival to NOT be an ascent - if the first attempt is aborted, it's still when they arrive
states_complete %>% 
  group_by(tag_code) %>% 
  subset(pathway %in% c("BON (adult)", "BON_aborted")) %>% 
  subset(life_stage != "Juvenile") %>% 
  filter(date_time == min(date_time)) %>% 
  dplyr::select(tag_code, date_time) %>% 
  dplyr::rename(BON_arrival = date_time) -> BON_arrival_df

# merge that back in
states_complete %>% 
  left_join(., BON_arrival_df, by = "tag_code") -> states_complete

# make sure it's a date
# states_complete %>% 
#   mutate(BON_arrival = ymd_hms(BON_arrival)) -> states_complete

# calculate time after adult arrival
# this line takes quite a while
states_complete %>% 
  mutate(time_after_arrival = date_time - as.POSIXct(BON_arrival)) -> states_complete





# One correction: if an implicit state follows a juvenile movement, it also has to be a juvenile movement. This way, we will remove
# the implicit mouth to BON state that we get when we see a juvenile in the adult ladder and return as an adult. That state
# visit in between shouldn't be modeled
states_complete %>% 
  mutate(life_stage = ifelse(pathway == "implicit" & lag(life_stage) == "Juvenile", "Juvenile", life_stage)) -> states_complete

checkpoint_states_complete <- states_complete

# Another correction: if mouth to BON follows a juvenile movement, it also has to be a juvenile movement. We want all fish starting when they reach BON (adult)
# Edit 2022-08-11: exclude mouth to BON that shows up because of an aborted BON attempt
# Second fix - make sure there are no NAs, so keep all juveniles as juveniles here
states_complete %>% 
  mutate(life_stage = ifelse(life_stage == "Juvenile", "Juvenile", 
                             # Make a fix here, so that any detections prior to arrival date at BON are juvenile. This will fix individuals seen at BHL
                             ifelse(time_after_arrival < 0, "Juvenile",
                             ifelse(state == "mainstem, mouth to BON" & lag(life_stage) == "Juvenile" & pathway != "BON_aborted", "Juvenile", life_stage)))) -> states_complete



# Let's look for juveniles
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "Juvenile")) -> juv_movements

# Let's look at the adult histories; do they all start with BON (adult)? They should
states_complete %>% 
  subset(life_stage == "Adult") %>% 
  group_by(tag_code) %>% 
  mutate(order_2 = row_number()) %>% 
  subset(order_2 == 1) -> first_adult_states

table(first_adult_states$pathway)
# This looks good! They're all starting either with BON (adult) or BON (aborted)

##### Identify kelt movements #####

# First, we need to identify the end of spawning movements.

# Identify which sites are tributaries and which are mainstem sites
# so we never even use this in this script!
# We jsut say that everything that isn't mainstem is a tributary. Which makes sense. So this isn't affected by division into upstream and river mouth interrogation sites.
tributary_mainstem <- data.frame(tributary = c("Asotin Creek", 
                                               "Clearwater River",
                                               "Deschutes River", 
                                               "Entiat River", 
                                               "Fifteenmile Creek", 
                                               "Grande Ronde River", 
                                               "Hood River",
                                               "Imnaha River",
                                               "John Day River", 
                                               "Methow River", 
                                               "Okanogan River", 
                                               "Salmon River", 
                                               "Tucannon River", 
                                               "Umatilla River",
                                               "Walla Walla River",
                                               "Wenatchee River", 
                                               "Yakima River",
                                               "Upstream LGR other tributaries",
                                               "ICH to LGR other tributaries",
                                               "BON to MCN other tributaries",
                                               "Upstream WEL other tributaries"
),
mainstem = c("mainstem, upstream of LGR",
             "mainstem, upstream of LGR",
             "mainstem, BON to MCN",
             "mainstem, RRE to WEL",
             "mainstem, BON to MCN",
             "mainstem, upstream of LGR",
             "mainstem, BON to MCN",
             "mainstem, upstream of LGR",
             "mainstem, BON to MCN",
             "mainstem, upstream of WEL",
             "mainstem, upstream of WEL",
             "mainstem, upstream of LGR",
             "mainstem, ICH to LGR",
             "mainstem, BON to MCN",
             "mainstem, MCN to ICH or PRA",
             "mainstem, RIS to RRE",
             "mainstem, MCN to ICH or PRA",
             "mainstem, upstream of LGR",
             "mainstem, ICH to LGR",
             "mainstem, BON to MCN",
             "mainstem, upstream of WEL"))

# Get the order of sites
# Order of sites (no tributaries), Columbia River
site_order_notrib_columbia <- c("mainstem, mouth to BON", "mainstem, BON to MCN", 
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


# Site order, no tributaries, Snake directly to upper Columbia
site_order_notrib_columbia_snake <- c("mainstem, upstream of WEL",
                                      "mainstem, RRE to WEL", 
                                      "mainstem, RIS to RRE",
                                      "mainstem, PRA to RIS",
                                      "mainstem, MCN to ICH or PRA",
                                      "mainstem, ICH to LGR",
                                      "mainstem, upstream of LGR")


# Site order, both:
site_order_combined <- data.frame(state = c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                                "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                "mainstem, upstream of WEL", "mainstem, ICH to LGR",
                                "mainstem, upstream of LGR"), state_order = c(1,2,3,4,5,6,7,4,5))


# Final spawning movements:
# Seen at BON adult (repeat spawner); then all downstream movement before that and after last being seen is kelt movement

states_complete %>% 
  group_by(tag_code) %>% 
  mutate(order = row_number()) %>% 
  mutate(max_order = max(order)) -> states_complete

# Determine if movement is occurring upstream or downstream
states_complete %>% 
  left_join(., site_order_combined, by = "state") %>% 
  # now add a line so that every tributary gets a 99, which makes it so that any movement 
  # out of a tributary is considered downstream movement and any movement into a tributary is considered an upstream movement
  mutate(state_order = ifelse(is.na(state_order), 99, state_order)) -> states_complete

# Now classify movements as upstream or downstream
states_complete %>% 
  group_by(tag_code) %>% 
  mutate(direction = ifelse(lag(state_order) > state_order, "downstream", "upstream")) %>% 
  mutate(direction = ifelse(order == 1, "upstream", direction))-> states_complete

# Find repeat spawners based on a sequence of downstream movements preceding a BON detection (which indicates that it has now become a repeat spawner)
# states_complete %>% 
#   mutate(life_stage = ifelse(tag_code == lag(tag_code) & # if it's the same fish
#                                pathway == "BON (adult)" & # if it was seen in the adult ladder at BON
#                               time_after_arrival > days(x = 180) &  # and a long time after arrival
#                                 event_month >= 6 & event_month <= 11 & # and it was seen at BON in a month consistent with adult run timing, based on https://www.cbr.washington.edu/dart/wrapper?type=php&fname=hrt_adult_1657840462_690.php
#                                lag(direction) == "downstream" & # & it was just seen moving downstream (kelt movement)
#                                 tag_code == lead(tag_code) & # if the tag code continues, then make sure that
#                                lead(direction) != "upstream", # we don't subsequently see it moving back upstream
#                             "repeat", # then call it a repeat spawner
#                             ifelse(tag_code == lag(tag_code) & # if it's the same fish
#                                      pathway == "BON (adult)" & # if it was seen in the adult ladder
#                                      time_after_arrival > days(x = 180) &  # and a long time after arrival
#                                      event_month >= 6 & event_month <= 11 & # and it was seen at BON in a month consistent with adult run timing, based on https://www.cbr.washington.edu/dart/wrapper?type=php&fname=hrt_adult_1657840462_690.php
#                                      lag(direction) == "downstream" & # & it was just seen moving downstream
#                                      tag_code != lead(tag_code), # if the tag code doesn't continue, then we're chillin
#                             "repeat", life_stage))) -> states_complete  # then call it a repeat spawner

# In the above code block, we're making sure that we are not calling fish repeat spawners
# that subsequently move right back upstream after going downstream to BON. BUT - 
# this is also not picking up on repeat spawners because they are going to 
# move upstream too. So Maybe we take that part out, and just increase the time constraint

# checkpoint
# states_complete_orig <- states_complete
# states_complete <- states_complete_orig

states_complete %>% 
  mutate(life_stage = ifelse(tag_code == lag(tag_code) & # if it's the same fish
                               pathway == "BON (adult)" & # if it was seen in the adult ladder at BON
                               time_after_arrival > days(x = 180) &  # and a long time after arrival
                               event_month >= 6 & event_month <= 11 & # and it was seen at BON in a month consistent with adult run timing, based on https://www.cbr.washington.edu/dart/wrapper?type=php&fname=hrt_adult_1657840462_690.php
                               lag(direction) == "downstream", # & it was just seen moving downstream (kelt movement)
                             "repeat",  life_stage)) -> states_complete  # then call it a repeat spawner

# Now, call all downstream movement prior to a return to BON kelt movement

# Make a new field - order of direction, based on tag code and direction
# every time it changes, change direction index
states_complete$direction_index <- c(0,cumsum(as.numeric(with(states_complete,direction[1:(length(direction)-1)] != direction[2:length(direction)]))))

# Now change the downstream movements prior to being a repeat spawner to kelt movement
# First, find the direction indices preceding repeat
states_complete %>% 
  subset(life_stage == "repeat") %>% 
  mutate(repeat_start = order) %>% 
  distinct(tag_code, .keep_all = TRUE)-> repeat_spawners

states_complete %>% 
  subset(tag_code == "3D9.1BF1CF4AA3") -> test_bad

# Get the direction indices of when individuals are repeat spawners
repeat_spawners$direction_index -> repeat_direction_indices

pre_repeat_movement <- repeat_direction_indices - 1

# Now call the downstream movement prior to repeat spawning kelt movement
states_complete %>% 
  mutate(life_stage = ifelse(direction_index %in% pre_repeat_movement & direction == "downstream", "kelt", life_stage)) -> states_complete

# Also find kelt movement that does not precede a return
# if a fish does not return, it's really hard to determine where it spawned - you can't really say
# for sure when fallback is downstream kelt movement and when it's just trying to get to trib

# steelhead are known to spawn in the spring, so kelts would be moving in spring/summer back downstream

# Let's look first at downstream movements at kelt timing
# states_complete %>% 
#   mutate(life_stage = ifelse(event_month %in% c(5,6,7) & # select months we might expect kelts
#                                direction == "downstream" & # make sure that they're going downstream and (next line) make sure it's a long time after
#                                time_after_arrival > days(180), "kelt", life_stage)) -> states_complete

# so by month doesn't work, because most of the kelt movement is implicit movements
# that aren't actually detected (due to bad detection in downstream movements)
# instead, we'll want to choose based on just time after arrival
# states_complete %>% 
#   mutate(life_stage = ifelse(direction == "downstream" & # make sure that they're going downstream and (next line) make sure it's a long time after
#                                time_after_arrival > days(180), "kelt", life_stage)) -> states_complete

# not going to work, again because of implicit move timing.
# Let's try lag(direction)
# also only allow it if it's been at least 90 days and movement is occurring when we'd expect it to (following spring spawning)
# states_complete_orig <- states_complete

states_complete %>%
  mutate(life_stage = ifelse(order == max_order & event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & direction == "downstream", "kelt",
                             ifelse(event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & direction == "downstream" & lag(direction) == "downstream" | # if two consecutive downstream movements, then they're possible kelt movements
                               event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & direction == "downstream" & lead(direction) == "downstream", "kelt", 
                              life_stage))) -> states_complete


# TEMPORARY CALL: Any terminal downstream movement is kelt movement

states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(is.na(life_stage))) -> NA_life_stages

# This sort of works, but also IDs lots of movement, especially between upper snake and upper columbia or vice versa, as kelt.
# We'll fix up the time interpolation, then subset by only spring/late winter (probably around March to July or something like that)
# Okay, so the problem with subsetting only spring/late winter downstream movement
# is that again, with implicit site movements, we don't have that time.
# SO: I think we need to identify repeat spawners, then identify all consecutive 
# downstream movement prior to it arriving back at BON as kelt movement

# Fix the misidentified kelt movements in the middle of adult histories
# Get the direction indices of when individuals kelts
kelts <- subset(states_complete, life_stage == "kelt")

# Put in a checkpoint in case we need to make edits
# states_complete_orig <- states_complete
# states_complete <- states_complete_orig

# Here we are figuring out what the next direction is - this is important so that we can ID when downstream movement is NOT kelt movement
# e.g., when we have fallback events in the middle of an adult history
states_complete %>% 
  ungroup() %>% 
  mutate(direction_index = direction_index - 1) %>% 
  subset(direction_index >= 0) %>%  # remove the first one
  # mutate(index = row_number()) %>% 
  dplyr::rename(next_direction = direction) %>% 
  dplyr::select(next_direction, direction_index, tag_code) %>% 
  # dplyr::rename(next_direction_tag_code = tag_code) %>% 
  distinct(., next_direction, direction_index) -> direction_index_lead_df

states_complete %>% 
  ungroup() %>% 
  # mutate(index = row_number()) %>% 
  left_join(., direction_index_lead_df, by = c("direction_index")) -> states_complete

# Okay, I think we can use group_by(tag_code) and max(direction index) to see if it still changes direction for the same tag code
# BUT - we need to make sure it's not a repeat spawner. But then theoretically this would exclude fallback by repeat spawners
# So it definitely does exclude fallback by repeat spawners, which is an issue, but it looks like only for tag code 3D9.1C2D7EBE8B, which seems to be the only affected individual
states_complete %>% 
  group_by(tag_code) %>% 
  mutate(max_direction_index = max(direction_index)) -> states_complete

# First, create a new field for if at any point it's a repeat spawner
states_complete %>% 
  mutate(eventual_repeat = ifelse(tag_code %in% repeat_spawners$tag_code, "eventual_repeat", "single")) %>% 
  left_join(., dplyr::select(repeat_spawners, c(tag_code, repeat_start)), by = "tag_code")-> states_complete

# Now, change every movement after repeat start to repeat
states_complete %>% 
  mutate(life_stage = ifelse(eventual_repeat == "eventual_repeat" & order >= repeat_start, "repeat", life_stage)) -> states_complete

# Now change those kelt life stage events that are not actually kelts (because they're followed by upstream movement) back to Adult
# edit 2022-08-12 - make the fix for tag code 3D9.1C2D7EBE8B, which is the only repeat spawner with fallback that's not kelt movement
states_complete %>% 
  mutate(life_stage = ifelse(life_stage == "kelt" & next_direction == "upstream" & direction_index != max_direction_index & eventual_repeat == "single",
                             "Adult", life_stage)) -> states_complete

# edit 2022-08-12 - make the fix for tag code 3D9.1C2D7EBE8B, which is the only repeat spawner with fallback that's not kelt movement
states_complete %>% 
  mutate(life_stage = ifelse(tag_code == "3D9.1C2D7EBE8B" & event_year == 2012 & event_month == 4 & event_day == 11, "Adult", life_stage)) -> states_complete

# states_complete_ckpt <- states_complete



# Let's check on our repeat spawners
states_complete %>% 
  subset(., eventual_repeat == "eventual_repeat") -> eventual_repeat_df


# Let's check on our kelts
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "kelt")) %>% 
  dplyr::select(tag_code, state, date_time, life_stage, direction, pathway, order, next_direction, max_direction_index) -> kelts


# Check some individuals
# check on our guy
states_complete %>% 
  subset(tag_code == "384.36F2B442DB") -> test_kelt

# check how the triple spawners are doing
subset(eventual_repeat_df, tag_code == "3D9.1C2C430C8D") -> test_triple

# find possible repeats
kelts %>% 
  group_by(tag_code) %>% 
  slice(tail(row_number(),1)) %>% 
  subset(life_stage != "kelt") -> non_terminal_kelts

non_terminal_kelts$tag_code -> non_terminal_kelts_tags

non_terminal_kelts_hist <- subset(kelts, tag_code %in% non_terminal_kelts_tags)


# This guy (384.1B796A8EE1) is weird - one random detection in BON adult. Likely same as with juveniles, where it was seen there but going downstream?
# Need to make it so that all movement following kelt detection is kelt, until we see them again as a repeat spawner
# 384.36F2B32373 also odd - definitely looks like kelt movement, but it would mean that it spawned by March. I guess that's not impossible
# 384.36F2B442DB: This is not kelt movement - some downstream movement in April, but then back upstream to Clearwater. 
# Need to not call things kelt movement if they then go back upstream after
  # 	384.36F2B48AEC - same story as above - two downstream movements in the middle, but back upstream after
# So kelt movement either has to be terminal in the detection history, or has to be followed by repeat spawner movement.
# we can make this fix I think
# not a kelt - 384.3B2397826C. Same as before (downstream in the middle, switching from upper columbia to snake)
# 384.3B239A8BD6 - individual where we're not identifying some of the kelt movements. Would be fixed by making the change to call everything prior
# to repeat spawning as kelt movements
# mis IDed fallback as kelt movement :	384.3B239CF5F8

# IDENTIFY REPEAT KELT MOVEMENT

# This guy definitely has kelt movement after being a repeat spawner and returning
# 3D9.1BF18C45D3

# Use the same code again, but only for repeat spawners

states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(is.na(life_stage))) -> NA_life_stages

states_complete %>% 
  mutate(life_stage = ifelse(life_stage == "repeat" & event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & 
                               direction == "downstream" & order == max_order, "repeat_kelt",
                             ifelse(life_stage == "repeat" & event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & 
                               direction == "downstream" & lag(direction) == "downstream" |
                               life_stage == "repeat" & event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & 
                               direction == "downstream" & lead(direction) == "downstream", "repeat_kelt", life_stage))) -> states_complete
  
# states_complete_ckpt <- states_complete

# Quick fix for individuals (e.g., 3D9.1C2CCDA88C) with an upstream movement that's clearly a kelt (sandwiched between downstream movements)
# This is introducing a lot of NAs because if it's the last line in the history, then there's no lead, and if it's the first, there's no lag
# fixed!
states_complete %>% 
  mutate(life_stage = ifelse(eventual_repeat == "eventual_repeat" & order == max_order | eventual_repeat == "eventual_repeat" & order == 1, life_stage,
                             ifelse(eventual_repeat == "eventual_repeat" & lag(life_stage == "kelt") & lead(life_stage == "kelt"), "kelt", life_stage))) -> states_complete

states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(is.na(life_stage))) -> NA_life_stages

# Look into our repeat kelts
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "repeat_kelt")) -> repeat_kelts

# Looks good!

##### Truncate and export final history, without juveniles and kelts #####

# Look into NA life stages
states_complete %>% 
  filter(any(is.na(life_stage))) -> NA_life_stages

# Export to share with rebecca
NA_life_stages %>% 
  dplyr::select(-c(release_year, release_month, release_day, event_year, event_month, event_day, order, max_order, state_order, direction_index)) %>% 
  mutate(life_stage = ifelse(is.na(life_stage), "possible_kelt", life_stage)) %>% 
  mutate(time_after_arrival = as.numeric(time_after_arrival, units="days")) %>% 
  dplyr::rename(days_after_arrival = time_after_arrival, days_after_release = time_after_release) %>% 
  mutate(days_after_arrival = round(days_after_arrival, 1), days_after_release = round(days_after_release, 1)) ->possible_kelts_for_export

# write.csv(possible_kelts_for_export, "possible_kelts_for_export.csv")

# We want to model adult movement, but not kelt movement, and note when an adult has become a repeat spawner.
# For that, we will create a new field based on tag code, with repeat appended for those individuals

# Edit 2022-08-12 - sometimes if there are repeat repeats, then we need to append a number.
# This might be the only fish that's a repeat repeat: 3D9.1BF27BDC8D
# So we can see if it was a repeat kelt, and if it was, then the next repeat after it is another repeat

states_complete %>% 
  group_by(tag_code) %>% 
  subset(life_stage == "repeat_kelt") %>% 
  mutate(repeat_kelt_end = max(order)) %>% 
  distinct(tag_code, .keep_all = TRUE)-> repeat_kelts

# states_complete %>% 
#   subset(tag_code %in% repeat_kelts$tag_code) -> repeat_kelts_test

states_complete %>% 
  left_join(., dplyr::select(repeat_kelts, c(tag_code, repeat_kelt_end)), by = "tag_code") %>% 
  # need to make this edit so that we aren't left with NA tag_code_2
  mutate(repeat_kelt_end = ifelse(is.na(repeat_kelt_end), 999, repeat_kelt_end)) %>% 
  mutate(tag_code_2 = ifelse(eventual_repeat == "eventual_repeat" & order > repeat_kelt_end, paste0(tag_code, "_repeat_2"),
    ifelse(eventual_repeat == "eventual_repeat" & order >= repeat_start, paste0(tag_code, "_repeat"), tag_code))) -> states_complete

# Check to make sure our new tag codes make sense
states_complete %>% 
  subset(., eventual_repeat == "eventual_repeat") -> eventual_repeat_df_2


# FILTER FINAL FILE FOR EXPORT
states_complete %>% 
  subset(., !(life_stage %in% c("kelt", "repeat_kelt", "Juvenile"))) -> adults_only_states

# Let's make sure we're starting upstream of BON
adults_only_states %>% 
  group_by(tag_code_2) %>% 
  filter(row_number() == 1) -> adults_first_states

length(unique(adults_first_states$tag_code_2))

# Some start downsteram of BON because they're aborted ascent attempts. I think that's okay? Or we could just remove it.
# BUT in all honesty, I don't think that movements from the mouth to BON state are all that interesting to us
table(adults_first_states$state)
subset(adults_first_states, state == "mainstem, mouth to BON") -> BON_abort_starts
subset(adults_first_states, pathway == "BON_fallback_arrays")
table(subset(adults_first_states, state != "mainstem, BON to MCN")$time_after_arrival)
subset(adults_first_states, time_after_arrival < 0) -> early_fish
# 7 fish: 5 are clearly juveniles, 2 were seen at BHL first


subset(adults_first_states, state != "mainstem, BON to MCN")$tag_code -> problem_tags
table(subset(adults_first_states, state != "mainstem, BON to MCN")$pathway)
table(subset(adults_first_states, state != "mainstem, BON to MCN")$life_stage)

# Literally are all juvenile movements that are misIDed because time after arrival is negative (meaning we know they happened before arrival)

table(subset(adults_first_states, state != "mainstem, BON to MCN")$pathway)

# Look into the ones we last saw at fallback arrays
subset(adults_first_states, pathway == "BON_fallback_arrays")$tag_code ->BFA_tags

subset(adults_only_states, tag_code %in% BFA_tags) -> BFA_firststates_adult_only
subset(states_complete, tag_code %in% BFA_tags) -> BFA_firststates_complete

# they're all BO1 descents
# descent at BO1 as juvenile: 3D9.1BF22DF252
# This guy gets called a descent, when really it's an aborted attempt: 3D9.1BF2350747
# let's just manually change this, because it's so unusual - first seen in the vertical slot detectors, then makes it to the very top, then slides all the way back down and out
# Exact same story for this tag: 3D9.1BF24A4EFF
# descent as juvenile: 3D9.1C2CCD36F8
# descent as juvenile, then misIDed aborted: 3D9.1C2DF09C30
# exact same story with aborted BO1: 3D9.257C61E76B
# descent as juvenile: 3D9.257C634620 
# descent as juvenile: 3DD.007743DC75

# BO1 misIDed as descent that are actually aborted have been fixed.
# Descents as juveniles

subset(adults_only_states, tag_code %in% problem_tags) -> adult_only_issues
setdiff(adults_first_states$tag_code_2, adults_only_states$tag_code_2)

##### Remove ghost tags - use criteria from DART (provided by Susannah Iltis) #####

# DART applies these criteria for identifying ghost/shed tags for our tributary queries:
# (1) detected 3 or more years after last detection for all species tagged as adults
# (2) detected 5 years after last detection for all stages and species
# (3) steelhead tagged as juvenile with first detection 4 or more years after tagging
# (4) salmon tagged as juvenile with first detection 2 or more years after tagging
# (5) all tagged as adults with first detection 2 or more years after tagging.

# for our script, we've already removed the entire juvenile history at this point. So we can apply criteria (1)
# and remove everything that's occurred >= 3 years after the last detection

# load previously identified bad tags
bad_tags <- read.csv(here::here("bad_tags.csv"))

adults_only_states %>% 
  # subset(tag_code %in% bad_tags$tag_code) %>%  # for testing
  mutate(date_time = ymd_hms(date_time)) %>% 
  mutate(time_between_detections = ifelse(tag_code == lag(tag_code), time_after_release - lag(time_after_release), NA)) %>%
  # mutate(time_between_detections = ifelse(tag_code == lag(tag_code), interval(date_time, lag(date_time)), NA)) %>% 
  mutate(ghost_tag = ifelse(time_between_detections >= 365*3, "ghost", 
                            ifelse(is.na(time_between_detections) | time_between_detections < 365*3, "not_ghost", NA))) -> adults_only_states

# now note that everything that occurs after when we first note it's a ghost is also ghost activity
eventual_ghost_tags <- unique(subset(adults_only_states, ghost_tag == "ghost")$tag_code)
adults_only_states %>% 
  mutate(ghost_start = ifelse(ghost_tag == "ghost", order, NA)) %>% 
  subset(ghost_tag == "ghost") %>% 
  dplyr::select(tag_code, ghost_start) -> ghost_tag_start

adults_only_states %>% 
  ungroup() %>% 
  left_join(., ghost_tag_start, by = "tag_code") %>% 
  # group_by(tag_code) %>% 
  mutate(ghost_tag = ifelse(tag_code %in% eventual_ghost_tags & order >= ghost_start, "ghost", "not_ghost")) -> adults_only_states

# check on our ghost tags
adults_only_states %>% 
  group_by(tag_code) %>% 
  filter(any(ghost_tag == "ghost")) -> ghost_tags

# remove all of the ghost tag activity
adults_only_states %>% 
  subset(ghost_tag != "ghost") -> adults_only_states


##### Now, export this file, removing irrelevant columns #####
adults_only_states %>% 
  dplyr::select(tag_code, state, date_time, pathway, life_stage, tag_code_2) -> adults_only_states_for_export
write.csv(adults_only_states_for_export, here::here("stan_actual", "adults_states_complete.csv"))
# adults_only_states_for_export <- read.csv(here::here("stan_actual", "adults_states_complete.csv"))


# Export four more files, by each ESU
# Snake: everything above ICH
# Upper Columbia: Everything above PRA
# Middle Columbia: Everything above BON, except Hood River

# load natal origins
natal_origin_table <- read.csv(here::here("covariate_data", "natal_origin_table.csv"))

# match natal origin to tag codes
tag_code_metadata %>% 
  left_join(natal_origin_table, by = "release_site_name") -> tag_code_metadata

tag_code_metadata %>% 
  mutate(ESU = ifelse(natal_origin %in% c("Tucannon_River", "Asotin_Creek", "Clearwater_River", "Salmon_River", "Grande_Ronde_River", "Imnaha_River"), "snake",
                      ifelse(natal_origin %in% c("Wenatchee_River", "Entiat_River", "Okanogan_River","Methow_River"), "upper_columbia",
                             ifelse(natal_origin %in% c("Deschutes_River", "Fifteenmile_Creek", "John_Day_River", "Umatilla_River", "Yakima_River", "Walla_Walla_River"), "middle_columbia",
                                    ifelse(natal_origin %in% c("Hood_River"), "lower_columbia", "ERROR"))))) -> tag_code_metadata

snake_tags <- subset(tag_code_metadata, ESU == "snake")$tag_code
upper_columbia_tags <- subset(tag_code_metadata, ESU == "upper_columbia")$tag_code
middle_columbia_tags <- subset(tag_code_metadata, ESU == "middle_columbia")$tag_code
lower_columbia_tags <- subset(tag_code_metadata, ESU == "lower_columbia")$tag_code

length(unique(adults_only_states_for_export$tag_code_2))
length(unique(snake_adults_only_states$tag_code_2))

snake_adults_only_states

subset(adults_only_states_for_export, tag_code %in% snake_tags) -> snake_adults_only_states

# export with metadata
snake_adults_only_states %>% 
  left_join(., dplyr::select(tag_code_metadata, tag_code, natal_origin), by = "tag_code") -> snake_adults_only_states_metadata

write.csv(snake_adults_only_states_metadata, here::here("stan_actual", "ESU_models", "snake", "snake_adults_states_complete_metadata.csv"))

# check distribution of snake adults by origin
snake_adults_only_states %>% 
  distinct(tag_code_2, .keep_all = TRUE) %>% 
  left_join(., tag_code_metadata, by = "tag_code") -> snake_tag_code_2_distinct

  table(snake_tag_code_2_distinct$natal_origin)


subset(adults_only_states_for_export, tag_code %in% upper_columbia_tags) -> upper_columbia_adults_only_states
subset(adults_only_states_for_export, tag_code %in% middle_columbia_tags) -> middle_columbia_adults_only_states
subset(adults_only_states_for_export, tag_code %in% lower_columbia_tags) -> lower_columbia_adults_only_states
write.csv(snake_adults_only_states, here::here("stan_actual", "ESU_models", "snake", "snake_adults_states_complete.csv"))
write.csv(upper_columbia_adults_only_states, here::here("stan_actual", "ESU_models", "upper_columbia", "upper_columbia_adults_states_complete.csv"))
write.csv(middle_columbia_adults_only_states, here::here("stan_actual", "ESU_models", "middle_columbia", "middle_columbia_adults_states_complete.csv"))
write.csv(lower_columbia_adults_only_states, here::here("stan_actual", "ESU_models", "lower_columbia", "lower_columbia_adults_states_complete.csv"))


# write them to new folders
write.csv(snake_adults_only_states, here::here("stan_actual", "deteff_ESU_models", "snake", "snake_adults_states_complete.csv"))
write.csv(upper_columbia_adults_only_states, here::here("stan_actual", "deteff_ESU_models", "upper_columbia", "upper_columbia_adults_states_complete.csv"))
write.csv(middle_columbia_adults_only_states, here::here("stan_actual", "deteff_ESU_models", "middle_columbia", "middle_columbia_adults_states_complete.csv"))
write.csv(lower_columbia_adults_only_states, here::here("stan_actual", "deteff_ESU_models", "lower_columbia", "lower_columbia_adults_states_complete.csv"))



# Some summary statistics:

# How many repeat spawners do we have?
length(unique(adults_only_states$tag_code_2)) - length(unique(adults_only_states$tag_code))
# only 120 according to our code
# about 0.2% of adults returning to BON are seen at BON again as repeat spawners

# How many individuals with kelt movement do we have?
states_complete %>% 
  filter(any(life_stage %in% c("kelt", "repeat_kelts"))) -> kelts_final

# 1,379
length(unique(kelts_final$tag_code))
# About 2% of adults (returns to BON) exhibit kelt movement

repeat_kelts %>% 
  dplyr::select(-c("release_year", "release_month",  "release_day", "event_year",    "event_month",  "event_day")) -> repeat_kelts_abbrev

# how often do we see RM AFTER upstream?
adults_only_states %>% 
  mutate(trib_det_type = ifelse(grepl("Upstream", state), "upstream",
                                ifelse(grepl("Mouth", state), "river_mouth",
                                       ifelse(grepl("mainstem", state), "mainstem",
                                       "neither")))) -> adults_only_states_2

adults_only_states_2 %>% 
  mutate(upstream_seq = ifelse(tag_code == lead(tag_code) & trib_det_type == "upstream" & lead(trib_det_type) == "river_mouth", "upstream to river mouth",
                               ifelse(tag_code != lead(tag_code) & trib_det_type == "upstream", "upstream to loss",
                                      ifelse(tag_code != lead(tag_code) & trib_det_type == "river_mouth", "RM to loss",
                                      ifelse(tag_code == lead(tag_code) & trib_det_type == "upstream" & lead(trib_det_type) == "mainstem", "upstream to mainstem",
                                             ifelse(tag_code == lead(tag_code) & trib_det_type == "river_mouth" & lead(trib_det_type) == "upstream", "RM to upstream",
                                                    ifelse(tag_code == lead(tag_code) & trib_det_type == "river_mouth" & lead(trib_det_type) == "mainstem", "RM to mainstem",
                           NA))))))) -> adults_only_states_2
adults_only_states_2 %>% 
  group_by(tag_code) %>% 
  filter(any(!(is.na(upstream_seq)))) -> trib_upstream_det

table(trib_upstream_det$upstream_seq)

trib_upstream_det %>% 
  group_by(tag_code) %>% 
  filter(any(upstream_seq == "upstream to river mouth")) -> upstoRM

# So we never see river mouth detection after upstream... makes sense? This should be removed by our script to ID kelts

# let's check the script where we don't filter out by life stage (includes kelts)
states_complete %>% 
  mutate(trib_det_type = ifelse(grepl("Upstream", state), "upstream",
                                ifelse(grepl("Mouth", state), "river_mouth",
                                       "neither"))) -> states_complete_2

states_complete_2 %>% 
  ungroup() %>% 
  mutate(upstream_seq = ifelse(tag_code == lead(tag_code) & trib_det_type == "upstream" & lead(trib_det_type) == "river_mouth", "upstream to river mouth",
                               ifelse(tag_code != lead(tag_code) & trib_det_type == "upstream", "upstream to loss",
                                      NA))) -> states_complete_2
states_complete_2 %>% 
  group_by(tag_code) %>% 
  filter(any(!(is.na(upstream_seq)))) -> trib_upstream_det

table(trib_upstream_det$upstream_seq)
# also none here

trib_upstream_det %>% 
  group_by(tag_code) %>% 
  filter(any(upstream_seq == "upstream to river mouth")) -> upstoRM




