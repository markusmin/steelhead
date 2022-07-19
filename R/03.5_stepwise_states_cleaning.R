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
read.csv(here::here("from_hyak_transfer", "2022-07-16-complete_det_hist", "states_complete.csv")) %>%
  dplyr::select(-X) -> states_complete

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

notime_tags_times <- subset(complete_det_hist, tag_code %in% notime_tags & event_site_name == "Bonneville Adult Fishways (combined)" & is.na(end_time))

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
  if (states_complete[i,"tag_code"] %in% notime_tags & states_complete[i,"pathway"] == "BON (adult)"){
    states_complete[i, "date_time"] <- subset(notime_tags_times, tag_code == states_complete[i,"tag_code"])$start_time
  }
  # end temporary fix

  else if (is.na(states_complete[i,"date_time"])){
    # Find the next known time
    # Truncate the missing date times to only ones after current state
    missing_date_time_subset <- missing_date_time[(i+1):nrow(states_complete)]
    next_time_index <- min(which(missing_date_time_subset %in% FALSE)) + i

    # Now, interpolate the missing time
    next_time <- ymd_hms(date_times[next_time_index])

    # Get the missing time - add the time difference divided by the number
    # of missing steps plus 1, multiply by which number step it is
    missing_time <- next_time

    # populate the missing date_time
    states_complete[i, "date_time"] <- missing_time

  }
  # If it's not, leave it alone (do nothing)
  else {

  }
}



# Put this checkpoint in so we don't have to re-run the for loop
# write.csv(states_complete, here::here("from_hyak_transfer", "2022-07-16-complete_det_hist", "states_complete_times_interpolated.csv"))

read.csv(here::here("from_hyak_transfer", "2022-07-16-complete_det_hist", "states_complete_times_interpolated.csv")) %>%
  dplyr::select(-X) -> states_complete

# Fix the bad tags
states_complete %>% 
  mutate(date_time = ymd_hms(date_time)) -> states_complete

# Figure out how long detections were after release time
states_complete %>% 
  mutate(time_after_release = date_time - as.POSIXct(release_date)) -> states_complete

# Read in the BON arrival data file, to tell us when adult migration started
BON_arrival_df <- read.csv(here::here("covariate_data", "complete_BON_arrival.csv"))

states_complete %>% 
  left_join(., BON_arrival_df, by = "tag_code") -> states_complete

states_complete %>% 
  mutate(BON_arrival = ymd_hms(BON_arrival)) -> states_complete

# Fix this file because somehow it got turned to GMT (off by 7 hours)
states_complete %>% 
  mutate(BON_arrival = BON_arrival - hours(x = 7)) -> states_complete


states_complete %>% 
  mutate(time_after_arrival = date_time - as.POSIXct(BON_arrival)) -> states_complete




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
  "Adult")) -> states_complete

# One correction: if an implicit state follows a juvenile movement, it also has to be a juvenile movement. This way, we will remove
# the implicit mouth to BON state that we get when we see a juvenile in the adult ladder and return as an adult. That state
# visit in between shouldn't be modeled

# One correction: if mouth to BON follows a juvenile movement, it also has to be a juvenile movement. We want all fish starting when they reach BON (adult)

states_complete %>% 
  mutate(life_stage = ifelse(state == "mainstem, mouth to BON" & lag(life_stage) == "Juvenile", "Juvenile", life_stage)) -> states_complete

states_complete %>% 
  mutate(life_stage = ifelse(pathway == "implicit" & lag(life_stage) == "Juvenile", "Juvenile", life_stage)) -> states_complete



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

subset(first_adult_states, pathway == "MCN (adult)")
# so they're not - that's bad. Let's look at all of those that aren't
subset(first_adult_states, pathway != "BON (adult)")$tag_code -> bad_adults

bad_adults_states <- subset(states_complete, tag_code %in% bad_adults)

# I'm pretty sure these are all the result of the previous script - will need to double check once that one is finished

##### Identify kelt movements #####

# First, we need to identify the end of spawning movements.

# Identify which sites are tributaries and which are mainstem sites
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
  mutate(order = row_number()) -> states_complete

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
states_complete %>% 
  mutate(life_stage = ifelse(tag_code == lag(tag_code) & # if it's the same fish
                               pathway == "BON (adult)" & # if it was seen in the adult ladder
                              time_after_arrival > days(x = 180) &  # and a long time after arrival
                                event_month >= 6 & event_month <= 11 & # and it was seen at BON in a month consistent with adult run timing, based on https://www.cbr.washington.edu/dart/wrapper?type=php&fname=hrt_adult_1657840462_690.php
                               lag(direction) == "downstream" & # & it was just seen moving downstream
                                tag_code == lead(tag_code) & # if the tag code continues, then make sure that
                               lead(direction) != "upstream", # we don't subsequently see it moving back upstream
                            "repeat", # then call it a repeat spawner
                            ifelse(tag_code == lag(tag_code) & # if it's the same fish
                                     pathway == "BON (adult)" & # if it was seen in the adult ladder
                                     time_after_arrival > days(x = 180) &  # and a long time after arrival
                                     event_month >= 6 & event_month <= 11 & # and it was seen at BON in a month consistent with adult run timing, based on https://www.cbr.washington.edu/dart/wrapper?type=php&fname=hrt_adult_1657840462_690.php
                                     lag(direction) == "downstream" & # & it was just seen moving downstream
                                     tag_code != lead(tag_code), # if the tag code doesn't continue, then we're chillin
                            "repeat", life_stage))) -> states_complete  # then call it a repeat spawner

# Also find kelt movement that does not precede a return
# if a fish does not return, it's really hard to determine where it spawned - you can't really say
# for sure when fallback is downstream kelt movement and when it's just trying to get to trib

# steelhead are known to spawn in the spring, so kelts would be moving in spring/summer back downstream

# Let's look first at downstream movements at kelt timing
# states_complete %>% 
#   mutate(life_stage = ifelse(event_month %in% c(5,6,7) & # select months we might expect kelts
#                                direction == "downstream" & # make sure that they're going downstream and (next line) make sure it's a long time after
#                                time_after_arrival > days(180), "possible_kelt", life_stage)) -> states_complete

# so by month doesn't work, because most of the kelt movement is implicit movements
# that aren't actually detected (due to bad detection in downstream movements)
# instead, we'll want to choose based on just time after arrival
# states_complete %>% 
#   mutate(life_stage = ifelse(direction == "downstream" & # make sure that they're going downstream and (next line) make sure it's a long time after
#                                time_after_arrival > days(180), "possible_kelt", life_stage)) -> states_complete

# not going to work, again because of implicit move timing.
# Let's try lag(direction)
# also only allow it if it's been at least 90 days and movement is occurring when we'd expect it to (following spring spawning)
# states_complete_orig <- states_complete

states_complete %>% 
  mutate(life_stage = ifelse(event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & direction == "downstream" & lag(direction) == "downstream" | # if two consecutive downstream movements, then they're possible kelt movements
                               event_month >= 3 & event_month <= 7 & time_after_arrival > days(x = 90) & direction == "downstream" & lead(direction) == "downstream", "possible_kelt", life_stage)) -> states_complete

# This sort of works, but also IDs lots of movement, especially between upper snake and upper columbia or vice versa, as kelt.
# We'll fix up the time interpolation, then subset by only spring/late winter (probably around March to July or something like that)

# To fix some of the issues:
states_complete %>% 
  group_by(tag_code) %>% 
  mutate(max_order = max(order)) -> states_complete

# Check these
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "possible_kelt")) %>% 
  dplyr::select(tag_code, state, date_time, life_stage, direction, pathway, order, max_order) -> possible_kelts

# for (i in 1:nrow(states_complete)){
#   if(states_complete[i, "life_stage"] == "possible_kelt"){
#     n_lead_states <- states_complete[i, "max_order"] - states_complete[i, "order"]
#     for (j in 1:n_lead_states){
#       direction_vector = vector(length = n_lead_states)
#       direction_vector[j] <- lead(states_complete$direction[i:(i + n_lead_states)], n = j)
#     }
#   }
#   
#   
# }

for (i in 1:nrow(non_terminal_kelts_hist)){
  if(non_terminal_kelts_hist[i, "life_stage"] == "possible_kelt"){
    n_lead_states <- (non_terminal_kelts_hist[i, "max_order"] - non_terminal_kelts_hist[i, "order"])[1,1]
    for (j in 1:n_lead_states){
      direction_vector = vector(length = n_lead_states)
      direction_vector[j] <- lead(non_terminal_kelts_hist$direction[i:(i + n_lead_states)], n = j)
    }
  }
  
  
}


# find possible repeats
possible_kelts %>% 
  group_by(tag_code) %>% 
  slice(tail(row_number(),1)) %>% 
  subset(life_stage != "possible_kelt") -> non_terminal_kelts

non_terminal_kelts$tag_code -> non_terminal_kelts_tags

non_terminal_kelts_hist <- subset(possible_kelts, tag_code %in% non_terminal_kelts_tags)

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

# check how we did
states_complete %>% 
  group_by(tag_code) %>% 
  filter(any(life_stage == "repeat")) %>% 
  dplyr::select(tag_code, state, date_time, life_stage, direction, pathway)-> repeat_spawners



# Now, classify all downstream movement preceding a repeat spawner arriving at BON as kelt movement
life_stages_vec <- states_complete$life_stage


# for (i in 1:nrow(states_complete)){
#   if()
#   
# }





