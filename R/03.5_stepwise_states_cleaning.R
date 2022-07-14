# stepwise_states_cleaning


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



library(here)
library(tidyverse)
library(lubridate)
library(janitor)

# Read in states complete
read.csv(here::here("from_hyak_transfer", "2022-05-25-complete_det_hist", "states_complete.csv")) %>%
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
  mutate(life_stage = ifelse(time_after_release <= days(x = 90) | # fish that were seen right after release
  # mutate(life_stage = ifelse(
                             # or fish that were seen the same year as release and before June 15
                             release_year == event_year & event_month <= 5 | 
                               release_year == event_year & event_month == 6 & event_day <= 15 |
                            # or fish that were released the year before (on or after July 1) and seen the subsequent year on or before June 15
                               release_year == event_year-1 & release_month >= 7 & event_month <= 5 | 
                               release_year == event_year-1 & release_month >= 7 & event_month == 6 & event_day <= 15, "Juvenile", "Adult")) -> states_complete



# Let's look for juveniles
juv_movements <- subset(states_complete, life_stage == "Juvenile")

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







