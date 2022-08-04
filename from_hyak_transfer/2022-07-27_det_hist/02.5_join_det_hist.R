# 02.5_join_det_hist

# Here we will be joining the different detection history files and joining them.

library(tidyverse)
library(lubridate)
library(janitor)
library(here)

# Load the files
det_hist1_4 <- read.csv(here::here("from_hyak_transfer", "2022-07-27_det_hist", "CTH1-4", "complete_det_hist_CTH1-4.csv"), row.names = 1)
det_hist5_8 <- read.csv(here::here("from_hyak_transfer", "2022-07-27_det_hist", "CTH5-8", "complete_det_hist_CTH5-8.csv"), row.names = 1)
det_hist9_11 <- read.csv(here::here("from_hyak_transfer", "2022-07-27_det_hist", "CTH9-11", "complete_det_hist_CTH9-11.csv"), row.names = 1)
det_hist12_14 <- read.csv(here::here("from_hyak_transfer", "2022-07-27_det_hist", "CTH12-14", "complete_det_hist_CTH12-14.csv"), row.names = 1)

# Fix the date/times (not sure why these got messed up, but of course they did)
det_hist1_4 %>% 
  mutate(start_time = ymd_hms(start_time)) %>% 
  mutate(end_time = ymd_hms(end_time)) -> det_hist1_4

det_hist5_8 %>% 
  mutate(start_time = ymd_hms(start_time)) %>% 
  mutate(end_time = ymd_hms(end_time)) -> det_hist5_8

det_hist9_11 %>% 
  mutate(start_time = ymd_hms(start_time)) %>% 
  mutate(end_time = ymd_hms(end_time)) -> det_hist9_11

det_hist12_14 %>% 
  mutate(start_time = ymd_hms(start_time)) %>% 
  mutate(end_time = ymd_hms(end_time)) -> det_hist12_14

# Join them:
det_hist1_4 %>% 
  bind_rows(., det_hist5_8) %>% 
  bind_rows(., det_hist9_11) %>% 
  bind_rows(., det_hist12_14) -> complete_det_hist

##### post-hoc fixes #####

# Fix 1: Fish that are seen in the exit antennas at BO2, but then skip BO4, are actually not aborting, but just missed detections at BO4
complete_det_hist %>% 
  mutate(non_ascent = ifelse(end_ant_group %in% c("BO2-WEIR 52", "BO2-WEIR 51", "BO2-UMT Entrance"), NA, non_ascent)) -> complete_det_hist

# Fix 2: Fish that are seen in the exit antennas at BO3, but then skip BO4, are actually not aborting, but just missed detections at BO4
# Use the exit coils, same as Susannah
complete_det_hist %>% 
  mutate(non_ascent = ifelse(end_ant_group %in% c("BO3-WEIR 59", "BO3-WEIR 58"), NA, non_ascent)) -> complete_det_hist

# Fix 3: Fish that are taking a very long time to pass through BO1
complete_det_hist %>% 
  mutate(to_remove = "no") -> complete_det_hist

# this for loop doesn't work yet

for (i in 1:(nrow(complete_det_hist)-1)){
  if (complete_det_hist[i, 'tag_code'] == complete_det_hist[i+1, 'tag_code'] & 
      complete_det_hist[i, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder" & 
      complete_det_hist[i+1, 'event_site_name'] == "BO1 - Bonneville Bradford Is. Ladder" &
      complete_det_hist[i, 'end_ant_group'] %in% c("A-BRANCH WEIR 50", "A-BRANCH WEIR 50",
                           "B-BRANCH WEIR 50", "B-BRANCH WEIR 50") &
      complete_det_hist[i+1, 'start_ant_group'] == "VERTICAL SLOT DETECTORS") {
    
    # Change non-ascent to not aborted
    complete_det_hist[i, "non_ascent"] <- NA
    
    # Move the end information from the subsequent row to the current row
    complete_det_hist[i, 'end_time'] <- complete_det_hist[i+1, 'end_time']
    complete_det_hist[i, 'end_antenna_id'] <- complete_det_hist[i+1, 'end_antenna_id']
    complete_det_hist[i, 'end_ant_group'] <- complete_det_hist[i+1, 'end_ant_group']
    
    # then, identify rows to remove - these are the rows following the current, from
    # which we have pulled the necessary information
    complete_det_hist[i+1, "to_remove"] <- "yes"
    
  }
}

# Now, remove those rows
complete_det_hist %>% 
  subset(to_remove == "no") -> complete_det_hist


# Export this file

write.csv(complete_det_hist, here::here("from_hyak_transfer", "2022-07-27_det_hist", "complete_det_hist.csv"))


# There are some clearly weird things happening here - these fixes are then incorporated above.

# Look into BO1 fish
complete_det_hist %>% 
  group_by(tag_code) %>% 
  filter(any(end_ant_group %in% c("A-BRANCH WEIR 50", "A-BRANCH WEIR 50",
                                  "B-BRANCH WEIR 50", "B-BRANCH WEIR 50"))) -> BO1_branch_ends

# so, majority of these fish are actually just hanging out in BO1 for a long time


# Look into detections at BO2-BO3-BO4 that end in BO2 or BO3. Based on the specific weirs,
# as well as where they were seen next, we can tell if these were really aborted ascension attempts or not

# These weirs are what Susannah calls the exit coils
complete_det_hist %>% 
  group_by(tag_code) %>% 
  filter(any(end_ant_group %in% c("BO3-WEIR 59", "BO3-WEIR 58"))) -> BO3_exit_ends

# These fish are almost all not aborting, just are being missed in BO4

# some are ambiguous, if they're next sen in BO4
BO3_exit_ends %>% 
  group_by(tag_code) %>% 
  filter(grepl("BO4",lead(end_ant_group)) & end_ant_group %in% c("BO3-WEIR 59", "BO3-WEIR 58") |
           grepl("BO4",end_ant_group) & lag(end_ant_group) %in% c("BO3-WEIR 59", "BO3-WEIR 58")) -> BO3_BO4s

BO3_BO4s %>% 
  group_by(tag_code) %>% 
  mutate(delay_time = ifelse(grepl("BO4", end_ant_group), start_time - lag(end_time), NA)) -> BO3_BO4_delay_time

ggplot(subset(BO3_BO4_delay_time, !(is.na(delay_time))), aes(x = delay_time)) +
  geom_histogram()

# I guess we just need to take the time component out of the BO2-BO3-BO4 complex entirely, since it seems like there are some fish that have 
# spent over 100 days in there

# We also need to allow fish to not be seen at BO4


complete_det_hist %>% 
  group_by(tag_code) %>% 
  filter(any(end_ant_group %in% c("BO2-WEIR 52", "BO2-WEIR 51", "BO2-UMT Entrance"))) -> BO2_aborts

# same story with BO2 really, sometimes they just take forever to get to BO4 (most of them), sometimes they just skip it (missed detections, I assume)

# check on PRA trap
complete_det_hist %>% 
  group_by(tag_code) %>% 
  filter(any(end_ant_group %in% c("LEFT [EAST] LADDER ADULT TRAP"))) -> PRA_trap_fish

# Yes, most of these fish are indeed clearly aborting an ascent attempt. However, some fish are continuing upstream. 
# It's impossible to say if they were just missed detections at the four antennas by the exit, or reascended and 
# missed the antennas again. It's a small number of  fish (single digits). the good thing is that it won't
# actually affect the number of ascents at PRA, there's just more uncertainty about the aborted attempts.


