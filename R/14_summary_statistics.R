# 14 - summary statistics

# This script will calculate sumamry statistics on overshoot, fallback, and other movements


# load libraries
library(here)
library(tidyverse)
library(lubridate)
library(janitor)

# Import data
# temporarily load old data to test script
states_complete <- read.csv(here::here("from_hyak_transfer", "2022-07-21-complete_det_hist", "states_complete.csv"), row.names = 1)

# load tag code metadata
tag_code_metadata <- read.csv(here::here("covariate_data", "tag_code_metadata.csv"))

# load natal origins
natal_origin_table <- read.csv(here::here("covariate_data", "natal_origin_table.csv"))

# match natal origin to tag codes
tag_code_metadata %>% 
  left_join(natal_origin_table, by = "release_site_name") -> tag_code_metadata



# Figure out what qualifies as overshoot for each natal origin

# What do we call Columbia fish showing up in the Snake or vice versa? This is sort of overshoot, but not exactly

# determine the first mainstem site that is overshoot for each natal origin
overshoot_mainstem_sites <- data.frame(natal_origin = sort(unique(tag_code_metadata$natal_origin)),
                                       overshoot_state = c(NA, # Asotin Creek
                                                           NA, # Clearwater
                                                           "mainstem, MCN to ICH or PRA", # Deschutes River
                                                           "mainstem, upstream of WEL", #Entiat
                                                           "mainstem, MCN to ICH or PRA", # Fifteenmile Creek
                                                           NA, # Grande Ronde River
                                                           "mainstem, MCN to ICH or PRA", # Hood River
                                                           NA, # Imnaha River
                                                           "mainstem, MCN to ICH or PRA", # John Day River
                                                           "mainstem, MCN to ICH or PRA",  # Klickitat River
                                                           NA, # Methow River
                                                           NA,  # Okanogan River
                                                           NA, # Salmon River
                                                           "mainstem, upstream of LGR", # Tucannon River
                                                           "mainstem, MCN to ICH or PRA", # Umatilla River
                                                           "mainstem, PRA to RIS", # Walla Walla River
                                                           "mainstem, RRE to WEL", # Wenatchee River
                                                           "mainstem, MCN to ICH or PRA",  # Wind River
                                                           "mainstem, PRA to RIS" # Yakima River
                                                           ))

# # add additional overshoot states for Yakima and Walla Walla, since they're at the branch
# overshoot_mainstem_sites %>% 
#   bind_rows(., data.frame(natal_origin = c("Yakima River", "Walla Walla River"),
#                           overshoot_state = c("mainstem, ICH to LGR", "mainstem, ICH to LGR"))) -> overshoot_mainstem_sites

# Now determine overshoot as being seen in any states upstream


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

states_complete %>% 
  left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin), by = "tag_code") %>% 
  left_join(., overshoot_mainstem_sites, by = "natal_origin") %>% 
  mutate(overshoot = ifelse(is.na(overshoot_state), "not_overshoot",
    ifelse(state == overshoot_state, "overshoot", "not_overshoot"))) -> states_complete_overshoot

# Let's count how many fish overshot at all
states_complete_overshoot %>% 
  group_by(tag_code) %>% 
  filter(any(overshoot == "overshoot")) -> overshooting_fish

# count how many fish per origin in our dataset

# first, get a count of how many fish we have per origin total
states_complete %>% 
  left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin), by = "tag_code") %>% 
  distinct(tag_code, .keep_all = TRUE) %>% 
  count(natal_origin) %>% 
  as.data.frame() %>% 
  dplyr::rename(total = n) -> natal_origin_fish_counts

# then, count how many overshooting fish
overshooting_fish %>% 
  ungroup() %>% 
  distinct(tag_code, .keep_all = TRUE) %>% 
  count(natal_origin) %>% 
  as.data.frame() %>% 
  dplyr::rename(n_overshot = n) -> natal_origin_overshooting_fish_counts


natal_origin_fish_counts %>% 
  left_join(., natal_origin_overshooting_fish_counts, by = "natal_origin") %>% 
  mutate(percent_overshoot = round(n_overshot/total,3)*100) -> overshoot_natal_origin_table
  
  
  
  
  
  

