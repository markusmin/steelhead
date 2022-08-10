# 14 - summary statistics

# This script will calculate sumamry statistics on overshoot, fallback, and other movements


# load libraries
library(here)
library(tidyverse)
library(lubridate)
library(janitor)

# Import data
# temporarily load old data to test script
# states_complete <- read.csv(here::here("from_hyak_transfer", "2022-07-21-complete_det_hist", "states_complete.csv"), row.names = 1)
states_complete <- read.csv(here::here("stan_actual", "adults_states_complete.csv"), row.names = 1)

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

# V2 - Add another column so that if a fish enters the other river OR overshoots by ascending dam upstream, then it's overshoot
# For yakima and Walla Walla, here we use the other river as the other_branch state
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
                                       ),
                                       other_branch_state = c("mainstem, PRA to RIS", # Asotin Creek
                                                           "mainstem, PRA to RIS", # Clearwater
                                                           NA, # Deschutes River
                                                           "mainstem, ICH to LGR", #Entiat
                                                           NA, # Fifteenmile Creek
                                                           "mainstem, PRA to RIS", # Grande Ronde River
                                                           NA, # Hood River
                                                           "mainstem, PRA to RIS", # Imnaha River
                                                           NA, # John Day River
                                                           NA,  # Klickitat River
                                                           "mainstem, ICH to LGR", # Methow River
                                                           "mainstem, ICH to LGR",  # Okanogan River
                                                           "mainstem, PRA to RIS", # Salmon River
                                                           "mainstem, PRA to RIS", # Tucannon River
                                                           NA, # Umatilla River
                                                           "mainstem, ICH to LGR", # Walla Walla River
                                                           "mainstem, ICH to LGR", # Wenatchee River
                                                           NA,  # Wind River
                                                           "mainstem, ICH to LGR" # Yakima River
                                       ))

# Need to add overshoot as going into the wrong river, as well

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
  # mutate(overshoot = ifelse(state %in% c(overshoot_state, other_branch_state), "overshoot", "not_overshoot")) -> states_complete_overshoot
  mutate(overshoot = ifelse(state == overshoot_state | state == other_branch_state, "overshoot", "not_overshoot")) %>% 
  mutate(overshoot = ifelse(is.na(overshoot), "not_overshoot", overshoot)) -> states_complete_overshoot

# Let's count how many fish overshot at all
states_complete_overshoot %>% 
  group_by(tag_code) %>% 
  filter(any(overshoot == "overshoot")) -> overshooting_fish

# count how many fish per origin in our dataset

# first, get a count of how many fish we have per origin total
states_complete %>% 
  left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin), by = "tag_code") %>% 
  distinct(tag_code, .keep_all = TRUE) %>% 
  ungroup() %>% 
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
  
# Count the number of fallback events, total
# first, create a df with the state order for mainstem sites

site_order_combined <- data.frame(state = c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                                            "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                            "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                            "mainstem, upstream of WEL", "mainstem, ICH to LGR",
                                            "mainstem, upstream of LGR"), state_order = c(1,2,3,4,5,6,7,4,5))
states_complete %>% 
  left_join(., site_order_combined, by = "state") %>% 
  mutate(fallback = ifelse(is.na(state_order), "tributary movement",
                           ifelse(tag_code == lag(tag_code) & state_order < lag(state_order),
                           "fallback", "ascent"))) -> fallback_movements

# count the total number of fallback movements
table(fallback_movements$fallback)

fallback_movements %>% 
  group_by(tag_code) %>% 
  count(fallback) %>% 
  subset(fallback == "fallback") -> fallback_counts_by_tag_code

table(fallback_counts_by_tag_code$n)

subset(fallback_counts_by_tag_code, n == 15)

subset(fallback_movements, tag_code == "3D9.1BF18CAC92") -> fallguy
# This is clearly a bug
# Actually it's not - this guy keeps switching ladders at Rock Island. bizarre

# 3D9.1BF2463473: did some crazy stuff around MCN
# Legitimately like 8 ascents + fallbacks at MCN

subset(fallback_counts_by_tag_code, n == 9)

subset(fallback_movements, tag_code == "3DD.003C02E3AA") -> fallguy2
# some of these are real, but there is some kelt movement in here as well that wasn't removed

# check on all of those that have >= 7 fallback movements
tag_codes_7ormore_fallbacks <- subset(fallback_counts_by_tag_code, n >= 7)$tag_code

subset(fallback_movements, tag_code %in% tag_codes_7ormore_fallbacks) -> fallguys

# dates for implicit states look off for this fish: 384.3B23AB763C
# So I think this is in the interpolating dates code. It's off by 8 hours, so I think it's a time zone issue.
# But it's also off by 7 hours in some cases

# Looks like a kelt: 384.3B23AB763C
# kelt: 3D9.1BF1896393
# kelt: 3D9.1BF207252D
# kelt: 3D9.1C2CB10926
# Kelt: 3D9.1C2D3D8D9D
# kelt: 3D9.1C2D756E9F
# kelt: 3D9.1C2D937BC6
# kelt: 3D9.257C62AE2D
# kelt: 3DD.003C02E3AA
# probably a kelt: 3DD.0077908190
# kelt: 3DD.00779FA9E8
# kelt: 3DD.0077A534EF
# kelt: 3DD.0077D68DE8
misIDedkelt_tag_codes <- c("384.3B23AB763C", "3D9.1BF1896393", "3D9.1BF207252D", "3D9.1C2CB10926", "3D9.1C2D3D8D9D")


# Distinguishing between post-overshoot fallback and delay fallback


  
  

