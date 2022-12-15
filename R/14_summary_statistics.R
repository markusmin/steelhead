# 14 - summary statistics

# This script will calculate sumamry statistics on overshoot, fallback, and other movements


# load libraries
library(here)
library(tidyverse)
library(lubridate)
library(janitor)


middle_columbia_origins <- c("Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River", "Walla Walla River", "Yakima River")
upper_columbia_origins <- c("Wenatchee River", "Entiat River",    "Okanogan River",  "Methow River")
snake_origins <- c("Tucannon River",   "Asotin Creek", "Clearwater River", "Salmon River",  "Grande Ronde River", "Imnaha River")

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

##### Fish x run year x origin #####

run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2022-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2023-05-31 23:59:59"), by = "years")
run_year_numeric = seq(4, 22, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)

# First, let's join the run year df with the states complete 
# check how long this takes:
# 15 minutes for the full dataset
print(Sys.time())
states_complete %>% 
  rowwise() %>% # this is apparently crucial to getting this to work
  dplyr::mutate(date_time = ymd_hms(date_time)) %>% 
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_time & run_year_end >= date_time)$run_year) -> states_complete
print(Sys.time())

states_complete %>% 
  left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin), by = "tag_code") %>% 
  distinct(tag_code_2, .keep_all = TRUE) %>% 
  ungroup() %>% 
  count(natal_origin, run_year) %>% 
  as.data.frame() %>% 
  dplyr::rename(total = n) -> natal_origin_fish_counts

### CREATE THE TABLE FOR REPORT - rows = natal origins, columns = natal origins, values = number of fish
natal_origin_fish_counts %>% 
  pivot_wider(., values_from = total, names_from = run_year) %>% 
  relocate("05/06", .before = "06/07") %>% 
  replace(is.na(.), 0) -> fish_year_origin_table

write.csv(fish_year_origin_table, here::here("CBR_report", "CBR_report_final", "tables", "fish_year_origin_table.csv"))
  

##### Overshoot summary statistics #####

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
                                       natal_origin_region = c("Snake", # Asotin Creek
                                                               "Snake", # Clearwater
                                                               "Middle Columbia", # Deschutes River
                                                               "mainstem, upstream of WEL", #Entiat
                                                               "Middle Columbia", # Fifteenmile Creek
                                                               "Snake", # Grande Ronde River
                                                               "Lower Columbia", # Hood River
                                                               "Snake", # Imnaha River
                                                               "Middle Columbia", # John Day River
                                                               "Middle Columbia",  # Klickitat River
                                                               "Upper Columbia", # Methow River
                                                               "Upper Columbia",  # Okanogan River
                                                               "Snake", # Salmon River
                                                               "Snake", # Tucannon River
                                                               "Middle Columbia", # Umatilla River
                                                               "Middle Columbia", # Walla Walla River
                                                               "Upper Columbia", # Wenatchee River
                                                               "Lower Columbia",  # Wind River
                                                               "Middle Columbia" # Yakima River
                                       ),
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

# Now, determine if the fish was in a post-overshoot mainstem site at that time

site_order_combined <- data.frame(state = c("mainstem, mouth to BON", "mainstem, BON to MCN", 
                                            "mainstem, MCN to ICH or PRA", "mainstem, PRA to RIS",
                                            "mainstem, RIS to RRE", "mainstem, RRE to WEL", 
                                            "mainstem, upstream of WEL", "mainstem, ICH to LGR",
                                            "mainstem, upstream of LGR"), state_order = c(1,2,3,4,5,6,7,4,5),
                                  state_region = c("Lower Columbia", "Lower Columbia", 
                                                   "Middle Columbia", "Upper Columbia",
                                                   "Upper Columbia", "Upper Columbia",
                                                   "Upper Columbia", "Snake",
                                                   "Snake"))

site_order_combined %>% 
  dplyr::rename(overshoot_state = state, overshoot_state_index = state_order) %>% 
  dplyr::select(-state_region)-> site_order_combined_overshoot

site_order_combined %>% 
  dplyr::rename(other_branch_state = state, other_branch_state_index = state_order) %>% 
  dplyr::select(-state_region) -> site_order_combined_other_branch



states_complete_overshoot %>% 
  left_join(., site_order_combined_overshoot, by = "overshoot_state") %>% 
  left_join(., site_order_combined_other_branch, by = "other_branch_state") -> states_complete_overshoot

# Now, figure out if the fish is currently in an overshoot state
states_complete_overshoot %>% 
  left_join(., site_order_combined, by = "state") %>% 
  mutate(overshoot_status = ifelse(natal_origin_region == "Upper Columbia" & state_region == "Snake", "overshoot_state",
                                          ifelse(natal_origin_region == "Snake" & state_region == "Upper Columbia", "overshoot_state",
                                                 ifelse(state_order >= overshoot_state_index, "overshoot_state",
                                                 "en_route_state")))) %>% 
  mutate(overshoot_status = ifelse(is.na(overshoot_status), "en_route_state", overshoot_status)) -> states_complete_overshoot

table(states_complete_overshoot$overshoot_status)



# Add a new column for directionality - this way we can count each overshoot movement
states_complete_overshoot %>% 
  group_by(tag_code) %>% 
  mutate(direction = ifelse(lag(state_order) > state_order, "downstream", "upstream")) -> states_complete_overshoot
 
# Let's count how many fish overshot at all
states_complete_overshoot %>% 
  group_by(tag_code) %>% 
  filter(any(overshoot == "overshoot")) -> overshooting_fish

# count how many fish per origin x run_year in our dataset

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

sum(overshoot_natal_origin_table$n_overshot)/sum(overshoot_natal_origin_table$total)

# calculate a standard error around this
sqrt(sum(overshoot_natal_origin_table$n_overshot)/sum(overshoot_natal_origin_table$total) * (1 - sum(overshoot_natal_origin_table$n_overshot)/sum(overshoot_natal_origin_table$total))/sum(overshoot_natal_origin_table$total))

# How many individuals overshot at least one dam?
length(unique(overshooting_fish$tag_code_2))
length(unique(states_complete$tag_code_2))

### CREATE THE TABLE FOR REPORT - columns = overshoot at each dam + total (fish that overshoot at any), rows are natal origins
overshooting_fish %>% 
  ungroup() %>% 
  subset(overshoot_status == "overshoot_state" & direction == "upstream") %>%
  # change all upstream movements that don't have the dam in the pathway to indicate the dam
  mutate(pathway = ifelse(pathway == "implicit" & state %in% c("mainstem, upstream of WEL", "mainstem, upstream of LGR"), gsub(".*\\s", "", state),
                          ifelse(pathway == "ICH_LGR_inriver", "ICH",
                                 ifelse(pathway == "implicit", gsub(" .*", "", gsub("mainstem, ", "", state)), pathway)))) %>% 
  mutate(pathway = gsub(" \\(adult\\)", "", pathway)) %>% 
  # Make sure you don't double count - only keep unique combinations of tag code + overshoot dam (pathway)
  # This eliminates about 12.5% of the overshoot movements. But I think this is what we want - proportion of fish, not absolute number of overshoot events
  filter(!(duplicated(paste0(tag_code_2, pathway)))) %>% 
  count(natal_origin, pathway) %>% 
  # add the total fish per origin, so we can calculate proportions
  left_join(natal_origin_fish_counts, by = "natal_origin") %>% 
  # mutate(value = paste0(formatC(n/total, digits = 2, format = "fg"), " (N = ", n, ")")) %>% 
  mutate(value = paste0(formatC(n/total * 100, digits = 2, format = "f"), "% (", n, "/", total, ")")) %>% 
  dplyr::select(-c(n, total)) %>% 
  pivot_wider(names_from = pathway, values_from = value) %>% 
  replace(is.na(.), "0") %>% 
  # reorder the dams
  relocate(MCN, .before = PRA)-> overshoot_table_for_report

# Replace 0 with NA for dams that are between the natal tributary and the ocean (since you can't overshoot those)
overshoot_table_for_report[1,c('MCN', 'ICH', 'LGR')] <- "NA" # Asotin Creek
overshoot_table_for_report[2,c('MCN', 'ICH', 'LGR')] <- "NA" # Clearwater River
# overshoot_table_for_report[3,c('MCN', 'ICH', 'LGR')] <- "NA" # Deschutes River
overshoot_table_for_report[4,c('MCN', 'PRA', 'RIS', 'RRE')] <- "NA" # Entiat River
# overshoot_table_for_report[5,c('MCN', 'ICH', 'LGR')] <- "NA" # Fifteenmile Creek
overshoot_table_for_report[6,c('MCN', 'ICH', 'LGR')] <- "NA" # Grande Ronde River
# overshoot_table_for_report[7,c('MCN', 'ICH', 'LGR')] <- "NA" # Hood River
overshoot_table_for_report[8,c('MCN', 'ICH', 'LGR')] <- "NA" # Imnaha River
# overshoot_table_for_report[9,c('MCN', 'ICH', 'LGR')] <- "NA" # John Day River
overshoot_table_for_report[10,c('MCN', 'PRA', 'RIS', 'RRE', 'WEL')] <- "NA" # Methow River
overshoot_table_for_report[11,c('MCN', 'PRA', 'RIS', 'RRE', 'WEL')] <- "NA" # Okanogan River
overshoot_table_for_report[12,c('MCN', 'ICH', 'LGR')] <- "NA" # Salmon River
overshoot_table_for_report[13,c('MCN', 'ICH')] <- "NA" # Tucannon River
# overshoot_table_for_report[14,c('MCN', 'ICH', 'LGR')] <- "NA" # Umatilla River
overshoot_table_for_report[15,c('MCN')] <- "NA" # Walla Walla River
overshoot_table_for_report[16,c('MCN', 'PRA', 'RIS')] <- "NA" # Wenatchee River
overshoot_table_for_report[17,c('MCN')] <- "NA" # Yakima River


write.csv(overshoot_table_for_report, here::here("CBR_report", "CBR_report_final", "tables", "overshoot_dam_origin_table.csv"))

##### How often did fish visit an ESU that they don't belong in? #####
# Snake in the upper Columbia
states_complete_overshoot %>% 
  subset(., natal_origin_region == "Snake") %>% 
  group_by(tag_code) %>% 
  filter(any(state_region == "Upper Columbia")) -> snake_origin_upper_columbia

length(unique(snake_origin_upper_columbia$tag_code_2))/length(unique(subset(states_complete_overshoot, natal_origin_region == "Snake")$tag_code_2))
460/33761

# Upper columbia in the Snake
states_complete_overshoot %>% 
  subset(., natal_origin_region == "Upper Columbia") %>% 
  group_by(tag_code) %>% 
  filter(any(state_region == "Snake")) -> upper_columbia_origin_snake

length(unique(upper_columbia_origin_snake$tag_code_2))
length(unique(subset(states_complete_overshoot, natal_origin_region == "Upper Columbia")$tag_code_2))
7/14093

# Middle Columbia in the Snake or upper Columbia
states_complete_overshoot %>% 
  subset(., natal_origin_region == "Middle Columbia") %>% 
  group_by(tag_code) %>% 
  filter(any(state_region %in% c("Snake", "Upper Columbia"))) -> middle_columbia_origin_snake_upper_columbia

length(unique(middle_columbia_origin_snake_upper_columbia$tag_code_2))
length(unique(subset(states_complete_overshoot, natal_origin_region == "Middle Columbia")$tag_code_2))
2163/9185

### Table for report - rows = natal origins, columns = how many individuals were observed outside of the ESU
snake_origin_upper_columbia %>% 
  ungroup() %>% 
  filter(!duplicated(tag_code_2)) %>% 
  count(natal_origin) -> snake_origin_upper_columbia_counts

upper_columbia_origin_snake %>% 
  ungroup() %>% 
  filter(!duplicated(tag_code_2)) %>% 
  count(natal_origin) -> upper_columbia_origin_snake_counts

middle_columbia_origin_snake_upper_columbia %>% 
  ungroup() %>% 
  filter(!duplicated(tag_code_2)) %>% 
  count(natal_origin) -> middle_columbia_origin_snake_upper_columbia_counts

middle_columbia_origin_snake_upper_columbia_counts %>% 
  bind_rows(., upper_columbia_origin_snake_counts) %>% 
  bind_rows(., snake_origin_upper_columbia_counts) -> out_of_ESU_counts

natal_origin_fish_counts %>% 
  group_by(natal_origin) %>% 
  summarise(total = sum(total)) -> natal_origin_total_fish_counts


out_of_ESU_counts %>% 
  left_join(., natal_origin_total_fish_counts, by = "natal_origin") %>% 
  mutate(percent = paste0(formatC(n/total * 100, digits = 2, format = "f"),"%")) %>% 
  relocate(total, .before = n) %>% 
  mutate(DPS = ifelse(natal_origin %in% gsub(" ", "_", middle_columbia_origins), "Middle Columbia",
                      ifelse(natal_origin %in% gsub(" ", "_", upper_columbia_origins), "Upper Columbia",
                             ifelse(natal_origin %in% gsub(" ", "_", snake_origins), "Snake River Basin", "ERROR")))) %>% 
  relocate(DPS, .before = total) -> out_of_ESU_counts_table

## Add in the missing ones
# no_out_of_ESU <- data.frame(natal_origin = c("Hood_River", "Entiat_River"), DPS = c("Lower Columbia", "Upper Columbia"), total = c(0,0), n = c(0,0), percent = c("0%","0%"))

out_of_ESU_counts_table %>% 
  add_row(natal_origin = c("Hood_River"), DPS = c("Lower Columbia"), total = c(subset(natal_origin_total_fish_counts, natal_origin == "Hood_River")$total), n = c(0), percent = c("0%"), .before = 1) %>% 
  add_row(natal_origin = c("Entiat_River"), DPS = c("Upper Columbia"), total = c(subset(natal_origin_total_fish_counts, natal_origin == "Entiat_River")$total), n = c(0), percent = c("0%"), .before = 8) -> out_of_ESU_counts_table

write.csv(out_of_ESU_counts_table, here::here("CBR_report", "CBR_report_final", "tables", "out_of_ESU_counts_table.csv"))

# Lower Columbia in the Snake or middle Colubmia or upper Columbia
states_complete_overshoot %>% 
  subset(., natal_origin_region == "Lower Columbia") %>% 
  group_by(tag_code) %>% 
  filter(any(state_region %in% c("Snake", "Upper Columbia", "Middle Columbia"))) -> lower_columbia_origin_snake_upper_columbia_middle_columbia

length(unique(lower_columbia_origin_snake_upper_columbia_middle_columbia$tag_code_2))
length(unique(subset(states_complete_overshoot, natal_origin_region == "Lower Columbia")$tag_code_2))
3/3007

##### Fallback summary statistics #####

# Count the number of fallback events, total
# first, create a df with the state order for mainstem sites

# states_complete %>% 
states_complete_overshoot %>% 
  # left_join(., site_order_combined, by = "state") %>% 
  mutate(fallback = ifelse(is.na(state_order), "tributary movement",
                           ifelse(tag_code_2 == lag(tag_code_2) & state_order < lag(state_order),
                           "fallback", "ascent"))) -> fallback_movements

# Determine at which dam fallback events are occurring
fallback_movements %>% 
  mutate(fallback_dam = ifelse(fallback != "fallback", "not fallback",
                               ifelse(state == "mainstem, mouth to BON" & lag(state) == "mainstem, BON to MCN", "BON",
                                      ifelse(state == "mainstem, BON to MCN" & lag(state) == "mainstem, MCN to ICH or PRA", "MCN",
                                             ifelse(state == "mainstem, MCN to ICH or PRA" & lag(state) =="mainstem, PRA to RIS", "PRA",
                                                    ifelse(state == "mainstem, PRA to RIS" & lag(state) == "mainstem, RIS to RRE", "RIS",
                                                           ifelse(state == "mainstem, RIS to RRE" & lag(state) == "mainstem, RRE to WEL", "RRE",
                                                                  ifelse(state == "mainstem, RRE to WEL" & lag(state) == "mainstem, upstream of WEL", "WEL",
                                                                         ifelse(state == "mainstem, MCN to ICH or PRA" & lag(state) == "mainstem, ICH to LGR", "ICH",
                                                                                ifelse(state == "mainstem, ICH to LGR" & lag(state) == "mainstem, upstream of LGR", "LGR", "error")))))))))) -> fallback_movements

# To provide context for interpreting results, also count the number of ascensions at each dam
fallback_movements %>% 
  mutate(pathway_detailed = ifelse(fallback == "fallback", paste0("fallback_", fallback_dam),
                               ifelse(lag(state) == "mainstem, mouth to BON" & state == "mainstem, BON to MCN", "ascent_BON",
                                      ifelse(lag(state) == "mainstem, BON to MCN" & state == "mainstem, MCN to ICH or PRA", "ascent_MCN",
                                             ifelse(lag(state) == "mainstem, MCN to ICH or PRA" & state =="mainstem, PRA to RIS", "ascent_PRA",
                                                    ifelse(lag(state) == "mainstem, PRA to RIS" & state == "mainstem, RIS to RRE", "ascent_RIS",
                                                           ifelse(lag(state) == "mainstem, RIS to RRE" & state == "mainstem, RRE to WEL", "ascent_RRE",
                                                                  ifelse(lag(state) == "mainstem, RRE to WEL" & state == "mainstem, upstream of WEL", "ascent_WEL",
                                                                         ifelse(lag(state) == "mainstem, MCN to ICH or PRA" & state == "mainstem, ICH to LGR", "ascent_ICH",
                                                                                ifelse(lag(state) == "mainstem, ICH to LGR" & state == "mainstem, upstream of LGR", "ascent_LGR", pathway)))))))))) -> fallback_movements

# get only the fallback movements themselves
fallback_movements %>% 
  subset(fallback == "fallback") -> fallback_only_movements

### CREATE THE TABLE FOR REPORT - columns = fallback at each dam + total (fish that overshoot at any), rows are natal origins
fallback_only_movements %>% 
  ungroup() %>% 
  # Make sure you don't double count - only keep unique combinations of tag code + fallback dam (pathway)
  filter(!(duplicated(paste0(tag_code_2, fallback_dam)))) %>% 
  count(natal_origin, fallback_dam) %>% 
  # add the total fish per origin, so we can calculate proportions
  left_join(natal_origin_fish_counts, by = "natal_origin") %>% 
  # mutate(value = paste0(formatC(n/total, digits = 2, format = "fg"), " (N = ", n, ")")) %>% 
  mutate(value = paste0(formatC(n/total * 100, digits = 2, format = "f"), "% (", n, "/", total, ")")) %>% 
  dplyr::select(-c(n, total)) %>% 
  pivot_wider(names_from = fallback_dam, values_from = value) %>% 
  replace(is.na(.), "0") %>% 
  # reorder the dams
  relocate(ICH, .after = WEL) %>% 
  relocate(LGR, .after = ICH)-> fallback_table_for_report

# Change 0 to NA if there is no opportunity for FB (i.e., no fish passed the dam in the first place)
fallback_table_for_report[3,c('PRA', 'RIS', 'RRE', 'WEL')] <- "NA" # Deschutes River
fallback_table_for_report[4,c('ICH', 'LGR')] <- "NA" # Entiat River
fallback_table_for_report[5,c('PRA', 'RIS', 'RRE', 'WEL')] <- "NA" # Fifteenmile Creek
fallback_table_for_report[7,c('PRA', 'RIS', 'RRE', 'WEL', 'ICH', 'LGR')] <- "NA" # Hood River


write.csv(fallback_table_for_report, here::here("CBR_report", "CBR_report_final", "tables", "fallback_dam_origin_table.csv"))

fallback_only_movements %>% 
  ungroup() %>% 
  # Make sure you don't double count - only keep unique combinations of tag code + fallback dam (pathway)
  filter(!(duplicated(paste0(tag_code_2, fallback_dam)))) %>% 
  count(natal_origin, fallback_dam) %>% 
  dplyr::select(-natal_origin) %>% 
  group_by(fallback_dam) %>% 
  summarise(total = sum(n)) -> total_fallback_by_dam


# 
# fallback_only_movements %>% 
#   left_join(dplyr::select(tag_code_metadata, tag_code, natal_origin), by = "tag_code") -> fallback_only_movements

# count the total number of fallback movements by natal origin
table(fallback_only_movements$natal_origin)

# count the total number of fallback movements by dam
table(fallback_only_movements$fallback_dam)

# compare to number of ascents by dam
table(fallback_movements$pathway_detailed) %>% 
  as.data.frame() %>% 
  dplyr::rename(pathway_detailed = Var1) -> pathway_detailed_DF

pathway_detailed_DF %>% 
  subset(grepl("ascent|fallback", pathway_detailed)) %>% 
  mutate(dam = str_replace(pathway_detailed, "ascent_|fallback_", "")) %>% 
  mutate(pathway_detailed = gsub("_.*", "", pathway_detailed)) %>% 
  pivot_wider(., names_from = dam, values_from = Freq) %>%
  column_to_rownames("pathway_detailed") %>% 
  t() -> ascent_fallback_dam_counts
  # pivot_longer(., cols = c("BON", "ICH", "LGR", "MCN", "PRA", "RIS", "RRE", "WEL"), names_to = "dam")

# Change the BON ascent counts to include all fish 
ascent_fallback_dam_counts[1,1] <- ascent_fallback_dam_counts[1,1] + length(unique(states_complete$tag_code_2))

# Calculate the proportion of fallback
ascent_fallback_dam_counts %>% 
  as.data.frame() %>% 
  mutate(percent_fallback = round(fallback/ascent * 100, 2)) -> ascent_fallback_dam_counts


# Distinguish post-overshoot fallback from delay (en route) fallback
fallback_movements %>% 
  # need to do lag overshoot_status, so we can see what the previous state was
  mutate(fallback_type = ifelse(fallback == "fallback" & lag(overshoot_status) == "overshoot_state", "post-overshoot fallback", 
                                ifelse(fallback == "fallback" & lag(overshoot_status) == "en_route_state", "en-route fallback",
                                       "NA"))) -> fallback_movements

table(fallback_movements$fallback_type)

# By only looking at fallback to natal tributary, we are ignore about half of all fallback that we can detect

# Inspect these
fallback_movements %>% 
  group_by(tag_code) %>% 
  dplyr::select("tag_code", "state", "date_time", "pathway", "tag_code_2", "natal_origin", "fallback", "fallback_type") %>% 
  filter(any(fallback == "fallback")) -> fallback_det_hist

fallback_movements %>% 
  group_by(tag_code) %>% 
  # dplyr::select("tag_code", "state", "date_time", "pathway", "tag_code_2", "natal_origin", "fallback", "fallback_type") %>% 
  filter(fallback == "fallback") -> fallback_movements_only

# check some individual natal origins
table(subset(fallback_movements_only, natal_origin == "John_Day_River")$fallback_type)

# check out all natal origins
fallback_movements_only %>% 
  group_by(natal_origin) %>% 
  count(fallback_type) %>% 
  as.data.frame() %>% 
  pivot_wider(., names_from = fallback_type, values_from = n) -> fallback_types_by_origin

# count how many unique fish fell back at least once
length(unique(fallback_movements_only$tag_code_2))
length(unique(states_complete$tag_code_2))



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



##### Final fates/success rate, by overshoot and fallback #####




  
  

