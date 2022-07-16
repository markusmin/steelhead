### 13 inspect distributions of transitions/data summaries

# Description: This R script summarizes the transitions that we see in our dataset. It allows
# us to 1) make sure the implicit site visit code is working, and 2) see how many transitions
# we have between each pair of states (and with which combination of covariates) to see
# if we have enough data to estimate parameters of interest


##### Load data, libraries #####
library(tidyverse)
library(lubridate)
library(janitor)
library(here)

# First, let's check the detection history script to make sure that the input file for states complete is good
# check the old file
read.csv(here::here("from_hyak_transfer", "2022-05-24_det_hist", "complete_det_hist.csv")) %>% 
  dplyr::select(-X) -> det_hist

det_hist %>% 
  group_by(tag_code) %>% 
  mutate(order = row_number()) %>% 
  subset(order == 1) -> first_states
# Yeah, so this one has a couple dozen errors where fish are starting way upstream


read.csv(here::here("from_hyak_transfer", "2022-07-16_det_hist", "complete_det_hist.csv")) %>% 
  dplyr::select(-X) -> det_hist

det_hist %>% 
  group_by(tag_code) %>% 
  mutate(order = row_number()) %>% 
  subset(order == 1) -> first_states
# This one looks good - bless

# Read in states complete
# read.csv(here::here("from_hyak_transfer", "2022-05-25-complete_det_hist", "states_complete.csv")) %>%
# read.csv(here::here("from_hyak_transfer", "2022-06-10-complete_det_hist", "states_complete.csv")) %>%

# check out the script is doing
read.csv(here::here("from_hyak_transfer", "2022-07-14-complete_det_hist", "states_complete.csv")) %>%
  dplyr::select(-X) -> states_complete

# Get rid of fake fish
states_complete <- subset(states_complete, !(tag_code == "dummy_fish"))

# Keep only individuals from the natal origins we're interested in

# Read in natal origin data
natal_origin_table <- read.csv(here::here("covariate_data", "natal_origin_table.csv"))
# Read in tag code metadata
tag_code_metadata <- read.csv(here::here("covariate_data", "tag_code_metadata.csv"))

tag_code_metadata %>% 
  dplyr::select(tag_code, release_site_name) %>% 
  left_join(., natal_origin_table, by = "release_site_name") -> tag_code_origins

# Remove fish from Klickitat and Wind River
KLIC_WIND_tag_codes <- subset(tag_code_origins, natal_origin %in% c("Klickitat_River", "Wind_River"))$tag_code
# Subset out those
states_complete %>% 
  subset(., !(tag_code %in% KLIC_WIND_tag_codes)) -> states_complete

# Keep fish from only years that arrays were active
tag_code_metadata %>% 
  left_join(., natal_origin_table, by = "release_site_name") %>% 
  dplyr::select(tag_code, run_year, natal_origin, rear_type_code) %>% 
  subset(., tag_code %in% states_complete$tag_code) -> origin_metadata

# Remove individuals from run years where arrays were not active
origin_metadata %>% 
  mutate(year_numeric = as.numeric(sub("\\/.*", "", run_year))) %>% 
  mutate(array_active = ifelse(natal_origin == "Methow_River" & year_numeric < 9, "inactive",
                               ifelse(natal_origin == "Klickitat_River" & year_numeric < 12, "inactive",
                                      ifelse(natal_origin == "Okanogan_River" & year_numeric < 14, "inactive",
                                             ifelse(natal_origin == "Wind_River" & year_numeric < 13, "inactive",
                                                    ifelse(natal_origin == "Asotin_Creek" & year_numeric < 12, "inactive", "active")))))) -> origin_metadata

##### Check implicit site visit code for any nonsensical transitions #####

states_complete %>% 
  mutate(prev_tag_code = lag(tag_code)) %>% 
  mutate(next_tag_code = lead(tag_code)) %>% 
  mutate(next_state = lead(state)) %>% 
  mutate(next_state = ifelse(next_tag_code != tag_code, "loss", next_state)) -> transitions

transitions %>% 
  dplyr::select(tag_code, state, next_state) -> transitions_2

transitions_2 %>% 
  group_by(state, next_state) %>% 
  summarise(count = n()) -> transitions_count
  
# Check for weird stuff
transitions_count %>% 
  # See if we're getting transitions to the same state - which should never happen
  subset(state == next_state) -> bad_transitions

# See if we're getting transitions between states that are not adjacent
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



transitions_count %>% 
  mutate(state_type = ifelse(state %in% tributary_mainstem$tributary, "tributary", "mainstem")) %>% 
  subset(state_type == "tributary" & next_state != "loss") -> tributary_transition_counts

# Issues:
# Fish going straight from Entiat to mainstem, RIS to RRE
# Fish going straight from Methow to mainstem, RRE to WEL

transitions_count %>% 
  mutate(state_type = ifelse(state %in% tributary_mainstem$tributary, "tributary", "mainstem")) %>% 
  mutate(next_state_type = ifelse(next_state %in% tributary_mainstem$tributary, "tributary", "mainstem")) %>% 
  subset(state_type == "mainstem" & next_state_type == "tributary") -> mainstem_tributary_transition_counts

# Issue: Hood River isn't listed in tributaries - FIXED
# The rest of these look okay

# Look at sites where there are too many implicit site visits, for example 3D9.1BF26D8CCC or 384.1B796A3D09 or 384.3B23987601 
subset(states_complete, tag_code == "3D9.1BF26D8CCC") # this one is clearly a repeat spawner

states_complete %>% 
  group_by(tag_code) %>% 
  subset(pathway == "implicit") %>% 
  summarise(count = n()) -> implicit_state_counts_by_tag_code

# First let's look at ones where we have a lot of implicit site visits
jumper_tag_codes <- subset(implicit_state_counts_by_tag_code, count >= 3)$tag_code

jumper_state_transitions <- subset(states_complete, tag_code %in% jumper_tag_codes)

# Most of these actually looks fine, there are just some where there are too many implicit site visits - let's make our subset more specific

states_complete %>% 
  subset(pathway == "implicit") %>% 
  group_by(tag_code, state) %>% 
  summarise(count = n()) -> implicit_state_counts_by_tag_code_state

# Look at tag codes
repeat_implicit_sites_tag_codes <- subset(implicit_state_counts_by_tag_code_state, count > 1)$tag_code

repeat_implicit_sites_state_transitions <- subset(states_complete, tag_code %in% repeat_implicit_sites_tag_codes)

##### Subset the tag codes that are bad #####
transitions %>% 
  subset(state == next_state) -> duplicate_state_transitions

repeat_state_tags <- unique(duplicate_state_transitions$tag_code)

subset(states_complete, tag_code %in% repeat_state_tags) -> repeat_sites_fish

transitions %>% 
  subset(state == "Entiat River" & next_state == "mainstem, RIS to RRE" | state == "Methow River" & next_state == "mainstem, RRE to WEL") -> skip_trib_transitions

skip_trib_transition_tags <- unique(skip_trib_transitions$tag_code)

states_complete %>% 
  subset(tag_code %in% skip_trib_transition_tags) -> skip_trib_fish

# let's take a quick look at known fallback
subset(transitions, pathway == "BON_fallback_arrays") %>% 
  mutate(date_time = ymd_hms(date_time)) %>% 
  mutate(month = month(date_time)) %>% 
  mutate(month = ifelse(month == 1, "Jan",
                        ifelse(month == 2, "Feb",
                               ifelse(month == 3, "Mar",
                                      ifelse(month == 4, "Apr",
                                             ifelse(month == 5, "May",
                                                    ifelse(month == 6, "Jun",
                                                           ifelse(month == 7, "Jul",
                                                                  ifelse(month == 8, "Aug",
                                                                         ifelse(month == 9, "Sep",
                                                                                ifelse(month == 10, "Oct",
                                                                                       ifelse(month == 11, "Nov",
                                                                                              ifelse(month == 12, "Dec", NA))))))))))))) %>% 
  mutate(month = factor(month, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))) -> BON_fallback_transitions_detected
ggplot(BON_fallback_transitions_detected, aes(x = date_time)) +
  geom_histogram()

ggplot(BON_fallback_transitions_detected, aes(x = month)) +
  geom_bar(stat = "count") +
  ggtitle("Fallback detections at BON (Juvenile Bypass + Corner Collector)")

# Look at which fish are being seeing in fallback at BON
BON_fallback_tag_codes <- unique(BON_fallback_transitions_detected$tag_code)

BON_fallback_fish <- subset(states_complete, tag_code %in% BON_fallback_tag_codes)



##### Create a table that contains every transition + the relevant metadata for that fish #####
transitions %>% 
  left_join(., origin_metadata, by = "tag_code") -> transitions_meta

write.csv(transitions_meta, here::here("figures", "transitions_table_2022-07-16.csv"))
# write.csv(transitions_meta, here::here("figures", "transitions_table_testpart1.csv"))

# Let's look at the frequency of going back to the mainstem vs. loss for tributaries

transitions %>% 
  subset(state %in% tributary_mainstem$tributary) %>% 
  group_by(state, next_state) %>% 
  summarise(count = n()) -> trib_transition_counts

# Manually remove the bad ones just for now
trib_transition_counts <- subset(trib_transition_counts, !(next_state %in% c("Umatilla River", "Deschutes River", "Hood River")))
trib_transition_counts <- subset(trib_transition_counts, !(state == "Entiat River" & next_state == "mainstem, RIS to RRE"))
trib_transition_counts <- subset(trib_transition_counts, !(state == "Methow River" & next_state == "mainstem, RRE to WEL"))
trib_transition_counts <- subset(trib_transition_counts, !(state == "Hood River" & next_state == "mainstem, ICH to LGR"))
trib_transition_counts <- subset(trib_transition_counts, !(state == "Hood River" & next_state == "mainstem, MCN to ICH or PRA"))
trib_transition_counts <- subset(trib_transition_counts, !(state == "Hood River" & next_state == "mainstem, mouth to BON"))
trib_transition_counts <- subset(trib_transition_counts, !(state == "Hood River" & next_state == "mainstem, PRA to RIS"))

# Just get rid of Hood entirely for now because I think there's a typo
trib_transition_counts <- subset(trib_transition_counts, !(state == "Hood River"))

trib_transition_counts %>% 
  group_by(state) %>% 
  summarise(total = sum(count)) -> total_trib_transitions


trib_transition_counts %>% 
  left_join(., total_trib_transitions, by = "state") %>% 
  mutate(prop = count/total) -> trib_transition_counts

trib_transition_counts %>% 
  subset(next_state == "loss") %>% 
  arrange(prop)


##### Create summaries for each state #####
full_transition_matrix <- read.csv(here::here("figures", "full_model_transition_matrix.csv")) 
colnames(full_transition_matrix) <- c("from", full_transition_matrix$X)


# Figure out which natal origins are missing for each mainstem state
origins_vec <- unique(natal_origin_table$natal_origin)

# Take out the Wind and Klickitat
origins_vec <- origins_vec[!(origins_vec %in% c("Wind_River", "Klickitat_River"))]


mainstem_sites <- c(unique(tributary_mainstem$mainstem), "mainstem, PRA to RIS", "mainstem, mouth to BON")
tributary_sites <- unique(tributary_mainstem$tributary)

transitions_meta %>% 
  subset(state %in% mainstem_sites) %>% 
  group_by(state, natal_origin) %>% 
  summarise(count = n()) -> mainstem_transitions_count_origin

for (i in 1:length(mainstem_sites)){
  subset(mainstem_transitions_count_origin, state == mainstem_sites[i]) -> mainstem_state_dat
  
  missing_origins <- setdiff(origins_vec, mainstem_state_dat$natal_origin)
  print(paste0("state: ", mainstem_sites[i]))
  print(paste0("missing origins: ", missing_origins))
}
  

##### Identify repeat spawners #####

# How are we going to identify repeat spawners?

# First let's look for individuals with really long det hists
states_complete %>% 
  mutate(date_time = ymd_hms(date_time)) %>% 
  group_by(tag_code) %>% 
  summarise(min(date_time, na.rm = TRUE)) -> tag_code_min_time

colnames(tag_code_min_time) <- c("tag_code", "min_time")

states_complete %>% 
  mutate(date_time = ymd_hms(date_time)) %>% 
  group_by(tag_code) %>% 
  summarise(max(date_time, na.rm = TRUE)) -> tag_code_max_time

colnames(tag_code_max_time) <- c("tag_code", "max_time")

tag_code_min_time %>% 
  left_join(., tag_code_max_time, by = "tag_code") %>% 
  mutate(total_time = max_time - min_time) -> tag_code_time

# Subset tags that are in the Columbia for over a year
tag_code_time %>% 
  subset(., total_time > days(x = 365)) -> year_tags

# Look for clear typos - three years
tag_code_time %>% 
  # subset(., total_time > days(x = 365*3))
  subset(., total_time > years(x = 3))-> three_year_tags

# Export the clearly wrong ones
bad_tags <- subset(three_year_tags, !(tag_code %in% c("3D9.1C2C430C8D", "3D9.1C2DA0CF28")))$tag_code
write.csv(bad_tags, file = here::here("bad_tags.csv"))

# Now inspect states complete for these tags
complete_det_hist <- read.csv(here::here("from_hyak_transfer", "2022-05-24_det_hist", "complete_det_hist.csv"), row.names = 1)
complete_det_hist %>% 
  subset(tag_code %in% year_tags$tag_code) -> year_tags_det_hist

complete_det_hist %>% 
  subset(tag_code %in% three_year_tags$tag_code) -> three_year_tags_det_hist


# Select individual fish

# Clear kelt 1 (WEL to BON) 3D9.1BF15F8910
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "3D9.1BF15F8910")

subset(states_complete, tag_code == "3D9.1BF15F8910") %>% dplyr::select(state, date_time)

# clear kelt 2 (Entiat/RRE juv to BON) 3D9.1C2D9303CC
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "3D9.1C2D9303CC")

subset(states_complete, tag_code == "3D9.1C2D9303CC") %>% dplyr::select(state, date_time)

# clear kelt 3 (Wells to BCC to BON) 3D9.1BF18BD51C
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "3D9.1BF18BD51C")

subset(states_complete, tag_code == "3D9.1BF18BD51C") %>% dplyr::select(state, date_time)

# Double BON 1
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "384.3B2396BB0B")

# Double BON 2 384.3B23AB34F1
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "384.3B23AB34F1")

# Double LGR - see CTH1
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "384.3B23A2FFB1")


# Triple spawner? 3D9.1C2C430C8D
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "3D9.1C2C430C8D")

# Triple spawner 2 3D9.1C2DA0CF29
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "3D9.1C2DA0CF28")


# Clear data errors: 3D9.1BF1104C52
complete_det_hist %>% 
  dplyr::select(tag_code, event_site_name, start_time) %>% 
  subset(., tag_code == "3D9.1BF1104C52")













##### UNUSED CODE #####

# Plot transitions function
state_trans_summary_plot <- function(from_state, data){
  # Subset state of interest
  state_dat <- subset(data, state == from_state)
  
  # Summarise transitions
  state_dat %>% 
    group_by(next_state) %>% 
    summarise(count = n()) -> state_counts
  
  state_counts %>% 
    mutate(index = row_number()) -> state_counts
  
  ggplot(state_counts, aes(y = index, label = next_state)) +
    geom_segment(aes(x = 1, y = (max(index)-1)/2, xend = 3, yend = index),
                 arrow = arrow(length = unit(0.5, "cm"))) +
    geom_text(aes(x = 3.1, y = index, hjust = 0)) +
    scale_x_continuous(limits = c(-2, 6)) +
    annotate(geom = "text", x = 0.9, hjust = 1, y = (max(state_counts$index)-1)/2, label = from_state)
}

