# 04.5-data-summaries

### Description ###
# This script creates summaries of the data itself, in order to generate some
# figures/tables for the paper, as well as to compare with the model results


#### Load libraries and data ####
library(tidyverse)
library(here)

# Load states complete (the complete detection history)
ASC <- read.csv(here::here("stan_actual", "adults_states_complete.csv"), row.names = 1)

# Load the metadata on each fish
tag_code_metadata <- read.csv(here::here("covariate_data", "tag_code_metadata.csv"))

# create a run_year_numeric field
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_numeric = seq(4, 22, 1)
run_year_df <- data.frame(run_year,run_year_numeric)

tag_code_metadata %>% 
  left_join(., run_year_df, by = "run_year") -> tag_code_metadata


# load natal origins
natal_origin_table <- read.csv(here::here("covariate_data", "natal_origin_table.csv"))

# update the natal origin table to drop the underscore (this will allow it to match state names)
natal_origin_table %>% 
  mutate(natal_origin = gsub("_", " ", natal_origin)) -> natal_origin_table

# match natal origin to tag codes
tag_code_metadata %>% 
  left_join(natal_origin_table, by = "release_site_name") -> tag_code_metadata

tag_code_metadata %>% 
  mutate(ESU = ifelse(natal_origin %in% c("Tucannon_River", "Asotin_Creek", "Clearwater_River", "Salmon_River", "Grande_Ronde_River", "Imnaha_River"), "snake",
                      ifelse(natal_origin %in% c("Wenatchee_River", "Entiat_River", "Okanogan_River","Methow_River"), "upper_columbia",
                             ifelse(natal_origin %in% c("Deschutes_River", "Fifteenmile_Creek", "John_Day_River", "Umatilla_River", "Yakima_River", "Walla_Walla_River"), "middle_columbia",
                                    ifelse(natal_origin %in% c("Hood_River"), "lower_columbia", "ERROR"))))) %>% 
  # correct rear type code so that unknowns are correctly identified as hatchery fish
  mutate(rear_type_code = ifelse(rear_type_code == "U", "H", rear_type_code))-> tag_code_metadata

ASC %>% 
  left_join(dplyr::select(tag_code_metadata, natal_origin, ESU, tag_code, 
                          run_year, run_year_numeric, rear_type_code), by = "tag_code") -> ASC


#### Summary data tables ####

ASC %>% 
  group_by(natal_origin, run_year, run_year_numeric, rear_type_code) %>% 
  summarise(total = n()) -> ASC_origin_year_counts




#### Final fates ####

# Keep only the last observation of each fish
ASC %>% 
  # combine upstream and mouth states
  mutate(state = gsub(" Upstream", "", state)) %>% 
  mutate(state = gsub(" Mouth", "", state)) %>% 
  group_by(tag_code_2) %>% 
  slice(n()) -> ASC_final_states

# for each run year and origin, get the final fates
ASC_final_states %>% 
  group_by(natal_origin, run_year, run_year_numeric, rear_type_code, state) %>% 
  summarise(count = n()) -> ASC_final_states_count

# Get the proportion of the tagged fish from that run year that ended up in each state
ASC_final_states_count %>% 
  left_join(ASC_origin_year_counts, by = c("natal_origin", "run_year", "run_year_numeric", "rear_type_code")) %>% 
  mutate(prop = count/total) -> ASC_final_states_count

# Get homing counts
ASC_final_states_count %>% 
  mutate(home = ifelse(natal_origin == state, TRUE, FALSE)) -> ASC_final_states_count

# Plot homing
plot_homing <- function(natal_origin_select, rear_type_code_select, data){
  plot_data <- filter(data, natal_origin == natal_origin_select & rear_type_code == rear_type_code_select)
  
  plot <- ggplot(plot_data, aes(x = run_year_numeric, y = count, fill = home)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0(natal_origin_select, " - ", rear_type_code_select))
  
  return(plot)
}

natal_origins <- unique(ASC_final_states_count$natal_origin)

homing_plot_list <- vector(mode = "list", length = length(natal_origins)*2)

for (i in 1:length(natal_origins)){
  homing_plot_list[[i]] <- plot_homing(natal_origin_select = natal_origins[i],
              rear_type_code_select = "W",
              data = ASC_final_states_count)
  
  homing_plot_list[[i*2]] <- plot_homing(natal_origin_select = natal_origins[i],
              rear_type_code_select = "H",
              data = ASC_final_states_count)
  
}

# Clearly there is a pattern related to tributary arrays - need to understand
# at what point we dropped the years of origins where we didn't have detection capacity





