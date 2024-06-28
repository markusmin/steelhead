# 04.5-data-summaries

### Description ###
# This script creates summaries of the data itself, in order to generate some
# figures/tables for the paper, as well as to compare with the model results


#### Load libraries and data ####
library(tidyverse)
library(here)
library(ggthemes)

# Load states complete (the complete detection history)
ASC <- read.csv(here::here("stan_actual", "adults_states_complete.csv"), row.names = 1)

# We currently have some fish that are categorized as "loss" because they were trapped and removed
# They are all hatchery fish
# For the purposes of data summary, we're going to reclassify these as "Trapped"
ASC %>% 
  mutate(state = ifelse(pathway == "WEL_trap_arrays", "Wells Trap", state)) -> ASC


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

# update the natal origin field to drop the underscore (this will allow it to match state names)
tag_code_metadata %>% 
  mutate(natal_origin = gsub("_", " ", natal_origin)) -> tag_code_metadata

ASC %>% 
  left_join(dplyr::select(tag_code_metadata, natal_origin, ESU, tag_code, 
                          run_year, run_year_numeric, rear_type_code), by = "tag_code") -> ASC


#### Summary data tables ####

ASC %>% 
  group_by(natal_origin, run_year, run_year_numeric, rear_type_code) %>% 
  filter(!(duplicated(tag_code_2))) %>% 
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




#### Plot homing ####
plot_homing <- function(natal_origin_select, rear_type_code_select, data){
  plot_data <- filter(data, natal_origin == natal_origin_select & rear_type_code == rear_type_code_select)
  
  # Determine which are stray tributaries
  all_tribs <- unique(ASC$natal_origin)
  stray_tribs <- all_tribs[!(all_tribs == natal_origin_select)]
  
  # Get homing counts
  plot_data %>% 
    mutate(final_fate = ifelse(state == natal_origin_select, "Home",
                         ifelse(state %in% stray_tribs, "Stray",
                                "Other"))) -> plot_data
  
  plot <- ggplot(plot_data, aes(x = run_year_numeric, y = count, fill = final_fate)) +
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
  
  homing_plot_list[[length(natal_origins)+i]] <- plot_homing(natal_origin_select = natal_origins[i],
              rear_type_code_select = "H",
              data = ASC_final_states_count)
  
}

# make a pdf of these
pdf(file = here::here("figures", "home_stray_from_data.pdf"))
homing_plot_list
dev.off()



# Plot final fate
plot_final_fate <- function(natal_origin_select, rear_type_code_select, data){
  plot_data <- filter(data, natal_origin == natal_origin_select & rear_type_code == rear_type_code_select)
  
  # Make it so that you only plot the top nine states, plus other
  plot_data %>% 
    group_by(state) %>% 
    summarise(total = sum(count)) -> state_visit_counts
  
  arrange(state_visit_counts, desc(total))$state[1:9] -> top9states
  
  plot_data %>% 
    mutate(final_fate = ifelse(state %in% top9states, state, "Other")) -> plot_data
  
  plot_data %>% 
    mutate(final_fate = factor(final_fate)) -> plot_data
  
  # make homing first category and other last category
  plot_data$final_fate <- forcats::fct_relevel(plot_data$final_fate, natal_origin_select)
  plot_data$final_fate <- forcats::fct_relevel(plot_data$final_fate, "Other", after = Inf)
  
  
  plot <- ggplot(plot_data, aes(x = run_year_numeric, y = count, fill = final_fate)) +
    geom_bar(stat = "identity") +
    ggtitle(paste0(natal_origin_select, " - ", rear_type_code_select)) +
    scale_fill_tableau(palette = "Tableau 10")
  
  return(plot)
}

natal_origins <- unique(ASC_final_states_count$natal_origin)

final_fate_plot_list <- vector(mode = "list", length = length(natal_origins)*2)

for (i in 1:length(natal_origins)){
  final_fate_plot_list[[i]] <- plot_final_fate(natal_origin_select = natal_origins[i],
                                       rear_type_code_select = "W",
                                       data = ASC_final_states_count)
  
  final_fate_plot_list[[length(natal_origins)+i]] <- plot_final_fate(natal_origin_select = natal_origins[i],
                                         rear_type_code_select = "H",
                                         data = ASC_final_states_count)
  
}

# make a pdf of these
pdf(file = here::here("figures", "final_fates_from_data.pdf"))
final_fate_plot_list
dev.off()


# Clearly there is a pattern related to tributary arrays - need to understand
# at what point we dropped the years of origins where we didn't have detection capacity


#### Investigate final fates across all years ####
# Let's drop all observations from non-DE years, to keep comparison most accurate
deschutes_river_trib_det_eff_capability <- data.frame(natal_origin = "Deschutes River",
                                                      run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                      DE = c(rep(0,9), rep(1,6), rep(0,3)))

john_day_river_trib_det_eff_capability <- data.frame(natal_origin = "John Day River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,8), rep(1,10)))

hood_river_trib_det_eff_capability <- data.frame(natal_origin = "Hood River",
                                                 run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                 DE = c(rep(0,9), rep(1,9)))

fifteenmile_creek_trib_det_eff_capability <- data.frame(natal_origin = "Fifteenmile Creek",
                                                        run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                        DE = c(rep(0,8), rep(1,10)))

umatilla_river_trib_det_eff_capability <- data.frame(natal_origin = "Umatilla River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,3), rep(1,15)))

yakima_river_trib_det_eff_capability <- data.frame(natal_origin = "Yakima River",     
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   # DE = c(rep(1,18)))
                                                   # For now - remove 04/05 run year
                                                   DE = c(0, rep(1,17)))

walla_walla_river_trib_det_eff_capability <- data.frame(natal_origin = "Walla Walla River",
                                                        run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                        DE = c(rep(0,1), rep(1,17)))

wenatchee_river_trib_det_eff_capability <- data.frame(natal_origin = "Wenatchee River",
                                                      run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                      DE = c(rep(0,7), rep(1,11)))

entiat_river_trib_det_eff_capability <- data.frame(natal_origin = "Entiat River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,4), rep(1,14)))

okanogan_river_trib_det_eff_capability <- data.frame(natal_origin = "Okanogan River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,9), rep(1,9)))

methow_river_trib_det_eff_capability <- data.frame(natal_origin = "Methow River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,5), rep(1,13)))

tucannon_river_trib_det_eff_capability <- data.frame(natal_origin = "Tucannon River",
                                                     run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                     DE = c(rep(0,7), rep(1,11)))

asotin_creek_trib_det_eff_capability <- data.frame(natal_origin = "Asotin Creek",      
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,7), rep(1,11)))

imnaha_river_trib_det_eff_capability <- data.frame(natal_origin = "Imnaha River",
                                                   run_year = run_year[1:18], # ignore 22/23 to keep consistent
                                                   DE = c(rep(0,7), rep(1,11)))


deschutes_river_trib_det_eff_capability %>% 
  bind_rows(., john_day_river_trib_det_eff_capability) %>% 
  bind_rows(., hood_river_trib_det_eff_capability) %>% 
  bind_rows(., fifteenmile_creek_trib_det_eff_capability) %>% 
  bind_rows(., umatilla_river_trib_det_eff_capability) %>% 
  bind_rows(., yakima_river_trib_det_eff_capability) %>% 
  bind_rows(., walla_walla_river_trib_det_eff_capability) %>% 
  bind_rows(., wenatchee_river_trib_det_eff_capability) %>% 
  bind_rows(., entiat_river_trib_det_eff_capability) %>% 
  bind_rows(., okanogan_river_trib_det_eff_capability) %>% 
  bind_rows(., methow_river_trib_det_eff_capability) %>% 
  bind_rows(., tucannon_river_trib_det_eff_capability) %>% 
  bind_rows(., asotin_creek_trib_det_eff_capability) %>% 
  bind_rows(., imnaha_river_trib_det_eff_capability) -> trib_det_eff_capability




plot_final_fate_proportions_DE <- function(natal_origin_select,  data){
  filter_data <- filter(data, natal_origin == natal_origin_select)
  
  # join with DE capability DF, drop non-DE years of data
  filter_data %>% 
    left_join(., trib_det_eff_capability, by = c("natal_origin", "run_year")) %>% 
    filter(DE == 1) -> filter_data
  
  
  # drop Wells Trap fish
  filter_data %>% 
    filter(state != "Wells Trap") -> filter_data
  
  # summarise across all years
  filter_data %>% 
    group_by(state, rear_type_code) %>% 
    summarise(count = sum(count)) -> final_fates_counts
  
  final_fates_counts %>% 
    ungroup() %>% 
    group_by(rear_type_code) %>% 
    summarise(total = sum(count)) -> rear_type_counts
  
  final_fates_counts %>% 
    left_join(., rear_type_counts, by = "rear_type_code") %>% 
    mutate(prop = count/total) -> final_fates_counts
      
      
  states_order_for_plot <- gsub(" Mouth| Upstream", "", model_states)
  states_order_for_plot <- states_order_for_plot[!(duplicated(states_order_for_plot))]
  # Make a couple of changes to make them be in the order from most downstream to most upstream
  states_order_for_plot[10] <- "Hood River"
  states_order_for_plot[11] <- "Fifteenmile Creek"
  states_order_for_plot[12] <- "Deschutes River"
  states_order_for_plot[13] <- "John Day River"
  states_order_for_plot[15] <- "Walla Walla River"
  states_order_for_plot[16] <- "Yakima River"
  states_order_for_plot[19] <- "Methow River"
  states_order_for_plot[20] <- "Okanogan River"
  
  states_order_for_plot[16:29] <- states_order_for_plot[15:28]
  states_order_for_plot[15] <- "BON to MCN other tributaries"
  states_order_for_plot[23:29] <- states_order_for_plot[22:28]
  states_order_for_plot[22] <- "Upstream WEL other tributaries"
  states_order_for_plot[29] <- "loss"
  
  states_order_for_plot[24] <- "Clearwater River"
  states_order_for_plot[25] <- "Asotin Creek"
  states_order_for_plot[26] <- "Grande Ronde River"
  states_order_for_plot[27] <- "Salmon River"
  
  rear_colors <- c(H = "#ff7f00", W = "#33a02c")
  
  
  # add all the unvisited states for plotting
  full_states <- data.frame(state = rep(states_order_for_plot,2),
                            rear_type_code = rep(c("H","W"), each = length(states_order_for_plot)),
                            prop = rep(0, length(states_order_for_plot)*2))
  
  final_fates_counts %>% 
    bind_rows(., full_states) %>% 
    distinct(state, rear_type_code, .keep_all = TRUE) -> final_fates_counts
  
  
  final_fates_counts$state <- fct_rev(factor(final_fates_counts$state, levels = states_order_for_plot))
  
  final_fates_counts_plot <- ggplot(final_fates_counts, aes(x = state, y = prop, color = rear_type_code)) +
    geom_point(size = 3.5, shape = 18, position=position_dodge(width=0.5)) +
    coord_flip() +
    ylab("Final Fates") +
    xlab("Model State") +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = rear_colors) +
    theme(plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    ggtitle(natal_origin_select)
  
  
  
}

final_fate_DE_data_plot_list <- vector(mode = "list", length = length(natal_origins))

for (i in 1:length(natal_origins)){
  final_fate_DE_data_plot_list[[i]] <- plot_final_fate_proportions_DE(natal_origin_select = natal_origins[i],
                                               data = ASC_final_states_count)
  
}

# make a pdf of these
pdf(file = here::here("figures", "final_fates_from_data_DE_years.pdf"))
final_fate_DE_data_plot_list
dev.off()

plot_final_fate_proportions_all_years <- function(natal_origin_select,  data){
  filter_data <- filter(data, natal_origin == natal_origin_select)
  
  
  # drop Wells Trap fish
  filter_data %>% 
    filter(state != "Wells Trap") -> filter_data
  
  # summarise across all years
  filter_data %>% 
    group_by(state, rear_type_code) %>% 
    summarise(count = sum(count)) -> final_fates_counts
  
  final_fates_counts %>% 
    ungroup() %>% 
    group_by(rear_type_code) %>% 
    summarise(total = sum(count)) -> rear_type_counts
  
  final_fates_counts %>% 
    left_join(., rear_type_counts, by = "rear_type_code") %>% 
    mutate(prop = count/total) -> final_fates_counts
  
  
  states_order_for_plot <- gsub(" Mouth| Upstream", "", model_states)
  states_order_for_plot <- states_order_for_plot[!(duplicated(states_order_for_plot))]
  # Make a couple of changes to make them be in the order from most downstream to most upstream
  states_order_for_plot[10] <- "Hood River"
  states_order_for_plot[11] <- "Fifteenmile Creek"
  states_order_for_plot[12] <- "Deschutes River"
  states_order_for_plot[13] <- "John Day River"
  states_order_for_plot[15] <- "Walla Walla River"
  states_order_for_plot[16] <- "Yakima River"
  states_order_for_plot[19] <- "Methow River"
  states_order_for_plot[20] <- "Okanogan River"
  
  states_order_for_plot[16:29] <- states_order_for_plot[15:28]
  states_order_for_plot[15] <- "BON to MCN other tributaries"
  states_order_for_plot[23:29] <- states_order_for_plot[22:28]
  states_order_for_plot[22] <- "Upstream WEL other tributaries"
  states_order_for_plot[29] <- "loss"
  
  states_order_for_plot[24] <- "Clearwater River"
  states_order_for_plot[25] <- "Asotin Creek"
  states_order_for_plot[26] <- "Grande Ronde River"
  states_order_for_plot[27] <- "Salmon River"
  
  rear_colors <- c(H = "#ff7f00", W = "#33a02c")
  
  
  # add all the unvisited states for plotting
  full_states <- data.frame(state = rep(states_order_for_plot,2),
                            rear_type_code = rep(c("H","W"), each = length(states_order_for_plot)),
                            prop = rep(0, length(states_order_for_plot)*2))
  
  final_fates_counts %>% 
    bind_rows(., full_states) %>% 
    distinct(state, rear_type_code, .keep_all = TRUE) -> final_fates_counts
  
  
  final_fates_counts$state <- fct_rev(factor(final_fates_counts$state, levels = states_order_for_plot))
  
  final_fates_counts_plot <- ggplot(final_fates_counts, aes(x = state, y = prop, color = rear_type_code)) +
    geom_point(size = 3.5, shape = 18, position=position_dodge(width=0.5)) +
    coord_flip() +
    ylab("Final Fates") +
    xlab("Model State") +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_color_manual(values = rear_colors) +
    theme(plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12)) +
    ggtitle(natal_origin_select)
  
  
  
}

final_fate_data_plot_list <- vector(mode = "list", length = length(natal_origins))

for (i in 1:length(natal_origins)){
  final_fate_data_plot_list[[i]] <- plot_final_fate_proportions_all_years(natal_origin_select = natal_origins[i],
                                                                      data = ASC_final_states_count)
  
}

# make a pdf of these
pdf(file = here::here("figures", "final_fates_from_data_all_years.pdf"))
final_fate_data_plot_list
dev.off()



