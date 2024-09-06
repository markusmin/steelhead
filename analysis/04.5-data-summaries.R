# 04.5-data-summaries

### Description ###
# This script creates summaries of the data itself, in order to generate some
# figures/tables for the paper, as well as to compare with the model results


#### Load libraries and data ####
library(tidyverse)
library(here)
library(ggthemes)

### DE information
# create a run_year_numeric field
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_numeric = seq(4, 22, 1)
run_year_df <- data.frame(run_year,run_year_numeric)

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



# Load states complete (the complete detection history)
ASC <- read.csv(here::here("stan_actual", "adults_states_complete.csv"), row.names = 1)

# We currently have some fish that are categorized as "loss" because they were trapped and removed
# They are all hatchery fish
# For the purposes of data summary, we're going to reclassify these as "Trapped"
ASC %>% 
  mutate(state = ifelse(pathway == "WEL_trap_arrays", "Wells Trap", state)) -> ASC


# Load the metadata on each fish
tag_code_metadata <- read.csv(here::here("covariate_data", "tag_code_metadata.csv"))


# add run year info
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

# New investigation: Look at movements for only DE years
# filter_ASC <- filter(ASC, natal_origin == natal_origin_select)
# 
# # join with DE capability DF, drop non-DE years of data
# filter_ASC %>%
#   left_join(., trib_det_eff_capability, by = c("natal_origin", "run_year")) %>%
#   filter(DE == 1) %>%
#   # drop the upstream transitions, since that happens in the modeling code for DE years
#   filter(!(grepl("Upstream", state))) -> filter_ASC
# 
# filter_ASC %>%
#   mutate(to_state = ifelse(tag_code_2 == lead(tag_code_2), lead(state), "loss")) -> filter_ASC
# 
# as.data.frame(table(subset(filter_ASC, state == "mainstem, BON to MCN")$to_state)) -> FIF_out2_transitions
# colnames(FIF_out2_transitions) <- c("state", "count")
# 
# FIF_out2_transitions %>%
#   mutate(prop = round(count/sum(count), 3)) -> FIF_out2_transitions


# New investigation: Let's look at amount of upstream vs. mouth detections at Deschutes

ASC %>%
  subset(grepl("Deschutes", state)) %>%
  # let's only look at DE years here, otherwise
  filter(run_year %in% c("13/14", "14/15", "15/16", "16/17", "17/18", "18/19")) %>%
  group_by(natal_origin, state) %>%
  summarise(count = n()) -> DES_trib_sites_count_by_origin

ASC %>%
  subset(grepl("Deschutes", state)) %>%
  # let's only look at DE years here, otherwise
  filter(run_year %in% c("13/14", "14/15", "15/16", "16/17", "17/18", "18/19")) %>%
  group_by(natal_origin) %>%
  summarise(total = n()) -> DES_trib_sites_count_total_by_origin

DES_trib_sites_count_by_origin %>%
  left_join(., DES_trib_sites_count_total_by_origin, by = "natal_origin") %>%
  mutate(prop = count/total) -> DES_trib_sites_count_by_origin

DES_det_upstream_mouth_plot <- ggplot(subset(DES_trib_sites_count_by_origin, state == "Deschutes River Mouth"), aes(x = fct_rev(as.factor(natal_origin)), y = prop, size = total)) +
  geom_point() +
  # theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(lim = c(0,1)) +
  xlab("Natal Origin") +
  ylab("Proportion of tributary detections at the mouth site, Deschutes River") +
  coord_flip()

ggsave(here::here("docs", "site_figures", "data_plots", "DES_det_upstream_mouth_plot.png"), DES_det_upstream_mouth_plot, height = 8, width = 8)


# Keep only the last observation of each fish
ASC %>% 
  # # combine upstream and mouth states
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

natal_origins <- unique(ASC_final_states_count$natal_origin)


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


plot_final_fate_proportions_DE <- function(natal_origin_select,  data, remove_upstream_trib_sites = F){
  if (remove_upstream_trib_sites == T){
    data %>% 
      # drop the upstream transitions, since that happens in the modeling code for DE years
        filter(!(grepl("Upstream", state))) %>% 
      # change the mouth detections to just say the tributary name
      mutate(state = gsub(" Mouth", "", state)) %>%
        group_by(tag_code_2) %>% 
        slice(n()) -> data_final_states
  } else {
    data %>% 
      # # combine upstream and mouth states
      mutate(state = gsub(" Upstream", "", state)) %>%
      mutate(state = gsub(" Mouth", "", state)) %>%
      group_by(tag_code_2) %>% 
      slice(n()) -> data_final_states
  }
  
  
  # for each run year and origin, get the final fates
  data_final_states %>% 
    group_by(natal_origin, run_year, run_year_numeric, rear_type_code, state) %>% 
    summarise(count = n()) -> data_final_states_count
  
  # Get the proportion of the tagged fish from that run year that ended up in each state
  data_final_states_count %>% 
    left_join(ASC_origin_year_counts, by = c("natal_origin", "run_year", "run_year_numeric", "rear_type_code")) %>% 
    mutate(prop = count/total) -> data_final_states_count
  
  
  filter_data <- filter(data_final_states_count, natal_origin == natal_origin_select)
  
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
  
  # Drop rear type if not modeled
  origin_rear_map <- data.frame(natal_origin = natal_origins,
                                rear = c("W", "both", "W",
                                         "W", "W", "both", "both",
                                         "both", "W", "both",
                                         "H", "both", "both",
                                         "both", "both", "both",
                                         "W"))
  
  subset(origin_rear_map, natal_origin == natal_origin_select)$rear -> rear
  
  if (rear == "W"){
    filter(final_fates_counts, rear_type_code == "W") -> final_fates_counts
  } else if (rear == "H"){
    filter(final_fates_counts, rear_type_code == "H") -> final_fates_counts
  } else{
    final_fates_counts -> final_fates_counts
  }
  
  
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
                                               data = ASC)
  
}

# make a pdf of these
pdf(file = here::here("figures", "final_fates_from_data_DE_years.pdf"))
final_fate_DE_data_plot_list
dev.off()

# save individually
for (i in 1:length(natal_origins)){
  ggsave(here::here("figures", "final_fates_from_data", paste0(gsub(" ", "_", natal_origins[i]), "_final_fates_from_data.png")),
         final_fate_DE_data_plot_list[[i]], height = 8, width = 8)
  
}

# Make a second version - now where we drop the upstream trib states, like we did in the model
final_fate_DE_no_upstream_trib_data_plot_list <- vector(mode = "list", length = length(natal_origins))

for (i in 1:length(natal_origins)){
  final_fate_DE_no_upstream_trib_data_plot_list[[i]] <- plot_final_fate_proportions_DE(natal_origin_select = natal_origins[i],
                                                                      data = ASC,
                                                                      remove_upstream_trib_sites = T)
  
}

# make a pdf of these
pdf(file = here::here("figures", "final_fates_from_data_DE_years_no_upstream_trib.pdf"))
final_fate_DE_no_upstream_trib_data_plot_list
dev.off()



# save individually
for (i in 1:length(natal_origins)){
   ggsave(here::here("figures", "final_fates_from_data_no_upstream_trib", paste0(gsub(" ", "_", natal_origins[i]), "_final_fates_from_data_no_upstream_trib.png")),
          final_fate_DE_no_upstream_trib_data_plot_list[[i]], height = 8, width = 8)
  
}

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

#### Final fates: fit to data ####
load(here::here("docs", "site_figures", "final_fates", "FF_comp_data.rda"))

bind_rows(FF_comp_data, .id = "natal_origin_abbrev") -> FF_sim_data
FF_sim_data$natal_origin_abbrev <- gsub("_ff_comp_median", "", FF_sim_data$natal_origin_abbrev)

natal_origin_abbrev <- data.frame(natal_origin = unique(subset(ASC_final_states_count, natal_origin != "Hood River")$natal_origin),
                                  natal_origin_abbrev = unique(FF_sim_data$natal_origin_abbrev)[order(unique(FF_sim_data$natal_origin_abbrev))])

FF_sim_data %>% 
  left_join(natal_origin_abbrev, by = "natal_origin_abbrev") -> FF_sim_data

### Prepare final fates data to join

# So, this is not easy to set up. If we join based on DE in natal tributary, we get
# more accurate estimates of homing. But then if there are large straying rates,
# we don't estimate those well, since we're not paying attention to which years 
# the stray tributaries have high DE rates

# 2024-08-28 edit: Don't join based on natal_origin, join based on the tributary!!!
model_states_df <- data.frame(state = 1:43, state_name = model_states)
trib_det_eff_capability %>% 
  mutate(tributary = paste0(natal_origin, " Mouth")) %>% 
  left_join(., model_states_df, join_by(tributary == state_name)) %>% 
  mutate(state = gsub(" Mouth", "", tributary)) %>% 
  dplyr::select(run_year, DE, state) -> trib_det_eff_capability_forjoin


# join with DE capability DF, drop non-DE years of data

### Option 1: this is the version where we look at years where the NATAL ORIGIN TRIBUTARY has 
# DE capability, for each origin
ASC_final_states_count %>% 
  left_join(., trib_det_eff_capability, by = c("natal_origin", "run_year")) %>% 
  filter(DE == 1) -> DE_years_final_states_count


# drop Wells Trap fish
DE_years_final_states_count %>% 
  filter(state != "Wells Trap") -> DE_years_final_states_count
# summarise across all years
DE_years_final_states_count %>% 
  group_by(state, rear_type_code, natal_origin) %>% 
  summarise(count = sum(count)) -> final_fates_counts
final_fates_counts %>% 
  ungroup() %>% 
  group_by(rear_type_code, natal_origin) %>% 
  summarise(total = sum(count)) -> origin_rear_type_counts
final_fates_counts %>% 
  left_join(., origin_rear_type_counts, by = c("rear_type_code", "natal_origin")) %>% 
  mutate(prop = count/total) %>% 
  mutate(rear_type = ifelse(rear_type_code == "W", "wild", "hatchery")) -> final_fates_counts


FF_sim_data %>% 
  left_join(final_fates_counts, by = c("rear_type", "natal_origin", "state")) %>% 
  # deal with NAs
  # mutate(rear_type = ifelse(rear_type_code == "W", "wild", "hatchery")) %>% 
  dplyr::select(-c(rear_type_code, total)) %>% 
  mutate(count = ifelse(is.na(count), 0, count)) %>% 
  mutate(prop = ifelse(is.na(prop), 0, prop)) -> FF_sim_data_comp

### Option 2: DESCHUTES FOCUS: Drop all non-DE years for Deschutes
DE_years_final_states_count %>% 
  filter(run_year %in% subset(trib_det_eff_capability_forjoin, state == "Deschutes River" & DE == 1)$run_year) -> DE_years_final_states_counts_Deschutes_DE

# drop Wells Trap fish
DE_years_final_states_counts_Deschutes_DE %>% 
  filter(state != "Wells Trap") -> DE_years_final_states_counts_Deschutes_DE
# summarise across all years
DE_years_final_states_counts_Deschutes_DE %>% 
  group_by(state, rear_type_code, natal_origin) %>% 
  summarise(count = sum(count)) -> final_fates_counts_Deschutes_DE
final_fates_counts_Deschutes_DE %>% 
  ungroup() %>% 
  group_by(rear_type_code, natal_origin) %>% 
  summarise(total = sum(count)) -> origin_rear_type_counts_Deschutes_DE
final_fates_counts_Deschutes_DE %>% 
  left_join(., origin_rear_type_counts_Deschutes_DE, by = c("rear_type_code", "natal_origin")) %>% 
  mutate(prop = count/total) %>% 
  mutate(rear_type = ifelse(rear_type_code == "W", "wild", "hatchery")) -> final_fates_counts_Deschutes_DE


FF_sim_data %>% 
  left_join(final_fates_counts_Deschutes_DE, by = c("rear_type", "natal_origin", "state")) %>% 
  # deal with NAs
  # mutate(rear_type = ifelse(rear_type_code == "W", "wild", "hatchery")) %>% 
  dplyr::select(-c(rear_type_code, total)) %>% 
  mutate(count = ifelse(is.na(count), 0, count)) %>% 
  mutate(prop = ifelse(is.na(prop), 0, prop)) -> FF_sim_data_comp_Deschutes_DE


#### Final fates fit to data plot ####

# Create the plot
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

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")

FF_sim_data_comp$state <- fct_rev(factor(FF_sim_data_comp$state, levels = states_order_for_plot))

# drop the three snake river origins without DE capability
FF_sim_data_comp_DE <- subset(FF_sim_data_comp, !(natal_origin %in% c("Clearwater River", "Grande Ronde River", "Salmon River")))

# show which is home vs. not home
# show which is Deschutes
FF_sim_data_comp_DE %>% 
  mutate(state_type = ifelse(state == natal_origin, "home",
                             ifelse(state == "Deschutes River", "Deschutes",
                                    ifelse(grepl("mainstem", state), "mainstem", "tributary")))) -> FF_sim_data_comp_DE

FF_sim_data_comp_plot <- ggplot(FF_sim_data_comp_DE, aes(x = prop, y = `0.5`,  ymin = `0.025`, ymax = `0.975`, color = rear_type, shape = state_type)) +
  annotate(geom = "segment", x = 0, y = 0, xend = 1, yend =1, lty = 2) +
  geom_point() +
  # geom_linerange(position=position_dodge(width=0.5)) +
  facet_wrap(~natal_origin) +
  scale_color_manual(values = rear_colors) +
  ylab("Simulated final fates") +
  xlab("Data Final Fates") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14))

ggsave(here::here("docs", "site_figures", "final_fates", "FF_sim_data_comp_plot.png"), FF_sim_data_comp_plot, height = 10, width = 10)

# create a second version, where we only show those with >1%
FF_sim_data_comp_DE_trimmed <- subset(FF_sim_data_comp_DE, prop > 0.01)
FF_sim_data_comp_plot_trimmed <- ggplot(FF_sim_data_comp_DE_trimmed, aes(x = prop, y = `0.5`,  ymin = `0.025`, ymax = `0.975`, color = rear_type, shape = state_type)) +
  annotate(geom = "segment", x = 0, y = 0, xend = 1, yend =1, lty = 2) +
  geom_point() +
  # geom_linerange(position=position_dodge(width=0.5)) +
  facet_wrap(~natal_origin) +
  scale_color_manual(values = rear_colors) +
  ylab("Simulated final fates") +
  xlab("Data Final Fates") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14))

# New plot: Show percent difference by state
FF_sim_data_comp_DE %>% 
  mutate(ratio = (`0.5`/prop)) %>% 
  mutate(ratio = ifelse(ratio == Inf, NaN, ratio)) -> FF_sim_data_comp_DE

# FF_sim_data_perc_comp_DE_plot <- ggplot(FF_sim_data_comp_DE, aes(x = state, y = ratio, color = rear_type, shape = state_type, size = count)) +
#   # annotate(geom = "segment", x = 0, y = 1, xend = 1, yend =1, lty = 2) +
#   geom_hline(yintercept = 1, linetype= "dashed") + 
#   geom_point() +
#   # geom_linerange(position=position_dodge(width=0.5)) +
#   facet_wrap(~natal_origin) +
#   scale_color_manual(values = rear_colors) +
#   ylab("Final Fates ratio for model:data") +
#   xlab("State") +
#   # ggtitle(" ") +
#   # Create a scale common to all
#   theme(plot.title = element_text(size = 12),
#         # axis.text.y = element_text(color = rev(state_significance_colors)),
#         axis.text.x = element_text(angle = 90),
#         axis.title = element_text(size = 14))

FF_sim_data_perc_comp_DE_plot <- ggplot(FF_sim_data_comp_DE, aes(x = count, y = ratio, color = rear_type, shape = state_type)) +
  # annotate(geom = "segment", x = 0, y = 1, xend = 1, yend =1, lty = 2) +
  geom_hline(yintercept = 1, linetype= "dashed") + 
  geom_point() +
  # geom_linerange(position=position_dodge(width=0.5)) +
  facet_wrap(~natal_origin) +
  scale_color_manual(values = rear_colors) +
  ylab("Final Fates ratio for model:data") +
  xlab("Count of final fates (data)") +
  # ggtitle(" ") +
  # Create a scale common to all
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 14))

ggsave(here::here("docs", "site_figures", "final_fates", "FF_sim_data_perc_comp_DE_plot.png"), FF_sim_data_perc_comp_DE_plot, height = 10, width = 10)

subset(FF_sim_data_comp_DE, state_type == "tributary" & ratio > 2)
subset(FF_sim_data_comp_DE, state_type == "Deschutes" & ratio > 2)
subset(FF_sim_data_comp_DE, ratio > 3 & `0.5`>0.02)

#### PLOT VERSION 2: Deschutes Focus ####

# Create the plot
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

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")

FF_sim_data_comp_Deschutes_DE$state <- fct_rev(factor(FF_sim_data_comp_Deschutes_DE$state, levels = states_order_for_plot))

# drop the three snake river origins without DE capability
FF_sim_data_comp_Deschutes_DE_DE <- subset(FF_sim_data_comp_Deschutes_DE, !(natal_origin %in% c("Clearwater River", "Grande Ronde River", "Salmon River")))

# show which is home vs. not home
# show which is Deschutes
FF_sim_data_comp_Deschutes_DE_DE %>% 
  mutate(state_type = ifelse(state == natal_origin, "home",
                             ifelse(state == "Deschutes River", "Deschutes",
                                    ifelse(grepl("mainstem", state), "mainstem", "tributary")))) -> FF_sim_data_comp_Deschutes_DE_DE

FF_sim_data_comp_Deschutes_DE_plot <- ggplot(FF_sim_data_comp_Deschutes_DE_DE, aes(x = prop, y = `0.5`,  ymin = `0.025`, ymax = `0.975`, color = rear_type, shape = state_type)) +
  annotate(geom = "segment", x = 0, y = 0, xend = 1, yend =1, lty = 2) +
  geom_point() +
  # geom_linerange(position=position_dodge(width=0.5)) +
  facet_wrap(~natal_origin) +
  scale_color_manual(values = rear_colors) +
  ylab("Simulated final fates") +
  xlab("Data Final Fates") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14))

ggsave(here::here("docs", "site_figures", "final_fates", "FF_sim_data_comp_Deschutes_DE_plot.png"), FF_sim_data_comp_Deschutes_DE_plot, height = 10, width = 10)

# create a second version, where we only show those with >1%
FF_sim_data_comp_Deschutes_DE_DE_trimmed <- subset(FF_sim_data_comp_Deschutes_DE_DE, prop > 0.01)
FF_sim_data_comp_Deschutes_DE_plot_trimmed <- ggplot(FF_sim_data_comp_Deschutes_DE_DE_trimmed, aes(x = prop, y = `0.5`,  ymin = `0.025`, ymax = `0.975`, color = rear_type, shape = state_type)) +
  annotate(geom = "segment", x = 0, y = 0, xend = 1, yend =1, lty = 2) +
  geom_point() +
  # geom_linerange(position=position_dodge(width=0.5)) +
  facet_wrap(~natal_origin) +
  scale_color_manual(values = rear_colors) +
  ylab("Simulated final fates") +
  xlab("Data Final Fates") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_x_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14))

# New plot: Show percent difference by state
FF_sim_data_comp_Deschutes_DE_DE %>% 
  mutate(ratio = (`0.5`/prop)) %>% 
  mutate(ratio = ifelse(ratio == Inf, NaN, ratio)) -> FF_sim_data_comp_Deschutes_DE_DE

# FF_sim_data_perc_comp_DE_plot <- ggplot(FF_sim_data_comp_Deschutes_DE_DE, aes(x = state, y = ratio, color = rear_type, shape = state_type, size = count)) +
#   # annotate(geom = "segment", x = 0, y = 1, xend = 1, yend =1, lty = 2) +
#   geom_hline(yintercept = 1, linetype= "dashed") + 
#   geom_point() +
#   # geom_linerange(position=position_dodge(width=0.5)) +
#   facet_wrap(~natal_origin) +
#   scale_color_manual(values = rear_colors) +
#   ylab("Final Fates ratio for model:data") +
#   xlab("State") +
#   # ggtitle(" ") +
#   # Create a scale common to all
#   theme(plot.title = element_text(size = 12),
#         # axis.text.y = element_text(color = rev(state_significance_colors)),
#         axis.text.x = element_text(angle = 90),
#         axis.title = element_text(size = 14))

FF_sim_data_perc_comp_DE_plot_Deschutes_DE <- ggplot(FF_sim_data_comp_Deschutes_DE_DE, aes(x = count, y = ratio, color = rear_type, shape = state_type)) +
  # annotate(geom = "segment", x = 0, y = 1, xend = 1, yend =1, lty = 2) +
  geom_hline(yintercept = 1, linetype= "dashed") + 
  geom_point() +
  # geom_linerange(position=position_dodge(width=0.5)) +
  facet_wrap(~natal_origin) +
  scale_color_manual(values = rear_colors) +
  ylab("Final Fates ratio for model:data") +
  xlab("Count of final fates (data)") +
  # ggtitle(" ") +
  # Create a scale common to all
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.text.x = element_text(angle = 90),
        axis.title = element_text(size = 14))

ggsave(here::here("docs", "site_figures", "final_fates", "FF_sim_data_perc_comp_DE_plot_Deschutes_DE.png"), FF_sim_data_perc_comp_DE_plot_Deschutes_DE, height = 10, width = 10)







# Lots of Methow fish in the Okanogan? Is it also a similar issue as the Deschutes?
# See plot below - it looks to be a similar issue, but not the same magnitude, so it doesn't entirely explain it

ASC %>% 
  subset(grepl("Okanogan", state)) %>% 
  # let's only look at DE years here, otherwise 
  filter(run_year %in% c("13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21", "21/22")) %>% 
  group_by(natal_origin, state) %>% 
  summarise(count = n()) -> OKA_trib_sites_count_by_origin

ASC %>% 
  subset(grepl("Okanogan", state)) %>% 
  # let's only look at DE years here, otherwise 
  filter(run_year %in% c("13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21", "21/22")) %>% 
  group_by(natal_origin) %>% 
  summarise(total = n()) -> OKA_trib_sites_count_total_by_origin

OKA_trib_sites_count_by_origin %>% 
  left_join(., OKA_trib_sites_count_total_by_origin, by = "natal_origin") %>% 
  mutate(prop = count/total) -> OKA_trib_sites_count_by_origin

ggplot(subset(OKA_trib_sites_count_by_origin, state == "Okanogan River Mouth"), aes(x = natal_origin, y = prop, size = total)) +
  geom_point() +
  # theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(lim = c(0,1)) +
  ylab("Proportion detections at mouth vs. upstream sites, Okanogan River") +
  coord_flip()

#### Umatilla deep dive ####
subset(FF_sim_data_comp_DE, natal_origin == "Umatilla River" & `0.5` > 0.01)




