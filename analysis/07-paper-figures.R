### Paper figures finalization ###

# This script creates the figures that aren't created by scripts run on hyak.

library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(rlang)
library(tibble)
library(forcats)
library(lubridate)
library(ggrepel)

#### Table 1: Sample sizes by run year and population ####

# Step 1: Load the model_data files to ensure that we are pulling the exact same
# data as what is going into the model
load(here::here("stan_actual", "reparameterization_v3", "upper_columbia_wild", "model_data.rda"))
UCW_data <- data
load(here::here("stan_actual", "reparameterization_v3", "upper_columbia_hatchery", "model_data.rda"))
UCH_data <- data
load(here::here("stan_actual", "reparameterization_v3", "middle_columbia_wild", "model_data.rda"))
MCW_data <- data
load(here::here("stan_actual", "reparameterization_v3", "middle_columbia_hatchery", "model_data.rda"))
MCH_data <- data
load(here::here("stan_actual", "reparameterization_v3", "snake_river_wild", "model_data.rda"))
SRW_data <- data
load(here::here("stan_actual", "reparameterization_v3", "snake_river_hatchery", "model_data.rda"))
SRH_data <- data

# load some useful info
# get the model states into a df, to help with interpretation
model_states = c(
  # Mainstem states (9)
  "mainstem, mouth to BON",
  "mainstem, BON to MCN",
  "mainstem, MCN to ICH or PRA",
  "mainstem, PRA to RIS",
  "mainstem, RIS to RRE",
  "mainstem, RRE to WEL",
  "mainstem, upstream of WEL",
  "mainstem, ICH to LGR",
  "mainstem, upstream of LGR",
  
  # Tributary states ()
  # With detection efficiencies in the model, we now have more tributary states,
  # since we have an upstream and a river mouth state
  
  # "Deschutes River", 
  "Deschutes River Mouth", "Deschutes River Upstream",
  # "John Day River", 
  "John Day River Mouth", "John Day River Upstream",
  # "Hood River",
  "Hood River Mouth", "Hood River Upstream",
  # "Fifteenmile Creek", 
  "Fifteenmile Creek Mouth", "Fifteenmile Creek Upstream",
  # "Umatilla River",
  "Umatilla River Mouth", "Umatilla River Upstream",
  # "Yakima River",
  "Yakima River Mouth", "Yakima River Upstream",
  # "Walla Walla River",
  "Walla Walla River Mouth", "Walla Walla River Upstream",
  # "Wenatchee River", 
  "Wenatchee River Mouth", "Wenatchee River Upstream",
  # "Entiat River", 
  "Entiat River Mouth", "Entiat River Upstream",
  # "Okanogan River", 
  "Okanogan River Mouth", "Okanogan River Upstream",
  # "Methow River", 
  "Methow River Mouth", "Methow River Upstream",
  # "Tucannon River",
  "Tucannon River Mouth", "Tucannon River Upstream",
  # "Asotin Creek", 
  "Asotin Creek Mouth", "Asotin Creek Upstream",
  "Clearwater River",
  "Salmon River",
  "Grande Ronde River",
  # "Imnaha River",
  "Imnaha River Mouth", "Imnaha River Upstream",
  "BON to MCN other tributaries",
  "Upstream WEL other tributaries",
  
  # Loss
  "loss"
)
natal_origins <- gsub(" Mouth| Upstream", "", model_states)
natal_origins <- natal_origins[!(duplicated(natal_origins))]
natal_origins <- natal_origins[!(grepl("mainstem", natal_origins))]
natal_origins <- natal_origins[!(grepl("other tributaries", natal_origins))]
natal_origins <- natal_origins[!(natal_origins == "loss")]

# Use the parameter map to index the right effects
origin_param_map <- data.frame(
  natal_origin = natal_origins,
  hatchery = c(NA, NA, NA, NA, 1, NA, 2, # MC
               1,NA,2,3, # UC
               5,NA,1,4,2,3), # SR,
  wild = c(1,3,NA,2,4,6,5, # MC
           1,2,NA,3, # UC
           6,1,2,5,3,4)) # SR

# Add the DPS
DPS_map <- data.frame(
  natal_origin = c("Middle Columbia", "Upper Columbia", "Snake River"),
  hatchery = rep(999, 3),
  wild = rep(999, 3)
)
origin_param_map %>% 
  bind_rows(DPS_map) -> origin_param_map




# Step 2: Write a function that takes a dataset and returns a table
# that contains each fish by origin and run year

extract_run_year_origin <- function(data, rear, DPS, origin_param_map = origin_param_map){
  origin_vector <- vector(length = nrow(data$cat_X_mat))
  for(i in 1:nrow(data$cat_X_mat)){
    origin_vector[i] <- which(data$cat_X_mat[i,]==1)
  }
  
  
  # for spill days - include the winter post-overshoot vector, which contains
  # info on whether they could have experienced winter spill conditions or not
  # add a new fish_ID column, which is not the same as tag code but will allow us to differentiate between fish
  pop_states_run_years <- data.frame(fish_ID = rep(1:length(origin_vector), each = ncol(data$y)),
                                     state = as.vector(t(data$y)),
                                     origin = rep(origin_vector, each = ncol(data$y)))
  
  # drop observations in the loss state and with index 0
  pop_states_run_years %>% 
    filter(!(state %in% c(0,43))) -> pop_states_run_years
  
  # add the run year column
  pop_states_run_years$run_year <- data$transition_run_years
  
  # keep only the first observation of each fish (this is when they are detected
  # at BON, and therefore is their run year)
  pop_states_run_years %>% 
    filter(!(duplicated(fish_ID))) -> fish_state_origin_year
  
  # join with key to translate numeric origin into actual origin
  if (DPS == "Middle Columbia") {
    origin_param_map_DPS <- subset(origin_param_map, natal_origin %in% origin_param_map$natal_origin[1:7])
  } else if (DPS == "Upper Columbia"){
    origin_param_map_DPS <- subset(origin_param_map, natal_origin %in% origin_param_map$natal_origin[8:11])
  } else{
    origin_param_map_DPS <- subset(origin_param_map, natal_origin %in% origin_param_map$natal_origin[12:17])
  }
  
  
  if (rear == "hatchery"){
    fish_state_origin_year %>% 
      left_join(., dplyr::select(origin_param_map_DPS, natal_origin, hatchery), 
                by = join_by(origin == hatchery)) -> fish_state_origin_year
    
    # Add a rear column
    fish_state_origin_year %>% 
      mutate(rear = "hatchery") -> fish_state_origin_year
  } else {
    fish_state_origin_year %>% 
      left_join(., dplyr::select(origin_param_map_DPS, natal_origin, wild), 
                by = join_by(origin == wild)) -> fish_state_origin_year
    
    # Add a rear column
    fish_state_origin_year %>% 
      mutate(rear = "natural") -> fish_state_origin_year
  }
  
  
  
  return(fish_state_origin_year)
  
}

UCW_origin_run_years <- extract_run_year_origin(data = UCW_data, rear = "natural", 
                                                DPS = "Upper Columbia", origin_param_map = origin_param_map)
UCH_origin_run_years <- extract_run_year_origin(data = UCH_data, rear = "hatchery", 
                                                DPS = "Upper Columbia", origin_param_map = origin_param_map)
MCW_origin_run_years <- extract_run_year_origin(data = MCW_data, rear = "natural", 
                                                DPS = "Middle Columbia", origin_param_map = origin_param_map)
MCH_origin_run_years <- extract_run_year_origin(data = MCH_data, rear = "hatchery", 
                                                DPS = "Middle Columbia", origin_param_map = origin_param_map)
SRW_origin_run_years <- extract_run_year_origin(data = SRW_data, rear = "natural", 
                                                DPS = "Snake River", origin_param_map = origin_param_map)
SRH_origin_run_years <- extract_run_year_origin(data = SRH_data, rear = "hatchery", 
                                                DPS = "Snake River", origin_param_map = origin_param_map)

UCW_origin_run_years %>% 
  bind_rows(., UCH_origin_run_years) %>% 
  bind_rows(., MCW_origin_run_years) %>% 
  bind_rows(., MCH_origin_run_years) %>% 
  bind_rows(., SRW_origin_run_years) %>% 
  bind_rows(., SRH_origin_run_years) -> full_origin_run_years

full_origin_run_years %>% 
  group_by(natal_origin, run_year, rear) %>% 
  summarise(count = n()) -> sample_size_long

# Convert numeric run years to actual
# first create the run year df
# NOTE: There are two ways that we index run year numeric (not great.)
# In the model, year 1 is 04/05; in other places where we are plotting run year,
# we have a numeric run year column where 4 = 04/05, 5 = 05/06, etc.
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2022-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2023-05-31 23:59:59"), by = "years")
run_year_index = 1:19

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_index)

sample_size_long %>% 
  dplyr::rename(run_year_index = run_year) %>% 
  left_join(run_year_df, by = "run_year_index") %>% 
  dplyr::select(-c(run_year_start, run_year_end)) -> sample_size_long

# change order to be from most downstream to most upstream
tributary_order <- c("Fifteenmile Creek", "Deschutes River", "John Day River",
                     "Umatilla River", "Walla Walla River", "Yakima River", "Wenatchee River",
                     "Entiat River", "Methow River", "Okanogan River",
                     "Tucannon River", "Clearwater River", "Asotin Creek", "Grande Ronde River",
                     "Salmon River", "Imnaha River")


# reformat so that we can export this is a table
sample_size_long %>% 
  ungroup() %>% 
  mutate(run_year = factor(run_year, levels = c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23"))) %>% 
  mutate(natal_origin = factor(natal_origin, levels = tributary_order)) %>% 
  dplyr::select(-c(run_year_index)) %>% 
  arrange(natal_origin, rear, run_year) %>% 
  mutate(population = paste0(natal_origin, " (", ifelse(rear == "natural", "N", "H"), ")")) %>% 
  dplyr::select(-c(natal_origin, rear)) %>% 
  pivot_wider(names_from = run_year, values_from = count) %>% 
  replace(is.na(.), 0) %>%
  relocate(`05/06`, .after = population) %>% 
  relocate(`06/07`, .after = `05/06`) %>% 
  relocate(`07/08`, .after = `06/07`) %>% 
  bind_rows(summarise(.,
                      across(where(is.numeric), sum),
                      across(where(is.character), ~"Total"))) %>% 
  mutate(Total = rowSums(across(where(is.numeric))))  -> sample_size_table

write.csv(sample_size_table, here::here("figures", "paper_figures", "Table1_sample_size.csv"))



#### Figure 7: Final fates vs. covariates ####

# Load the files
load(here::here("figures", "paper_figures", "final_fates_covariates", "UMA_MCN_winterspill_homing.rda"))
load(here::here("figures", "paper_figures", "final_fates_covariates", "ENT_WEL_winterspill_homing.rda"))
load(here::here("figures", "paper_figures", "final_fates_covariates", "JDR_MCN_winterspill_homing.rda"))
load(here::here("figures", "paper_figures", "final_fates_covariates", "TUC_LGR_winterspill_homing.rda"))
load(here::here("figures", "paper_figures", "final_fates_covariates", "WAWA_ICH_winterspill_homing.rda"))
load(here::here("figures", "paper_figures", "final_fates_covariates", "WEN_RRE_winterspill_homing.rda"))
load(here::here("figures", "paper_figures", "final_fates_covariates", "YAK_PRA_winterspill_homing.rda"))


plot_ff_cov <- function(data, dam_name, dam_spill_column){
  rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
  rear_shapes <- c(hatchery = 24, natural = 21)
  condition_colors <- c("coldest" = "#2c7bb6", "average" = "goldenrod2", "warmest" = "#d7191c")
  
  
  # fix the condition levels
  data %>% 
    mutate(conditions = ifelse(conditions == "cool", "coldest",
                               ifelse(conditions == "warm", "warmest", conditions))) -> data
  
  # fix the rear levels
  data %>% 
    mutate(rear_type = ifelse(rear_type == "wild", "natural", rear_type)) -> data
  
  # add a condition x rear
  data %>% 
    mutate(rear_condition = paste(rear_type, conditions)) -> data
  
  rear_condition_fills <- c("natural coldest" =  "#2c7bb6",  "hatchery coldest" = "white",
                      "natural average" = "goldenrod2",  "hatchery average" = "white",
                      "natural warmest" = "#d7191c", "hatchery warmest" = "white")
  
  data$conditions <- factor(data$conditions, levels = c("coldest", "average", "warmest"))
  
  plot <- ggplot(data, aes(x = eval(parse(text = dam_spill_column)), y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = conditions,
                fill = rear_condition, shape = rear_type)) +
    geom_point(aes(size = rear_type), position=position_dodge(width=3)) +
    geom_linerange(position=position_dodge(width=3)) +
    ylab("Homing Probability") +
    xlab(paste0("Days of Winter Spill at ", dam_name)) +
    # ggtitle(" ") +
    # Create a scale common to all
    scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_x_continuous(lim = c(-3, 95), breaks = c(0, 30, 60, 90)) +
    scale_fill_manual(values = rear_condition_fills, guide = "none") +
    scale_size_manual(values = c("hatchery" = 4, "natural" = 3.5), guide = "none") +
    scale_shape_manual(values = rear_shapes) +
    scale_color_manual(values = condition_colors) +
    theme(plot.title = element_text(size = 12),
          # axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = "black"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.4),
          legend.position = "none",
          plot.margin = unit(c(1.2, 0.2, 0.2, 0.2),"cm"))
  
  return(plot)
}


# John Day
JDR_FF_cov_plot <- plot_ff_cov(data = JDR_MCN_winterspill_homing, dam_name = "McNary Dam", dam_spill_column = "MCN_winterspill_actual")
# Umatilla
UMA_FF_cov_plot <- plot_ff_cov(data = UMA_MCN_winterspill_homing, dam_name = "McNary Dam", dam_spill_column = "MCN_winterspill_actual")
# Walla Walla
WAWA_FF_cov_plot <- plot_ff_cov(data = WAWA_ICH_winterspill_homing, dam_name = "Ice Harbor Dam", dam_spill_column = "ICH_winterspill_actual")
# Wenatchee
WEN_FF_cov_plot <- plot_ff_cov(data = WEN_RRE_winterspill_homing, dam_name = "Rocky Reach Dam", dam_spill_column = "RRE_winterspill_actual")
# Tucannon
TUC_FF_cov_plot <- plot_ff_cov(data = TUC_LGR_winterspill_homing, dam_name = "Lower Granite Dam", dam_spill_column = "LGR_winterspill_actual")
# Entiat
ENT_FF_cov_plot <- plot_ff_cov(data = ENT_WEL_winterspill_homing, dam_name = "Wells Dam", dam_spill_column = "WEL_winterspill_actual")
# Yakima
YAK_FF_cov_plot <- plot_ff_cov(data = YAK_PRA_winterspill_homing, dam_name = "Priest Rapids Dam", dam_spill_column = "PRA_winterspill_actual")

# Create the legend figure by itself
rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
rear_shapes <- c(hatchery = 24, natural = 21)
condition_colors <- c("coldest" = "#2c7bb6", "average" = "goldenrod2", "warmest" = "#d7191c")


# fix the condition levels
WAWA_ICH_winterspill_homing %>% 
  mutate(conditions = ifelse(conditions == "cool", "coldest",
                             ifelse(conditions == "warm", "warmest", conditions))) -> WAWA_ICH_winterspill_homing

# fix the rear levels
WAWA_ICH_winterspill_homing %>% 
  mutate(rear_type = ifelse(rear_type == "wild", "natural", rear_type)) -> WAWA_ICH_winterspill_homing

# add a condition x rear
WAWA_ICH_winterspill_homing %>% 
  mutate(rear_condition = paste(rear_type, conditions)) -> WAWA_ICH_winterspill_homing

rear_condition_fills <- c("natural coldest" =  "#2c7bb6",  "hatchery coldest" = "white",
                          "natural average" = "goldenrod2",  "hatchery average" = "white",
                          "natural warmest" = "#d7191c", "hatchery warmest" = "white")

WAWA_ICH_winterspill_homing$conditions <- factor(WAWA_ICH_winterspill_homing$conditions, levels = c("coldest", "average", "warmest"))

plot_for_legend <- ggplot(WAWA_ICH_winterspill_homing, aes(x = ICH_winterspill_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = conditions,
                         fill = rear_condition, shape = rear_type)) +
  geom_point(size = 8, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3), show.legend = FALSE) +
  ylab("Homing Probability") +
  xlab(paste0("Days of Winter Spill at Ice Harbor Dam")) +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_fill_manual(values = rear_condition_fills, guide = "none") +
  scale_size_manual(values = c("hatchery" = 10, "natural" = 8), guide = "none") +
  scale_shape_manual(values = rear_shapes) +
  # scale_shape_manual(values = c(hatchery = 24, natural = 19)) +
  scale_color_manual(values = condition_colors) +
  guides(shape=guide_legend(title="Rearing Type", override.aes = list(fill = c("white", "black"),
                                                                      size = c(8, 10))),
         color = guide_legend(title = "Basin Conditions")) +
  # theme(plot.title = element_text(size = 12),
  #       # axis.text.y = element_text(color = rev(state_significance_colors)),
  #       axis.title = element_text(size = 12),
  #       axis.text.x = element_text(size = 12))
  theme(panel.grid.major = element_line(color = "gray90"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
        legend.key.height = unit(1.25, "cm"),
        legend.key.width = unit(1.25, "cm"),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 15),
        # these plot margins are to leave space for the population name on the big figure
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))

FF_cov_legend <- ggpubr::get_legend(plot_for_legend)
FF_cov_legend_gg <- as_ggplot(FF_cov_legend) + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))

FF_cov_combined_plot <- cowplot::ggdraw(ggarrange(JDR_FF_cov_plot, UMA_FF_cov_plot, WAWA_FF_cov_plot,
                                  WEN_FF_cov_plot, TUC_FF_cov_plot, ENT_FF_cov_plot,
                                  YAK_FF_cov_plot, FF_cov_legend_gg, nrow = 2, ncol = 4,
                                          labels = c("(A) JDR", " (B) UMA", "(C) WAWA", 
                                                     "(D) WEN", "(E) TUC", "(F) ENT", "(G) YAK"),
                                  label.x = 0.05, label.y = 0.95,  font.label = list(size = 14, face = "plain"),
                                          hjust = 0, vjust = 0)) + theme(plot.background = element_rect(fill="white", color = NA))


ggsave(here::here("figures", "paper_figures", "Fig7_FF_cov_combined_plot.png"), FF_cov_combined_plot, height = 12, width = 16)



