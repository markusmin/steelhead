# 19_stan_actual_deteff_results


# This script will load and investigate the outputs from our stan models.
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(tidyverse)
library(here)
library(ggpubr)


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

# Get info about state names and numbers
from_state_number_names <- data.frame(from = seq(1,43,1), from_name = model_states)
to_state_number_names <- data.frame(to = seq(1,43,1), to_name = model_states)



# Read in 200iter runs
snake_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "snake", "200iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))
middle_columbia_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "middle_columbia", "200iter_parallel_middle_columbia_stan_actual_int_origin_stan_fit.rds"))
upper_columbia_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "upper_columbia", "200iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))

# Get each of the summary objects
snake_summary <- snake_fit$summary()
middle_columbia_summary <- middle_columbia_fit$summary()
upper_columbia_summary <- upper_columbia_fit$summary()

snake_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat))  %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable)))  %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Tucannon_River",
                         ifelse(origin == "origin2", "Asotin_Creek",
                                ifelse(origin == "origin3",  "Clearwater_River",
                                       ifelse(origin == "origin4", "Salmon_River",
                                              ifelse(origin == "origin5", "Grande_Ronde",
                                                     ifelse(origin == "origin6", "Imnaha_River", "error"))))))) -> snake_derived_probabilities 

middle_columbia_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat))  %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable)))  %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Deschutes_River",
                         ifelse(origin == "origin2", "Fifteenmile_Creek",
                                ifelse(origin == "origin3",  "John_Day_River",
                                       ifelse(origin == "origin4", "Umatilla_River",
                                              ifelse(origin == "origin5", "Yakima_River",
                                                     ifelse(origin == "origin6", "Walla_Walla_River", "error"))))))) -> middle_columbia_derived_probabilities 

upper_columbia_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat))  %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable)))  %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Wenatchee_River",
                         ifelse(origin == "origin2", "Entiat_River",
                                ifelse(origin == "origin3",  "Okanogan_River",
                                       ifelse(origin == "origin4", "Methow_River", "error"))))) -> upper_columbia_derived_probabilities 

# Split this up into DE and NDE probabilities
snake_derived_probabilities %>% 
  filter(.,grepl("_NDE", variable)) -> snake_NDE_derived_probabilities
snake_derived_probabilities %>% 
  filter(.,grepl("_DE", variable)) -> snake_DE_derived_probabilities

middle_columbia_derived_probabilities %>% 
  filter(.,grepl("_NDE", variable)) -> middle_columbia_NDE_derived_probabilities
middle_columbia_derived_probabilities %>% 
  filter(.,grepl("_DE", variable)) -> middle_columbia_DE_derived_probabilities

upper_columbia_derived_probabilities %>% 
  filter(.,grepl("_NDE", variable)) -> upper_columbia_NDE_derived_probabilities
upper_columbia_derived_probabilities %>% 
  filter(.,grepl("_DE", variable)) -> upper_columbia_DE_derived_probabilities


##### SNAKE - Export tables of parameters for report #####

# DE CORRECTED
snake_DE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> snake_DE_prob_table

# subset the probabilities shared across the DPS
snake_DE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                     "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River Mouth",             
                     "John Day River Mouth",  "Hood River Mouth",  "Fifteenmile Creek Mouth",  "Umatilla River Mouth", "Yakima River Mouth",
                     "Walla Walla River Mouth", "Wenatchee River Mouth", "Entiat River Mouth",     "Okanogan River Mouth",    "Methow River Mouth",
                     "Deschutes River Upstream",             
                     "John Day River Upstream",  "Hood River Upstream",  "Fifteenmile Creek Upstream",  "Umatilla River Upstream", "Yakima River Upstream",
                     "Walla Walla River Upstream", "Wenatchee River Upstream", "Entiat River Upstream",     "Okanogan River Upstream",    "Methow River Upstream",
                     "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> snake_DE_DPS_probs

snake_DE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                       "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River Mouth",             
                       "John Day River Mouth",  "Hood River Mouth",  "Fifteenmile Creek Mouth",  "Umatilla River Mouth", "Yakima River Mouth",
                       "Walla Walla River Mouth", "Wenatchee River Mouth", "Entiat River Mouth",     "Okanogan River Mouth",    "Methow River Mouth",
                       "Deschutes River Upstream",             
                       "John Day River Upstream",  "Hood River Upstream",  "Fifteenmile Creek Upstream",  "Umatilla River Upstream", "Yakima River Upstream",
                       "Walla Walla River Upstream", "Wenatchee River Upstream", "Entiat River Upstream",     "Okanogan River Upstream",    "Methow River Upstream",
                       "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> snake_DE_origin_probs

# export both tables

# Change output to be mean (5% - 95% CI), rather than three separate columns
snake_DE_DPS_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability) -> snake_DE_DPS_probs

snake_DE_origin_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability, origin) -> snake_DE_origin_probs

write.csv(snake_DE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "snake_DE_DPS_probs.csv"))
write.csv(snake_DE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "snake_DE_origin_probs.csv"))



# NOT DE CORRECTED
snake_NDE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> snake_NDE_prob_table

# subset the probabilities shared across the DPS
snake_NDE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                     "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River Mouth",             
                     "John Day River Mouth",  "Hood River Mouth",  "Fifteenmile Creek Mouth",  "Umatilla River Mouth", "Yakima River Mouth",
                     "Walla Walla River Mouth", "Wenatchee River Mouth", "Entiat River Mouth",     "Okanogan River Mouth",    "Methow River Mouth",
                     "Deschutes River Upstream",             
                     "John Day River Upstream",  "Hood River Upstream",  "Fifteenmile Creek Upstream",  "Umatilla River Upstream", "Yakima River Upstream",
                     "Walla Walla River Upstream", "Wenatchee River Upstream", "Entiat River Upstream",     "Okanogan River Upstream",    "Methow River Upstream",
                     "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> snake_NDE_DPS_probs

snake_NDE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON", "mainstem, BON to MCN", "mainstem, PRA to RIS",  "mainstem, RIS to RRE",
                       "mainstem, RRE to WEL", "mainstem, upstream of WEL",  "Deschutes River Mouth",             
                       "John Day River Mouth",  "Hood River Mouth",  "Fifteenmile Creek Mouth",  "Umatilla River Mouth", "Yakima River Mouth",
                       "Walla Walla River Mouth", "Wenatchee River Mouth", "Entiat River Mouth",     "Okanogan River Mouth",    "Methow River Mouth",
                       "Deschutes River Upstream",             
                       "John Day River Upstream",  "Hood River Upstream",  "Fifteenmile Creek Upstream",  "Umatilla River Upstream", "Yakima River Upstream",
                       "Walla Walla River Upstream", "Wenatchee River Upstream", "Entiat River Upstream",     "Okanogan River Upstream",    "Methow River Upstream",
                       "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> snake_NDE_origin_probs

# Change output to be mean (5% - 95% CI), rather than three separate columns
snake_NDE_DPS_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability) -> snake_NDE_DPS_probs

snake_NDE_origin_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability, origin) -> snake_NDE_origin_probs

# export both tables
write.csv(snake_NDE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "snake_NDE_DPS_probs.csv"))
write.csv(snake_NDE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "snake_NDE_origin_probs.csv"))

##### MIDDLE COLUMBIA - Export tables of parameters for report #####

# DE CORRECTED
middle_columbia_DE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> middle_columbia_DE_prob_table

# subset the probabilities shared across the DPS
middle_columbia_DE_prob_table %>% 
  subset(from %in% c("mainstem, RIS to RRE", "mainstem, RRE to WEL", "mainstem, upstream of WEL",
                     "mainstem, upstream of LGR", "Wenatchee River Mouth",  "Entiat River Mouth",
                     "Okanogan River Mouth", "Methow River Mouth",  "Tucannon River Mouth",  "Asotin Creek Mouth",
                     "Clearwater River", "Salmon River", "Grande Ronde River", "Imnaha River Mouth",
                     "Wenatchee River Upstream",  "Entiat River Upstream","Okanogan River Upstream", "Methow River Upstream",  
                     "Tucannon River Upstream",  "Asotin Creek Upstream", "Imnaha River Upstream")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> middle_columbia_DE_DPS_probs

middle_columbia_DE_prob_table %>% 
  subset(!(from %in% c("mainstem, RIS to RRE", "mainstem, RRE to WEL", "mainstem, upstream of WEL",
                       "mainstem, upstream of LGR", "Wenatchee River Mouth",  "Entiat River Mouth",
                       "Okanogan River Mouth", "Methow River Mouth",  "Tucannon River Mouth",  "Asotin Creek Mouth",
                       "Clearwater River", "Salmon River", "Grande Ronde River", "Imnaha River Mouth",
                       "Wenatchee River Upstream",  "Entiat River Upstream","Okanogan River Upstream", "Methow River Upstream",  
                       "Tucannon River Upstream",  "Asotin Creek Upstream", "Imnaha River Upstream"))) -> middle_columbia_DE_origin_probs

# Change output to be mean (5% - 95% CI), rather than three separate columns
middle_columbia_DE_DPS_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability) -> middle_columbia_DE_DPS_probs

middle_columbia_DE_origin_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability, origin) -> middle_columbia_DE_origin_probs


# export both tables
write.csv(middle_columbia_DE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "middle_columbia_DE_DPS_probs.csv"))
write.csv(middle_columbia_DE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "middle_columbia_DE_origin_probs.csv"))



# NOT DE CORRECTED
middle_columbia_NDE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> middle_columbia_NDE_prob_table

# subset the probabilities shared across the DPS
middle_columbia_NDE_prob_table %>% 
  subset(from %in% c("mainstem, RIS to RRE", "mainstem, RRE to WEL", "mainstem, upstream of WEL",
                     "mainstem, upstream of LGR", "Wenatchee River Mouth",  "Entiat River Mouth",
                     "Okanogan River Mouth", "Methow River Mouth",  "Tucannon River Mouth",  "Asotin Creek Mouth",
                     "Clearwater River", "Salmon River", "Grande Ronde River", "Imnaha River Mouth",
                     "Wenatchee River Upstream",  "Entiat River Upstream","Okanogan River Upstream", "Methow River Upstream",  
                     "Tucannon River Upstream",  "Asotin Creek Upstream", "Imnaha River Upstream")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> middle_columbia_NDE_DPS_probs

middle_columbia_NDE_prob_table %>% 
  subset(!(from %in% c("mainstem, RIS to RRE", "mainstem, RRE to WEL", "mainstem, upstream of WEL",
                       "mainstem, upstream of LGR", "Wenatchee River Mouth",  "Entiat River Mouth",
                       "Okanogan River Mouth", "Methow River Mouth",  "Tucannon River Mouth",  "Asotin Creek Mouth",
                       "Clearwater River", "Salmon River", "Grande Ronde River", "Imnaha River Mouth",
                       "Wenatchee River Upstream",  "Entiat River Upstream","Okanogan River Upstream", "Methow River Upstream",  
                       "Tucannon River Upstream",  "Asotin Creek Upstream", "Imnaha River Upstream"))) -> middle_columbia_NDE_origin_probs

# Change output to be mean (5% - 95% CI), rather than three separate columns
middle_columbia_NDE_DPS_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability) -> middle_columbia_NDE_DPS_probs

middle_columbia_NDE_origin_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability, origin) -> middle_columbia_NDE_origin_probs

# export both tables
write.csv(middle_columbia_NDE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "middle_columbia_NDE_DPS_probs.csv"))
write.csv(middle_columbia_NDE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "middle_columbia_NDE_origin_probs.csv"))

##### UPPER COLUMBIA - Export tables of parameters for report #####

# DE CORRECTED
upper_columbia_DE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> upper_columbia_DE_prob_table

# subset the probabilities shared across the DPS
upper_columbia_DE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON",  "mainstem, BON to MCN", "mainstem, ICH to LGR",           "mainstem, upstream of LGR",      "Deschutes River Mouth",
                     "Deschutes River Upstream", "John Day River Mouth", "Hood River Mouth",  "Fifteenmile Creek Mouth",       
                     "Umatilla River Mouth", "Umatilla River Upstream",  "Yakima River Mouth",    "Walla Walla River Mouth",            
                     "Tucannon River Mouth", "Asotin Creek Mouth", "Asotin Creek Upstream", "Clearwater River",              
                     "Salmon River", "Grande Ronde River", "Imnaha River Mouth", "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> upper_columbia_DE_DPS_probs

upper_columbia_DE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON",  "mainstem, BON to MCN", "mainstem, ICH to LGR",           "mainstem, upstream of LGR",      "Deschutes River Mouth",
                       "Deschutes River Upstream", "John Day River Mouth", "Hood River Mouth",  "Fifteenmile Creek Mouth",       
                       "Umatilla River Mouth", "Umatilla River Upstream",  "Yakima River Mouth",    "Walla Walla River Mouth",            
                       "Tucannon River Mouth", "Asotin Creek Mouth", "Asotin Creek Upstream", "Clearwater River",              
                       "Salmon River", "Grande Ronde River", "Imnaha River Mouth", "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> upper_columbia_DE_origin_probs

# Change output to be mean (5% - 95% CI), rather than three separate columns
upper_columbia_DE_DPS_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability) -> upper_columbia_DE_DPS_probs

upper_columbia_DE_origin_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability, origin) -> upper_columbia_DE_origin_probs
# export both tables
write.csv(upper_columbia_DE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "upper_columbia_DE_DPS_probs.csv"))
write.csv(upper_columbia_DE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "upper_columbia_DE_origin_probs.csv"))



# NOT DE CORRECTED
upper_columbia_NDE_derived_probabilities %>% 
  left_join(., from_state_number_names, by = "from") %>% 
  left_join(., to_state_number_names, by = "to") %>% 
  dplyr::select(mean, q5, q95, origin, from_name, to_name) %>% 
  relocate(to_name) %>% 
  relocate(from_name) %>% 
  dplyr::rename(from = from_name, to = to_name) -> upper_columbia_NDE_prob_table

# subset the probabilities shared across the DPS
upper_columbia_NDE_prob_table %>% 
  subset(from %in% c("mainstem, mouth to BON",  "mainstem, BON to MCN", "mainstem, ICH to LGR",           "mainstem, upstream of LGR",      "Deschutes River Mouth",
                     "Deschutes River Upstream", "John Day River Mouth", "Hood River Mouth",  "Fifteenmile Creek Mouth",       
                     "Umatilla River Mouth", "Umatilla River Upstream",  "Yakima River Mouth",    "Walla Walla River Mouth",            
                     "Tucannon River Mouth", "Asotin Creek Mouth", "Asotin Creek Upstream", "Clearwater River",              
                     "Salmon River", "Grande Ronde River", "Imnaha River Mouth", "BON to MCN other tributaries", "Upstream WEL other tributaries")) %>% 
  distinct(from, to, .keep_all = TRUE) %>% 
  dplyr::select(-origin) -> upper_columbia_NDE_DPS_probs

upper_columbia_NDE_prob_table %>% 
  subset(!(from %in% c("mainstem, mouth to BON",  "mainstem, BON to MCN", "mainstem, ICH to LGR",           "mainstem, upstream of LGR",      "Deschutes River Mouth",
                       "Deschutes River Upstream", "John Day River Mouth", "Hood River Mouth",  "Fifteenmile Creek Mouth",       
                       "Umatilla River Mouth", "Umatilla River Upstream",  "Yakima River Mouth",    "Walla Walla River Mouth",            
                       "Tucannon River Mouth", "Asotin Creek Mouth", "Asotin Creek Upstream", "Clearwater River",              
                       "Salmon River", "Grande Ronde River", "Imnaha River Mouth", "BON to MCN other tributaries", "Upstream WEL other tributaries"))) -> upper_columbia_NDE_origin_probs

# Change output to be mean (5% - 95% CI), rather than three separate columns
upper_columbia_NDE_DPS_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability) -> upper_columbia_NDE_DPS_probs

upper_columbia_NDE_origin_probs %>% 
  mutate(probability = paste0(formatC(mean, digits = 2, format = "fg"), " (", formatC(q5, digits = 2, format = "fg"), 
                              " - ", formatC(q95, digits = 2, format = "fg"), ")")) %>% 
  dplyr::select(from, to, probability, origin) -> upper_columbia_NDE_origin_probs

# export both tables
write.csv(upper_columbia_NDE_DPS_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "upper_columbia_NDE_DPS_probs.csv"))
write.csv(upper_columbia_NDE_origin_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "upper_columbia_NDE_origin_probs.csv"))

##### CONDITINONAL PROBABILITIES #####

# First - load functions
# 1) Function that reformats movement probabilities (probabilities) into a transition matrix
transition_matrix_reformat <- function(movement_probs, variable_prefix){
  movement_probs %>% 
    mutate(from = gsub(",.*$", "", variable)) %>% 
    mutate(from = gsub(paste0(variable_prefix, "\\["), "", from)) %>% 
    mutate(to = gsub(paste0(variable_prefix, "\\[.*,"), "", variable)) %>% 
    mutate(to = gsub("\\]", "", to)) -> movement_probs
  
  # distribute into a matrix
  # populate with means
  movement_probs %>% 
    dplyr::select(mean, from, to) %>% 
    pivot_wider(., names_from = to, values_from = mean) %>% 
    column_to_rownames("from") %>% 
    as.matrix() -> transition_matrix
  
  return(transition_matrix)
  
}

# 2) final fates function (if fish start in a certain state, what is their final distribution?)
non_upstream_states <- model_states[-grep("Creek Upstream|River Upstream", model_states)]

final_fates_simulation <- function(nsim, transition_matrix, start_state = 2){
  trans <- transition_matrix
  
  # Create an empty state matrix
  state_matrix <- matrix(rep(0, 29), nrow = 29, ncol = 1)
  
  transition_states <- non_upstream_states
  
  rownames(state_matrix) <- transition_states
  
  state_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("state") %>% 
    dplyr::rename(t1 = V1) -> state_matrix
  
  
  # start all individuals above Bonneville (release site in multistate model)
  state_matrix[start_state,2] <- nsim
  
  # create a final fate matrix
  final_fate_matrix <- matrix(rep(0, 29), nrow = 1, ncol = 29)
  colnames(final_fate_matrix) <- transition_states
  
  # run while loop until all fish have been lost
  i <- 2
  condition <- "not finished"
  while (condition != "finished"){
    # First, append a row for the next states
    state_matrix %>% 
      bind_cols(., t = rep(0, 29)) %>% 
      # dplyr::rename(paste0("t", i) = "newcol")
      dplyr::rename_with(t, .fn = ~paste0(., i), cols = c("t")) -> state_matrix
    
    # Now, calculate the probabilities from each starting state
    movement_matrix <- matrix(rep(0, 29*29), ncol = 29, nrow = 29)
    rownames(movement_matrix) <- transition_states
    colnames(movement_matrix) <- transition_states
    for (j in 1:nrow(state_matrix)){
      movement_matrix[,j] <- rmultinom(n = 1, size = state_matrix[j,i], prob = trans[j,])
    }
    # sum across rows, store in the state matrix
    state_matrix[,i+1] <- rowSums(movement_matrix)
    
    # add loss row to final fate
    final_fate_matrix %>% 
      rbind(., movement_matrix['loss',]) -> final_fate_matrix
    
    # increase counter
    i <- i + 1
    
    if (state_matrix[29, ncol(state_matrix)] == nsim){
      condition <- "finished"
    }
  }
  
  data.frame(colSums(final_fate_matrix[,1:28])) %>% 
    dplyr::rename(count = colSums.final_fate_matrix...1.28..) %>% 
    mutate(prop = count/nsim)-> final_fates
  
  outputs <- list(state_matrix, final_fates)
  # return(final_fates)
  return(outputs)
  
}

### Function for getting uncertainty using raw MCMC draws
# first, load upstream indices
upstream_indices <- grep("Creek Upstream|River Upstream", model_states)

transition_matrix_reformat_draws <- function(draws_object, variable_prefix, niter){
  draws_object %>% 
    dplyr::select(colnames(draws_object)[grepl(variable_prefix, colnames(draws_object))]) %>%  # keep the actual probabilities
    t() %>% 
    as.data.frame() %>% 
    rownames_to_column("variable") %>% 
    mutate(from = gsub(",.*$", "", variable)) %>% 
    mutate(from = gsub(paste0("1.", variable_prefix, "\\["), "", from)) %>% 
    mutate(to = gsub(paste0("1.", variable_prefix, "\\[.*,"), "", variable)) %>% 
    mutate(to = gsub("\\]", "", to)) -> movement_probs
  
  rownames(movement_probs) <- NULL
  
  # Then, turn each of these into a transition matrix
  # First, make a list to store these
  transition_matrix_list <- list()
  
  for (i in 1:niter){
    movement_probs %>% 
      as_tibble() %>% 
      dplyr::select(i+1, from, to) %>% 
      dplyr::rename(values = paste0(i)) %>% 
      pivot_wider(., names_from = to, values_from = values) %>% 
      column_to_rownames("from") %>% 
      as.matrix() -> transition_matrix
    
    # Remove the upstream states from each
    transition_matrix <- transition_matrix[-upstream_indices, -upstream_indices]
    
    # Manually change each so that transition probability from loss to loss is 1
    transition_matrix[29,29] <- 1
    
    # Fix 2022-12-07: Sometimes there are small rounding errors that lead to very very small negative probabilities (e.g., -2 * 10^-16).
    # Having negative probabilities leads to errors when using rmultinom.
    # Solution: Round any negative values up to zero
    transition_matrix[transition_matrix < 0] <- 0
    
    transition_matrix_list[[i]] <- transition_matrix
    
  }
  
  
  
  return(transition_matrix_list)
}

### Function to loop through all of these
final_fates_uncertainty <- function(transition_matrix_list, niter, start_state = 2){
  # For the first run, need to initialize data frames to store everything
  final_fates_summary <- final_fates_simulation(nsim = 1000000, start_state = start_state, transition_matrix = transition_matrix_list[[1]])[2][[1]]
  
  final_fates_summary %>% 
    rownames_to_column("state") %>% 
    dplyr::select(-count) %>% 
    dplyr::rename(!!quo_name(1) := prop) -> final_fates_summary
  
  # Now, loop through the remaining draws
  for (i in 2:niter){
    final_fates <- final_fates_simulation(nsim = 1000000, start_state = start_state, transition_matrix = transition_matrix_list[[i]])[2][[1]]
    
    # Bind this to the previous runs
    final_fates %>% 
      rownames_to_column("state") %>% 
      dplyr::select(-count) %>% 
      dplyr::rename(!!quo_name(i) := prop) -> final_fates
    
    final_fates_summary %>% 
      left_join(final_fates, by = "state") -> final_fates_summary
  }
  
  
  # Now, evaluate uncertainty
  final_fates_summary %>% 
    pivot_longer(., cols = colnames(final_fates_summary)[2:(niter + 1)], names_to = "draw") %>% 
  # Figure out the 95% CI
    group_by(state) %>% 
    summarise(prob = quantile(value, c(0.025, 0.5, 0.975)), q = c(0.025, 0.5, 0.975)) %>% 
    # mutate(percent = round(prob * 100,2)) %>% 
    # mutate(probability = format(prob, digits = 2)) %>% 
    # mutate(probability = prob) %>% 
    # dplyr::select(-prob) %>% 
    pivot_wider(names_from = q, values_from = prob) -> final_fates_quantiles
  
  return(final_fates_quantiles)
}

### Function to reformat data for plot
# ORDER THE STATES
states_order_for_plot <- gsub(" Mouth", "", model_states[-upstream_indices])
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



final_fate_plot_data_reformat <- function(final_fates_quantiles){
  final_fates_quantiles %>% 
    # pivot_longer(., cols = c(`0.025`, `0.5`, `0.975`))
    mutate(state = gsub(" Mouth", "", state)) %>% 
    mutate(state = fct_rev(factor(state, levels = states_order_for_plot))) %>% 
    arrange(desc(state)) %>% # reorder by factor levels
    # mutate(label = paste0(`0.5`, "% (", `0.025`, "% - ", `0.975`, "%)")) %>% 
    mutate(label = paste0(format(`0.5`, digits = 2, scientific = FALSE), " (", format(`0.025`, digits = 2, scientific = FALSE), 
                          " - ", format(`0.975`, digits = 2, scientific = FALSE), ")")) %>% 
    mutate(label.pos = 0) -> final_fates_quantiles_forplot
  
  return(final_fates_quantiles_forplot)
}


### Function for combined table + plot figure
final_fates_plot <- function(final_fates_quantiles_forplot, homing_state, downstream_tributaries,
                             upstream_tributaries){
  # Generate a vector of colors for different states and their meanings
  state_significance_colors <- ifelse(final_fates_quantiles_forplot$state == homing_state, "#33a02c", # homing
                                      ifelse(final_fates_quantiles_forplot$state %in% downstream_tributaries, "#ff7f00", # Tributaries downstream of home
                                             ifelse(final_fates_quantiles_forplot$state %in% upstream_tributaries, "#e31a1c", # Tributaries upstream of home
                                                    ifelse(grepl("River|Creek|tributaries", final_fates_quantiles_forplot$state), "#6a3d9a", # Tributaries out of the basin (this should be everything that's not upstream or downstream)
                                                        "black")))) # mainstem loss
  
  # Create the quantiles plot
  quantiles_plot <- ggplot(final_fates_quantiles_forplot, aes(x = state, y = `0.5`, ymin = `0.025`, ymax = `0.975`)) +
    geom_point() +
    geom_linerange() +
    coord_flip() +
    ylab("Final Distribution Probability") +
    xlab("Model State") +
    # ggtitle(" ") +
    theme(plot.title = element_text(size = 12),
          axis.text.y = element_text(color = rev(state_significance_colors)),
          axis.title = element_text(size = 14),
          axis.text.x = element_text(size = 12))
  
  
  
  # Create quantiles table
  quantiles_table <- ggplot(final_fates_quantiles_forplot, aes(x = state, y = label.pos, label = label)) +
    geom_text(hjust = 0) +
    coord_flip() +
    ylab(" ") +
    xlab(" ") +
    scale_y_continuous(lim = c(0, 0.05)) +
    # ggtitle("Percent") +
    theme(
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.title.y = element_blank(),
      axis.title.x = element_text(color = "white", size = 14),
      axis.text.x = element_text(color = "white", size = 12),
      axis.text.y = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.title = element_text(size = 12))
  
  # create a legend
  quantiles_legend <- data.frame(loss_type = c("Mainstem", "Downstream Tributaries", "Upstream Tributaries", "Out-of-Basin Tributaries", "Home"), x = c(1,2,3.3,4.6,5.5), y = 0)
  
  
  ggplot(quantiles_legend, aes(label = loss_type, x = x, y = y, color = loss_type)) +
    geom_text(size = 5) +
    scale_color_manual(values = c("Mainstem" = "black", "Downstream Tributaries" = "#ff7f00", "Upstream Tributaries" = "#e31a1c", 
                                  "Out-of-Basin Tributaries" = "#6a3d9a", "Home" = "#33a02c")) +
    scale_x_continuous(lim = c(0.75,6)) +
    theme(panel.background = element_rect(fill = "white"),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          plot.title = element_text(size = 16)) +
    ggtitle("State Legend:") -> quantiles_legend_plot
  
  # Arrange quantiles plot, table, and legend in one plot
  ggarrange(ggarrange(quantiles_plot, quantiles_table, nrow = 1, ncol = 2, widths = c(6, 2)),
            quantiles_legend_plot, nrow = 2, ncol = 1, heights = c(8, 1)) -> plot_for_report
  
  # return the final plot
  return(plot_for_report)
  
  
}






# Get the draws object for each of the three ESUs
as.data.frame(snake_fit$draws()) -> snake_draws
as.data.frame(middle_columbia_fit$draws()) -> middle_columbia_draws
as.data.frame(upper_columbia_fit$draws()) -> upper_columbia_draws

##### MIDDLE COLUMBIA, UPPER COLUMBIA, AND SNAKE - prepare data for plotting ####
print(Sys.time())
### Middle Columbia
# deschutes River
deschutes_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = middle_columbia_draws, variable_prefix = "origin1_probs_DE", niter = 200)
deschutes_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = deschutes_DE_transition_matrix_list, niter = 200)
deschutes_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = deschutes_final_fates_quantiles)

# john_day River
john_day_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = middle_columbia_draws, variable_prefix = "origin3_probs_DE", niter = 200)
john_day_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = john_day_DE_transition_matrix_list, niter = 200)
john_day_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = john_day_final_fates_quantiles)

# fifteenmile creek
fifteenmile_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = middle_columbia_draws, variable_prefix = "origin2_probs_DE", niter = 200)
fifteenmile_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = fifteenmile_DE_transition_matrix_list, niter = 200)
fifteenmile_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = fifteenmile_final_fates_quantiles)

# umatilla River
umatilla_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = middle_columbia_draws, variable_prefix = "origin4_probs_DE", niter = 200)
umatilla_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = umatilla_DE_transition_matrix_list, niter = 200)
umatilla_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = umatilla_final_fates_quantiles)

# yakima River
yakima_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = middle_columbia_draws, variable_prefix = "origin5_probs_DE", niter = 200)
yakima_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = yakima_DE_transition_matrix_list, niter = 200)
yakima_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = yakima_final_fates_quantiles)

# walla_walla River
walla_walla_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = middle_columbia_draws, variable_prefix = "origin6_probs_DE", niter = 200)
walla_walla_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = walla_walla_DE_transition_matrix_list, niter = 200)
walla_walla_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = walla_walla_final_fates_quantiles)

### Upper Columbia

# wenatchee River
wenatchee_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = upper_columbia_draws, variable_prefix = "origin1_probs_DE", niter = 200)
wenatchee_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = wenatchee_DE_transition_matrix_list, niter = 200)
wenatchee_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = wenatchee_final_fates_quantiles)

# entiat River
entiat_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = upper_columbia_draws, variable_prefix = "origin2_probs_DE", niter = 200)
entiat_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = entiat_DE_transition_matrix_list, niter = 200)
entiat_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = entiat_final_fates_quantiles)

# okanogan River
okanogan_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = upper_columbia_draws, variable_prefix = "origin3_probs_DE", niter = 200)
okanogan_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = okanogan_DE_transition_matrix_list, niter = 200)
okanogan_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = okanogan_final_fates_quantiles)

# methow River
methow_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = upper_columbia_draws, variable_prefix = "origin4_probs_DE", niter = 200)
methow_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = methow_DE_transition_matrix_list, niter = 200)
methow_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = methow_final_fates_quantiles)


### Snake River
# Tucannon River
tucannon_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = snake_draws, variable_prefix = "origin1_probs_DE", niter = 200)
tucannon_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = tucannon_DE_transition_matrix_list, niter = 200)
tucannon_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = tucannon_final_fates_quantiles)

# asotin Creek
asotin_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = snake_draws, variable_prefix = "origin2_probs_DE", niter = 200)
asotin_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = asotin_DE_transition_matrix_list, niter = 200)
asotin_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = asotin_final_fates_quantiles)

# clearwater River
clearwater_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = snake_draws, variable_prefix = "origin3_probs_DE", niter = 200)
clearwater_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = clearwater_DE_transition_matrix_list, niter = 200)
clearwater_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = clearwater_final_fates_quantiles)

# imnaha River
imnaha_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = snake_draws, variable_prefix = "origin6_probs_DE", niter = 200)
imnaha_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = imnaha_DE_transition_matrix_list, niter = 200)
imnaha_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = imnaha_final_fates_quantiles)

# salmon River
salmon_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = snake_draws, variable_prefix = "origin4_probs_DE", niter = 200)
salmon_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = salmon_DE_transition_matrix_list, niter = 200)
salmon_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = salmon_final_fates_quantiles)

# grande_ronde River
grande_ronde_DE_transition_matrix_list <- transition_matrix_reformat_draws(draws_object = snake_draws, variable_prefix = "origin5_probs_DE", niter = 200)
grande_ronde_final_fates_quantiles <- final_fates_uncertainty(transition_matrix_list = grande_ronde_DE_transition_matrix_list, niter = 200)
grande_ronde_final_fates_quantiles_forplot <- final_fate_plot_data_reformat(final_fates_quantiles = grande_ronde_final_fates_quantiles)

print(Sys.time())




##### MIDDLE COLUMBIA, UPPER COLUMBIA, AND SNAKE - plot final fates ####

### Middle Columbia
deschutes_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = deschutes_final_fates_quantiles_forplot, homing_state = "Deschutes River", 
                                             downstream_tributaries = c("Hood River", "Fifteenmile Creek"),
                                             upstream_tributaries = c("John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                      "Yakima River", "Walla Walla River",
                                                                      "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                      "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "deschutes_fates_plot.pdf"), deschutes_plot_for_report, height = 8, width = 10)
ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "deschutes_fates_plot.pdf"), deschutes_plot_for_report, height = 8, width = 10)

john_day_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = john_day_final_fates_quantiles_forplot, homing_state = "John Day River", 
                                             downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River"),
                                             upstream_tributaries = c("BON to MCN other tributaries", "Umatilla River",
                                                                      "Yakima River", "Walla Walla River",
                                                                      "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                      "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "john_day_fates_plot.pdf"), john_day_plot_for_report, height = 8, width = 10)

fifteenmile_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = fifteenmile_final_fates_quantiles_forplot, homing_state = "Fifteenmile Creek", 
                                                downstream_tributaries = c("Hood River"),
                                                upstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                         "Yakima River", "Walla Walla River",
                                                                         "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                         "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "fifteenmile_fates_plot.pdf"), fifteenmile_plot_for_report, height = 8, width = 10)

umatilla_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = umatilla_final_fates_quantiles_forplot, homing_state = "Umatilla River", 
                                             downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River"),
                                             upstream_tributaries = c("BON to MCN other tributaries", 
                                                                      "Yakima River", "Walla Walla River",
                                                                      "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                      "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "umatilla_fates_plot.pdf"), umatilla_plot_for_report, height = 8, width = 10)

walla_walla_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = walla_walla_final_fates_quantiles_forplot, homing_state = "Walla Walla River", 
                                                downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River", "BON to MCN other tributaries"),
                                                upstream_tributaries = c("Yakima River", "Walla Walla River",
                                                                         "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                         "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "walla_walla_fates_plot.pdf"), walla_walla_plot_for_report, height = 8, width = 10)

yakima_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = yakima_final_fates_quantiles_forplot, homing_state = "Yakima River", 
                                           downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River", "BON to MCN other tributaries", "Walla Walla River"),
                                           upstream_tributaries = c("Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                    "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "yakima_fates_plot.pdf"), yakima_plot_for_report, height = 8, width = 10)




### Upper Columbia
wenatchee_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = wenatchee_final_fates_quantiles_forplot, homing_state = "Wenatchee River", 
                                             downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                        "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Yakima River"),
                                             upstream_tributaries = c("Entiat River", "Methow River", "Okanogan River", "Upstream WEL other tributaries"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "wenatchee_fates_plot.pdf"), wenatchee_plot_for_report, height = 8, width = 10)

entiat_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = entiat_final_fates_quantiles_forplot, homing_state = "Entiat River", 
                                           downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                      "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Yakima River",
                                                                      "Wenatchee River"),
                                           upstream_tributaries = c("Methow River", "Okanogan River", "Upstream WEL other tributaries"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "entiat_fates_plot.pdf"), entiat_plot_for_report, height = 8, width = 10)

methow_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = methow_final_fates_quantiles_forplot, homing_state = "Methow River", 
                                           downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                      "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Yakima River",
                                                                      "Methow River", "Wenatchee River", "Entiat River"),
                                           upstream_tributaries = c("Okanogan River", "Upstream WEL other tributaries"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "methow_fates_plot.pdf"), methow_plot_for_report, height = 8, width = 10)

okanogan_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = okanogan_final_fates_quantiles_forplot, homing_state = "Okanogan River", 
                                             downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                        "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Yakima River",
                                                                        "Methow River", "Wenatchee River", "Methow River", "Entiat River"),
                                             upstream_tributaries = c("Upstream WEL other tributaries"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "okanogan_fates_plot.pdf"), okanogan_plot_for_report, height = 8, width = 10)


### Snake
tucannon_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = tucannon_final_fates_quantiles_forplot, homing_state = "Tucannon River", 
                                             downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                        "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River"),
                                             upstream_tributaries = c("Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "tucannon_fates_plot.pdf"), tucannon_plot_for_report, height = 8, width = 10)

asotin_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = asotin_final_fates_quantiles_forplot, homing_state = "Asotin Creek", 
                                             downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                        "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Tucannon River",
                                                                        "Clearwater River"),
                                             upstream_tributaries = c("Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "asotin_fates_plot.pdf"), asotin_plot_for_report, height = 8, width = 10)

clearwater_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = clearwater_final_fates_quantiles_forplot, homing_state = "Clearwater River", 
                                             downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                         "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Tucannon River"),
                                             upstream_tributaries = c("Imnaha River", "Grande Ronde River", "Salmon River", "Asotin Creek"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "clearwater_fates_plot.pdf"), clearwater_plot_for_report, height = 8, width = 10)

salmon_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = salmon_final_fates_quantiles_forplot, homing_state = "Salmon River", 
                                           downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                      "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Tucannon River",
                                                                      "Asotin Creek", "Clearwater River", "Grande Ronde River"),
                                           upstream_tributaries = c("Imnaha River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "salmon_fates_plot.pdf"), salmon_plot_for_report, height = 8, width = 10)

grande_ronde_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = grande_ronde_final_fates_quantiles_forplot, homing_state = "Grande Ronde River", 
                                                 downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                            "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Tucannon River",
                                                                            "Asotin Creek", "Clearwater River"),
                                                 upstream_tributaries = c("Imnaha River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "grande_ronde_fates_plot.pdf"), grande_ronde_plot_for_report, height = 8, width = 10)

imnaha_plot_for_report <- final_fates_plot(final_fates_quantiles_forplot = imnaha_final_fates_quantiles_forplot, homing_state = "Imnaha River", 
                                           downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                      "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Tucannon River",
                                                                      "Asotin Creek", "Clearwater River", "Salmon River", "Grande Ronde River"),
                                           upstream_tributaries = NA)

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_fates_plots", "imnaha_fates_plot.pdf"), imnaha_plot_for_report, height = 8, width = 10)

transition_matrix_reformat <- function(movement_probs, variable_prefix){
    movement_probs %>% 
      mutate(from = gsub(",.*$", "", variable)) %>% 
      mutate(from = gsub(paste0(variable_prefix, "\\["), "", from)) %>% 
      mutate(to = gsub(paste0(variable_prefix, "\\[.*,"), "", variable)) %>% 
      mutate(to = gsub("\\]", "", to)) -> movement_probs
    
    # distribute into a matrix
    # populate with means
    movement_probs %>% 
      dplyr::select(mean, from, to) %>% 
      pivot_wider(., names_from = to, values_from = mean) %>% 
      column_to_rownames("from") %>% 
      as.matrix() -> transition_matrix
    
    
    return(transition_matrix)
    
}

  
subset(snake_summary, grepl("origin1_probs_DE", variable)) -> tucannon_DE_movement_probs
subset(snake_summary, grepl("origin2_probs_DE", variable)) -> asotin_DE_movement_probs
subset(snake_summary, grepl("origin3_probs_DE", variable)) -> clearwater_DE_movement_probs
subset(snake_summary, grepl("origin4_probs_DE", variable)) -> salmon_DE_movement_probs
subset(snake_summary, grepl("origin5_probs_DE", variable)) -> grande_ronde_DE_movement_probs
subset(snake_summary, grepl("origin6_probs_DE", variable)) -> imnaha_DE_movement_probs

transition_matrix_reformat(movement_probs = tucannon_DE_movement_probs, variable_prefix = "origin1_probs_DE") -> tucannon_DE_transition_matrix
transition_matrix_reformat(movement_probs = asotin_DE_movement_probs, variable_prefix = "origin2_probs_DE") -> asotin_DE_transition_matrix
transition_matrix_reformat(movement_probs = clearwater_DE_movement_probs, variable_prefix = "origin3_probs_DE") -> clearwater_DE_transition_matrix
transition_matrix_reformat(movement_probs = salmon_DE_movement_probs, variable_prefix = "origin4_probs_DE") -> salmon_DE_transition_matrix
transition_matrix_reformat(movement_probs = grande_ronde_DE_movement_probs, variable_prefix = "origin5_probs_DE") -> grande_ronde_DE_transition_matrix
transition_matrix_reformat(movement_probs = imnaha_DE_movement_probs, variable_prefix = "origin6_probs_DE") -> imnaha_DE_transition_matrix

# Take out the upstream states from each of these
upstream_indices <- grep("Creek Upstream|River Upstream", model_states)

tucannon_DE_transition_matrix <- tucannon_DE_transition_matrix[-upstream_indices, -upstream_indices]
asotin_DE_transition_matrix <- asotin_DE_transition_matrix[-upstream_indices, -upstream_indices]
clearwater_DE_transition_matrix <- clearwater_DE_transition_matrix[-upstream_indices, -upstream_indices]
salmon_DE_transition_matrix <- salmon_DE_transition_matrix[-upstream_indices, -upstream_indices]
grande_ronde_DE_transition_matrix <- grande_ronde_DE_transition_matrix[-upstream_indices, -upstream_indices]
imnaha_DE_transition_matrix <- imnaha_DE_transition_matrix[-upstream_indices, -upstream_indices]

# Manually change each so that transition probability from loss to loss is 1
tucannon_DE_transition_matrix[29,29] <- 1
asotin_DE_transition_matrix[29,29] <- 1
clearwater_DE_transition_matrix[29,29] <- 1
salmon_DE_transition_matrix[29,29] <- 1
grande_ronde_DE_transition_matrix[29,29] <- 1
imnaha_DE_transition_matrix[29,29] <- 1

# Get final fates of fish for each origin, starting at BON

tucannon_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = tucannon_DE_transition_matrix)[2][[1]]
asotin_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = asotin_DE_transition_matrix)[2][[1]]
clearwater_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = clearwater_DE_transition_matrix)[2][[1]]
salmon_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = salmon_DE_transition_matrix)[2][[1]]
grande_ronde_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = grande_ronde_DE_transition_matrix)[2][[1]]
imnaha_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = imnaha_DE_transition_matrix)[2][[1]]

##### MIDDLE COLUMBIA - Calculate final fates for different origins ####

subset(middle_columbia_summary, grepl("origin1_probs_DE", variable)) -> deschutes_DE_movement_probs
subset(middle_columbia_summary, grepl("origin2_probs_DE", variable)) -> fifteenmile_DE_movement_probs
subset(middle_columbia_summary, grepl("origin3_probs_DE", variable)) -> john_day_DE_movement_probs
subset(middle_columbia_summary, grepl("origin4_probs_DE", variable)) -> umatilla_DE_movement_probs
subset(middle_columbia_summary, grepl("origin5_probs_DE", variable)) -> yakima_DE_movement_probs
subset(middle_columbia_summary, grepl("origin6_probs_DE", variable)) -> walla_walla_DE_movement_probs

transition_matrix_reformat(movement_probs = deschutes_DE_movement_probs, variable_prefix = "origin1_probs_DE") -> deschutes_DE_transition_matrix
transition_matrix_reformat(movement_probs = fifteenmile_DE_movement_probs, variable_prefix = "origin2_probs_DE") -> fifteenmile_DE_transition_matrix
transition_matrix_reformat(movement_probs = john_day_DE_movement_probs, variable_prefix = "origin3_probs_DE") -> john_day_DE_transition_matrix
transition_matrix_reformat(movement_probs = umatilla_DE_movement_probs, variable_prefix = "origin4_probs_DE") -> umatilla_DE_transition_matrix
transition_matrix_reformat(movement_probs = yakima_DE_movement_probs, variable_prefix = "origin5_probs_DE") -> yakima_DE_transition_matrix
transition_matrix_reformat(movement_probs = walla_walla_DE_movement_probs, variable_prefix = "origin6_probs_DE") -> walla_walla_DE_transition_matrix

# Take out the upstream states from each of these
upstream_indices <- grep("Creek Upstream|River Upstream", model_states)

deschutes_DE_transition_matrix <- deschutes_DE_transition_matrix[-upstream_indices, -upstream_indices]
fifteenmile_DE_transition_matrix <- fifteenmile_DE_transition_matrix[-upstream_indices, -upstream_indices]
john_day_DE_transition_matrix <- john_day_DE_transition_matrix[-upstream_indices, -upstream_indices]
umatilla_DE_transition_matrix <- umatilla_DE_transition_matrix[-upstream_indices, -upstream_indices]
yakima_DE_transition_matrix <- yakima_DE_transition_matrix[-upstream_indices, -upstream_indices]
walla_walla_DE_transition_matrix <- walla_walla_DE_transition_matrix[-upstream_indices, -upstream_indices]

# Manually change each so that transition probability from loss to loss is 1
deschutes_DE_transition_matrix[29,29] <- 1
fifteenmile_DE_transition_matrix[29,29] <- 1
john_day_DE_transition_matrix[29,29] <- 1
umatilla_DE_transition_matrix[29,29] <- 1
yakima_DE_transition_matrix[29,29] <- 1
walla_walla_DE_transition_matrix[29,29] <- 1

# Get final fates of fish for each origin, starting at BON

deschutes_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = deschutes_DE_transition_matrix)[2][[1]]
fifteenmile_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = fifteenmile_DE_transition_matrix)[2][[1]]
john_day_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = john_day_DE_transition_matrix)[2][[1]]
umatilla_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = umatilla_DE_transition_matrix)[2][[1]]
yakima_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = yakima_DE_transition_matrix)[2][[1]]
walla_walla_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = walla_walla_DE_transition_matrix)[2][[1]]

##### UPPER COLUMBIA - Calculate final fates for different origins ####

subset(upper_columbia_summary, grepl("origin1_probs_DE", variable)) -> wenatchee_DE_movement_probs
subset(upper_columbia_summary, grepl("origin2_probs_DE", variable)) -> entiat_DE_movement_probs
subset(upper_columbia_summary, grepl("origin3_probs_DE", variable)) -> okanogan_DE_movement_probs
subset(upper_columbia_summary, grepl("origin4_probs_DE", variable)) -> methow_DE_movement_probs

transition_matrix_reformat(movement_probs = wenatchee_DE_movement_probs, variable_prefix = "origin1_probs_DE") -> wenatchee_DE_transition_matrix
transition_matrix_reformat(movement_probs = entiat_DE_movement_probs, variable_prefix = "origin2_probs_DE") -> entiat_DE_transition_matrix
transition_matrix_reformat(movement_probs = okanogan_DE_movement_probs, variable_prefix = "origin3_probs_DE") -> okanogan_DE_transition_matrix
transition_matrix_reformat(movement_probs = methow_DE_movement_probs, variable_prefix = "origin4_probs_DE") -> methow_DE_transition_matrix

# Take out the upstream states from each of these
upstream_indices <- grep("Creek Upstream|River Upstream", model_states)

wenatchee_DE_transition_matrix <- wenatchee_DE_transition_matrix[-upstream_indices, -upstream_indices]
entiat_DE_transition_matrix <- entiat_DE_transition_matrix[-upstream_indices, -upstream_indices]
okanogan_DE_transition_matrix <- okanogan_DE_transition_matrix[-upstream_indices, -upstream_indices]
methow_DE_transition_matrix <- methow_DE_transition_matrix[-upstream_indices, -upstream_indices]

# Manually change each so that transition probability from loss to loss is 1
wenatchee_DE_transition_matrix[29,29] <- 1
entiat_DE_transition_matrix[29,29] <- 1
okanogan_DE_transition_matrix[29,29] <- 1
methow_DE_transition_matrix[29,29] <- 1

# Get final fates of fish for each origin, starting at BON

wenatchee_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = wenatchee_DE_transition_matrix)[2][[1]]
entiat_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = entiat_DE_transition_matrix)[2][[1]]
okanogan_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = okanogan_DE_transition_matrix)[2][[1]]
methow_final_fates <- final_fates_simulation(nsim = 1000000, transition_matrix = methow_DE_transition_matrix)[2][[1]]


##### FINAL FATES/DISTRIBUTION - where do fish end up? Export for report #####
final_fates_reformat <- function(final_fates, trib_name){
  final_fates %>% 
    rownames_to_column("state") %>% 
    mutate(percent = round(prop*100,2)) %>% 
    dplyr::select(state,percent) %>% 
    dplyr::rename(!!quo_name(trib_name) := percent) -> final_fates_table
  
  return(final_fates_table)
}


# Middle Columbia
final_fates_reformat(final_fates = deschutes_final_fates, trib_name = "Deschutes River") -> deschutes_final_distribution
final_fates_reformat(final_fates = fifteenmile_final_fates, trib_name = "Fifteenmile Creek") -> fifteenmile_final_distribution
final_fates_reformat(final_fates = john_day_final_fates, trib_name = "John Day River") -> john_day_final_distribution
final_fates_reformat(final_fates = umatilla_final_fates, trib_name = "Umatila River") -> umatilla_final_distribution
final_fates_reformat(final_fates = yakima_final_fates, trib_name = "Yakima River") -> yakima_final_distribution
final_fates_reformat(final_fates = walla_walla_final_fates, trib_name = "Walla Walla River") -> walla_walla_final_distribution

# Upper Columbia
final_fates_reformat(final_fates = wenatchee_final_fates, trib_name = "Wenatchee River") -> wenatchee_final_distribution
final_fates_reformat(final_fates = entiat_final_fates, trib_name = "Entiat River") -> entiat_final_distribution
final_fates_reformat(final_fates = okanogan_final_fates, trib_name = "Okanogan River") -> okanogan_final_distribution
final_fates_reformat(final_fates = methow_final_fates, trib_name = "Methow River") -> methow_final_distribution

# Snake River
final_fates_reformat(final_fates = tucannon_final_fates, trib_name = "Tucannon River") -> tucannon_final_distribution
final_fates_reformat(final_fates = asotin_final_fates, trib_name = "Asotin Creek") -> asotin_final_distribution
final_fates_reformat(final_fates = clearwater_final_fates, trib_name = "Clearwater River") -> clearwater_final_distribution
final_fates_reformat(final_fates = salmon_final_fates, trib_name = "Salmon River") -> salmon_final_distribution
final_fates_reformat(final_fates = grande_ronde_final_fates, trib_name = "Grande Ronde River") -> grande_ronde_final_distribution
final_fates_reformat(final_fates = imnaha_final_fates, trib_name = "Imnaha River") -> imnaha_final_distribution

# Join these all together for export
deschutes_final_distribution %>% 
  left_join(., fifteenmile_final_distribution, by = "state") %>% 
  left_join(., john_day_final_distribution, by = "state") %>% 
  left_join(., umatilla_final_distribution, by = "state") %>% 
  left_join(., yakima_final_distribution, by = "state") %>% 
  left_join(., walla_walla_final_distribution, by = "state") %>% 
  left_join(., wenatchee_final_distribution, by = "state") %>% 
  left_join(., entiat_final_distribution, by = "state") %>% 
  left_join(., okanogan_final_distribution, by = "state") %>% 
  left_join(., methow_final_distribution, by = "state") %>% 
  left_join(., tucannon_final_distribution, by = "state") %>% 
  left_join(., asotin_final_distribution, by = "state") %>% 
  left_join(., clearwater_final_distribution, by = "state") %>% 
  left_join(., salmon_final_distribution, by = "state") %>% 
  left_join(., grande_ronde_final_distribution, by = "state") %>% 
  left_join(., imnaha_final_distribution, by = "state") %>% 
  mutate(state = gsub(" Mouth", "", state)) -> all_origins_final_distribution

write.csv(all_origins_final_distribution, here::here("stan_actual", "deteff_ESU_models", "output_tables", "final_distribution.csv"))



##### Get final fates of fish - conditional on overshoot #####

# Populations that we can measure overshoot for: 
# 1) AT MCN: Deschutes River, John Day River, Hood River, Fifteenmile Creek, Umatilla River
# 2) At PRA or ICH: Yakima River and Walla Walla River
# 3) At RRE: Wenatchee River
# 4) At WEL: Entiat River
# 5) At LGR: Tucannon River

# Populations that we can't measure overshoot for:
# Okanogan River, Methow River, Asotin Creek, Clearwater River, Salmon River, Grande Ronde River

# In each case, we will compare the final fate of fish that are in the immediate overshoot state to the final fate
# that are in the correct part of the mainstem. This second probability for comparison does necessarily include
# some fish that will overshoot, but it's our best comparison. Just need to explain in the text.

### Middle Columbia

deschutes_final_fates_quantiles_overshoot <- final_fates_uncertainty(transition_matrix_list = deschutes_DE_transition_matrix_list, niter = 200, start_state = 3)
deschutes_final_fates_quantiles_forplot_overshoot <- final_fate_plot_data_reformat(final_fates_quantiles = deschutes_final_fates_quantiles_overshoot)
deschutes_overshoot_fates <- final_fates_plot(final_fates_quantiles_forplot = deschutes_final_fates_quantiles_forplot_overshoot, homing_state = "Deschutes River", 
                                              downstream_tributaries = c("Hood River", "Fifteenmile Creek"),
                                              upstream_tributaries = c("John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                       "Yakima River", "Walla Walla River",
                                                                       "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                       "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "deschutes_overshoot_fates.pdf"), deschutes_overshoot_fates, height = 8, width = 10)

fifteenmile_final_fates_quantiles_overshoot <- final_fates_uncertainty(transition_matrix_list = fifteenmile_DE_transition_matrix_list, niter = 200, start_state = 3)
fifteenmile_final_fates_quantiles_forplot_overshoot <- final_fate_plot_data_reformat(final_fates_quantiles = fifteenmile_final_fates_quantiles_overshoot)
fifteenmile_overshoot_fates <- final_fates_plot(final_fates_quantiles_forplot = fifteenmile_final_fates_quantiles_forplot_overshoot, homing_state = "Fifteenmile Creek", 
                                              downstream_tributaries = c("Hood River"),
                                              upstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                       "Yakima River", "Walla Walla River",
                                                                       "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                       "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "fifteenmile_overshoot_fates.pdf"), fifteenmile_overshoot_fates, height = 8, width = 10)

john_day_final_fates_quantiles_overshoot <- final_fates_uncertainty(transition_matrix_list = john_day_DE_transition_matrix_list, niter = 200, start_state = 3)
john_day_final_fates_quantiles_forplot_overshoot <- final_fate_plot_data_reformat(final_fates_quantiles = john_day_final_fates_quantiles_overshoot)
john_day_overshoot_fates <- final_fates_plot(final_fates_quantiles_forplot = john_day_final_fates_quantiles_forplot_overshoot, homing_state = "John Day River", 
                                              downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River"),
                                              upstream_tributaries = c("BON to MCN other tributaries", "Umatilla River",
                                                                       "Yakima River", "Walla Walla River",
                                                                       "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                       "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "john_day_overshoot_fates.pdf"), john_day_overshoot_fates, height = 8, width = 10)

umatilla_final_fates_quantiles_overshoot <- final_fates_uncertainty(transition_matrix_list = umatilla_DE_transition_matrix_list, niter = 200, start_state = 3)
umatilla_final_fates_quantiles_forplot_overshoot <- final_fate_plot_data_reformat(final_fates_quantiles = umatilla_final_fates_quantiles_overshoot)
umatilla_overshoot_fates <- final_fates_plot(final_fates_quantiles_forplot = umatilla_final_fates_quantiles_forplot_overshoot, homing_state = "Umatilla River", 
                                              downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River"),
                                              upstream_tributaries = c("BON to MCN other tributaries", "Umatilla River",
                                                                       "Yakima River", "Walla Walla River",
                                                                       "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                       "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "umatilla_overshoot_fates.pdf"), umatilla_overshoot_fates, height = 8, width = 10)

# Note - Walla Walla can overshoot at two places (PRA and ICH). Create two plots for each
walla_walla_final_fates_quantiles_overshoot_PRA <- final_fates_uncertainty(transition_matrix_list = walla_walla_DE_transition_matrix_list, niter = 200, start_state = 4)
walla_walla_final_fates_quantiles_forplot_overshoot_PRA <- final_fate_plot_data_reformat(final_fates_quantiles = walla_walla_final_fates_quantiles_overshoot_PRA)
walla_walla_overshoot_PRA_fates <- final_fates_plot(final_fates_quantiles_forplot = walla_walla_final_fates_quantiles_forplot_overshoot_PRA, homing_state = "Walla Walla River", 
                                              downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River"),
                                              upstream_tributaries = c("BON to MCN other tributaries",
                                                                       "Yakima River",
                                                                       "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                       "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "walla_walla_overshoot_PRA_fates.pdf"), walla_walla_overshoot_PRA_fates, height = 8, width = 10)

walla_walla_final_fates_quantiles_overshoot_ICH <- final_fates_uncertainty(transition_matrix_list = walla_walla_DE_transition_matrix_list, niter = 200, start_state = 8)
walla_walla_final_fates_quantiles_forplot_overshoot_ICH <- final_fate_plot_data_reformat(final_fates_quantiles = walla_walla_final_fates_quantiles_overshoot_ICH)
walla_walla_overshoot_ICH_fates <- final_fates_plot(final_fates_quantiles_forplot = walla_walla_final_fates_quantiles_forplot_overshoot_ICH, homing_state = "Walla Walla River", 
                                                    downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River"),
                                                    upstream_tributaries = c("BON to MCN other tributaries",
                                                                             "Yakima River",
                                                                             "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                             "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "walla_walla_overshoot_ICH_fates.pdf"), walla_walla_overshoot_ICH_fates, height = 8, width = 10)

yakima_final_fates_quantiles_overshoot_PRA <- final_fates_uncertainty(transition_matrix_list = yakima_DE_transition_matrix_list, niter = 200, start_state = 4)
yakima_final_fates_quantiles_forplot_overshoot_PRA <- final_fate_plot_data_reformat(final_fates_quantiles = yakima_final_fates_quantiles_overshoot_PRA)
yakima_overshoot_PRA_fates <- final_fates_plot(final_fates_quantiles_forplot = yakima_final_fates_quantiles_forplot_overshoot_PRA, homing_state = "Yakima River", 
                                           downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River", "Walla Wala River"),
                                           upstream_tributaries = c("BON to MCN other tributaries",
                                                                    "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                    "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "yakima_overshoot_PRA_fates.pdf"), yakima_overshoot_PRA_fates, height = 8, width = 10)

yakima_final_fates_quantiles_overshoot_ICH <- final_fates_uncertainty(transition_matrix_list = yakima_DE_transition_matrix_list, niter = 200, start_state = 8)
yakima_final_fates_quantiles_forplot_overshoot_ICH <- final_fate_plot_data_reformat(final_fates_quantiles = yakima_final_fates_quantiles_overshoot_ICH)
yakima_overshoot_ICH_fates <- final_fates_plot(final_fates_quantiles_forplot = yakima_final_fates_quantiles_forplot_overshoot_ICH, homing_state = "Yakima River", 
                                               downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River", "Walla Wala River"),
                                               upstream_tributaries = c("BON to MCN other tributaries",
                                                                        "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                        "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "yakima_overshoot_ICH_fates.pdf"), yakima_overshoot_ICH_fates, height = 8, width = 10)

### Upper Columbia
wenatchee_final_fates_quantiles_overshoot <- final_fates_uncertainty(transition_matrix_list = wenatchee_DE_transition_matrix_list, niter = 200, start_state = 6)
wenatchee_final_fates_quantiles_forplot_overshoot <- final_fate_plot_data_reformat(final_fates_quantiles = wenatchee_final_fates_quantiles_overshoot)
wenatchee_overshoot_fates <- final_fates_plot(final_fates_quantiles_forplot = wenatchee_final_fates_quantiles_forplot_overshoot, homing_state = "Deschutes River", 
                                              downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                         "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Yakima River"),
                                              upstream_tributaries = c("Entiat River", "Methow River", "Okanogan River", "Upstream WEL other tributaries"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "wenatchee_overshoot_fates_plot.pdf"), wenatchee_overshoot_fates, height = 8, width = 10)

entiat_final_fates_quantiles_overshoot <- final_fates_uncertainty(transition_matrix_list = entiat_DE_transition_matrix_list, niter = 200, start_state = 7)
entiat_final_fates_quantiles_forplot_overshoot <- final_fate_plot_data_reformat(final_fates_quantiles = entiat_final_fates_quantiles_overshoot)
entiat_overshoot_fates <- final_fates_plot(final_fates_quantiles_forplot = entiat_final_fates_quantiles_forplot_overshoot, homing_state = "Deschutes River", 
                                           downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                      "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Yakima River",
                                                                      "Wenatchee River"),
                                           upstream_tributaries = c("Methow River", "Okanogan River", "Upstream WEL other tributaries"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "entiat_overshoot_fates_plot.pdf"), entiat_overshoot_fates, height = 8, width = 10)

### Snake

tucannon_final_fates_quantiles_overshoot <- final_fates_uncertainty(transition_matrix_list = tucannon_DE_transition_matrix_list, niter = 200, start_state = 9)
tucannon_final_fates_quantiles_forplot_overshoot <- final_fate_plot_data_reformat(final_fates_quantiles = tucannon_final_fates_quantiles_overshoot)
tucannon_overshoot_fates <- final_fates_plot(final_fates_quantiles_forplot = tucannon_final_fates_quantiles_forplot_overshoot, homing_state = "Deschutes River", 
                                             downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                        "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River"),
                                             upstream_tributaries = c("Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_fates_plots", "tucannon_overshoot_fates_plot.pdf"), tucannon_overshoot_fates, height = 8, width = 10)

##### Get final fates of fish - conditional on being in the correct mainstem state #####

# These values will be compared to fish in the overshoot state

### Middle Columbia

deschutes_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = deschutes_DE_transition_matrix_list, niter = 200, start_state = 2)
deschutes_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = deschutes_final_fates_quantiles_correct_mainstem)
deschutes_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = deschutes_final_fates_quantiles_forplot_correct_mainstem, homing_state = "Deschutes River", 
                                                     downstream_tributaries = c("Hood River", "Fifteenmile Creek"),
                                                     upstream_tributaries = c("John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                              "Yakima River", "Walla Walla River",
                                                                              "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                              "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "deschutes_correct_mainstem_fates.pdf"), deschutes_correct_mainstem_fates, height = 8, width = 10)

fifteenmile_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = fifteenmile_DE_transition_matrix_list, niter = 200, start_state = 2)
fifteenmile_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = fifteenmile_final_fates_quantiles_correct_mainstem)
fifteenmile_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = fifteenmile_final_fates_quantiles_forplot_correct_mainstem, homing_state = "Fifteenmile Creek", 
                                                       downstream_tributaries = c("Hood River"),
                                                       upstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                                "Yakima River", "Walla Walla River",
                                                                                "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                                "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "fifteenmile_correct_mainstem_fates.pdf"), fifteenmile_correct_mainstem_fates, height = 8, width = 10)

john_day_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = john_day_DE_transition_matrix_list, niter = 200, start_state = 2)
john_day_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = john_day_final_fates_quantiles_correct_mainstem)
john_day_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = john_day_final_fates_quantiles_forplot_correct_mainstem, homing_state = "John Day River", 
                                                    downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River"),
                                                    upstream_tributaries = c("BON to MCN other tributaries", "Umatilla River",
                                                                             "Yakima River", "Walla Walla River",
                                                                             "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                             "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "john_day_correct_mainstem_fates.pdf"), john_day_correct_mainstem_fates, height = 8, width = 10)

umatilla_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = umatilla_DE_transition_matrix_list, niter = 200, start_state = 2)
umatilla_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = umatilla_final_fates_quantiles_correct_mainstem)
umatilla_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = umatilla_final_fates_quantiles_forplot_correct_mainstem, homing_state = "Umatilla River", 
                                                    downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River"),
                                                    upstream_tributaries = c("BON to MCN other tributaries", "Umatilla River",
                                                                             "Yakima River", "Walla Walla River",
                                                                             "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                             "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "umatilla_correct_mainstem_fates.pdf"), umatilla_correct_mainstem_fates, height = 8, width = 10)

walla_walla_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = walla_walla_DE_transition_matrix_list, niter = 200, start_state = 3)
walla_walla_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = walla_walla_final_fates_quantiles_correct_mainstem)
walla_walla_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = walla_walla_final_fates_quantiles_forplot_correct_mainstem, homing_state = "Walla Walla River", 
                                                       downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River"),
                                                       upstream_tributaries = c("BON to MCN other tributaries",
                                                                                "Yakima River",
                                                                                "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                                "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "walla_walla_correct_mainstem_fates.pdf"), walla_walla_correct_mainstem_fates, height = 8, width = 10)


yakima_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = yakima_DE_transition_matrix_list, niter = 200, start_state = 3)
yakima_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = yakima_final_fates_quantiles_correct_mainstem)
yakima_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = yakima_final_fates_quantiles_forplot_correct_mainstem, homing_state = "Yakima River", 
                                                  downstream_tributaries = c("Hood River", "Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River", "Walla Wala River"),
                                                  upstream_tributaries = c("BON to MCN other tributaries",
                                                                           "Wenatchee River", "Entiat River", "Okanogan River", "Methow River", "Upstream WEL other tributaries",
                                                                           "Tucannon River", "Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "yakima_correct_mainstem_fates.pdf"), yakima_correct_mainstem_fates, height = 8, width = 10)


### Upper Columbia
wenatchee_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = wenatchee_DE_transition_matrix_list, niter = 200, start_state = 5)
wenatchee_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = wenatchee_final_fates_quantiles_correct_mainstem)
wenatchee_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = wenatchee_final_fates_quantiles_forplot_correct_mainstem, homing_state = "Deschutes River", 
                                                     downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                                "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Yakima River"),
                                                     upstream_tributaries = c("Entiat River", "Methow River", "Okanogan River", "Upstream WEL other tributaries"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "wenatchee_correct_mainstem_fates_plot.pdf"), entiat_correct_mainstem_fates, height = 8, width = 10)

entiat_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = entiat_DE_transition_matrix_list, niter = 200, start_state = 6)
entiat_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = entiat_final_fates_quantiles_correct_mainstem)
entiat_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = entiat_final_fates_quantiles_forplot_correct_mainstem, homing_state = "Deschutes River", 
                                                  downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                             "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River", "Yakima River",
                                                                             "Wenatchee River"),
                                                  upstream_tributaries = c("Methow River", "Okanogan River", "Upstream WEL other tributaries"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "entiat_correct_mainstem_fates_plot.pdf"), entiat_correct_mainstem_fates, height = 8, width = 10)

### Snake

tucannon_final_fates_quantiles_correct_mainstem <- final_fates_uncertainty(transition_matrix_list = tucannon_DE_transition_matrix_list, niter = 200, start_state = 8)
tucannon_final_fates_quantiles_forplot_correct_mainstem <- final_fate_plot_data_reformat(final_fates_quantiles = tucannon_final_fates_quantiles_correct_mainstem)
tucannon_correct_mainstem_fates <- final_fates_plot(final_fates_quantiles_forplot = tucannon_final_fates_quantiles_forplot_correct_mainstem, homing_state = "Deschutes River", 
                                                    downstream_tributaries = c("Deschutes River", "John Day River", "BON to MCN other tributaries", "Umatilla River",
                                                                               "Hood River", "Fifteenmile Creek", "Yakima River", "Walla Walla River"),
                                                    upstream_tributaries = c("Asotin Creek", "Clearwater River", "Imnaha River", "Grande Ronde River", "Salmon River"))

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "correct_mainstem_fates_plots", "tucannon_correct_mainstem_fates_plot.pdf"), tucannon_correct_mainstem_fates, height = 8, width = 10)


##### Compare final fates of fish, overshoot vs. being in the correct mainstem #####

overshoot_homing_probabilities <- data.frame(natal_origin = c("Fifteenmile Creek", "Deschutes River", "John Day River",
                                                              "Umatilla River", "Walla Walla River", "Yakima River",
                                                              "Wenatchee River", "Entiat River", "Tucannon River"),
                                             non_overshoot = c(subset(fifteenmile_final_fates_quantiles_forplot_correct_mainstem, state == "Fifteenmile Creek")$label,
                                                               subset(deschutes_final_fates_quantiles_forplot_correct_mainstem, state == "Deschutes River")$label,
                                                               subset(john_day_final_fates_quantiles_forplot_correct_mainstem, state == "John Day River")$label,
                                                               subset(umatilla_final_fates_quantiles_forplot_correct_mainstem, state == "Umatilla River")$label,
                                                               subset(walla_walla_final_fates_quantiles_forplot_correct_mainstem, state == "Walla Walla River")$label,
                                                               subset(yakima_final_fates_quantiles_forplot_correct_mainstem, state == "Yakima River")$label,
                                                               subset(wenatchee_final_fates_quantiles_forplot_correct_mainstem, state == "Wenatchee River")$label,
                                                               subset(entiat_final_fates_quantiles_forplot_correct_mainstem, state == "Entiat River")$label,
                                                               subset(tucannon_final_fates_quantiles_forplot_correct_mainstem, state == "Tucannon River")$label),
                                             overshoot = c(subset(fifteenmile_final_fates_quantiles_forplot_overshoot, state == "Fifteenmile Creek")$label,
                                                               subset(deschutes_final_fates_quantiles_forplot_overshoot, state == "Deschutes River")$label,
                                                               subset(john_day_final_fates_quantiles_forplot_overshoot, state == "John Day River")$label,
                                                               subset(umatilla_final_fates_quantiles_forplot_overshoot, state == "Umatilla River")$label,
                                                               paste0("PRA: ", subset(walla_walla_final_fates_quantiles_forplot_overshoot_PRA, state == "Walla Walla River")$label,
                                                                      " / ICH: ", subset(walla_walla_final_fates_quantiles_forplot_overshoot_ICH, state == "Walla Walla River")$label),
                                                           paste0("PRA: ", subset(yakima_final_fates_quantiles_forplot_overshoot_PRA, state == "Yakima River")$label,
                                                                  " / ICH: ", subset(yakima_final_fates_quantiles_forplot_overshoot_ICH, state == "Yakima River")$label),
                                                               subset(wenatchee_final_fates_quantiles_forplot_overshoot, state == "Wenatchee River")$label,
                                                               subset(entiat_final_fates_quantiles_forplot_overshoot, state == "Entiat River")$label,
                                                               subset(tucannon_final_fates_quantiles_forplot_overshoot, state == "Tucannon River")$label))

write.csv(overshoot_homing_probabilities, here::here("stan_actual", "deteff_ESU_models", "output_tables", "overshoot_homing_probabilities.csv"))


##### Fallback probabilities - dam + origin specific #####
middle_columbia_origins <- c("Fifteenmile Creek", "Deschutes River", "John Day River", "Umatilla River", "Walla Walla River", "Yakima River")
upper_columbia_origins <- c("Wenatchee River", "Entiat River",    "Okanogan River",  "Methow River")
snake_origins <- c("Tucannon River",   "Asotin Creek", "Clearwater River", "Salmon River",  "Grande Ronde", "Imnaha River")


### Bonneville Dam
# from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON"
# Origin specific for Middle Columbia only
BON_fallback_prob <- data.frame(population = c(middle_columbia_origins,
                                                upper_columbia_origins,
                                                snake_origins),
                                 probability = c(subset(middle_columbia_DE_origin_probs, from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON" & origin == "Fifteenmile_Creek")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON" & origin == "Deschutes_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON" & origin == "John_Day_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON" & origin == "Umatilla_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON" & origin == "Walla_Walla_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON" & origin == "Yakima_River")$probability,
                                                 rep(subset(upper_columbia_DE_DPS_probs, from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON")$probability, length(upper_columbia_origins)),
                                                 rep(subset(snake_DE_DPS_probs, from == "mainstem, BON to MCN" & to == "mainstem, mouth to BON")$probability, length(snake_origins))))


### McNary Dam
# from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN"
# For McNary, each probability is unique, because this is the transition state between all ESUs

MCN_fallback_prob <- data.frame(population = c(middle_columbia_origins,
                                                upper_columbia_origins,
                                                snake_origins),
                                 probability = c(subset(middle_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Fifteenmile_Creek")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Deschutes_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "John_Day_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Umatilla_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Walla_Walla_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Yakima_River")$probability,
                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Wenatchee_River")$probability,
                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Entiat_River")$probability,
                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Okanogan_River")$probability,
                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Methow_River")$probability,
                                 # rep(subset(upper_columbia_DE_DPS_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN")$probability, length(upper_columbia_origins)),
                                           subset(snake_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Tucannon_River")$probability,
                                           subset(snake_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Asotin_Creek")$probability,
                                           subset(snake_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Clearwater_River")$probability,
                                           subset(snake_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Salmon_River")$probability,
                                           subset(snake_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Grande_Ronde")$probability,
                                           subset(snake_DE_origin_probs, from == "mainstem, MCN to ICH or PRA" & to == "mainstem, BON to MCN" & origin == "Imnaha_River")$probability))



### Priest Rapids Dam
# from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA"
# Specific for middle and upper columbia
PRA_fallback_prob <- data.frame(population = c(middle_columbia_origins,
                                               upper_columbia_origins,
                                               snake_origins),
                                 probability = c(subset(middle_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Fifteenmile_Creek")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Deschutes_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "John_Day_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Umatilla_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Walla_Walla_River")$probability,
                                                 subset(middle_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Yakima_River")$probability,
                                                 # rep(subset(upper_columbia_DE_DPS_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA")$probability, length(upper_columbia_origins)),
                                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Wenatchee_River")$probability,
                                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Entiat_River")$probability,
                                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Okanogan_River")$probability,
                                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA" & origin == "Methow_River")$probability,
                                                 rep(subset(snake_DE_DPS_probs, from == "mainstem, PRA to RIS" & to == "mainstem, MCN to ICH or PRA")$probability, length(snake_origins))))



### Rock Island Dam
# Specific for Upper columbia only
# from == "mainstem, RIS to RRE" & to == "mainstem, PRA to RIS"
RIS_fallback_prob <- data.frame(population = c(middle_columbia_origins,
                                               upper_columbia_origins,
                                               snake_origins),
                                 probability = c(rep(subset(middle_columbia_DE_DPS_probs, from == "mainstem, RIS to RRE" & to == "mainstem, PRA to RIS")$probability, length(middle_columbia_origins)),
                                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, RIS to RRE" & to == "mainstem, PRA to RIS" & origin == "Wenatchee_River")$probability,
                                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, RIS to RRE" & to == "mainstem, PRA to RIS" & origin == "Entiat_River")$probability,
                                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, RIS to RRE" & to == "mainstem, PRA to RIS" & origin == "Okanogan_River")$probability,
                                                 subset(upper_columbia_DE_origin_probs, from == "mainstem, RIS to RRE" & to == "mainstem, PRA to RIS" & origin == "Methow_River")$probability,
                                                 rep(subset(snake_DE_DPS_probs, from == "mainstem, RIS to RRE" & to == "mainstem, PRA to RIS")$probability, length(snake_origins))))



### Rocky Reach Dam
# Specific for Upper columbia only
# from == "mainstem, RRE to WEL" & to == "mainstem, RIS to RRE"
RRE_fallback_prob <- data.frame(population = c(middle_columbia_origins,
                                               upper_columbia_origins,
                                               snake_origins),
                                probability = c(rep(subset(middle_columbia_DE_DPS_probs, from == "mainstem, RRE to WEL" & to == "mainstem, RIS to RRE")$probability, length(middle_columbia_origins)),
                                                subset(upper_columbia_DE_origin_probs, from == "mainstem, RRE to WEL" & to == "mainstem, RIS to RRE" & origin == "Wenatchee_River")$probability,
                                                subset(upper_columbia_DE_origin_probs, from == "mainstem, RRE to WEL" & to == "mainstem, RIS to RRE" & origin == "Entiat_River")$probability,
                                                subset(upper_columbia_DE_origin_probs, from == "mainstem, RRE to WEL" & to == "mainstem, RIS to RRE" & origin == "Okanogan_River")$probability,
                                                subset(upper_columbia_DE_origin_probs, from == "mainstem, RRE to WEL" & to == "mainstem, RIS to RRE" & origin == "Methow_River")$probability,
                                                rep(subset(snake_DE_DPS_probs, from == "mainstem, RRE to WEL" & to == "mainstem, RIS to RRE")$probability, length(snake_origins))))



### Wells Dam
# Specific for Upper columbia only
# from == "mainstem, upstream of WEL" & to == "mainstem, RRE to WEL"
WEL_fallback_prob <- data.frame(population = c(middle_columbia_origins,
                                               upper_columbia_origins,
                                               snake_origins),
                                probability = c(rep(subset(middle_columbia_DE_DPS_probs, from == "mainstem, upstream of WEL" & to == "mainstem, RRE to WEL")$probability, length(middle_columbia_origins)),
                                                subset(upper_columbia_DE_origin_probs, from == "mainstem, upstream of WEL" & to == "mainstem, RRE to WEL" & origin == "Wenatchee_River")$probability,
                                                subset(upper_columbia_DE_origin_probs, from == "mainstem, upstream of WEL" & to == "mainstem, RRE to WEL" & origin == "Entiat_River")$probability,
                                                subset(upper_columbia_DE_origin_probs, from == "mainstem, upstream of WEL" & to == "mainstem, RRE to WEL" & origin == "Okanogan_River")$probability,
                                                subset(upper_columbia_DE_origin_probs, from == "mainstem, upstream of WEL" & to == "mainstem, RRE to WEL" & origin == "Methow_River")$probability,
                                                rep(subset(snake_DE_DPS_probs, from == "mainstem, upstream of WEL" & to == "mainstem, RRE to WEL")$probability, length(snake_origins))))



### Ice Harbor Dam
# Specific for Snake and Middle Columbia
# from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA"
ICH_fallback_prob <- data.frame(population = c(middle_columbia_origins,
                                               upper_columbia_origins,
                                               snake_origins),
                                probability = c(subset(middle_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Fifteenmile_Creek")$probability,
                                                subset(middle_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Deschutes_River")$probability,
                                                subset(middle_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "John_Day_River")$probability,
                                                subset(middle_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Umatilla_River")$probability,
                                                subset(middle_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Walla_Walla_River")$probability,
                                                subset(middle_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Yakima_River")$probability,
                                                # subset(upper_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Wenatchee_River")$probability,
                                                # subset(upper_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Entiat_River")$probability,
                                                # subset(upper_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Okanogan_River")$probability,
                                                # subset(upper_columbia_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Methow_River")$probability,
                                                rep(subset(upper_columbia_DE_DPS_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA")$probability, length(upper_columbia_origins)),
                                                subset(snake_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Tucannon_River")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Asotin_Creek")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Clearwater_River")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Salmon_River")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Grande_Ronde")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, ICH to LGR" & to == "mainstem, MCN to ICH or PRA" & origin == "Imnaha_River")$probability))


### Lower Granite Dam
# Specific only for Snake
# from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR"
LGR_fallback_prob <- data.frame(population = c(middle_columbia_origins,
                                               upper_columbia_origins,
                                               snake_origins),
                                probability = c(rep(subset(middle_columbia_DE_DPS_probs, from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR")$probability, length(middle_columbia_origins)),
                                                rep(subset(upper_columbia_DE_DPS_probs, from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR")$probability, length(upper_columbia_origins)),
                                                subset(snake_DE_origin_probs, from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR" & origin == "Tucannon_River")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR" & origin == "Asotin_Creek")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR" & origin == "Clearwater_River")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR" & origin == "Salmon_River")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR" & origin == "Grande_Ronde")$probability,
                                                subset(snake_DE_origin_probs, from == "mainstem, upstream of LGR" & to == "mainstem, ICH to LGR" & origin == "Imnaha_River")$probability))

### Join all together
dplyr::rename(BON_fallback_prob, BON = probability) %>% 
  left_join(., dplyr::rename(MCN_fallback_prob, MCN = probability), by = "population") %>% 
  left_join(., dplyr::rename(PRA_fallback_prob, PRA = probability), by = "population") %>% 
  left_join(., dplyr::rename(RIS_fallback_prob, RIS = probability), by = "population") %>% 
  left_join(., dplyr::rename(RRE_fallback_prob, RRE = probability), by = "population") %>% 
  left_join(., dplyr::rename(WEL_fallback_prob, WEL = probability), by = "population") %>% 
  left_join(., dplyr::rename(ICH_fallback_prob, ICH = probability), by = "population") %>% 
  left_join(., dplyr::rename(LGR_fallback_prob, LGR = probability), by = "population") %>% 
  dplyr::rename(origin = population) -> dam_specific_fallback_probs


write.csv(dam_specific_fallback_probs, here::here("stan_actual", "deteff_ESU_models", "output_tables", "dam_specific_fallback_probs.csv"), row.names = FALSE)





##### Plot detection efficiency #####

# Version 2 - using draws
deteff_draws <- function(draws_object, alpha, beta, niter){
  
  # Make this work even if we don't have discharge values - if beta = NA
  if (is.na(beta)){
    draws_object %>% 
      dplyr::select(colnames(draws_object)[grepl(alpha, colnames(draws_object))]) %>%  # keep the actual probabilities
      as.data.frame() -> DE_param_draws
    DE_param_draws$beta <- 0
  } else {
    draws_object %>% 
      dplyr::select(c(colnames(draws_object)[grepl(alpha, colnames(draws_object))],
                      colnames(draws_object)[grepl(beta, colnames(draws_object))])) %>%  # keep the actual probabilities
      as.data.frame() -> DE_param_draws
    
  }
  

  

  
  
  # create a sequence of Z-scores
  seq(-2,2, 0.01) -> zscore_discharge_values
  
  DE_predicted <- matrix(nrow = length(zscore_discharge_values), ncol = niter)
  # Then predict detection probability at each draws of alpha and beta for each
  for (i in 1:length(zscore_discharge_values)){
    for (j in 1:niter){
      DE_predicted[i,j] <- plogis(DE_param_draws[j,1] + DE_param_draws[j,2] * zscore_discharge_values[i])
    }
  }
  
  # Then get quantiles at each zscore
  DE_predicted_quantiles <- data.frame(discharge = zscore_discharge_values,
                                       median = NA,
                                       q2.5 = NA,
                                       q97.5 = NA)
  
  for (i in 1:nrow(DE_predicted)){
    DE_predicted_quantiles$median[i] <- quantile(DE_predicted[i,], 0.5)
    DE_predicted_quantiles$q2.5[i] <- quantile(DE_predicted[i,], 0.025)
    DE_predicted_quantiles$q97.5[i] <- quantile(DE_predicted[i,], 0.975)
  }
  
  DE_predicted_quantiles %>% 
    pivot_longer(., names_to = "name", cols = c("q2.5", "median", "q97.5")) -> DE_predicted_quantiles_long
  
  DE_predicted_quantiles_long$name <- factor(DE_predicted_quantiles_long$name, levels = c("q2.5", "median", "q97.5"))
  
  return(DE_predicted_quantiles_long)
}

# Step 1: Extract correct parameters
det_eff_params <- c("asotin_alpha1", "asotin_alpha2", "deschutes_alpha1", "entiat_alpha1"           ,
                    "fifteenmile_alpha1", "hood_alpha1", "imnaha_alpha1", "john_day_alpha1", "methow_alpha1", "methow_alpha2"           ,
                    "okanogan_alpha1", "tucannon_alpha1", "tucannon_alpha2", "umatilla_alpha1", "umatilla_alpha2", "walla_walla_alpha1",
                    "walla_walla_alpha2", "walla_walla_alpha3", "wenatchee_alpha1", "yakima_alpha1", "asotin_beta", "deschutes_beta"         ,
                    "entiat_beta", "hood_beta", "john_day_beta", "methow_beta"  ,
                    "okanogan_beta", "tucannon_beta", "umatilla_beta", "walla_walla_beta", "wenatchee_beta", "yakima_beta")

# Step 3: Make a function to plot them
DE_plot_function <- function(df, trib_name){
  plot <- ggplot(df, aes(x = discharge, y = value, color = name, linetype = name)) +
    geom_line(size = 1) +
    xlab("Discharge (Z-scored)") +
    ylab("Detection Efficiency") +
    scale_color_manual(values = c("grey50", "black", "grey50")) +
    scale_linetype_manual(values = c(2,1,2)) +
    ggtitle(trib_name) +
    theme(legend.position = "none",
          axis.text = element_text(size = 10),
          axis.title = element_text(size = 14)) +
    scale_y_continuous(limits = c(0,1), expand = c(0,0))
  
  return(plot)
}


# Step 4: Plot them
# Middle Columbia Tributaries
#  deteff_draws <- function(draws_object, alpha, beta, niter)
deteff_draws(alpha = "umatilla_alpha1", beta = "umatilla_beta", draws_object = middle_columbia_draws, niter = 200) -> umatilla_predicted_DE1
umatilla_DE_plot1 <- DE_plot_function(df = umatilla_predicted_DE1, trib_name = "Umatilla River (07/08-13/14)")

deteff_draws(alpha = "umatilla_alpha2", beta = "umatilla_beta", draws_object = middle_columbia_draws, niter = 200) -> umatilla_predicted_DE2
umatilla_DE_plot2 <- DE_plot_function(df = umatilla_predicted_DE2, trib_name = "Umatilla River (14/15-21/22)")

deteff_draws(alpha = "walla_walla_alpha1", beta = "walla_walla_beta", draws_object = middle_columbia_draws, niter = 200) -> walla_walla_predicted_DE1
walla_walla_DE_plot1 <- DE_plot_function(df = walla_walla_predicted_DE1, trib_name = "Walla Walla River (05/06-11/12)")

deteff_draws(alpha = "walla_walla_alpha2", beta = "walla_walla_beta", draws_object = middle_columbia_draws, niter = 200) -> walla_walla_predicted_DE2
walla_walla_DE_plot2 <- DE_plot_function(df = walla_walla_predicted_DE2, trib_name = "Walla Walla River (12/13-18/19)")

deteff_draws(alpha = "walla_walla_alpha3", beta = "walla_walla_beta", draws_object = middle_columbia_draws, niter = 200) -> walla_walla_predicted_DE3
walla_walla_DE_plot3 <- DE_plot_function(df = walla_walla_predicted_DE3, trib_name = "Walla Walla River (19/20-21/22)")

deteff_draws(alpha = "yakima_alpha1", beta = "yakima_beta", draws_object = middle_columbia_draws, niter = 200) -> yakima_predicted_DE1
yakima_DE_plot1 <- DE_plot_function(df = yakima_predicted_DE1, trib_name = "Yakima River (04/05-21/22)")

deteff_draws(alpha = "deschutes_alpha1", beta = "deschutes_beta", draws_object = middle_columbia_draws, niter = 200) -> deschutes_predicted_DE1
deschutes_DE_plot1 <- DE_plot_function(df = deschutes_predicted_DE1, trib_name = "Deschutes River (12/13-19/20)")

deteff_draws(alpha = "john_day_alpha1", beta = "john_day_beta", draws_object = middle_columbia_draws, niter = 200) -> john_day_predicted_DE1
john_day_DE_plot1 <- DE_plot_function(df = john_day_predicted_DE1, trib_name = "John Day River (12/13-21/22)")

deteff_draws(alpha = "fifteenmile_alpha1", beta = NA, draws_object = middle_columbia_draws, niter = 200) -> fifteenmile_predicted_DE1
fifteenmile_DE_plot1 <- DE_plot_function(df = fifteenmile_predicted_DE1, trib_name = "Fifteenmile Creek (12/13-18/19)")

# Hood River is technically lower Columbia
deteff_draws(alpha = "hood_alpha1", beta = "hood_beta", draws_object = middle_columbia_draws, niter = 200) -> hood_predicted_DE1
hood_DE_plot1 <- DE_plot_function(df = hood_predicted_DE1, trib_name = "Hood River (15/16-21/22)")


# Upper Columbia Tributaries
deteff_draws(alpha = "okanogan_alpha1", beta = "okanogan_beta", draws_object = upper_columbia_draws, niter = 200) -> okanogan_predicted_DE
okanogan_DE_plot <- DE_plot_function(df = okanogan_predicted_DE, trib_name = "Okanogan (12/13-21/22)")

deteff_draws(alpha = "wenatchee_alpha1", beta = "wenatchee_beta", draws_object = upper_columbia_draws, niter = 200) -> wenatchee_predicted_DE1
wenatchee_DE_plot1 <- DE_plot_function(df = wenatchee_predicted_DE1, trib_name = "Wenatchee River (10/11-21/22)")

deteff_draws(alpha = "entiat_alpha1", beta = "entiat_beta", draws_object = upper_columbia_draws, niter = 200) -> entiat_predicted_DE1
entiat_DE_plot1 <- DE_plot_function(df = entiat_predicted_DE1, trib_name = "Entiat River (07/08-21/22)")

deteff_draws(alpha = "methow_alpha1", beta = "methow_beta", draws_object = upper_columbia_draws, niter = 200) -> methow_predicted_DE1
methow_DE_plot1 <- DE_plot_function(df = methow_predicted_DE1, trib_name = "Methow River (09/10-16/17)")

deteff_draws(alpha = "methow_alpha2", beta = "methow_beta", draws_object = upper_columbia_draws, niter = 200) -> methow_predicted_DE2
methow_DE_plot2 <- DE_plot_function(df = methow_predicted_DE2, trib_name = "Methow River (17/18-21/22)")


# Snake River Tributaries
deteff_draws(alpha = "asotin_alpha1", beta = "asotin_beta", draws_object = snake_draws, niter = 200) -> asotin_predicted_DE1
asotin_DE_plot1 <- DE_plot_function(df = asotin_predicted_DE1, trib_name = "Asotin Creek (11/12-17/18)")

deteff_draws(alpha = "asotin_alpha2", beta = "asotin_beta", draws_object = snake_draws, niter = 200) -> asotin_predicted_DE2
asotin_DE_plot2 <- DE_plot_function(df = asotin_predicted_DE2, trib_name = "Asotin Creek (18/19-21/22)")

deteff_draws(alpha = "tucannon_alpha1", beta = "tucannon_beta", draws_object = snake_draws, niter = 200) -> tucannon_predicted_DE1
tucannon_DE_plot1 <- DE_plot_function(df = tucannon_predicted_DE1, trib_name = "Tucannon River (10/11-19/20)")

deteff_draws(alpha = "tucannon_alpha2", beta = "tucannon_beta", draws_object = snake_draws, niter = 200) -> tucannon_predicted_DE2
tucannon_DE_plot2 <- DE_plot_function(df = tucannon_predicted_DE2, trib_name = "Tucannon River (20/21-21/22)")

deteff_draws(alpha = "imnaha_alpha1", beta = NA, draws_object = snake_draws, niter = 200) -> imnaha_predicted_DE1
imnaha_DE_plot1 <- DE_plot_function(df = imnaha_predicted_DE1, trib_name = "Imnaha River (12/13-21/22)")


# ARRANGE ALL DE PLOTS

ggarrange(umatilla_DE_plot1, 
          umatilla_DE_plot2,
                       walla_walla_DE_plot1,
                       walla_walla_DE_plot2,
                       walla_walla_DE_plot3,
                       yakima_DE_plot1,
                       deschutes_DE_plot1,
                       john_day_DE_plot1,
                       hood_DE_plot1,
                       okanogan_DE_plot,
                       wenatchee_DE_plot1,
                       entiat_DE_plot1,
                       methow_DE_plot1,
                       methow_DE_plot2,
                       asotin_DE_plot1,
                       asotin_DE_plot2,
                       tucannon_DE_plot1,
                       tucannon_DE_plot2,
          fifteenmile_DE_plot1,
                       imnaha_DE_plot1, ncol = 4, nrow = 5) -> DE_plots

ggsave(here::here("stan_actual", "deteff_ESU_models", "output_tables", "DE_plots.pdf"), DE_plots, height = 15, width = 20)




