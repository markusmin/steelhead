# Visualize the rear, temp, and year effect model outputs from stan



# This script will load and investigate the outputs from our stan models.
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(tidyverse)
library(here)
library(ggpubr)
library(stringr)

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

# Get info about state names and numbers
from_state_number_names <- data.frame(from = seq(1,43,1), from_name = model_states)
to_state_number_names <- data.frame(to = seq(1,43,1), to_name = model_states)

# Read in the 100iter test run
# note that this file is mis-named
upper_columbia_fit <- readRDS(here::here("stan_actual", "year_only", "upper_columbia", "seed101_500iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))

# this takes a long time to load with RE year effect
upper_columbia_summary <- upper_columbia_fit$summary()

# some quick traceplots to look at chains
mcmc_trace(upper_columbia_fit$draws(), pars = c("b0_matrix_1_2"))

sigma_year_matrix_1_2_trace <- mcmc_trace(upper_columbia_fit$draws(), pars = c("sigma_year_matrix_1_2"))

# tracking all of these is what's causing this to be huge
byear_raw_vector_9_36_8_trace <- mcmc_trace(upper_columbia_fit$draws(), pars = c("byear_raw_vector_9_38[6]"))
sigma_year_matrix_1_2_trace
byear_raw_vector_9_36_8_trace
# should also look at sigma_year parameters
upper_columbia_summary %>% 
  filter(grepl("sigma_year_matrix_[[:digit:]]", variable)) -> sigma_year_params

# what do these look like?
ggplot(sigma_year_params, aes(x = variable, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  coord_flip() +
  ylab("Year Effect Scaling")

# Inspect the distribution of byear parameters
upper_columbia_summary %>% 
  filter(grepl("byear_raw_vector", variable)) %>% 
  mutate(from = as.numeric(sub("[^_]*_[^_]*_[^_]*_", "", str_extract(variable, "[^_]*_[^_]*_[^_]*_[^_]*")))) %>%
  # mutate(to = as.numeric(sub("\D*", "", sub("[^_]*_[^_]*_[^_]*_[^_]*_", "", variable))))  %>% 
  mutate(to = as.numeric(str_extract(sub("[^_]*_[^_]*_[^_]*_[^_]*_", "", variable), "\\d+")))  %>% 
  left_join(from_state_number_names, by = "from") %>% 
  left_join(to_state_number_names, by = "to") %>% 
  arrange(from, to) -> byear_params

# are any of these different from zero?
subset(byear_params, q5 > 0 | q95 < 0)
# zero

# some plots
byear_2_3_plot <- ggplot(subset(byear_params, from  == 2 & to == 3), aes(x = variable, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  coord_flip() +
  ylab("Year Effect")

# what about deschutes movement? that should show up with year based on warm vs.  cold years
byear_2_10_plot <- ggplot(subset(byear_params, from  == 2 & to == 10), aes(x = variable, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  coord_flip() +
  ylab("Year Effect")


##### investigate two temperature parameters #####

# keep only temperatuer parameters that don't have brackets (these are the actual parameters)
temp_variables <- upper_columbia_summary$variable[grep("btemp", upper_columbia_summary$variable)][1:82]

# 2temp model - look at upstream of WEL
upper_columbia_summary %>% 
  filter(grepl("btemp", variable)) %>% 
  filter(!(grepl("\\[", variable))) %>% 
  mutate(from = as.numeric(sub("[^_]*_[^_]*_", "", str_extract(variable, "[^_]*_[^_]*_[^_]*")))) %>%
  mutate(to = as.numeric(sub("_.*", "", sub("[^_]*_[^_]*_[^_]*_", "", variable))))  %>% 
  left_join(from_state_number_names, by = "from") %>% 
  left_join(to_state_number_names, by = "to") %>% 
  mutate(origin = ifelse(grepl("btemp_", variable), "DPS",
                         ifelse(grepl("origin1", variable), "Wenatchee",
                                ifelse(grepl("origin2", variable), "Entiat",
                                       ifelse(grepl("origin3", variable), "Methow", NA))))) %>% 
  mutate(to_name = sub(" Mouth", "", to_name)) %>% 
  # for now, drop all NDE versions of parameters 
  filter(!(grepl("NDE", variable))) %>% 
  mutate(movement = paste0(from_name, " -> ", to_name)) -> temp_param_summary


# Let's just look at which are "significant"
subset(temp_param_summary, q5 > 0 |  q95 < 0) -> significant_temp_variables_2

# Movement into the Deschutes - 2 -> 10 - there's a huge positive temperature effect, regardless of time of year
# Entiat River - overshoot when it's hot


# Check Entiat, since it had a statistically significant effect in Shelby's

# overshoot is 6 -> 7 for Entiat
# Entiat are origin 2
subset(upper_columbia_summary, variable %in% c("btemp0xorigin2_matrix_6_7", "btemp1xorigin2_matrix_6_7"))
# yes, these are both significant




# compare to the old temperature model
  
write.csv(significant_temp_variables_2, here::here("stan_actual", "rear_temp", "processed", "significant_temp_variable_2.csv"))

significant_temp_variables_2 %>% 
  mutate(movement = paste0(from, "_", to)) %>% 
  dplyr::select(variable, median, q5, q95, movement, from, to) -> significant_temp_variables_2_for_join


# read in the other one
significant_temp_variables_1 <- read.csv(here::here("stan_actual", "rear_temp", "processed", "significant_temp_variable_1.csv"), row.names = 1)
  

significant_temp_variables_1 %>% 
  # for now, drop all NDE versions of parameters 
  filter(!(grepl("NDE", variable))) %>% 
  dplyr::select(variable, median, q5, q95) %>% 
  # dplyr::rename(median_old = median, q5_old = q5, q95_old = q95) %>% 
  mutate(from = as.numeric(sub("[^_]*_[^_]*_", "", str_extract(variable, "[^_]*_[^_]*_[^_]*")))) %>%
  mutate(to = as.numeric(sub("_.*", "", sub("[^_]*_[^_]*_[^_]*_", "", variable)))) %>% 
  mutate(movement = paste0(from, "_", to)) -> significant_temp_variables_1

significant_temp_variables_2_for_join %>% 
  bind_rows(., significant_temp_variables_1) %>% 
  mutate(parameter = ifelse(grepl("temp0", variable), "winter_spring_temp",
                            ifelse(grepl("temp1", variable), "summer_fall_temp", "single_temp"))) %>% 
  mutate(origin_type = ifelse(grepl("origin1", variable), "Wenatchee",
                              ifelse(grepl("origin2", variable), "Entiat",
                                     ifelse(grepl("origin3", variable), "Methow", "DPS-wide")))) %>% 
  left_join(from_state_number_names, by = "from") %>% 
  left_join(to_state_number_names, by = "to") %>% 
  mutate(movement = paste0(from_name, "->", to_name)) -> temp_comp

temp_comp_plot <- ggplot(temp_comp, aes(x = movement, y = median, color = parameter)) +
  geom_point(position = position_dodge(width = 0.75)) +
  geom_errorbar(aes(ymax = q95, ymin = q5), position = position_dodge(width = 0.75)) +
  coord_flip() +
  facet_wrap(~origin_type) +
  geom_hline(yintercept = 0, lty = 2) +
  ylab("Temperature Parameter Magnitude")

ggsave(here::here("stan_actual", "rear_temp_year", "processed", "temp_comparison_plot.pdf"), plot = temp_comp_plot, height = 6, width = 10)








