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
# upper_columbia_fit <- readRDS(here::here("stan_actual", "year_only", "upper_columbia", "seed101_500iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))


# Read in the 200iter run from stf
# this file is still misnamed as 300 iter - it's actually 100 warmup and 100 sampling
# well apparently I used 200 warmup and 50 sampling, which is very dumb
upper_columbia_fit <- readRDS(here::here("stan_actual", "year_only_v2", "upper_columbia", "seed101_300iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))

# this takes a long time to load with RE year effect
upper_columbia_summary <- upper_columbia_fit$summary()

# some quick traceplots to look at chains
mcmc_trace(upper_columbia_fit$draws(), pars = c("b0_matrix_1_2"))
mcmc_trace(upper_columbia_fit$draws(), pars = c("sigma_yearxorigin1_matrix_4_5"))
# this still looks spiky, but that's to be expected with a cauchy prior (and won't be a big deal once we 
# have a credible interval that basically drops these)
# tracking all of these is what's causing this to be huge
mcmc_trace(upper_columbia_fit$draws(), pars = c("byearxorigin2_raw_vector_6_7[17]"))


# should also look at sigma_year parameters
upper_columbia_summary %>% 
  filter(grepl("sigma_yearxorigin[[:digit:]]_matrix_[[:digit:]]", variable)) -> sigma_year_params

# what do these look like?
ggplot(sigma_year_params, aes(x = variable, y = median)) +
  geom_point() +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  coord_flip() +
  ylab("Year Effect Scaling")



# first create the run year df
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2022-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2023-05-31 23:59:59"), by = "years")
run_year_numeric = seq(4, 22, 1)
year_index = seq(1,19, 1)

run_year_df <- data.frame(year_index, run_year, run_year_start, run_year_end, run_year_numeric)

# Inspect the distribution of byear parameters
upper_columbia_summary %>% 
  filter(grepl("byearxorigin[[:digit:]]_raw_vector", variable)) %>% 
  mutate(from = as.numeric(sub("[^_]*_[^_]*_[^_]*_", "", str_extract(variable, "[^_]*_[^_]*_[^_]*_[^_]*")))) %>%
  # mutate(to = as.numeric(sub("\D*", "", sub("[^_]*_[^_]*_[^_]*_[^_]*_", "", variable))))  %>% 
  mutate(to = as.numeric(str_extract(sub("[^_]*_[^_]*_[^_]*_[^_]*_", "", variable), "\\d+")))  %>% 
  left_join(from_state_number_names, by = "from") %>% 
  left_join(to_state_number_names, by = "to") %>% 
  # add the origin field
  mutate(origin = ifelse(grepl("origin1", variable), "Wenatchee",
                         ifelse(grepl("origin2", variable), "Entiat",
                                ifelse(grepl("origin3", variable), "Methow", "ERROR")))) %>% 
  mutate(year_index = as.numeric(gsub("\\[|\\]", "", str_extract(variable, "\\[[[:digit:]]*]")))) %>% 
  arrange(from, to) %>% 
  # add the run year
  left_join(., dplyr::select(run_year_df, year_index, run_year), by = "year_index")-> byear_params

# are any of these different from zero?
subset(byear_params, q5 > 0 | q95 < 0)
# there are now a decent amount that are different from zero

# some plots
byear_5_24_NDE_Wenatchee_plot <- ggplot(subset(byear_params, from  == 5 & to == 24 & origin == "Wenatchee" & grepl("_NDE", variable))[1:7,], 
                                       # aes(x = fct_reorder(variable, desc(year)), y = median)) +
                                       aes(x = fct_rev(run_year), y = median)) +
  geom_point() +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  coord_flip() +
  ylab("Year Effect") +
  xlab("Run Year") +
  ggtitle("Movement from the mainstem into the Wenatchee River")

byear_5_24_NDE_Wenatchee_plot

byear_5_24_NDE_Wenatchee_plot <- ggplot(subset(byear_params, from  == 5 & to == 24 & origin == "Wenatchee" & grepl("_NDE", variable)), 
                                        # aes(x = fct_reorder(variable, desc(year)), y = median)) +
                                        aes(x = fct_rev(run_year), y = median)) +
  geom_point() +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  coord_flip() +
  ylab("Year Effect") +
  xlab("Run Year") +
  ggtitle("Movement from the mainstem into the Wenatchee River")

byear_5_24_NDE_Wenatchee_plot


# but how does this compare with what years were DE vs. NDE?
# for movement into the Wenatchee:
# run_year_DE_array[5,24,7:18] <- 1
# so this is weird, why is there a trend in [5] and [6]? Those are years leading
# up to it becoming DE; corresponds to run years 08/09 and 09/10
# Tons of arrays were installed in 2008 and 2009. This might explain this.
# this is probably a good thing that year effects can account for differences in DE?

# Overshoot for Wenatchee

byear_5_6_Wenatchee_plot <- ggplot(subset(byear_params, from  == 5 & to == 6 & origin == "Wenatchee" & run_year != "04/05"), 
                                        # aes(x = fct_reorder(variable, desc(year)), y = median)) +
                                        aes(x = fct_rev(run_year), y = median)) +
  geom_point() +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  coord_flip() +
  ylab("Year Effect") +
  xlab("Run Year") +
  ggtitle("Overshoot for Wenatchee River Steelhead")

byear_5_6_Wenatchee_plot

# what about temperature in these years?
read.csv(here::here("covariate_data", "for_model", "temp", "temp_mod_est.csv")) %>% 
  mutate(date = ymd(date)) %>% 
  rowwise() %>%
  dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year) %>% 
  dplyr::mutate(run_year_numeric = subset(run_year_df, run_year_start <= date & run_year_end >= date)$run_year_numeric) %>% 
  # keep only run years up to 21/22 (last year of full data)
  filter(run_year_numeric <= 21) -> temp_data

# summarise into an annual average for RIS
temp_data %>% 
  group_by(run_year) %>% 
  summarise(avg_temp = mean(RIS)) %>% 
  ungroup() -> RIS_annual_temp

# summarise into an annual average for RRE
temp_data %>% 
  group_by(run_year) %>% 
  summarise(avg_temp = mean(RRE)) %>% 
  ungroup() -> RRE_annual_temp


# add a line for temperature
RIS_temp_plot <- ggplot(RIS_annual_temp, aes(x = fct_rev(run_year), y = avg_temp)) +
  geom_point() +
  geom_line(group = 1) +
  coord_flip() +
  ylab("Average Temperature") +
  xlab("Run Year") +
  ggtitle("Mean annual temperature at RIS")

RIS_temp_plot

# what if we combine the temp and param dfs and then run a lm?
subset(byear_params, from  == 5 & to == 6 & origin == "Wenatchee" & run_year != "04/05") %>% 
  left_join(RIS_annual_temp, by = "run_year") -> wen_5_6_temp

library(ggrepel)
ggplot(wen_5_6_temp, aes(x = avg_temp, y = median)) +
  geom_point() +
  ylab("Year Effect Parameter Median") +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label = run_year))

year_temp_wen <- lm(median~avg_temp, data = wen_5_6_temp)

summary(year_temp_wen)


# Overshoot for Entiat

byear_6_7_Entiat_plot <- ggplot(subset(byear_params, from  == 6 & to == 7 & origin == "Entiat" & run_year != "04/05"), 
                                   # aes(x = fct_reorder(variable, desc(year)), y = median)) +
                                   aes(x = fct_rev(run_year), y = median)) +
  geom_point() +
  geom_errorbar(aes(ymax = q95, ymin = q5)) +
  coord_flip() +
  ylab("Year Effect") +
  xlab("Run Year") +
  ggtitle("Overshoot for Entiat River Steelhead")

byear_6_7_Entiat_plot


# what if we combine the temp and param dfs and then run a lm?
subset(byear_params, from  == 6 & to == 7 & origin == "Entiat" & run_year != "04/05") %>% 
  left_join(RRE_annual_temp, by = "run_year") -> ent_6_7_temp

library(ggrepel)
ggplot(ent_6_7_temp, aes(x = avg_temp, y = median)) +
  geom_point() +
  ylab("Year Effect Parameter Median") +
  geom_smooth(method='lm') +
  geom_text_repel(aes(label = run_year))

year_temp_ent <- lm(median~avg_temp, data = ent_6_7_temp)

summary(year_temp_ent)


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








