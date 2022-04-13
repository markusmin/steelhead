# 07_covariate_data_plots

# Exploratory plots of various covariates

# Load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)
library(ggthemes)

temp_cov_df <- read.csv(here::here("covariate_data", "temperature_by_state.csv"), row.names = 1)
spill_cov_df <- read.csv(here::here("covariate_data", "spill_by_state.csv"), row.names = 1)
flow_cov_df <- read.csv(here::here("covariate_data", "flow_by_state.csv"), row.names = 1)

# Create DF of state + dam
dam_state <- data.frame(state = c("mouth.to.BON",
                                  "MCN.to.ICH.or.PRA..ICH.",
                                  "ICH.to.LGR",
                                  "BON.to.MCN",           
                                  "MCN.to.ICH.or.PRA..PRA.",
                                  "PRA.to.RIS", 
                                  "RIS.to.RRE", 
                                  "RRE.to.WEL"),
                        dam = c("BON",
                                "ICH",
                                "LGR",
                                "MCN",
                                "PRA",
                                "RIS",
                                "RRE",
                                "WEL"))


##### Temperature #####

# Pivot longer
temp_cov_df %>% 
  pivot_longer(., cols = colnames(temp_cov_df[2:ncol(temp_cov_df)])) %>% 
  dplyr::rename(state = name, temp = value) %>% 
  mutate(date = ymd(date)) %>% 
  left_join(., dam_state, by = "state")-> temp_long

# Remove anomalously high values - above 25 C
temp_long %>% 
  subset(temp <= 25) -> temp_long

ggplot(temp_long, aes(x = date, y = temp, color = dam)) +
         geom_line() +
  scale_color_tableau(palette = "Tableau 10")
       

##### Spill #####

# Pivot longer
spill_cov_df %>% 
  pivot_longer(., cols = colnames(spill_cov_df[2:ncol(spill_cov_df)])) %>% 
  dplyr::rename(state = name, spill = value) %>% 
  mutate(date = ymd(date)) %>% 
  left_join(., dam_state, by = "state") -> spill_long

ggplot(spill_long, aes(x = date, y = spill, color = dam)) +
  geom_line() +
  scale_color_tableau(palette = "Tableau 10")


##### Flow #####

# Pivot longer
flow_cov_df %>% 
  pivot_longer(., cols = colnames(flow_cov_df[2:ncol(flow_cov_df)])) %>% 
  dplyr::rename(state = name, flow = value) %>% 
  mutate(date = ymd(date)) %>% 
  left_join(., dam_state, by = "state") -> flow_long

ggplot(flow_long, aes(x = date, y = flow, color = dam)) +
  geom_line() +
  scale_color_tableau(palette = "Tableau 10")



##### Inspect MCN #####
MCN_temp_long <- subset(temp_long, dam == "MCN" & date >= ymd("2010-01-01"))
MCN_spill_long <- subset(spill_long, dam == "MCN" & date >= ymd("2010-01-01"))
MCN_flow_long <- subset(flow_long, dam == "MCN" & date >= ymd("2010-01-01"))

# Let's look at flow at just MCN, 2010 to present
ggplot(MCN_temp_long, aes(x = date, y = temp, color = dam)) +
  geom_line() +
  scale_color_tableau(palette = "Tableau 10")

ggplot(MCN_spill_long, aes(x = date, y = spill, color = dam)) +
  geom_line() +
  scale_color_tableau(palette = "Tableau 10")

ggplot(MCN_flow_long, aes(x = date, y = flow, color = dam)) +
  geom_line() +
  scale_color_tableau(palette = "Tableau 10")





