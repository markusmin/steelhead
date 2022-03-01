### 07 - covariates

# This script will reformat covariate data for inclusion in the multistate model


# Load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)

# Load data from McNary

# Tailrace temperature
MCN_t_tailrace <- clean_names(read.csv(here::here("covariate_data","MCN","basin_tempc_tailrace_MCN_allyears.csv")))

# fix column names
colnames(MCN_t_tailrace) <- gsub("x", "", colnames(MCN_t_tailrace))

# Pivot longer
MCN_t_tailrace %>% 
  pivot_longer(., cols = colnames(MCN_t_tailrace)[2:ncol(MCN_t_tailrace)]) %>% 
  dplyr::rename(year = name, temp = value) %>% 
  mutate(day = as.numeric(day), year = as.factor(year)) %>% 
  arrange(year, day) -> MCN_t_tailrace_long

# Plot it

ggplot(MCN_t_tailrace_long, aes(x = day, y = temp, color = year)) +
  geom_line()

# Add a data filtering step to remove what appears to be instrument error

# Data filtering function
data_filter <- function(data, limit, var){
  # Flag values that are more than the limit different from the previous
  data %>% 
    mutate(drop = ifelse(is.na(eval(parse(text = var))), FALSE,
                         ifelse(abs(eval(parse(text = var)) - lag(eval(parse(text = var)))) > limit, TRUE, FALSE))) -> data
  
  clean_data <- list(clean = subset(data, drop == FALSE), dropped = subset(data, drop == TRUE))
  
  return(clean_data)
}

subset(data, drop == TRUE)

# this isn't quite working right and we don't need it, so we will just manually subset


# Store clean data
MCN_t_tailrace_long_clean <- subset(MCN_t_tailrace_long, temp < 28)

ggplot(MCN_t_tailrace_long_clean, aes(x = day, y = temp, color = year)) +
  geom_line()

#### Reformat into mean temperature for certain portions of year

MCN_t_tailrace_long_clean %>% 
  mutate(window = ifelse(day <= 90 & day >= 0, "winter",
                         ifelse(day <= 250 & day >= 160, "summer", "none"))) -> MCN_t_tailrace_long_clean

MCN_t_tailrace_long_clean %>% 
  group_by(year, window) %>% 
  summarise(mean(temp)) %>% 
  dplyr::rename(mean_temp = `mean(temp)`)-> MCN_t_tailrace_window_means
  
##### Look at distribution of run timing #####

adult_returns <- clean_names(read.csv(here::here("PTAGIS_queries", "intermediate_files", "all_adult_returns.csv")))

adult_returns %>% 
  mutate(event_date = mdy(first_obs_date_max))-> adult_returns

# Get the date of arrival at BON - take min in case it fell back over Bonneville
adult_returns %>% 
  group_by(tag_code) %>% 
  subset(site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                                "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF", 
                                "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) %>% 
  filter(event_date == min(event_date)) %>% 
  dplyr::select(tag_code, event_date) %>% 
  dplyr::rename(BON_arrival = event_date) -> adult_BON_arrival

# Add run year info

# Sort into run years
run_year <- c("05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21", "21/22")
run_year_start <- seq(ymd("2005-06-01"), ymd("2021-06-01"), by = "years")
run_year_end <- seq(ymd("2006-05-31"), ymd("2022-05-31"), by = "years")

run_year_df <- data.frame(run_year, run_year_start, run_year_end)

adult_BON_arrival %>% 
  mutate(dummy =TRUE) %>% 
  left_join(run_year_df %>% mutate(dummy=TRUE), by = "dummy") %>% 
  filter(run_year_start <= BON_arrival, run_year_end >= BON_arrival) %>% 
  select(-c(dummy, run_year_start, run_year_end)) -> adult_BON_arrival

# Get date timing
adult_BON_arrival %>% 
  mutate(DayMonth = format(as.Date(BON_arrival), "%m-%d")) %>% 
  mutate(DayMonth = yday(BON_arrival)) -> adult_BON_arrival


# Plot date
ggplot(data = adult_BON_arrival, aes(x = BON_arrival)) +
  geom_bar(stat = "count")

ggplot(data = adult_BON_arrival, aes(x = DayMonth)) +
  geom_bar(stat = "count") +
  xlab("Day of year") +
  ggtitle("Arrival at Bonneville, all adult Steelhead, 2005-2022")


# Plot John Day River steelhead, arrival at natal tributary
JDR_CTH_1 <- clean_names(read.csv(here::here("PTAGIS_queries", "complete_tag_histories", "john_day_river", "2022-01-14-john_day_river_CTH_2005_2015_1.csv")))
JDR_CTH_2 <- clean_names(read.csv(here::here("PTAGIS_queries", "complete_tag_histories", "john_day_river", "2022-01-14-john_day_river_CTH_2005_2015_2.csv")))

### Combine files, fix some column data types
JDR_CTH_1 %>% 
  bind_rows(., JDR_CTH_2) %>% 
  mutate(event_date_time_value = mdy_hms(event_date_time_value))-> JDR_CTH

# Get names of natal tributary sites
JDR_det_hist <- read.csv(file = here::here("model_files", "JDR_det_hist.csv"))
JDR_det_hist %>% 
  group_by(event_site_name) %>% 
  summarise(n()) %>% 
  dplyr::rename(count = `n()`)-> JDR_event_det_counts

# Get the metadata
JDR_det_hist %>% 
  dplyr::select(-c(tag_code, start_time, end_time)) %>% 
  distinct(event_site_name, .keep_all = TRUE) -> JDR_event_site_metadata

BON_MCN_natal_sites <- JDR_event_site_metadata$event_site_name[grep("John Day", JDR_event_site_metadata$event_site_basin_name)]

# Get the date of arrival at BON - take min in case it fell back over Bonneville
JDR_CTH %>% 
  group_by(tag_code) %>% 
  subset(event_site_name %in% c("BHL - Adult Fishway at BONH", "BO1 - Bonneville Bradford Is. Ladder", 
                                "BO2 - Bonneville Cascades Is. Ladder", "BO3 - Bonneville WA Shore Ladder/AFF", 
                                "BO4 - Bonneville WA Ladder Slots", "BONAFF - BON - Adult Fish Facility")) %>% 
  filter(event_date_time_value == min(event_date_time_value)) %>% 
  dplyr::select(tag_code, event_date_time_value) %>% 
  dplyr::rename(BON_arrival = event_date_time_value) -> JDR_BON_arrival

JDR_CTH %>% 
  group_by(tag_code) %>% 
  subset(event_site_name %in% BON_MCN_natal_sites) %>% 
  filter(event_date_time_value == max(event_date_time_value)) %>% # take max, since min would be release
  dplyr::select(tag_code, event_date_time_value) %>% 
  dplyr::rename(nat_trib_arrival = event_date_time_value) -> JDR_nat_trib_arrival

# Join the two, keep only those where the nat_trib_arrival is greater than the BON arrival

JDR_nat_trib_arrival %>% 
  left_join(., JDR_BON_arrival, by = "tag_code") %>% 
  filter(nat_trib_arrival > BON_arrival) %>% 
  # Get day and month
  mutate(DayMonth = yday(nat_trib_arrival)) -> JDR_nat_trib_final

ggplot(data = JDR_nat_trib_final, aes(x = DayMonth)) +
  geom_bar(stat = "count") +
  xlab("Day of year") +
  ggtitle("Arrival at natal tributary, JDR steelhead 05/06-14/15")

##### Window of time for in-stream covariates #####
# 

