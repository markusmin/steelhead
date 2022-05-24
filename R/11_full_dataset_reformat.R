# 11 Full dataset - reformat for JAGS model

# In this R script, we will reformat all of the detection histories so that we can fit the JAGS model

### Load libraries
library(here)
library(tidyverse)
library(lubridate)
library(janitor)

### Load data
states_times <- read.csv(here::here("model_files", "states.csv"))

unique(states_times$tag_code)

sort(unique(states_times$state))


## Load data from hyak run
states_complete <- read.csv(here::here("from_hyak_transfer", "2022-05-23-complete_det_hist", "states_complete.csv"))

length(unique(states_complete$tag_code))

run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22")
run_year_start <- seq(ymd("2004-06-01"), ymd("2021-06-01"), by = "years")
run_year_end <- seq(ymd("2005-05-31"), ymd("2022-05-31"), by = "years")

run_year_df <- data.frame(run_year, run_year_start, run_year_end)

det_hist %>%
  group_by(tag_code) %>%
  filter(row_number() == 1) %>%
  mutate(run_year = subset(run_year_df, run_year_start < start_time & run_year_end > start_time)$run_year) %>%
  dplyr::select(tag_code, run_year) -> tag_codes_run_year