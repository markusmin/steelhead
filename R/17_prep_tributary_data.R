# 17 - prep tributary data from USGS for stan model 


# Load libraries
library(plyr)
library(tidyverse)
library(here)
library(lubridate)
library(janitor)
library(kableExtra)
library(ggrepel)

# Function to load USGS data from file

# first create run year df
run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2022-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2023-05-31 23:59:59"), by = "years")

run_year_df <- data.frame(run_year, run_year_start, run_year_end)

# function to load data
usgs_data_load <- function(variables, tributary_name, usgs_file_name){
  if(variables == "discharge only") {
    usgs_headers <- read.table(here::here("tributary_data", usgs_file_name), skip = 28, header = F, nrows = 1, as.is = T, sep = "\t")
    usgs <- read.table(here::here("tributary_data", usgs_file_name), skip = 30, header = F, sep = "\t")
    colnames(usgs) <- usgs_headers
    
    if (tributary_name == "Imnaha River"){
      # drop the time zone column
      usgs %>% 
        dplyr::select(-tz_cd) -> usgs
    }
    
    
    # edit column names and drop irrelevant ones
    colnames(usgs) <- c("agency", "site_no", "date_time", "discharge_cfs", "discharge_cfs_approved")
    
    # If it's the imnaha, it's reported every 15 minutes. We want just the daily; remove the time
    if (tributary_name == "Imnaha River"){
      usgs %>% 
        mutate(date_time =substr(date_time, start = 1, stop = 10)) -> usgs
    }
    
    usgs %>% 
      dplyr::select("date_time", "discharge_cfs") %>% 
      mutate(discharge_cfs = as.numeric(discharge_cfs)) %>% 
      dplyr::mutate(date_time = ymd(date_time)) %>% 
      dplyr::mutate(month = month(date_time)) %>% 
      dplyr::mutate(year = year(date_time))-> usgs
    
    # summarise into monthly averages
    usgs %>% 
      group_by(year,month) %>% 
      dplyr::summarise(mean_discharge_cfs = mean(discharge_cfs, na.rm = TRUE)) %>% 
      dplyr::mutate(year_month = ymd(paste(year, month,"01",sep="-")))-> monthly_discharge
    
    # calculate annual (run year) averages
    usgs %>% 
      rowwise() %>%
      dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_time & run_year_end >= date_time)$run_year) %>% 
      subset(!(run_year %in% c("04/05", "22/23"))) %>% 
      group_by(run_year) %>% 
      dplyr::summarise(mean_discharge_cfs = mean(discharge_cfs, na.rm = TRUE)) %>% 
      dplyr::mutate(year_month = ymd(paste(substr(run_year, start = 4, stop = 5),"01-01",sep="-"))) %>% 
      # add a line to paste tributary name
      mutate(tributary = tributary_name) %>% 
      # drop the year_month column
      dplyr::select(-year_month) -> annual_discharge
    
    # return data for stan model
    
    return(annual_discharge)
    
    
  } else{
    usgs_headers <- read.table(here::here("tributary_data", usgs_file_name), skip = 30, header = F, nrows = 1, as.is = T, sep = "\t")
    usgs <- read.table(here::here("tributary_data", usgs_file_name), skip = 32, header = F, sep = "\t")
    colnames(usgs) <- usgs_headers
    
    
    # edit column names and drop irrelevant ones
    colnames(usgs) <- c("agency", "site_no", "date_time", "discharge_cfs", "discharge_cfs_approved", "gage_height_feet", "gage_height_feet_approved")
    
    usgs %>% 
      dplyr::select("date_time", "discharge_cfs", "gage_height_feet") %>% 
      mutate(discharge_cfs = as.numeric(discharge_cfs)) %>% 
      mutate(gage_height_feet = as.numeric(gage_height_feet)) %>% 
      dplyr::mutate(date_time = ymd(date_time)) %>% 
      dplyr::mutate(month = month(date_time)) %>% 
      dplyr::mutate(year = year(date_time)) -> usgs
    
    # summarise into monthly averages
    usgs %>% 
      group_by(year,month) %>% 
      dplyr::summarise(mean_discharge_cfs = mean(discharge_cfs, na.rm = TRUE)) %>% 
      dplyr::mutate(year_month = ymd(paste(year, month,"01",sep="-")))-> monthly_discharge
    
    usgs %>% 
      group_by(year,month) %>% 
      dplyr::summarise(mean_gage_height_ft = mean(gage_height_feet, na.rm = TRUE)) %>% 
      dplyr::mutate(year_month = ymd(paste(year, month,"01",sep="-")))-> monthly_gage_height
    
    # calculate annual (run year) averages
    usgs %>% 
      rowwise() %>%
      dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_time & run_year_end >= date_time)$run_year) %>% 
      subset(!(run_year %in% c("04/05", "22/23"))) %>% 
      group_by(run_year) %>% 
      dplyr::summarise(mean_discharge_cfs = mean(discharge_cfs, na.rm = TRUE)) %>% 
      dplyr::mutate(year_month = ymd(paste(substr(run_year, start = 4, stop = 5),"01-01",sep="-"))) %>% 
      # add a line to paste tributary name
      mutate(tributary = tributary_name) %>% 
      # drop the year_month column
      dplyr::select(-year_month) -> annual_discharge
    
    usgs %>% 
      rowwise() %>%
      dplyr::mutate(run_year = subset(run_year_df, run_year_start <= date_time & run_year_end >= date_time)$run_year) %>% 
      subset(!(run_year %in% c("04/05", "22/23"))) %>% 
      group_by(run_year) %>% 
      dplyr::summarise(mean_gage_height_ft = mean(gage_height_feet, na.rm = TRUE)) %>% 
      dplyr::mutate(year_month = ymd(paste(substr(run_year, start = 4, stop = 5),"01-01",sep="-"))) %>% 
      # add a line to paste tributary name
      mutate(tributary = tributary_name) %>% 
      # drop the year_month column
      dplyr::select(-year_month) -> annual_gage_height
    
    
    return(list(annual_discharge, annual_gage_height))
  }
  
}

# Get data for each of the tributaries
asotin_discharge <- usgs_data_load(variables = "discharge and gage height", tributary_name = "Asotin Creek", usgs_file_name = "asotin_creek_13335050.txt")[[1]]

deschutes_discharge <- usgs_data_load(variables = "discharge only", tributary_name = "Deschutes River", usgs_file_name = "deschutes_river_14103000.txt")

entiat_discharge <- usgs_data_load(variables = "discharge and gage height", tributary_name = "Entiat River", usgs_file_name = "entiat_river_12452990.txt")[[1]]

hood_discharge <- usgs_data_load(variables = "discharge only", tributary_name = "Hood River", usgs_file_name = "hood_river_14120000.txt")

imnaha_discharge <- usgs_data_load(variables = "discharge only", tributary_name = "Imnaha River", usgs_file_name = "imnaha_river_13292000.txt")

john_day_discharge <- usgs_data_load(variables = "discharge only", tributary_name = "John Day River", usgs_file_name = "john_day_river_14048000.txt")

methow_discharge <- usgs_data_load(variables = "discharge and gage height", tributary_name = "Methow River", usgs_file_name = "methow_river_12449950.txt")[[1]]

okanogan_discharge <- usgs_data_load(variables = "discharge and gage height", tributary_name = "Okanogan River", usgs_file_name = "okanogan_12447200.txt")[[1]]

tucannon_discharge <- usgs_data_load(variables = "discharge and gage height", tributary_name = "Tucannon River", usgs_file_name = "tucannon_river_13344500.txt")[[1]]

umatilla_discharge <- usgs_data_load(variables = "discharge only", tributary_name = "Umatilla River", usgs_file_name = "umatilla_river_14033500.txt")

walla_walla_discharge <- usgs_data_load(variables = "discharge and gage height", tributary_name = "Walla Walla River", usgs_file_name = "walla_walla_14018500.txt")[[1]]

wenatchee_discharge <- usgs_data_load(variables = "discharge and gage height", tributary_name = "Wenatchee River", usgs_file_name = "wenatchee_river_12462500.txt")[[1]]

yakima_discharge <- usgs_data_load(variables = "discharge and gage height", tributary_name = "Yakima River", usgs_file_name = "yakima_river_12510500.txt")[[1]]

# combine these all into one df for export

asotin_discharge %>% 
  bind_rows(., deschutes_discharge) %>% 
  bind_rows(., entiat_discharge) %>% 
  bind_rows(., hood_discharge) %>% 
  bind_rows(., imnaha_discharge) %>% 
  bind_rows(., john_day_discharge) %>% 
  bind_rows(., methow_discharge) %>% 
  bind_rows(., okanogan_discharge) %>% 
  bind_rows(., tucannon_discharge) %>% 
  bind_rows(., umatilla_discharge) %>% 
  bind_rows(., walla_walla_discharge) %>% 
  bind_rows(., wenatchee_discharge) %>% 
  bind_rows(., yakima_discharge) -> all_tribs_discharge

# Export this file
write.csv(all_tribs_discharge, here::here("covariate_data", "tributary_discharge_data.csv"), row.names = FALSE)

# Export a second file - z-score the values (by each of the tributaries), then export
all_tribs_discharge %>% 
  group_by(tributary) %>% 
  mutate(mean_discharge_zscore = (mean_discharge_cfs - mean(mean_discharge_cfs))/sd(mean_discharge_cfs)) -> all_tribs_discharge_zscore

write.csv(all_tribs_discharge_zscore, here::here("covariate_data", "tributary_discharge_data_zscore.csv"), row.names = FALSE)


