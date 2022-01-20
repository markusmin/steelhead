### 03 - Fit multistate model

# NOTE: To fit this model, you will first need to run 02_detection_histories.R
# for each of your river systems, with each output named "<river_system>_det_hist"


### Load libraries
library(here)
library(tidyverse)
library(lubridate)
library(janitor)

# Look at what sites
JDR_site_count <- as.data.frame(table(JDR_det_hist$event_site_name))

# Get unique
JDR_event_site_metadata

# Look at when fish were seen at what site
# BON corner collector - looks like pretty much every year
table(subset(JDR_det_hist, event_site_name == "BCC - BON PH2 Corner Collector")$start_time)

# Dams in order
dams <- c("Bonneville Adult Fishways (combined)", "McNary Adult Fishways (Combined)", 
          "ICH - Ice Harbor Dam (Combined)", "PRA - Priest Rapids Adult",
          "RIA - Rock Island Adult", "WEA - Wells Dam, DCPUD Adult Ladders",
          "GRA - Lower Granite Dam Adult")

# Natal tributary detection sites
dam1_dam2_natal_sites <- JDR_event_site_metadata$event_site_name[grep("John Day", JDR_event_site_metadata$event_site_basin_name)]

# Straying sites between BON and MCN
dam1_dam2_stray_sites <- c("WHITEC - White Creek, Klickitat River Basin",
                           "SHERFT - Sherars Falls Fishway Trap, Deschutes River",
                           "TROU2C - Trout Creek, Deschutes River Watershed",
                           )

# Straying detection sites




