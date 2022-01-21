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

# Dams in order - separate Columbia and Snake
columbia_dams <- c("Bonneville Adult Fishways (combined)", "McNary Adult Fishways (Combined)", 
          "PRA - Priest Rapids Adult","RIA - Rock Island Adult", 
          "RRF - Rocky Reach Fishway", "WEA - Wells Dam, DCPUD Adult Ladders")

snake_dams <- c("ICH - Ice Harbor Dam (Combined)",  "GRA - Lower Granite Dam Adult")

# Natal tributary detection sites
BON_MCN_natal_sites <- JDR_event_site_metadata$event_site_name[grep("John Day", JDR_event_site_metadata$event_site_basin_name)]

# Straying sites between BON and MCN:
# Basins: Descuhtes River
# Subbasins: Klickitat, Umatilla, Fifteenmile Creek
# grep: Hood in event_site_name, since "Middle Columbia-Hood" contains mainstem sites
# Note: When using grep, you need to take out the river mouth arrays
BON_MCN_stray_sites <- c(JDR_event_site_metadata$event_site_name[grep("Hood|Fifteenmile|Deschutes",
                                                                      JDR_event_site_metadata$event_site_basin_name)],
                           JDR_event_site_metadata$event_site_name[grep("Klickitat|Umatilla",
                                                                        JDR_event_site_metadata$event_site_subbasin_name)],
                         JDR_event_site_metadata$event_site_name[grep("Hood",
                                                                      JDR_event_site_metadata$event_site_name)]
                         )
# Take out river mouth arrays
BON_MCN_stray_sites <- BON_MCN_stray_sites[!grepl("mouth|Mouth", BON_MCN_stray_sites)]
 

# Straying sites between MCN and branches (PRA or ICH):
# Basins: Yakima
# Subbasins: Walla Walla
MCN_PRA_ICH_stray_sites <- c(JDR_event_site_metadata$event_site_name[grep("Yakima",
                                                                          JDR_event_site_metadata$event_site_basin_name)],
                             JDR_event_site_metadata$event_site_name[grep("Walla Walla",
                                                                          JDR_event_site_metadata$event_site_subbasin_name)])
# Take out river mouth arrays
MCN_PRA_ICH_stray_sites <- MCN_PRA_ICH_stray_sites[!grepl("mouth|Mouth", MCN_PRA_ICH_stray_sites)]

# Straying sites between ICH and LGR:
# Basins: NA
# Subbasins: NA
# Here, you have to grep the site name, because Little Goose and Lower Monumental
# are considered part of the "Lower Snake-Tucannon" subbasin
ICH_LGR_stray_sites <- c(JDR_event_site_metadata$event_site_name[grep("Tucannon",
                                                                      JDR_event_site_metadata$event_site_name)])

# Take out river mouth arrays
# Here I'm not sure if the Lower Tucannon River Array is a similar distance as a
# river mouth array
ICH_LGR_stray_sites <- ICH_LGR_stray_sites[!grepl("mouth|Mouth", ICH_LGR_stray_sites)]


# Straying sites upstream of LGR:
# Basins: Salmon
# Subbasins: Upper Grande Ronde, Clearwater, Lower Grande Ronde
LGR_upstream_stray_sites <- c(JDR_event_site_metadata$event_site_name[grep("Salmon",
                                                                          JDR_event_site_metadata$event_site_basin_name)],
                             JDR_event_site_metadata$event_site_name[grep("Upper Grande Ronde|Clearwater|Lower Grande Ronde",
                                                                          JDR_event_site_metadata$event_site_subbasin_name)])

# Take out river mouth arrays
# Don't do this for these sites - because there are creek mouths in here,
# i.e., "LAP - Lapwai Creek, near its mouth" and "SWT - Sweetwater Cr. near its mouth"  
# LGR_upstream_stray_sites <- LGR_upstream_stray_sites[!grepl("mouth|Mouth", LGR_upstream_stray_sites)]


# Let's see which sites aren't captured in our sites
setdiff(JDR_event_site_metadata$event_site_name, c(LGR_upstream_stray_sites,
                                                   ICH_LGR_stray_sites,
                                                   MCN_PRA_ICH_stray_sites,
                                                   BON_MCN_stray_sites,
                                                   BON_MCN_natal_sites))
