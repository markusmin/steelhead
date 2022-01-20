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

# Look at when fish were seen at what site
# BON corner collector - looks like pretty much every year
table(subset(JDR_det_hist, event_site_name == "BCC - BON PH2 Corner Collector")$start_time)






