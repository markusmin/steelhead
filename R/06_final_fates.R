### 06 - final fates

# Description: This R script uses a transition matrix to convert movement probabilities into final fates.

library(here)
library(tidyverse)

param_table <- read.csv(here::here("model_files", "JDR_param_estimates.csv"), row.names = 1)



