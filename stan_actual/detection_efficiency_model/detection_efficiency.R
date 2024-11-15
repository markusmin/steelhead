# detection_efficiency_script
# this script applies to all tributaries, such that only one of these models must be run 
# this R script will reformat data for inputs into the accompanying stan script, // detection_efficiency

# load libraries
library(plyr)
library(tidyverse)
library(here)
library(lubridate)
library(janitor)
library(cmdstanr)

# load tributary data
# tributary_data <- read.csv(here::here("covariate_data", "tributary_discharge_data.csv"))
tributary_data <- read.csv(here::here("covariate_data", "tributary_discharge_data_zscore.csv"))

# split the data up by tributary

subset(tributary_data, tributary == "Asotin Creek") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> asotin_discharge
subset(tributary_data, tributary == "Deschutes River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> deschutes_discharge
subset(tributary_data, tributary == "Entiat River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> entiat_discharge
subset(tributary_data, tributary == "Fifteenmile Creek") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> fifteenmile_discharge
# Create fake data for fifteenmile, with all zeros for discharge - allows fitting of intercept-only model
fifteenmile_discharge <- data.frame(run_year = asotin_discharge$run_year, mean_discharge_zscore = 0)

subset(tributary_data, tributary == "Hood River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> hood_discharge
subset(tributary_data, tributary == "Imnaha River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> imnaha_discharge
# Create fake data for imnaha, with all zeros for discharge - allows fitting of intercept-only model
imnaha_discharge <- data.frame(run_year = asotin_discharge$run_year, mean_discharge_zscore = 0)

subset(tributary_data, tributary == "John Day River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> john_day_discharge
subset(tributary_data, tributary == "Methow River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> methow_discharge
subset(tributary_data, tributary == "Okanogan River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> okanogan_discharge
subset(tributary_data, tributary == "Tucannon River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> tucannon_discharge
subset(tributary_data, tributary == "Umatilla River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> umatilla_discharge
subset(tributary_data, tributary == "Walla Walla River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> walla_walla_discharge
subset(tributary_data, tributary == "Wenatchee River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> wenatchee_discharge
subset(tributary_data, tributary == "Yakima River") %>% 
  dplyr::select(-c(tributary, mean_discharge_cfs)) -> yakima_discharge

# load fish detections
states_complete <- read.csv(here::here("stan_actual", "adults_states_complete.csv"), row.names = 1)

# load tag code metadata (contains run year info)
tag_code_metadata <- read.csv(here::here("covariate_data", "tag_code_metadata.csv"))

# join fish detections with run year data
states_complete %>% 
  left_join(., dplyr::select(tag_code_metadata, tag_code, run_year), by = "tag_code") -> states_complete




# Now use these states to subset the states_complete df and reformat data for stan model

# write a function:
det_eff_data_prep <- function(rm_state, ups_state){
  subset(states_complete, state %in% c(rm_state, ups_state)) -> deteff_data
  deteff_data_ind <- data.frame(tag_code = unique(deteff_data$tag_code), upstream = NA, river_mouth = NA)
  for (i in 1:nrow(deteff_data_ind)){
    data <- subset(deteff_data, tag_code == deteff_data_ind$tag_code[i])
    if(rm_state %in% data$state){
      deteff_data_ind$river_mouth[i] <- 1
    } else {
      deteff_data_ind$river_mouth[i] <- 0
    }
    if(ups_state %in% data$state){
      deteff_data_ind$upstream[i] <- 1
    } else {
      deteff_data_ind$upstream[i] <- 0
    }
  }
  
  # add run year
  deteff_data_ind %>% 
    left_join(., dplyr::select(tag_code_metadata, tag_code, run_year), by = "tag_code") -> deteff_data_ind
  
  # Now reformat; new column for detected; 0 if not, 1 if detected
  
  deteff_data_ind %>% 
    # we can only look at fish that are seen at the upstream site
    subset(upstream == 1) %>% 
    mutate(detected = ifelse(river_mouth == 1, 1, 0)) %>% 
    dplyr::select(tag_code, run_year, detected) -> deteff_data_ind
  
  return(deteff_data_ind)
}

# get data for each of the snake tribs

asotin_deteff_data <- det_eff_data_prep(rm_state = "Asotin Creek Mouth", ups_state = "Asotin Creek Upstream")
deschutes_deteff_data <- det_eff_data_prep(rm_state = "Deschutes River Mouth", ups_state = "Deschutes River Upstream")
entiat_deteff_data <- det_eff_data_prep(rm_state = "Entiat River Mouth", ups_state = "Entiat River Upstream")
fifteenmile_deteff_data <- det_eff_data_prep(rm_state = "Fifteenmile Creek Mouth", ups_state = "Fifteenmile Creek Upstream")
hood_deteff_data <- det_eff_data_prep(rm_state = "Hood River Mouth", ups_state = "Hood River Upstream")
imnaha_deteff_data <- det_eff_data_prep(rm_state = "Imnaha River Mouth", ups_state = "Imnaha River Upstream")
john_day_deteff_data <- det_eff_data_prep(rm_state = "John Day River Mouth", ups_state = "John Day River Upstream")
methow_deteff_data <- det_eff_data_prep(rm_state = "Methow River Mouth", ups_state = "Methow River Upstream")
okanogan_deteff_data <- det_eff_data_prep(rm_state = "Okanogan River Mouth", ups_state = "Okanogan River Upstream")
tucannon_deteff_data <- det_eff_data_prep(rm_state = "Tucannon River Mouth", ups_state = "Tucannon River Upstream")
umatilla_deteff_data <- det_eff_data_prep(rm_state = "Umatilla River Mouth", ups_state = "Umatilla River Upstream")
walla_walla_deteff_data <- det_eff_data_prep(rm_state = "Walla Walla River Mouth", ups_state = "Walla Walla River Upstream")
wenatchee_deteff_data <- det_eff_data_prep(rm_state = "Wenatchee River Mouth", ups_state = "Wenatchee River Upstream")
yakima_deteff_data <- det_eff_data_prep(rm_state = "Yakima River Mouth", ups_state = "Yakima River Upstream")


run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2022-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2023-05-31 23:59:59"), by = "years")
run_year_numeric = seq(4, 22, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)

### Record categorical eras
# first create a template that you can modify for each tributary
run_year_eras <- data.frame(run_year = run_year,
                       run_year_numeric = run_year_numeric,
                       era = 0)



run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 11 & run_year_numeric <= 17, 1,
                      ifelse(run_year_numeric >= 18, 2, 0))) %>% 
  dplyr::select(-run_year_numeric) -> asotin_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 13 & run_year_numeric <= 18, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> deschutes_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 7, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> entiat_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 11 & run_year_numeric <= 18, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> fifteenmile_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 12, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> hood_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 10, 1, 0)) %>% # note that Imnaha is a special case because it doesn't have discharge data
  dplyr::select(-run_year_numeric) -> imnaha_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 12, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> john_day_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 9 & run_year_numeric <= 16, 1,
                      ifelse(run_year_numeric >= 17, 2, 0))) %>% 
  dplyr::select(-run_year_numeric) -> methow_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 13, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> okanogan_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 10 & run_year_numeric <= 19, 1,
                      ifelse(run_year_numeric >= 20, 2, 0))) %>% 
  dplyr::select(-run_year_numeric) -> tucannon_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 6 & run_year_numeric <= 13, 1,
                      ifelse(run_year_numeric >= 14, 2, 0))) %>% 
  dplyr::select(-run_year_numeric) -> umatilla_eras

# edit 2023-05-08: Add one more era to WAWA, when ORB and PRV sites overlapped
# for the purposes of this analysis, we will consider this one RM site, since they are only 1 KM apart
# so we will make the assumption that fish traveled past both
# This is for run years 12/13-14/15
run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 5 & run_year_numeric <= 11, 1,
                      ifelse(run_year_numeric >= 12 & run_year_numeric <= 14, 2,
                             ifelse(run_year_numeric >= 15 & run_year_numeric <= 18, 3,
                                  ifelse(run_year_numeric >= 19, 4, 0))))) %>% 
  dplyr::select(-run_year_numeric) -> walla_walla_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 10, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> wenatchee_eras

run_year_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 4, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> yakima_eras


# Collate all info into design matrix

asotin_deteff_data %>% 
  left_join(., asotin_eras, by = "run_year") %>% 
  left_join(., asotin_discharge, by = "run_year") %>% 
  mutate(tributary = "Asotin")-> asotin_deteff_data

deschutes_deteff_data %>% 
  left_join(., deschutes_eras, by = "run_year") %>% 
  left_join(., deschutes_discharge, by = "run_year") %>% 
  mutate(tributary = "Deschutes")-> deschutes_deteff_data

entiat_deteff_data %>% 
  left_join(., entiat_eras, by = "run_year") %>% 
  left_join(., entiat_discharge, by = "run_year") %>% 
  mutate(tributary = "Entiat")-> entiat_deteff_data

fifteenmile_deteff_data %>% 
  left_join(., fifteenmile_eras, by = "run_year") %>% 
  left_join(., fifteenmile_discharge, by = "run_year") %>% 
  mutate(tributary = "Fifteenmile")-> fifteenmile_deteff_data

hood_deteff_data %>% 
  left_join(., hood_eras, by = "run_year") %>% 
  left_join(., hood_discharge, by = "run_year") %>% 
  mutate(tributary = "Hood")-> hood_deteff_data

imnaha_deteff_data %>% 
  left_join(., imnaha_eras, by = "run_year") %>% 
  left_join(., imnaha_discharge, by = "run_year") %>% 
  mutate(tributary = "Imnaha")-> imnaha_deteff_data

john_day_deteff_data %>% 
  left_join(., john_day_eras, by = "run_year") %>% 
  left_join(., john_day_discharge, by = "run_year") %>% 
  mutate(tributary = "John Day")-> john_day_deteff_data

methow_deteff_data %>% 
  left_join(., methow_eras, by = "run_year") %>% 
  left_join(., methow_discharge, by = "run_year") %>% 
  mutate(tributary = "Methow")-> methow_deteff_data

okanogan_deteff_data %>% 
  left_join(., okanogan_eras, by = "run_year") %>% 
  left_join(., okanogan_discharge, by = "run_year") %>% 
  mutate(tributary = "Okanogan")-> okanogan_deteff_data

tucannon_deteff_data %>% 
  left_join(., tucannon_eras, by = "run_year") %>% 
  left_join(., tucannon_discharge, by = "run_year") %>% 
  mutate(tributary = "Tucannon")-> tucannon_deteff_data

umatilla_deteff_data %>% 
  left_join(., umatilla_eras, by = "run_year") %>% 
  left_join(., umatilla_discharge, by = "run_year") %>% 
  mutate(tributary = "Umatilla")-> umatilla_deteff_data

walla_walla_deteff_data %>% 
  left_join(., walla_walla_eras, by = "run_year") %>% 
  left_join(., walla_walla_discharge, by = "run_year") %>% 
  mutate(tributary = "Walla Walla")-> walla_walla_deteff_data

wenatchee_deteff_data %>% 
  left_join(., wenatchee_eras, by = "run_year") %>% 
  left_join(., wenatchee_discharge, by = "run_year") %>% 
  mutate(tributary = "Wenatchee")-> wenatchee_deteff_data

yakima_deteff_data %>% 
  left_join(., yakima_eras, by = "run_year") %>% 
  left_join(., yakima_discharge, by = "run_year") %>% 
  mutate(tributary = "Yakima")-> yakima_deteff_data

# Bind together
asotin_deteff_data %>% 
  bind_rows(., deschutes_deteff_data) %>% 
  bind_rows(., entiat_deteff_data) %>% 
  bind_rows(., fifteenmile_deteff_data) %>% 
  bind_rows(., hood_deteff_data) %>% 
  bind_rows(., imnaha_deteff_data) %>% 
  bind_rows(., john_day_deteff_data) %>% 
  bind_rows(., methow_deteff_data) %>% 
  bind_rows(., okanogan_deteff_data) %>% 
  bind_rows(., tucannon_deteff_data) %>% 
  bind_rows(., umatilla_deteff_data) %>% 
  bind_rows(., walla_walla_deteff_data) %>% 
  bind_rows(., wenatchee_deteff_data) %>% 
  bind_rows(., yakima_deteff_data) -> deteff_data

# convert into a design matrix
design_matrix <- matrix(nrow = nrow(deteff_data), ncol = 35)

# First 21 columns: alphas: aso1, aso2, des1, ent1, fif1, hood1, imn1, jdr1, met1, met2, oka1, tuc1, tuc2, uma1, uma2, wawa1, wawa2, wawa3, wawa4, wen1, yak1
# Columns 21-34: betas: beta_aso, beta_des, beta_ent, beta_fif, beta_hood, beta_imn, beta_jdr, beta_met, beta_oka, beta_tuc, beta_uma, beta_wawa, beta_wen, beta_yak
for (i in 1:nrow(deteff_data)){
  # alphas
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Asotin") {
    design_matrix[i,1] <- 1
  } else {
    design_matrix[i,1] <- 0
  }
  if (deteff_data$era[i] == 2 & deteff_data$tributary[i] == "Asotin") {
    design_matrix[i,2] <- 1
  } else {
    design_matrix[i,2] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Deschutes") {
    design_matrix[i,3] <- 1
  } else {
    design_matrix[i,3] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Entiat") {
    design_matrix[i,4] <- 1
  } else {
    design_matrix[i,4] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Fifteenmile") {
    design_matrix[i,5] <- 1
  } else {
    design_matrix[i,5] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Hood") {
    design_matrix[i,6] <- 1
  } else {
    design_matrix[i,6] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Imnaha") {
    design_matrix[i,7] <- 1
  } else {
    design_matrix[i,7] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "John Day") {
    design_matrix[i,8] <- 1
  } else {
    design_matrix[i,8] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Methow") {
    design_matrix[i,9] <- 1
  } else {
    design_matrix[i,9] <- 0
  }
  if (deteff_data$era[i] == 2 & deteff_data$tributary[i] == "Methow") {
    design_matrix[i,10] <- 1
  } else {
    design_matrix[i,10] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Okanogan") {
    design_matrix[i,11] <- 1
  } else {
    design_matrix[i,11] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Tucannon") {
    design_matrix[i,12] <- 1
  } else {
    design_matrix[i,12] <- 0
  }
  if (deteff_data$era[i] == 2 & deteff_data$tributary[i] == "Tucannon") {
    design_matrix[i,13] <- 1
  } else {
    design_matrix[i,13] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Umatilla") {
    design_matrix[i,14] <- 1
  } else {
    design_matrix[i,14] <- 0
  }
  if (deteff_data$era[i] == 2 & deteff_data$tributary[i] == "Umatilla") {
    design_matrix[i,15] <- 1
  } else {
    design_matrix[i,15] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Walla Walla") {
    design_matrix[i,16] <- 1
  } else {
    design_matrix[i,16] <- 0
  }
  if (deteff_data$era[i] == 2 & deteff_data$tributary[i] == "Walla Walla") {
    design_matrix[i,17] <- 1
  } else {
    design_matrix[i,17] <- 0
  }
  if (deteff_data$era[i] == 3 & deteff_data$tributary[i] == "Walla Walla") {
    design_matrix[i,18] <- 1
  } else {
    design_matrix[i,18] <- 0
  }
  if (deteff_data$era[i] == 4 & deteff_data$tributary[i] == "Walla Walla") {
    design_matrix[i,19] <- 1
  } else {
    design_matrix[i,19] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Wenatchee") {
    design_matrix[i,20] <- 1
  } else {
    design_matrix[i,20] <- 0
  }
  if (deteff_data$era[i] == 1 & deteff_data$tributary[i] == "Yakima") {
    design_matrix[i,21] <- 1
  } else {
    design_matrix[i,21] <- 0
  }
  # Betas
  if (deteff_data$tributary[i] == "Asotin") {
    design_matrix[i,22] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,22] <- 0
  }
  if (deteff_data$tributary[i] == "Deschutes") {
    design_matrix[i,23] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,23] <- 0
  }
  if (deteff_data$tributary[i] == "Entiat") {
    design_matrix[i,24] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,24] <- 0
  }
  if (deteff_data$tributary[i] == "Fifteenmile") {
    design_matrix[i,25] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,25] <- 0
  }
  if (deteff_data$tributary[i] == "Hood") {
    design_matrix[i,26] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,26] <- 0
  }
  if (deteff_data$tributary[i] == "Imnaha") {
    design_matrix[i,27] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,27] <- 0
  }
  if (deteff_data$tributary[i] == "John Day") {
    design_matrix[i,28] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,28] <- 0
  }
  if (deteff_data$tributary[i] == "Methow") {
    design_matrix[i,29] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,29] <- 0
  }
  if (deteff_data$tributary[i] == "Okanogan") {
    design_matrix[i,30] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,30] <- 0
  }
  if (deteff_data$tributary[i] == "Tucannon") {
    design_matrix[i,31] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,31] <- 0
  }
  if (deteff_data$tributary[i] == "Umatilla") {
    design_matrix[i,32] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,32] <- 0
  }
  if (deteff_data$tributary[i] == "Walla Walla") {
    design_matrix[i,33] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,33] <- 0
  }
  if (deteff_data$tributary[i] == "Wenatchee") {
    design_matrix[i,34] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,34] <- 0
  }
  if (deteff_data$tributary[i] == "Yakima") {
    design_matrix[i,35] <- deteff_data$mean_discharge_zscore[i]
  } else {
    design_matrix[i,35] <- 0
  }
}

# Remove all rows from design that have missing data (is there a way around this 
# in the future to still estimate an intercept but no slope?)
# na.omit(design_matrix) -> design_matrix2


##### Fit stan model #####

# step 0: data in a list #
data <- list(N = nrow(design_matrix),
             # y = subset(deteff_data, !(is.na(mean_discharge_cfs)))$detected,
             y = deteff_data$detected,
             J = 21,
             K = 14,
             X = design_matrix)

# Step 1: load the model
# mod <- cmdstan_model("01_stan_sim_int_only.stan", compile = FALSE)
mod <- cmdstan_model(here::here("stan_actual", "deteff_ESU_models", "detection_efficiency.stan"), compile = FALSE)

# Step 2: Compile the model, set up to run in parallel
mod$compile(cpp_options = list(stan_threads = TRUE))

# Step 3: Run MCMC (HMC)
fit <- mod$sample(
  data = data, 
  # seed = 123, 
  chains = 3,
  # parallel_chains = 1,
  # parallel_chains = 3,
  refresh = 10, # print update every iter
  # iter_sampling = 1000,
  # iter_warmup = 1000,
  iter_warmup = 5000,
  iter_sampling = 5000,
  # threads_per_chain = 28,
  threads_per_chain = 7,
  init = 1,
  thin = 10
)

mcmc_trace(fit$draws(), pars = c("alpha[1]"))


# Extract parameter estimates
fit$save_object(file = here::here("stan_actual", "detection_efficiency_model", "detection_efficiency_stan_fit_new.rds"))

# fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "detection_efficiency_stan_fit_new.rds"))

# Calculate detection probability from parameter estimates
deteff_stan_fit_summary <- fit$summary()

# simple diagnostics
acf(fit$draws()[,2,3])
# mcmc_trace(fit$draws(), pars = c("alpha[1]"))

deteff_stan_fit_summary[1:10,]

##### EXPORT PARAMETER ESTIMATES TO USE AS PRIORS IN PRIMARY MODEL #####
# create matrix to store values
det_eff_param_posteriors <- matrix(nrow = 35, ncol = 2)
det_eff_param_posteriors[,1] <- deteff_stan_fit_summary$mean[2:36]
det_eff_param_posteriors[,2] <- deteff_stan_fit_summary$sd[2:36]

write.csv(det_eff_param_posteriors, here::here("stan_actual", "detection_efficiency_model", "det_eff_param_posteriors_new.csv"))





# Get the parameters in a vector, then multiply by the matrix of discharges and calculate p for each run year/trib combination
snake_deteff_param_est <- snake_deteff_stan_fit_summary$median[2:9]

# Create the matrix of run year/era/tributry
snake_run_year_trib_deteff <- matrix(nrow = length(seq(5,21,1)), ncol = 3)
rownames(snake_run_year_trib_deteff) <- tucannon_discharge$run_year
colnames(snake_run_year_trib_deteff) <- c("Tucannon River", "Asotin Creek", "Imnaha River")

# mini design matrix for tucannon
# Ignore 22/23 run year, since we don't have full data for this year
# nrow = 12, because this is the number of years we can calculate a detection efficiency for (i.e., we have a river mouth site)
# ncol = 3; two for different alphas, one for beta
tucannon_design_matrix <- matrix(nrow = 12, ncol = 3)
tucannon_design_matrix[,1] <- c(rep(1, 10), rep(0, 2))
tucannon_design_matrix[,2] <- c(rep(0, 10), rep(1, 2))
tucannon_design_matrix[,3] <- tucannon_discharge$mean_discharge_cfs[6:nrow(tucannon_discharge)]
# for ease of interpretation (and for sanity check), add run years:
rownames(tucannon_design_matrix) <- tucannon_discharge$run_year[6:nrow(tucannon_discharge)]

# Get the linear predictor eta
# design_matrix2 %*% as.matrix(snake_deteff_param_est, nrow = 8, ncol = 1) -> snake_eta_vec
tucannon_design_matrix %*% as.matrix(snake_deteff_param_est[c(1,2,6)], nrow = 8, ncol = 1) -> tucannon_eta_vec

# Get the estimated detection efficiency using inv logit
snake_run_year_trib_deteff[6:nrow(snake_run_year_trib_deteff), 1] <- exp(tucannon_eta_vec)/(1 + exp(tucannon_eta_vec))



# Now, get the Imnaha detection efficiency without any covariates
subset(deteff_data, is.na(mean_discharge_cfs) & tributary == "Imnaha") -> imn_nodischarge

# Get p, plus 90% CI
imn_nodisp <- sum(imn_nodischarge$detected)/nrow(imn_nodischarge)

alpha = 0.10 # 90% confidence interval
nd = sum(imn_nodischarge$detected)
n = nrow(imn_nodischarge)

## Define functions for upper and lower limits
## for a 90 % confidence interval.
fl = function(p){pbinom(nd-1,n,p) - (1-alpha/2)}
fu = function(p){pbinom(nd,n,p) - alpha/2}
# all in one line:
lower = uniroot(function(p){pbinom(nd-1,n,p) - (1-alpha/2)}, c(0.01, 0.99))$root
upper = uniroot(function(p){pbinom(nd,n,p) - alpha/2}, c(0.01, 0.99))$root
print(paste0("mean: ", round(imn_nodisp, 3), "; q5: ", round(lower, 3), "; q95: ", round(upper, 3)))




