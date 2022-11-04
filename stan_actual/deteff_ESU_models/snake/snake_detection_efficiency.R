# 03_snake_detection_efficiency - R companion script
# this R script will reformat data for inputs into the accompanying stan script, // 03_snake_detection_efficiency

# load libraries
library(plyr)
library(tidyverse)
library(here)
library(lubridate)
library(janitor)

# load tributary data
tributary_data <- read.csv("tributary_discharge_data.csv")

# subset the tributaries of interest
subset(tributary_data, tributary == "Tucannon River") %>% 
  dplyr::select(-tributary) -> tuc_discharge
subset(tributary_data, tributary == "Asotin Creek") %>% 
  dplyr::select(-tributary) -> asotin_discharge
subset(tributary_data, tributary == "Imnaha River") %>% 
  dplyr::select(-tributary) -> imnaha_discharge

# load fish detections
states_complete <- read.csv("snake_adults_states_complete.csv", row.names = 1)

# load tag code metadata (contains run year info)
tag_code_metadata <- read.csv("tag_code_metadata.csv")

# join fish detections with run year data
states_complete %>% 
  left_join(., dplyr::select(tag_code_metadata, tag_code, run_year), by = "tag_code") -> states_complete




# Note which tributaries are in the snake
# "Tucannon_River", "Asotin_Creek", "Clearwater_River", "Salmon_River", "Grande_Ronde_River", "Imnaha_River"
# Of these, we can only calculate detection efficiency for 3: Tucannon, Asotin Creek, and Imnaha

snake_rm_states <- c("Tucannon River Mouth", "Asotin Creek Mouth", "Imnaha River Mouth")
snake_ups_states <- c("Tucannon River Upstream", "Asotin Creek Upstream", "Imnaha River Upstream")

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
tuc_deteff_data <- det_eff_data_prep(rm_state = "Tucannon River Mouth", ups_state = "Tucannon River Upstream")
asotin_deteff_data <- det_eff_data_prep(rm_state = "Asotin Creek Mouth", ups_state = "Asotin Creek Upstream")
imnaha_deteff_data <- det_eff_data_prep(rm_state = "Imnaha River Mouth", ups_state = "Imnaha River Upstream")


run_year <- c("04/05", "05/06", "06/07", "07/08", "08/09", "09/10", "10/11", "11/12", "12/13", "13/14", "14/15", "15/16", "16/17", "17/18", "18/19", "19/20", "20/21","21/22", "22/23")
run_year_start <- seq(ymd_hms("2004-06-01 00:00:00"), ymd_hms("2022-06-01 00:00:00"), by = "years")
run_year_end <- seq(ymd_hms("2005-05-31 23:59:59"), ymd_hms("2023-05-31 23:59:59"), by = "years")
run_year_numeric = seq(4, 22, 1)

run_year_df <- data.frame(run_year, run_year_start, run_year_end, run_year_numeric)

### Record categorical eras
tuc_eras <- data.frame(run_year = run_year,
                       run_year_numeric = run_year_numeric,
                       era = 0)

tuc_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 10 & run_year_numeric <= 19, 1,
                      ifelse(run_year_numeric >= 20, 2, 0))) %>% 
  dplyr::select(-run_year_numeric) -> tuc_eras

asotin_eras <- data.frame(run_year = run_year,
                       run_year_numeric = run_year_numeric,
                       era = 0)

asotin_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 11 & run_year_numeric <= 17, 1,
                      ifelse(run_year_numeric >= 18, 2, 0))) %>% 
  dplyr::select(-run_year_numeric) -> asotin_eras

imnaha_eras <- data.frame(run_year = run_year,
                          run_year_numeric = run_year_numeric,
                          era = 0)

imnaha_eras %>% 
  mutate(era = ifelse(run_year_numeric >= 12, 1, 0)) %>% 
  dplyr::select(-run_year_numeric) -> imnaha_eras


# Collate all info into design matrix
tuc_deteff_data %>% 
  left_join(., tuc_eras, by = "run_year") %>% 
  left_join(., tuc_discharge, by = "run_year") %>% 
  mutate(tributary = "Tucannon")-> tuc_deteff_data
asotin_deteff_data %>% 
  left_join(., asotin_eras, by = "run_year") %>% 
  left_join(., asotin_discharge, by = "run_year") %>% 
  mutate(tributary = "Asotin")-> asotin_deteff_data
imnaha_deteff_data %>% 
  left_join(., imnaha_eras, by = "run_year") %>% 
  left_join(., imnaha_discharge, by = "run_year") %>% 
  mutate(tributary = "Imnaha")-> imnaha_deteff_data

# Bind together
tuc_deteff_data %>% 
  bind_rows(., asotin_deteff_data) %>% 
  bind_rows(., imnaha_deteff_data) -> snake_deteff_data

# convert into a design matrix
design_matrix <- matrix(nrow = nrow(snake_deteff_data), ncol = 8)

# First five columns: alphas: tuc1, tuc2, aso1, aso2, imn1
# Columns 6-8: betas: beta_tuc, beta_aso, beta_imn
for (i in 1:nrow(snake_deteff_data)){
  
  if (snake_deteff_data$era[i] == 1 & snake_deteff_data$tributary[i] == "Tucannon") {
    design_matrix[i,1] <- 1
  } else {
    design_matrix[i,1] <- 0
  }
  if (snake_deteff_data$era[i] == 2 & snake_deteff_data$tributary[i] == "Tucannon") {
    design_matrix[i,2] <- 1
  } else {
    design_matrix[i,2] <- 0
  }
  if (snake_deteff_data$era[i] == 1 & snake_deteff_data$tributary[i] == "Asotin") {
    design_matrix[i,3] <- 1
  } else {
    design_matrix[i,3] <- 0
  }
  if (snake_deteff_data$era[i] == 2 & snake_deteff_data$tributary[i] == "Asotin") {
    design_matrix[i,4] <- 1
  } else {
    design_matrix[i,4] <- 0
  }
  if (snake_deteff_data$era[i] == 1 & snake_deteff_data$tributary[i] == "Imnaha") {
    design_matrix[i,5] <- 1
  } else {
    design_matrix[i,5] <- 0
  }
  if (snake_deteff_data$tributary[i] == "Tucannon") {
    design_matrix[i,6] <- snake_deteff_data$mean_discharge_cfs[i]
  } else {
    design_matrix[i,6] <- 0
  }
  if (snake_deteff_data$tributary[i] == "Asotin") {
    design_matrix[i,7] <- snake_deteff_data$mean_discharge_cfs[i]
  } else {
    design_matrix[i,7] <- 0
  }
  if (snake_deteff_data$tributary[i] == "Imnaha") {
    design_matrix[i,8] <- snake_deteff_data$mean_discharge_cfs[i]
  } else {
    design_matrix[i,8] <- 0
  }
}

# Remove all rows from design that have missing data (is there a way around this 
# in the future to still estimate an intercept but no slope?)
na.omit(design_matrix) -> design_matrix2


##### Fit stan model #####

# step 0: data in a list #
data <- list(N = nrow(design_matrix2),
             y = subset(snake_deteff_data, !(is.na(mean_discharge_cfs)))$detected,
             J = 5,
             K = 3,
             X = design_matrix2)

# Step 1: load the model
# mod <- cmdstan_model("01_stan_sim_int_only.stan", compile = FALSE)
mod <- cmdstan_model("snake_detection_efficiency.stan", compile = FALSE)

# Step 2: Compile the model, set up to run in parallel
mod$compile(cpp_options = list(stan_threads = TRUE))

# Step 3: Run MCMC (HMC)
fit <- mod$sample(
  data = data, 
  # seed = 123, # this seed gets stuck around 22-24, goes really fast and then at that iteration it slows way down
  # seed = 456,
  # chains = 3, 
  chains = 3,
  # parallel_chains = 1,
  # parallel_chains = 3,
  refresh = 10, # print update every iter
  # iter_sampling = 1000,
  # iter_warmup = 1000,
  iter_warmup = 1000,
  iter_sampling = 1000,
  # threads_per_chain = 28,
  threads_per_chain = 1,
  init = 1
)


# Extract parameter estimates
fit$save_object(file = "snake_detection_efficiency_stan_fit.rds")

# Calculate detection probability from parameter estimates
snake_deteff_stan_fit_summary <- fit$summary()

# Get the parameters in a vector, then multiply by the matrix of discharges and calculate p for each run year/trib combination
snake_deteff_param_est <- snake_deteff_stan_fit_summary$median[2:9]

# Create the matrix of run year/era/tributry
snake_run_year_trib_deteff <- matrix(nrow = length(seq(5,21,1)), ncol = 3)
rownames(snake_run_year_trib_deteff) <- tuc_discharge$run_year
colnames(snake_run_year_trib_deteff) <- c("Tucannon River", "Asotin Creek", "Imnaha River")

# mini design matrix for tucannon
# Ignore 22/23 run year, since we don't have full data for this year
# nrow = 12, because this is the number of years we can calculate a detection efficiency for (i.e., we have a river mouth site)
# ncol = 3; two for different alphas, one for beta
tuc_design_matrix <- matrix(nrow = 12, ncol = 3)
tuc_design_matrix[,1] <- c(rep(1, 10), rep(0, 2))
tuc_design_matrix[,2] <- c(rep(0, 10), rep(1, 2))
tuc_design_matrix[,3] <- tuc_discharge$mean_discharge_cfs[6:nrow(tuc_discharge)]
# for ease of interpretation (and for sanity check), add run years:
rownames(tuc_design_matrix) <- tuc_discharge$run_year[6:nrow(tuc_discharge)]

# Get the linear predictor eta
# design_matrix2 %*% as.matrix(snake_deteff_param_est, nrow = 8, ncol = 1) -> snake_eta_vec
tuc_design_matrix %*% as.matrix(snake_deteff_param_est[c(1,2,6)], nrow = 8, ncol = 1) -> tuc_eta_vec

# Get the estimated detection efficiency using inv logit
snake_run_year_trib_deteff[6:nrow(snake_run_year_trib_deteff), 1] <- exp(tuc_eta_vec)/(1 + exp(tuc_eta_vec))



# Now, get the Imnaha detection efficiency without any covariates
subset(snake_deteff_data, is.na(mean_discharge_cfs) & tributary == "Imnaha") -> imn_nodischarge

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




