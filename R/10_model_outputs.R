# 10 Plot model results

# In this script, we will be loading the JAGS objects that we saved in script 
# 09_fit_simulation to generate diagnostic plots and examine our model outputs.

# Load libraries
library(tidyverse)
library(lubridate)
library(janitor)
library(rjags)
library(jagsUI)
library(R2jags)
library(coda)

# Load states in the simulation
sim_states = c("mainstem, mouth to BON", 
               "mainstem, BON to MCN", 
               "mainstem, MCN to ICH or PRA", 
               "mainstem, PRA to RIS",
               "mainstem, ICH to LGR",
               "Deschutes River", 
               "John Day River", 
               "Tucannon River", 
               "Yakima River",
               "loss")

JAGS_obj <- readRDS(here::here("simulation", "JAGS_nocov_3chains_15kiter_5kburnin_4.rds"))
# Note: in the 'bugs' object (JAGS_obj$BUGSoutput), [7] are the chains, [10] is the summary
JAGS_obj$BUGSoutput[7]
mod_mcmc <- as.mcmc(JAGS_obj)
# plot(mod_mcmc[[1]])
plot(mod_mcmc)

# Some of these don't look good - estimates seem far too low
# [4,3], [6,2], [7,2], [8,5], [9,3]
# I think these are all states with only one "to" option. State 1 ([1,2]) is the only one where there's only one "to" state that looks okay

param_est <- as.data.frame(JAGS_obj$BUGSoutput[10])

colnames(param_est) <- gsub("summary.", "", colnames(param_est))

param_est %>% 
  rownames_to_column("parameter") %>% 
  mutate(in_sim = ifelse(parameter %in% c("b0_matrix[2,1]",
"b0_matrix[1,2]",
"b0_matrix[3,2]",
"b0_matrix[6,2]",
"b0_matrix[7,2]",
"b0_matrix[2,3]",
"b0_matrix[4,3]",
"b0_matrix[5,3]",
"b0_matrix[9,3]",
"b0_matrix[3,4]",
"b0_matrix[3,5]",
"b0_matrix[8,5]",
"b0_matrix[2,6]",
"b0_matrix[2,7]",
"b0_matrix[5,8]",
"b0_matrix[3,9]"), "yes", "no")) -> param_est

params_in_model <- c("b0_matrix[2,1]",
                     "b0_matrix[1,2]",
                     "b0_matrix[3,2]",
                     "b0_matrix[6,2]",
                     "b0_matrix[7,2]",
                     "b0_matrix[2,3]",
                     "b0_matrix[4,3]",
                     "b0_matrix[5,3]",
                     "b0_matrix[9,3]",
                     "b0_matrix[3,4]",
                     "b0_matrix[3,5]",
                     "b0_matrix[8,5]",
                     "b0_matrix[2,6]",
                     "b0_matrix[2,7]",
                     "b0_matrix[5,8]",
                     "b0_matrix[3,9]")

param_est %>% 
  subset(in_sim == "yes") -> param_est_subset

# Organize parameters back into a matrix
b0_matrix_mean <- matrix(param_est$mean[1:90], nrow = 10, ncol = 9)
colnames(b0_matrix_mean) <- sim_states[1:9]
rownames(b0_matrix_mean) <- sim_states[1:10]

# Create a matrix to store movement probabilities plus loss
movement_probs <- matrix(nrow = 10, ncol = 10)

rownames(movement_probs) <- sim_states
colnames(movement_probs) <- sim_states

##### Evaluate parameters in mlogit ####
for (i in 1:10){
  movement_probs[i,] <- c(exp(b0_matrix_mean[i,])/(1 + sum(exp(b0_matrix_mean[i,]))), 1 - sum(exp(b0_matrix_mean[i,])/(1 + sum(exp(b0_matrix_mean[i,])))))
}

# Check out these movement probabilities
round(movement_probs[1,],2)
round(movement_probs[2,],2)
round(movement_probs[3,],2)
round(movement_probs[4,],2)
round(movement_probs[5,],2)
round(movement_probs[6,],2)
round(movement_probs[7,],2)
round(movement_probs[8,],2)
round(movement_probs[9,],2)
round(movement_probs[10,],2)

# All of these are estimating that the probability of loss is zero?
# Which means that for any of the states where the only options are one other state or loss, it's 100% going
# to that state. 1, 4, 6, 7, 8, 9
# these are actually the states where the traces look the best.
# In every other state, the traces all show the same patterns, which makes sense because
# they're all supposed to be equal - so they track each other
# Okay - so if we can sort out why we're getting 0% chance of loss, we should be able
# to correct all of these traces and get the correct estimates

# So I think the problem is that we're allowing dmulti to have non-zero values
# for states that aren't possible - back to the exp(0) issue. So I think that what
# JAGS is doing is increasing the b0 values really high to make sure that the relative
# probability of going to any of the impossible states is close to zero. But the
# associated problem is that loss also goes to zero. So I need to fix that issue 
# to get these chains to converge

# This is what we should be getting for movement_probs[3,] (4 possible next states + loss:)
exp(1)/ (1 + sum(exp(rep(1,4))))
exp(1)/ (1 + sum(exp(rep(1,1))))
exp(1)/ (1 + sum(exp(rep(1,2))))


##### Plot simulation results - no covariates #####
simulation_plots_nocov <- function(JAGS_list){
  JAGS_runs_comp <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  for (i in 1:length(JAGS_list)){
    as.data.frame(JAGS_list[[i]]$BUGSoutput$summary) %>% 
      rownames_to_column("parameter") %>% 
      subset(parameter %in% params_in_model) %>% 
      dplyr::rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% 
      dplyr::select(parameter, mean, q2.5, q97.5) %>% 
      mutate(run = paste0(i)) -> summary
    
    JAGS_runs_comp %>% 
      bind_rows(., summary) -> JAGS_runs_comp
  }
  
  
  # Create the plot
  JAGS_runs_comp <- subset(JAGS_runs_comp, !(is.na(parameter))) 
  
  
  plot <- ggplot(JAGS_runs_comp, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(0.5, 1, 1.5), limits = c(0,2)) +
    xlab("Run")
  
  return(plot)
}

JAGS_600_list <- readRDS(here::here("simulation", "JAGS_nocov_600_list.rds"))
sim600_plots <- simulation_plots_nocov(JAGS_list = JAGS_600_list)
ggsave(here::here("simulation", "figures", "sim600_plots.png"), height = 6, width = 10, sim600_plots)

JAGS_1200_list <- readRDS(here::here("simulation", "JAGS_nocov_1200_list.rds"))
sim1200_plots <- simulation_plots_nocov(JAGS_list = JAGS_1200_list)
ggsave(here::here("simulation", "figures", "sim1200_plots.png"), height = 6, width = 10, sim1200_plots)

# JAGS_3000_list <- readRDS(here::here("simulation", "JAGS_nocov_3000_list.rds"))
# sim3000_plots <- simulation_plots_nocov(JAGS_list = JAGS_3000_list)
# ggsave(here::here("simulation", "figures", "sim3000_plots.png"), height = 6, width = 10, sim3000_plots)

##### Plot simulation results - covariates #####
JAGS_600_cov_list <- readRDS(here::here("simulation", "JAGS_cov_600_list.rds"))

# Return 5 facet-wrapped plots, one for each beta matrix
# The hlines used as referernce will have to be changed for each
simulation_plots_cov <- function(JAGS_list){
  JAGS_runs_comp <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  for (i in 1:length(JAGS_list)){
    as.data.frame(JAGS_list[[i]]$BUGSoutput$summary) %>% 
      rownames_to_column("parameter") %>% 
      subset(parameter %in% params_in_model) %>% 
      dplyr::rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% 
      dplyr::select(parameter, mean, q2.5, q97.5) %>% 
      mutate(run = paste0(i)) -> summary
    
    # Split the summary into five, one for each parameter
    summary %>% 
      subset(contains("b0"), parameter) -> b0_summary
    
    summary %>% 
      subset(contains("bflow"), parameter) -> bflow_summary
    
    summary %>% 
      subset(contains("btemp"), parameter) -> btemp_summary
    
    summary %>% 
      subset(contains("borigin"), parameter) -> borigin_summary
    
    summary %>% 
      subset(contains("brear"), parameter) -> brear_summary
      
    
    JAGS_runs_comp %>% 
      bind_rows(., summary) -> JAGS_runs_comp
  }
  
  
  # Create the plot
  JAGS_runs_comp <- subset(JAGS_runs_comp, !(is.na(parameter))) 
  
  
  plot <- ggplot(JAGS_runs_comp, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(0.5, 1, 1.5), limits = c(0,2)) +
    xlab("Run")
  
  return(plot)
}



