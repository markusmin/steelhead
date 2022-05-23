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

# JAGS_obj <- readRDS(here::here("simulation", "JAGS_nocov_3chains_15kiter_5kburnin_4.rds"))
# JAGS_obj <- out.jags
# Note: in the 'bugs' object (JAGS_obj$BUGSoutput), [7] are the chains, [10] is the summary
# JAGS_obj$BUGSoutput[7]
# mod_mcmc <- as.mcmc(JAGS_obj)
# plot(mod_mcmc)

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
exp(-50)/ (1 + sum(exp(-50)))
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
JAGS_obj <- JAGS_600_cov_list[[1]]
chains <- JAGS_obj$BUGSoutput[7]
mod_mcmc <- as.mcmc(JAGS_obj)
plot(mod_mcmc)

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



##### Plot results from covariate model, just one run #####
# Plot the traces
JAGS_obj <- out.jags
# Note: in the 'bugs' object (JAGS_obj$BUGSoutput), you can $ again to access various things
mod_mcmc <- as.mcmc(JAGS_obj)
mod_mcmc["b0_matrix"]
plot(mod_mcmc[[1]])
# plot(mod_mcmc)
mod_mcmc[[1]][,39] #b0, 4,3
mod_mcmc[[1]][,(180 + 39*3)] #borigin, 4,3,3
mod_mcmc[[1]][,(90 + 90 + 270 + 39*2)] # brear, 4,3,2

b_43 <- data.frame(b0 = mod_mcmc[[1]][,39] , borigin = mod_mcmc[[1]][,(180 + 39*3)], brear = mod_mcmc[[1]][,(90 + 90 + 270 + 39*2)], trace = seq(1,2000,1))
colnames(b_43) <- c("b0", "brear", "borigin", "iteration")
b_43 %>% 
  mutate(total = b0 + brear + borigin) %>% 
  pivot_longer(cols = c("b0", "brear", "borigin", "total")) -> b43_long

ggplot(b43_long, aes( x = iteration, y = value, color = name)) +
  geom_line() +
  geom_hline(yintercept = 3, lty = 2) +
  theme(legend.text = element_text(size = 12))



# b0, flow, origin, rear
plot(y = mod_mcmc[[1]][,2], x = seq(5001, 24991, by = 10))

b0_matrix_chains <- JAGS_obj$BUGSoutput$sims.list$b0_matrix
b0_43 <- data.frame(trace = b0_matrix_chains[,4,3][seq(1, length(b0_matrix_chains[,4,3]), by = 10)], sample = seq(1, length(b0_matrix_chains[,4,3])/10))
ggplot(b0_43, aes( x = sample, y = trace)) +
  geom_line()

brear_matrix_chains <- JAGS_obj$BUGSoutput$sims.list$brear_array
borigin_matrix_chains <- JAGS_obj$BUGSoutput$sims.list$borigin_matrix


# out.jags <- readRDS(here::here("simulation", "JAGS_cov_1200_onerun.rds"))
as.data.frame(out.jags$BUGSoutput[10]) %>% 
  rownames_to_column("parameter") -> param_est_1200

# Let's look at the parameters for moving from state 4
param_est_1200 %>% 
  filter(grepl("4,3", parameter)) -> from4
# This is from PRA to RIS state, which can only go back to MCN to ICH or PRA or loss

from4_means <- from4$summary.mean

##### Plot results from temperature model, just one run #####
# Plot the traces
JAGS_obj <- out.jags
# Note: in the 'bugs' object (JAGS_obj$BUGSoutput), you can $ again to access various things
mod_mcmc <- as.mcmc(JAGS_obj)
plot(mod_mcmc[[1]])

##### Plot simulation results - temperature model #####
JAGS_600_temp_list <- readRDS( here::here("simulation", "JAGS_cov_temp_600_list.rds"))
JAGS_obj <- JAGS_600_temp_list[[1]]
chains <- JAGS_obj$BUGSoutput[7]
mod_mcmc <- as.mcmc(JAGS_obj)
plot(mod_mcmc)


# Here, remove all of the tributary temperature effects - there is no temp data for the tributaries
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
                     "b0_matrix[3,9]",
                     "btemp_matrix[2,1]",
                     "btemp_matrix[1,2]",
                     "btemp_matrix[3,2]",
                     # "btemp_matrix[6,2]",
                     # "btemp_matrix[7,2]",
                     "btemp_matrix[2,3]",
                     "btemp_matrix[4,3]",
                     "btemp_matrix[5,3]",
                     # "btemp_matrix[9,3]",
                     "btemp_matrix[3,4]",
                     "btemp_matrix[3,5]",
                     # "btemp_matrix[8,5]",
                     "btemp_matrix[2,6]",
                     "btemp_matrix[2,7]",
                     "btemp_matrix[5,8]",
                     "btemp_matrix[3,9]")

# Return 2 facet-wrapped plots, one for each beta matrix
# The hlines used as referernce will have to be changed for each
simulation_plots_temp <- function(JAGS_list, params_in_model){
  JAGS_runs_comp_b0 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_btemp <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  for (i in 1:length(JAGS_list)){
    as.data.frame(JAGS_list[[i]]$BUGSoutput$summary) %>% 
      rownames_to_column("parameter") %>% 
      subset(parameter %in% params_in_model) %>% 
      dplyr::rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% 
      dplyr::select(parameter, mean, q2.5, q97.5) %>% 
      mutate(run = paste0(i)) -> summary
    
    # Change run to a factor for plotting
    summary %>% 
      mutate(run = factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) -> summary
    
    # Split the summary into five, one for each parameter
    summary %>% 
      filter(grepl("b0", parameter)) -> b0_summary
    summary %>% 
      filter(grepl("btemp", parameter))  -> btemp_summary
    
    JAGS_runs_comp_b0 %>% 
      bind_rows(., b0_summary) -> JAGS_runs_comp_b0
    
    JAGS_runs_comp_btemp %>% 
      bind_rows(., btemp_summary) -> JAGS_runs_comp_btemp
  }
  
  
  # Create the plot
  JAGS_runs_comp_b0 <- subset(JAGS_runs_comp_b0, !(is.na(parameter))) 
  JAGS_runs_comp_btemp <- subset(JAGS_runs_comp_btemp, !(is.na(parameter))) 
  
  b0_plot <- ggplot(JAGS_runs_comp_b0, aes(y = mean, x = run)) +
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
  
  btemp_plot <- ggplot(JAGS_runs_comp_btemp, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  # Annotations
  actual_values <- data.frame(parameter = c("btemp_matrix[1,2]",
                                            "btemp_matrix[2,1]",
                                            "btemp_matrix[2,3]",
                                            "btemp_matrix[2,6]",
                                            "btemp_matrix[2,7]",
                                              "btemp_matrix[3,2]",
                                            "btemp_matrix[3,4]",
                                            "btemp_matrix[3,5]",
                                            "btemp_matrix[3,9]",
                                              "btemp_matrix[4,3]",
                                              "btemp_matrix[5,3]",
                                              "btemp_matrix[5,8]"
                                              ),
                              yint = c(0,0,0.5,0,0,-0.5,0.3,1,0,0,0.2,0))
  
  btemp_plot +
    geom_hline(data = actual_values, aes(yintercept = yint), lty = 2) -> btemp_plot
  
  return(list(b0_plot, btemp_plot))
}

JAGS_600_temp_list <- readRDS( here::here("simulation", "JAGS_cov_temp_600_list.rds"))
sim600_temp_plots <- simulation_plots_temp(JAGS_list = JAGS_600_temp_list)
# sim600_temp_plots[[1]]
# sim600_temp_plots[[2]]
ggsave(here::here("simulation", "figures", "sim600_temp_plots_b0.png"), height = 6, width = 10, sim600_temp_plots[[1]])
ggsave(here::here("simulation", "figures", "sim600_temp_plots_btemp.png"), height = 6, width = 10, sim600_temp_plots[[2]])
ggsave(here::here("simulation", "figures", "sim600_temp_plots_btemp_bad.png"), height = 6, width = 10, sim600_temp_plots[[3]])








##### Plot simulation results - temperature + flow model #####
# Here, remove all of the tributary temperature effects - there is no temp data for the tributaries
params_in_temp_flow_model <- c("b0_matrix[2,1]",
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
                     "b0_matrix[3,9]",
                     "btemp_matrix[2,1]",
                     "btemp_matrix[1,2]",
                     "btemp_matrix[3,2]",
                     # "btemp_matrix[6,2]",
                     # "btemp_matrix[7,2]",
                     "btemp_matrix[2,3]",
                     "btemp_matrix[4,3]",
                     "btemp_matrix[5,3]",
                     # "btemp_matrix[9,3]",
                     "btemp_matrix[3,4]",
                     "btemp_matrix[3,5]",
                     # "btemp_matrix[8,5]",
                     "btemp_matrix[2,6]",
                     "btemp_matrix[2,7]",
                     "btemp_matrix[5,8]",
                     "btemp_matrix[3,9]",
                     "bflow_matrix[2,1]",
                     "bflow_matrix[1,2]",
                     "bflow_matrix[3,2]",
                     # "bflow_matrix[6,2]",
                     # "bflow_matrix[7,2]",
                     "bflow_matrix[2,3]",
                     "bflow_matrix[4,3]",
                     "bflow_matrix[5,3]",
                     # "bflow_matrix[9,3]",
                     "bflow_matrix[3,4]",
                     "bflow_matrix[3,5]",
                     # "bflow_matrix[8,5]",
                     "bflow_matrix[2,6]",
                     "bflow_matrix[2,7]",
                     "bflow_matrix[5,8]",
                     "bflow_matrix[3,9]")

# Return 2 facet-wrapped plots, one for each beta matrix
# The hlines used as referernce will have to be changed for each
simulation_plots_temp_flow <- function(JAGS_list, params_in_model){
  JAGS_runs_comp_b0 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_btemp <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_bflow <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  for (i in 1:length(JAGS_list)){
    as.data.frame(JAGS_list[[i]]$BUGSoutput$summary) %>% 
      rownames_to_column("parameter") %>% 
      subset(parameter %in% params_in_model) %>% 
      dplyr::rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% 
      dplyr::select(parameter, mean, q2.5, q97.5) %>% 
      mutate(run = paste0(i)) -> summary
    
    # Change run to a factor for plotting
    summary %>% 
      mutate(run = factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) -> summary
    
    # Split the summary into five, one for each parameter
    summary %>% 
      filter(grepl("b0", parameter)) -> b0_summary
    summary %>% 
      filter(grepl("btemp", parameter))  -> btemp_summary
    summary %>% 
      filter(grepl("bflow", parameter))  -> bflow_summary
    
    JAGS_runs_comp_b0 %>% 
      bind_rows(., b0_summary) -> JAGS_runs_comp_b0
    
    JAGS_runs_comp_btemp %>% 
      bind_rows(., btemp_summary) -> JAGS_runs_comp_btemp
    
    JAGS_runs_comp_bflow %>% 
      bind_rows(., bflow_summary) -> JAGS_runs_comp_bflow
  }
  
  
  # Create the plot
  JAGS_runs_comp_b0 <- subset(JAGS_runs_comp_b0, !(is.na(parameter))) 
  JAGS_runs_comp_btemp <- subset(JAGS_runs_comp_btemp, !(is.na(parameter))) 
  JAGS_runs_comp_bflow <- subset(JAGS_runs_comp_bflow, !(is.na(parameter))) 
  
  b0_plot <- ggplot(JAGS_runs_comp_b0, aes(y = mean, x = run)) +
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
  
  btemp_plot <- ggplot(JAGS_runs_comp_btemp, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  bflow_plot <- ggplot(JAGS_runs_comp_bflow, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  # Annotations
  actual_values_temp <- data.frame(parameter = c("btemp_matrix[1,2]",
                                            "btemp_matrix[2,1]",
                                            "btemp_matrix[2,3]",
                                            "btemp_matrix[2,6]",
                                            "btemp_matrix[2,7]",
                                            "btemp_matrix[3,2]",
                                            "btemp_matrix[3,4]",
                                            "btemp_matrix[3,5]",
                                            "btemp_matrix[3,9]",
                                            "btemp_matrix[4,3]",
                                            "btemp_matrix[5,3]",
                                            "btemp_matrix[5,8]"
  ),
  yint = c(0,0,0.5,0,0,-0.5,0.3,0,0,0,0.2,0))
  
  actual_values_flow <- data.frame(parameter = c("bflow_matrix[1,2]",
                                                 "bflow_matrix[2,1]",
                                                 "bflow_matrix[2,3]",
                                                 "bflow_matrix[2,6]",
                                                 "bflow_matrix[2,7]",
                                                 "bflow_matrix[3,2]",
                                                 "bflow_matrix[3,4]",
                                                 "bflow_matrix[3,5]",
                                                 "bflow_matrix[3,9]",
                                                 "bflow_matrix[4,3]",
                                                 "bflow_matrix[5,3]",
                                                 "bflow_matrix[5,8]"
  ),
  yint = c(0,0.5,0,0,0,0.05,0,-0.3,0,0,0,0))
  
  btemp_plot +
    geom_hline(data = actual_values_temp, aes(yintercept = yint), lty = 2) -> btemp_plot
  
  bflow_plot +
    geom_hline(data = actual_values_flow, aes(yintercept = yint), lty = 2) -> bflow_plot
  
  return(list(b0_plot, btemp_plot, bflow_plot))
}






##### Plot simulation results - ORIGIN + REAR model #####

params_in_origin_rear_model <- c("b0_matrix[2,1]",
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
                     "b0_matrix[3,9]",
                     
                     "brear_matrix[2,1]",
                     "brear_matrix[1,2]",
                     "brear_matrix[3,2]",
                     "brear_matrix[6,2]",
                     "brear_matrix[7,2]",
                     "brear_matrix[2,3]",
                     "brear_matrix[4,3]",
                     "brear_matrix[5,3]",
                     "brear_matrix[9,3]",
                     "brear_matrix[3,4]",
                     "brear_matrix[3,5]",
                     "brear_matrix[8,5]",
                     "brear_matrix[2,6]",
                     "brear_matrix[2,7]",
                     "brear_matrix[5,8]",
                     "brear_matrix[3,9]",
                     
                     "borigin1_matrix[2,1]",
                     "borigin1_matrix[1,2]",
                     "borigin1_matrix[3,2]",
                     "borigin1_matrix[6,2]",
                     "borigin1_matrix[7,2]",
                     "borigin1_matrix[2,3]",
                     "borigin1_matrix[4,3]",
                     "borigin1_matrix[5,3]",
                     "borigin1_matrix[9,3]",
                     "borigin1_matrix[3,4]",
                     "borigin1_matrix[3,5]",
                     "borigin1_matrix[8,5]",
                     "borigin1_matrix[2,6]",
                     "borigin1_matrix[2,7]",
                     "borigin1_matrix[5,8]",
                     "borigin1_matrix[3,9]",
                     
                     "borigin2_matrix[2,1]",
                     "borigin2_matrix[1,2]",
                     "borigin2_matrix[3,2]",
                     "borigin2_matrix[6,2]",
                     "borigin2_matrix[7,2]",
                     "borigin2_matrix[2,3]",
                     "borigin2_matrix[4,3]",
                     "borigin2_matrix[5,3]",
                     "borigin2_matrix[9,3]",
                     "borigin2_matrix[3,4]",
                     "borigin2_matrix[3,5]",
                     "borigin2_matrix[8,5]",
                     "borigin2_matrix[2,6]",
                     "borigin2_matrix[2,7]",
                     "borigin2_matrix[5,8]",
                     "borigin2_matrix[3,9]")

# Return 4 facet-wrapped plots, one for each beta matrix
# The hlines used as referernce will have to be changed for each
simulation_plots_origin_rear <- function(JAGS_list, params_in_model){
  JAGS_runs_comp_b0 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_brear <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_borigin1 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_borigin2 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  
  for (i in 1:length(JAGS_list)){
    as.data.frame(JAGS_list[[i]]$BUGSoutput$summary) %>% 
      rownames_to_column("parameter") %>% 
      subset(parameter %in% params_in_model) %>% 
      dplyr::rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% 
      dplyr::select(parameter, mean, q2.5, q97.5) %>% 
      mutate(run = paste0(i)) -> summary
    
    # Change run to a factor for plotting
    summary %>% 
      mutate(run = factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) -> summary
    
    # Split the summary into five, one for each parameter
    summary %>% 
      filter(grepl("b0", parameter)) -> b0_summary
    summary %>% 
      filter(grepl("brear", parameter))  -> brear_summary
    summary %>% 
      filter(grepl("borigin1", parameter))  -> borigin1_summary
    summary %>% 
      filter(grepl("borigin2", parameter))  -> borigin2_summary
    
    JAGS_runs_comp_b0 %>% 
      bind_rows(., b0_summary) -> JAGS_runs_comp_b0
    
    JAGS_runs_comp_brear %>% 
      bind_rows(., brear_summary) -> JAGS_runs_comp_brear
    
    JAGS_runs_comp_borigin1 %>% 
      bind_rows(., borigin1_summary) -> JAGS_runs_comp_borigin1
    
    JAGS_runs_comp_borigin2 %>% 
      bind_rows(., borigin2_summary) -> JAGS_runs_comp_borigin2
  }
  
  
  # Create the plot
  JAGS_runs_comp_b0 <- subset(JAGS_runs_comp_b0, !(is.na(parameter))) 
  JAGS_runs_comp_brear <- subset(JAGS_runs_comp_brear, !(is.na(parameter))) 
  JAGS_runs_comp_borigin1 <- subset(JAGS_runs_comp_borigin1, !(is.na(parameter))) 
  JAGS_runs_comp_borigin2 <- subset(JAGS_runs_comp_borigin2, !(is.na(parameter))) 
  
  
  b0_plot <- ggplot(JAGS_runs_comp_b0, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(0, 1, 2), limits = c(-1,3)) +
    xlab("Run")
  
  brear_plot <- ggplot(JAGS_runs_comp_brear, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  borigin1_plot <- ggplot(JAGS_runs_comp_borigin1, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  borigin2_plot <- ggplot(JAGS_runs_comp_borigin2, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  # Annotations
  actual_values_brear <- data.frame(parameter = c("brear_matrix[1,2]",
                                            "brear_matrix[2,1]",
                                            "brear_matrix[2,3]",
                                            "brear_matrix[2,6]",
                                            "brear_matrix[2,7]",
                                            "brear_matrix[3,2]",
                                            "brear_matrix[3,4]",
                                            "brear_matrix[3,5]",
                                            "brear_matrix[3,9]",
                                            "brear_matrix[4,3]",
                                            "brear_matrix[5,3]",
                                            "brear_matrix[5,8]",
                                            "brear_matrix[6,2]",
                                            "brear_matrix[7,2]",
                                            "brear_matrix[8,5]",
                                            "brear_matrix[9,3]"
                                            
  ),
  yint = c(0,0,0,0.5,0,0,0,0,-0.2,0,0,0,-0.5,0,0,0))
  
  actual_values_borigin1 <- data.frame(parameter = c("borigin1_matrix[1,2]",
                                                  "borigin1_matrix[2,1]",
                                                  "borigin1_matrix[2,3]",
                                                  "borigin1_matrix[2,6]",
                                                  "borigin1_matrix[2,7]",
                                                  "borigin1_matrix[3,2]",
                                                  "borigin1_matrix[3,4]",
                                                  "borigin1_matrix[3,5]",
                                                  "borigin1_matrix[3,9]",
                                                  "borigin1_matrix[4,3]",
                                                  "borigin1_matrix[5,3]",
                                                  "borigin1_matrix[5,8]",
                                                  "borigin1_matrix[6,2]",
                                                  "borigin1_matrix[7,2]",
                                                  "borigin1_matrix[8,5]",
                                                  "borigin1_matrix[9,3]"
  ),
  yint = c(0,0,-0.5,0.25,0.5,0.5,0,-0.5,-0.25,0,0.25,-0.25,0,-0.5,0.25,0.25))
  
  actual_values_borigin2 <- data.frame(parameter = c("borigin2_matrix[1,2]",
                                                  "borigin2_matrix[2,1]",
                                                  "borigin2_matrix[2,3]",
                                                  "borigin2_matrix[2,6]",
                                                  "borigin2_matrix[2,7]",
                                                  "borigin2_matrix[3,2]",
                                                  "borigin2_matrix[3,4]",
                                                  "borigin2_matrix[3,5]",
                                                  "borigin2_matrix[3,9]",
                                                  "borigin2_matrix[4,3]",
                                                  "borigin2_matrix[5,3]",
                                                  "borigin2_matrix[5,8]",
                                                  "borigin2_matrix[6,2]",
                                                  "borigin2_matrix[7,2]",
                                                  "borigin2_matrix[8,5]",
                                                  "borigin2_matrix[9,3]"
  ),
  yint = c(0,0,0.25,-0.125,-0.25,-0.25,0,0,0.5,0,0.25,-0.25,0,0.25,0.25,-0.5))
  
  brear_plot +
    geom_hline(data = actual_values_brear, aes(yintercept = yint), lty = 2) -> brear_plot
  
  borigin1_plot +
    geom_hline(data = actual_values_borigin1, aes(yintercept = yint), lty = 2) -> borigin1_plot
  
  borigin2_plot +
    geom_hline(data = actual_values_borigin2, aes(yintercept = yint), lty = 2) -> borigin2_plot
  
  
  return(list(b0_plot, brear_plot, borigin1_plot, borigin2_plot))
}

##### Plot simulation results - origin only model #####
origin_params <- c("b0_matrix[2,1]",
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
                     "b0_matrix[3,9]",
                     "borigin1_matrix[2,1]",
                     "borigin1_matrix[1,2]",
                     "borigin1_matrix[3,2]",
                     "borigin1_matrix[6,2]",
                     "borigin1_matrix[7,2]",
                     "borigin1_matrix[2,3]",
                     "borigin1_matrix[4,3]",
                     "borigin1_matrix[5,3]",
                     "borigin1_matrix[9,3]",
                     "borigin1_matrix[3,4]",
                     "borigin1_matrix[3,5]",
                     "borigin1_matrix[8,5]",
                     "borigin1_matrix[2,6]",
                     "borigin1_matrix[2,7]",
                     "borigin1_matrix[5,8]",
                     "borigin1_matrix[3,9]",
                     "borigin2_matrix[2,1]",
                     "borigin2_matrix[1,2]",
                     "borigin2_matrix[3,2]",
                     "borigin2_matrix[6,2]",
                     "borigin2_matrix[7,2]",
                     "borigin2_matrix[2,3]",
                     "borigin2_matrix[4,3]",
                     "borigin2_matrix[5,3]",
                     "borigin2_matrix[9,3]",
                     "borigin2_matrix[3,4]",
                     "borigin2_matrix[3,5]",
                     "borigin2_matrix[8,5]",
                     "borigin2_matrix[2,6]",
                     "borigin2_matrix[2,7]",
                     "borigin2_matrix[5,8]",
                     "borigin2_matrix[3,9]")

# Return 2 facet-wrapped plots, one for each beta matrix
# The hlines used as referernce will have to be changed for each
simulation_plots_origin <- function(JAGS_list, params_in_model = origin_params){
  JAGS_runs_comp_b0 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_borigin1 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_borigin2 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  for (i in 1:length(JAGS_list)){
    as.data.frame(JAGS_list[[i]]$BUGSoutput$summary) %>% 
      rownames_to_column("parameter") %>% 
      subset(parameter %in% params_in_model) %>% 
      dplyr::rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% 
      dplyr::select(parameter, mean, q2.5, q97.5) %>% 
      mutate(run = paste0(i)) -> summary
    
    # Change run to a factor for plotting
    summary %>% 
      mutate(run = factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) -> summary
    
    # Split the summary into five, one for each parameter
    summary %>% 
      filter(grepl("b0", parameter)) -> b0_summary
    summary %>% 
      filter(grepl("borigin1", parameter))  -> borigin1_summary
    summary %>% 
      filter(grepl("borigin2", parameter))  -> borigin2_summary
    
    JAGS_runs_comp_b0 %>% 
      bind_rows(., b0_summary) -> JAGS_runs_comp_b0
    
    JAGS_runs_comp_borigin1 %>% 
      bind_rows(., borigin1_summary) -> JAGS_runs_comp_borigin1
    
    JAGS_runs_comp_borigin2 %>% 
      bind_rows(., borigin2_summary) -> JAGS_runs_comp_borigin2
  }
  
  
  # Create the plot
  JAGS_runs_comp_b0 <- subset(JAGS_runs_comp_b0, !(is.na(parameter))) 
  JAGS_runs_comp_borigin1 <- subset(JAGS_runs_comp_borigin1, !(is.na(parameter))) 
  JAGS_runs_comp_borigin2 <- subset(JAGS_runs_comp_borigin2, !(is.na(parameter))) 
  
  
  b0_plot <- ggplot(JAGS_runs_comp_b0, aes(y = mean, x = run)) +
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
  
  borigin1_plot <- ggplot(JAGS_runs_comp_borigin1, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0 ,1), limits = c(-2,2)) +
    xlab("Run")
  
  borigin2_plot <- ggplot(JAGS_runs_comp_borigin2, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0 ,1), limits = c(-2,2)) +
    xlab("Run")
  
  # Annotations
  
  ### start old values
  borigin1_values <- data.frame(parameter = c("borigin1_matrix[2,1]",
                                              "borigin1_matrix[1,2]",
                                              "borigin1_matrix[3,2]",
                                              "borigin1_matrix[6,2]",
                                              "borigin1_matrix[7,2]",
                                              "borigin1_matrix[2,3]",
                                              "borigin1_matrix[4,3]",
                                              "borigin1_matrix[5,3]",
                                              "borigin1_matrix[9,3]",
                                              "borigin1_matrix[3,4]",
                                              "borigin1_matrix[3,5]",
                                              "borigin1_matrix[8,5]",
                                              "borigin1_matrix[2,6]",
                                              "borigin1_matrix[2,7]",
                                              "borigin1_matrix[5,8]",
                                              "borigin1_matrix[3,9]"
  ),
  yint = c(0, 0, 1, 0, -1, -1, 0, 0.5, 0, 0, -1, 0.5, 0.25, 1, -0.5, 0))
  
  borigin2_values <- data.frame(parameter = c("borigin2_matrix[2,1]",
                                              "borigin2_matrix[1,2]",
                                              "borigin2_matrix[3,2]",
                                              "borigin2_matrix[6,2]",
                                              "borigin2_matrix[7,2]",
                                              "borigin2_matrix[2,3]",
                                              "borigin2_matrix[4,3]",
                                              "borigin2_matrix[5,3]",
                                              "borigin2_matrix[9,3]",
                                              "borigin2_matrix[3,4]",
                                              "borigin2_matrix[3,5]",
                                              "borigin2_matrix[8,5]",
                                              "borigin2_matrix[2,6]",
                                              "borigin2_matrix[2,7]",
                                              "borigin2_matrix[5,8]",
                                              "borigin2_matrix[3,9]"
  ),
  yint = c(0, 0, 0, 0, 0, 0, 0, 0.5, -1, 0, -0.5, 0.5, 0, 0, -0.5, 1))
  
  ### end old values
  
  
  actual_values_borigin1 <- data.frame(parameter = c("borigin1_matrix[1,2]",
                                                     "borigin1_matrix[2,1]",
                                                     "borigin1_matrix[2,3]",
                                                     "borigin1_matrix[2,6]",
                                                     "borigin1_matrix[2,7]",
                                                     "borigin1_matrix[3,2]",
                                                     "borigin1_matrix[3,4]",
                                                     "borigin1_matrix[3,5]",
                                                     "borigin1_matrix[3,9]",
                                                     "borigin1_matrix[4,3]",
                                                     "borigin1_matrix[5,3]",
                                                     "borigin1_matrix[5,8]",
                                                     "borigin1_matrix[6,2]",
                                                     "borigin1_matrix[7,2]",
                                                     "borigin1_matrix[8,5]",
                                                     "borigin1_matrix[9,3]"
  ),
  yint = c(0,0,-0.5,0.25,0.5,0.5,0,-0.5,-0.25,0,0.25,-0.25,0,-0.5,0.25,0.25))
  
  actual_values_borigin2 <- data.frame(parameter = c("borigin2_matrix[1,2]",
                                                     "borigin2_matrix[2,1]",
                                                     "borigin2_matrix[2,3]",
                                                     "borigin2_matrix[2,6]",
                                                     "borigin2_matrix[2,7]",
                                                     "borigin2_matrix[3,2]",
                                                     "borigin2_matrix[3,4]",
                                                     "borigin2_matrix[3,5]",
                                                     "borigin2_matrix[3,9]",
                                                     "borigin2_matrix[4,3]",
                                                     "borigin2_matrix[5,3]",
                                                     "borigin2_matrix[5,8]",
                                                     "borigin2_matrix[6,2]",
                                                     "borigin2_matrix[7,2]",
                                                     "borigin2_matrix[8,5]",
                                                     "borigin2_matrix[9,3]"
  ),
  yint = c(0,0,0.25,-0.125,-0.25,-0.25,0,0,0.5,0,0.25,-0.25,0,0.25,0.25,-0.5))
  
  borigin1_plot +
    geom_hline(data = actual_values_borigin1, aes(yintercept = yint), lty = 2) -> borigin1_plot
  
  borigin2_plot +
    geom_hline(data = actual_values_borigin2, aes(yintercept = yint), lty = 2) -> borigin2_plot
  
  
  return(list(b0_plot, borigin1_plot, borigin2_plot))
}



##### Plot simulation results - all 4 covariates #####

# Here, remove all of the tributary temperature effects - there is no temp data for the tributaries
params_in_all4_model <- c("b0_matrix[2,1]",
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
                               "b0_matrix[3,9]",
                               "btemp_matrix[2,1]",
                               "btemp_matrix[1,2]",
                               "btemp_matrix[3,2]",
                               # "btemp_matrix[6,2]",
                               # "btemp_matrix[7,2]",
                               "btemp_matrix[2,3]",
                               "btemp_matrix[4,3]",
                               "btemp_matrix[5,3]",
                               # "btemp_matrix[9,3]",
                               "btemp_matrix[3,4]",
                               "btemp_matrix[3,5]",
                               # "btemp_matrix[8,5]",
                               "btemp_matrix[2,6]",
                               "btemp_matrix[2,7]",
                               "btemp_matrix[5,8]",
                               "btemp_matrix[3,9]",
                               "bflow_matrix[2,1]",
                               "bflow_matrix[1,2]",
                               "bflow_matrix[3,2]",
                               # "bflow_matrix[6,2]",
                               # "bflow_matrix[7,2]",
                               "bflow_matrix[2,3]",
                               "bflow_matrix[4,3]",
                               "bflow_matrix[5,3]",
                               # "bflow_matrix[9,3]",
                               "bflow_matrix[3,4]",
                               "bflow_matrix[3,5]",
                               # "bflow_matrix[8,5]",
                               "bflow_matrix[2,6]",
                               "bflow_matrix[2,7]",
                               "bflow_matrix[5,8]",
                               "bflow_matrix[3,9]",
                          
                          
                          "brear_matrix[2,1]",
                          "brear_matrix[1,2]",
                          "brear_matrix[3,2]",
                          "brear_matrix[6,2]",
                          "brear_matrix[7,2]",
                          "brear_matrix[2,3]",
                          "brear_matrix[4,3]",
                          "brear_matrix[5,3]",
                          "brear_matrix[9,3]",
                          "brear_matrix[3,4]",
                          "brear_matrix[3,5]",
                          "brear_matrix[8,5]",
                          "brear_matrix[2,6]",
                          "brear_matrix[2,7]",
                          "brear_matrix[5,8]",
                          "brear_matrix[3,9]",
                          
                          "borigin1_matrix[2,1]",
                          "borigin1_matrix[1,2]",
                          "borigin1_matrix[3,2]",
                          "borigin1_matrix[6,2]",
                          "borigin1_matrix[7,2]",
                          "borigin1_matrix[2,3]",
                          "borigin1_matrix[4,3]",
                          "borigin1_matrix[5,3]",
                          "borigin1_matrix[9,3]",
                          "borigin1_matrix[3,4]",
                          "borigin1_matrix[3,5]",
                          "borigin1_matrix[8,5]",
                          "borigin1_matrix[2,6]",
                          "borigin1_matrix[2,7]",
                          "borigin1_matrix[5,8]",
                          "borigin1_matrix[3,9]",
                          
                          "borigin2_matrix[2,1]",
                          "borigin2_matrix[1,2]",
                          "borigin2_matrix[3,2]",
                          "borigin2_matrix[6,2]",
                          "borigin2_matrix[7,2]",
                          "borigin2_matrix[2,3]",
                          "borigin2_matrix[4,3]",
                          "borigin2_matrix[5,3]",
                          "borigin2_matrix[9,3]",
                          "borigin2_matrix[3,4]",
                          "borigin2_matrix[3,5]",
                          "borigin2_matrix[8,5]",
                          "borigin2_matrix[2,6]",
                          "borigin2_matrix[2,7]",
                          "borigin2_matrix[5,8]",
                          "borigin2_matrix[3,9]")

# Return 2 facet-wrapped plots, one for each beta matrix
# The hlines used as referernce will have to be changed for each
simulation_plots_all4 <- function(JAGS_list, params_in_model){
  JAGS_runs_comp_b0 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_btemp <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_bflow <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_brear <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_borigin1 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  JAGS_runs_comp_borigin2 <- data.frame("parameter" = NA, "mean" = NA, "q2.5" = NA, "q97.5" = NA)
  
  for (i in 1:length(JAGS_list)){
    as.data.frame(JAGS_list[[i]]$BUGSoutput$summary) %>% 
      rownames_to_column("parameter") %>% 
      subset(parameter %in% params_in_model) %>% 
      dplyr::rename(q2.5 = `2.5%`, q97.5 = `97.5%`) %>% 
      dplyr::select(parameter, mean, q2.5, q97.5) %>% 
      mutate(run = paste0(i)) -> summary
    
    # Change run to a factor for plotting
    summary %>% 
      mutate(run = factor(run, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))) -> summary
    
    # Split the summary into six, one for each parameter
    summary %>% 
      filter(grepl("b0", parameter)) -> b0_summary
    summary %>% 
      filter(grepl("brear", parameter))  -> brear_summary
    summary %>% 
      filter(grepl("borigin1", parameter))  -> borigin1_summary
    summary %>% 
      filter(grepl("borigin2", parameter))  -> borigin2_summary
    summary %>% 
      filter(grepl("btemp", parameter))  -> btemp_summary
    summary %>% 
      filter(grepl("bflow", parameter))  -> bflow_summary
    
    JAGS_runs_comp_b0 %>% 
      bind_rows(., b0_summary) -> JAGS_runs_comp_b0
    
    JAGS_runs_comp_brear %>% 
      bind_rows(., brear_summary) -> JAGS_runs_comp_brear
    
    JAGS_runs_comp_borigin1 %>% 
      bind_rows(., borigin1_summary) -> JAGS_runs_comp_borigin1
    
    JAGS_runs_comp_borigin2 %>% 
      bind_rows(., borigin2_summary) -> JAGS_runs_comp_borigin2
    
    JAGS_runs_comp_btemp %>% 
      bind_rows(., btemp_summary) -> JAGS_runs_comp_btemp
    
    JAGS_runs_comp_bflow %>% 
      bind_rows(., bflow_summary) -> JAGS_runs_comp_bflow
  }
  
  JAGS_runs_comp_b0 <- subset(JAGS_runs_comp_b0, !(is.na(parameter))) 
  JAGS_runs_comp_btemp <- subset(JAGS_runs_comp_btemp, !(is.na(parameter))) 
  JAGS_runs_comp_bflow <- subset(JAGS_runs_comp_bflow, !(is.na(parameter))) 
  JAGS_runs_comp_brear <- subset(JAGS_runs_comp_brear, !(is.na(parameter))) 
  JAGS_runs_comp_borigin1 <- subset(JAGS_runs_comp_borigin1, !(is.na(parameter))) 
  JAGS_runs_comp_borigin2 <- subset(JAGS_runs_comp_borigin2, !(is.na(parameter))) 
  
  
  # Create the plots
  
  
  brear_plot <- ggplot(JAGS_runs_comp_brear, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  borigin1_plot <- ggplot(JAGS_runs_comp_borigin1, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  borigin2_plot <- ggplot(JAGS_runs_comp_borigin2, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  # Annotations
  actual_values_brear <- data.frame(parameter = c("brear_matrix[1,2]",
                                                  "brear_matrix[2,1]",
                                                  "brear_matrix[2,3]",
                                                  "brear_matrix[2,6]",
                                                  "brear_matrix[2,7]",
                                                  "brear_matrix[3,2]",
                                                  "brear_matrix[3,4]",
                                                  "brear_matrix[3,5]",
                                                  "brear_matrix[3,9]",
                                                  "brear_matrix[4,3]",
                                                  "brear_matrix[5,3]",
                                                  "brear_matrix[5,8]",
                                                  "brear_matrix[6,2]",
                                                  "brear_matrix[7,2]",
                                                  "brear_matrix[8,5]",
                                                  "brear_matrix[9,3]"
                                                  
  ),
  yint = c(0,0,0,0.5,0,0,0,0,-0.2,0,0,0,-0.5,0,0,0))
  
  actual_values_borigin1 <- data.frame(parameter = c("borigin1_matrix[1,2]",
                                                     "borigin1_matrix[2,1]",
                                                     "borigin1_matrix[2,3]",
                                                     "borigin1_matrix[2,6]",
                                                     "borigin1_matrix[2,7]",
                                                     "borigin1_matrix[3,2]",
                                                     "borigin1_matrix[3,4]",
                                                     "borigin1_matrix[3,5]",
                                                     "borigin1_matrix[3,9]",
                                                     "borigin1_matrix[4,3]",
                                                     "borigin1_matrix[5,3]",
                                                     "borigin1_matrix[5,8]",
                                                     "borigin1_matrix[6,2]",
                                                     "borigin1_matrix[7,2]",
                                                     "borigin1_matrix[8,5]",
                                                     "borigin1_matrix[9,3]"
  ),
  yint = c(0,0,-0.5,0.25,0.5,0.5,0,-0.5,-0.25,0,0.25,-0.25,0,-0.5,0.25,0.25))
  
  actual_values_borigin2 <- data.frame(parameter = c("borigin2_matrix[1,2]",
                                                     "borigin2_matrix[2,1]",
                                                     "borigin2_matrix[2,3]",
                                                     "borigin2_matrix[2,6]",
                                                     "borigin2_matrix[2,7]",
                                                     "borigin2_matrix[3,2]",
                                                     "borigin2_matrix[3,4]",
                                                     "borigin2_matrix[3,5]",
                                                     "borigin2_matrix[3,9]",
                                                     "borigin2_matrix[4,3]",
                                                     "borigin2_matrix[5,3]",
                                                     "borigin2_matrix[5,8]",
                                                     "borigin2_matrix[6,2]",
                                                     "borigin2_matrix[7,2]",
                                                     "borigin2_matrix[8,5]",
                                                     "borigin2_matrix[9,3]"
  ),
  yint = c(0,0,0.25,-0.125,-0.25,-0.25,0,0,0.5,0,0.25,-0.25,0,0.25,0.25,-0.5))
  
  brear_plot +
    geom_hline(data = actual_values_brear, aes(yintercept = yint), lty = 2) -> brear_plot
  
  borigin1_plot +
    geom_hline(data = actual_values_borigin1, aes(yintercept = yint), lty = 2) -> borigin1_plot
  
  borigin2_plot +
    geom_hline(data = actual_values_borigin2, aes(yintercept = yint), lty = 2) -> borigin2_plot
  
  
  b0_plot <- ggplot(JAGS_runs_comp_b0, aes(y = mean, x = run)) +
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
  
  btemp_plot <- ggplot(JAGS_runs_comp_btemp, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  bflow_plot <- ggplot(JAGS_runs_comp_bflow, aes(y = mean, x = run)) +
    geom_point() +
    geom_linerange(aes(ymin = q2.5, ymax = q97.5)) +
    # geom_hline(yintercept = 1, lty = 2) +
    facet_wrap(~parameter) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          strip.text.x = element_text(size = 12)) +
    ylab("Estimate") +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-2, 2)) +
    xlab("Run")
  
  # Annotations
  actual_values_temp <- data.frame(parameter = c("btemp_matrix[1,2]",
                                                 "btemp_matrix[2,1]",
                                                 "btemp_matrix[2,3]",
                                                 "btemp_matrix[2,6]",
                                                 "btemp_matrix[2,7]",
                                                 "btemp_matrix[3,2]",
                                                 "btemp_matrix[3,4]",
                                                 "btemp_matrix[3,5]",
                                                 "btemp_matrix[3,9]",
                                                 "btemp_matrix[4,3]",
                                                 "btemp_matrix[5,3]",
                                                 "btemp_matrix[5,8]"
  ),
  yint = c(0,0,0.5,0,0,-0.5,0.3,1,0,0,0.2,0))
  
  actual_values_flow <- data.frame(parameter = c("bflow_matrix[1,2]",
                                                 "bflow_matrix[2,1]",
                                                 "bflow_matrix[2,3]",
                                                 "bflow_matrix[2,6]",
                                                 "bflow_matrix[2,7]",
                                                 "bflow_matrix[3,2]",
                                                 "bflow_matrix[3,4]",
                                                 "bflow_matrix[3,5]",
                                                 "bflow_matrix[3,9]",
                                                 "bflow_matrix[4,3]",
                                                 "bflow_matrix[5,3]",
                                                 "bflow_matrix[5,8]"
  ),
  yint = c(0,0.5,0,0,0,0.05,0,-0.3,0,0,0,0))
  
  btemp_plot +
    geom_hline(data = actual_values_temp, aes(yintercept = yint), lty = 2) -> btemp_plot
  
  bflow_plot +
    geom_hline(data = actual_values_flow, aes(yintercept = yint), lty = 2) -> bflow_plot
  
  return(list(b0_plot, btemp_plot, bflow_plot, brear_plot, borigin1_plot, borigin2_plot))
}





##### 600 fish, origin only ######

JAGS_600_origin_list <- readRDS( here::here("simulation", "JAGS_cov_origin_600_list.rds"))
JAGS_obj <- JAGS_600_origin_list[[1]]
mod_mcmc <- as.mcmc(JAGS_obj)
plot(mod_mcmc)


JAGS_600_origin_list <- readRDS( here::here("simulation", "JAGS_cov_origin_600_list.rds"))
sim600_origin_plots <- simulation_plots_origin(JAGS_list = JAGS_600_origin_list)
sim600_origin_plots[[1]]
sim600_origin_plots[[2]]
sim600_origin_plots[[3]]
ggsave(here::here("simulation", "figures", "sim600_origin_plots_b0.png"), height = 6, width = 10, sim600_origin_plots[[1]])
ggsave(here::here("simulation", "figures", "sim600_origin_plots_borigin1.png"), height = 6, width = 10, sim600_origin_plots[[2]])
ggsave(here::here("simulation", "figures", "sim600_origin_plots_borigin2.png"), height = 6, width = 10, sim600_origin_plots[[3]])


##### 1200 fish, origin only ######

origin_1200_1 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_1200", "JAGS_cov_origin_1200_list_1.rds"))
origin_1200_2 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_1200", "JAGS_cov_origin_1200_list_2.rds"))
origin_1200_3 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_1200", "JAGS_cov_origin_1200_list_3.rds"))
origin_1200_4 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_1200", "JAGS_cov_origin_1200_list_4.rds"))
origin_1200_5 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_1200", "JAGS_cov_origin_1200_list_5.rds"))

# Combine these into a single list
JAGS_1200_origin_list <- c(list(origin_1200_1[[1]]), list(origin_1200_1[[2]]), list(origin_1200_2[[3]]), list(origin_1200_2[[4]]),
                    list(origin_1200_3[[5]]), list(origin_1200_3[[6]]), list(origin_1200_4[[7]]), list(origin_1200_4[[8]]),
                    list(origin_1200_5[[9]]), list(origin_1200_5[[10]]))

# JAGS_1200_origin_list <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_1200", "JAGS_cov_origin_1200_list.rds"))
# JAGS_obj <- JAGS_1200_origin_list[[3]]
mod_mcmc <- as.mcmc(JAGS_obj)
# summary(mod_mcmc)
# plot(mod_mcmc)

sim1200_origin_plots <- simulation_plots_origin(JAGS_list = JAGS_1200_origin_list)
# sim1200_origin_plots[[1]]
# sim1200_origin_plots[[2]]
# sim1200_origin_plots[[3]]
ggsave(here::here("simulation", "figures", "origin_only", "sim1200_origin_plots_b0.png"), height = 6, width = 10, sim1200_origin_plots[[1]])
ggsave(here::here("simulation", "figures", "origin_only", "sim1200_origin_plots_borigin1.png"), height = 6, width = 10, sim1200_origin_plots[[2]])
ggsave(here::here("simulation", "figures", "origin_only", "sim1200_origin_plots_borigin2.png"), height = 6, width = 10, sim1200_origin_plots[[3]])


##### 1200 fish, temperature and flow ######

JAGS_1200_cov_continuous_list <- readRDS( here::here("simulation", "JAGS_cov_continuous_1200_list.rds"))
JAGS_obj <- JAGS_1200_cov_continuous_list[[1]]
JAGS_obj <- readRDS(here::here("simulation", "JAGS_cov_continuous_1200_onerun.rds"))
mod_mcmc <- as.mcmc(JAGS_obj)
plot(mod_mcmc)

sim1200_origin_plots <- simulation_plots_origin(JAGS_list = JAGS_1200_origin_list)
sim1200_origin_plots[[1]]
sim1200_origin_plots[[2]]
sim1200_origin_plots[[3]]
ggsave(here::here("simulation", "figures", "sim1200_temp_plots_b0.png"), height = 6, width = 10, sim1200_origin_plots[[1]])
ggsave(here::here("simulation", "figures", "sim1200_temp_plots_borigin1.png"), height = 6, width = 10, sim1200_origin_plots[[2]])
ggsave(here::here("simulation", "figures", "sim1200_temp_plots_borigin2.png"), height = 6, width = 10, sim1200_origin_plots[[3]])





##### 1200 fish, origin + rear #####
OR_1200_1 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_rear_1200", "JAGS_origin_rear_1200_list_1.rds"))
OR_1200_2 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_rear_1200", "JAGS_origin_rear_1200_list_2.rds"))
OR_1200_3 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_rear_1200", "JAGS_origin_rear_1200_list_3.rds"))
OR_1200_4 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_rear_1200", "JAGS_origin_rear_1200_list_4.rds"))
OR_1200_5 <- readRDS(here::here("from_hyak_transfer", "2022-05-17_origin_rear_1200", "JAGS_origin_rear_1200_list_5.rds"))

# Combine these into a single list
OR_1200 <- c(list(OR_1200_1[[1]]), list(OR_1200_1[[2]]), list(OR_1200_2[[3]]), list(OR_1200_2[[4]]),
             list(OR_1200_3[[5]]), list(OR_1200_3[[6]]), list(OR_1200_4[[7]]), list(OR_1200_4[[8]]),
             list(OR_1200_5[[9]]), list(OR_1200_5[[10]]))


# Use plotting function
sim1200_origin_rear_plots <- simulation_plots_origin_rear(JAGS_list = OR_1200, params_in_model = params_in_origin_rear_model)
# sim600_origin_rear_plots[[1]]
# sim600_origin_rear_plots[[2]]
ggsave(here::here("simulation", "figures", "origin_rear", "sim1200_origin_rear_plots_b0.png"), height = 6, width = 10, sim1200_origin_rear_plots[[1]])
ggsave(here::here("simulation", "figures", "origin_rear", "sim1200_origin_rear_plots_brear.png"), height = 6, width = 10, sim1200_origin_rear_plots[[2]])
ggsave(here::here("simulation", "figures", "origin_rear", "sim1200_origin_rear_plots_borigin1.png"), height = 6, width = 10, sim1200_origin_rear_plots[[3]])
ggsave(here::here("simulation", "figures", "origin_rear", "sim1200_origin_rear_plots_borigin2.png"), height = 6, width = 10, sim1200_origin_rear_plots[[4]])

# Okay, so run 4 was terrible - let's look at the traces
JAGS_obj <- OR_1200[[4]]
mod_mcmc <- as.mcmc(JAGS_obj)
# plot(mod_mcmc)










##### 1200 fish, temperature only #####

params_in_temp_model <- c("b0_matrix[2,1]",
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
                     "b0_matrix[3,9]",
                     "btemp_matrix[2,1]",
                     "btemp_matrix[1,2]",
                     "btemp_matrix[3,2]",
                     # "btemp_matrix[6,2]",
                     # "btemp_matrix[7,2]",
                     "btemp_matrix[2,3]",
                     "btemp_matrix[4,3]",
                     "btemp_matrix[5,3]",
                     # "btemp_matrix[9,3]",
                     "btemp_matrix[3,4]",
                     "btemp_matrix[3,5]",
                     # "btemp_matrix[8,5]",
                     "btemp_matrix[2,6]",
                     "btemp_matrix[2,7]",
                     "btemp_matrix[5,8]",
                     "btemp_matrix[3,9]")

temp_1200_1 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_1200", "JAGS_temp_1200_list_1.rds"))
temp_1200_2 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_1200", "JAGS_temp_1200_list_2.rds"))
temp_1200_3 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_1200", "JAGS_temp_1200_list_3.rds"))
temp_1200_4 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_1200", "JAGS_temp_1200_list_4.rds"))
temp_1200_5 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_1200", "JAGS_temp_1200_list_5.rds"))

# Combine these into a single list
temp_1200 <- c(list(temp_1200_1[[1]]), list(temp_1200_1[[2]]), list(temp_1200_2[[3]]), list(temp_1200_2[[4]]),
             list(temp_1200_3[[5]]), list(temp_1200_3[[6]]), list(temp_1200_4[[7]]), list(temp_1200_4[[8]]),
             list(temp_1200_5[[9]]), list(temp_1200_5[[10]]))


# Use plotting function
sim1200_temp_plots <- simulation_plots_temp(JAGS_list = temp_1200, params_in_model = params_in_temp_model)
# sim600_temp_plots[[1]]
# sim600_temp_plots[[2]]
ggsave(here::here("simulation", "figures", "temperature_only", "sim1200_temp_plots_b0.png"), height = 6, width = 10, sim1200_temp_plots[[1]])
ggsave(here::here("simulation", "figures", "temperature_only", "sim1200_temp_plots_btemp.png"), height = 6, width = 10, sim1200_temp_plots[[2]])


##### 1200 fish, temperature + flow only #####


temp_flow_1200_1 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_flow_1200", "JAGS_temp_flow_1200_list_1.rds"))
temp_flow_1200_2 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_flow_1200", "JAGS_temp_flow_1200_list_2.rds"))
temp_flow_1200_3 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_flow_1200", "JAGS_temp_flow_1200_list_3.rds"))
temp_flow_1200_4 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_flow_1200", "JAGS_temp_flow_1200_list_4.rds"))
temp_flow_1200_5 <- readRDS(here::here("from_hyak_transfer", "2022-05-18_temp_flow_1200", "JAGS_temp_flow_1200_list_5.rds"))

# Combine these into a single list
temp_flow_1200 <- c(list(temp_flow_1200_1[[1]]), list(temp_flow_1200_1[[2]]), list(temp_flow_1200_2[[3]]), list(temp_flow_1200_2[[4]]),
               list(temp_flow_1200_3[[5]]), list(temp_flow_1200_3[[6]]), list(temp_flow_1200_4[[7]]), list(temp_flow_1200_4[[8]]),
               list(temp_flow_1200_5[[9]]), list(temp_flow_1200_5[[10]]))

# Inspect chains for one run
JAGS_obj <- temp_flow_1200[[1]]
mod_mcmc <- as.mcmc(JAGS_obj)
# plot(mod_mcmc)

# Use plotting function
sim1200_temp_flow_plots <- simulation_plots_temp_flow(JAGS_list = temp_flow_1200, params_in_model = params_in_temp_flow_model)
# sim600_temp_flow_plots[[1]]
# sim600_temp_flow_plots[[2]]
ggsave(here::here("simulation", "figures", "temperature_flow", "sim1200_temp_flow_plots_b0.png"), height = 6, width = 10, sim1200_temp_flow_plots[[1]])
ggsave(here::here("simulation", "figures", "temperature_flow", "sim1200_temp_flow_plots_btemp.png"), height = 6, width = 10, sim1200_temp_flow_plots[[2]])
ggsave(here::here("simulation", "figures", "temperature_flow", "sim1200_temp_flow_plots_bflow.png"), height = 6, width = 10, sim1200_temp_flow_plots[[3]])







##### 1200 fish, all 4 covariates #####

all4_1200_1 <- readRDS(here::here("from_hyak_transfer", "2022-05-19_all4_1200", "JAGS_all4_1200_list_1.rds"))
all4_1200_2 <- readRDS(here::here("from_hyak_transfer", "2022-05-19_all4_1200", "JAGS_all4_1200_list_2.rds"))
all4_1200_3 <- readRDS(here::here("from_hyak_transfer", "2022-05-19_all4_1200", "JAGS_all4_1200_list_3.rds"))
all4_1200_4 <- readRDS(here::here("from_hyak_transfer", "2022-05-19_all4_1200", "JAGS_all4_1200_list_4.rds"))
all4_1200_5 <- readRDS(here::here("from_hyak_transfer", "2022-05-19_all4_1200", "JAGS_all4_1200_list_5.rds"))

# Combine these into a single list
all4_1200 <- c(list(all4_1200_1[[1]]), list(all4_1200_1[[2]]), list(all4_1200_2[[3]]), list(all4_1200_2[[4]]),
                    list(all4_1200_3[[5]]), list(all4_1200_3[[6]]), list(all4_1200_4[[7]]), list(all4_1200_4[[8]]),
                    list(all4_1200_5[[9]]), list(all4_1200_5[[10]]))

# Inspect chains for one run
JAGS_obj <- all4_1200[[1]]
mod_mcmc <- as.mcmc(JAGS_obj)
# plot(mod_mcmc)

# Use plotting function
sim1200_all4_plots <- simulation_plots_all4(JAGS_list = all4_1200, params_in_model = params_in_all4_model)
# sim600_temp_flow_plots[[1]]
# sim600_temp_flow_plots[[2]]
ggsave(here::here("simulation", "figures", "all4_covariates", "sim1200_all4_plots_b0.png"), height = 6, width = 10, sim1200_all4_plots[[1]])
ggsave(here::here("simulation", "figures", "all4_covariates", "sim1200_all4_plots_btemp.png"), height = 6, width = 10, sim1200_all4_plots[[2]])
ggsave(here::here("simulation", "figures", "all4_covariates", "sim1200_all4_plots_bflow.png"), height = 6, width = 10, sim1200_all4_plots[[3]])
ggsave(here::here("simulation", "figures", "all4_covariates", "sim1200_all4_plots_brear.png"), height = 6, width = 10, sim1200_all4_plots[[4]])
ggsave(here::here("simulation", "figures", "all4_covariates", "sim1200_all4_plots_borigin1.png"), height = 6, width = 10, sim1200_all4_plots[[5]])
ggsave(here::here("simulation", "figures", "all4_covariates", "sim1200_all4_plots_borigin2.png"), height = 6, width = 10, sim1200_all4_plots[[6]])




