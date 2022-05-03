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


