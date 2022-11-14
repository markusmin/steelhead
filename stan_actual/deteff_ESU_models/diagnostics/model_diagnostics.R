# stan DE model diagnostics

library(cmdstanr)
library(posterior)
library(tidyverse)
library(lubridate)
library(bayesplot)
library(here)

# Load the model run
snake_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "diagnostics", "snake", "50iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))

snake_summary <- snake_fit$summary()
# Looking at this summary, it never moved. That makes sense since they were all divergent transitions.
# From stan reference manual:
# "The positions along the simulated trajectory after the Hamiltonian diverges will never be selected as the next draw of the MCMC algorithm, 
# potentially reducing Hamiltonian Monte Carlo to a simple random walk and biasing estimates by not being able to thoroughly explore the posterior distribution.

# plot some traces
# acf(fit$draws()[,2,3])
mcmc_trace(fit$draws(), pars = c("b0_matrix_1_2"))
