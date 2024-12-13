# 06-model-diagnostics

# This script takes the output from the stan model runs in 05-stan-runs and
# runs diagnostic checks

# Here, we are analyzing all of the output from the v2 runs, where we made
# the minor changes described on the github.io website

# Note that the package installation suddenly became much more complicated after that July 9 2024 Hyak update, not sure why

# If packages have been deleted from gscratch/scrubbed, run the following lines:
# cmdstanr and rstan both apparently require some additional dependencies, try this:
install.packages("checkmate")
install.packages("data.table")
install.packages("jsonlite")
install.packages("processx")
install.packages("R6")
install.packages("withr")
install.packages("rlang")
install.packages("backports")
install.packages("abind")
install.packages("distributional")
install.packages("vctrs")
install.packages("generics")
install.packages("tidyverse")
install.packages("lifecycle")
install.packages("tibble")
install.packages("magrittr")
install.packages("ggplot2", dependcies = TRUE)
install.packages("gtable")
install.packages("dplyr")
install.packages("tidyselect")
install.packages("rprojroot")
install.packages("scales")
install.packages("munsell")
install.packages("colorspace")
install.packages("ggsignif")
install.packages("StanHeaders")
install.packages("inline")
install.packages("gridExtra")
install.packages("Rcpp")
install.packages("RcppParallel")
install.packages("loo")
install.packages("pkgbuild")
install.packages("QuickJSR")
install.packages("RcppEigen")
install.packages("BH")
install.packages("purrr")
install.packages("rstatix")
install.packages("broom")
install.packages("tidyr")
install.packages("car")
install.packages("carData")
install.packages("stringr")
install.packages("stringi")
install.packages("forcats")
install.packages("cowplot")
install.packages("labeling")
install.packages("farver")
install.packages("reshape2")
install.packages("plyr")
install.packages("ggthemes")
install.packages("timechange")
install.packages("lubridate")
# # install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")), dependencies = TRUE)
# # alternatively: remotes::install_github("stan-dev/cmdstanr")
# alternatively: renv::install("stan-dev/cmdstanr") # this ended up working
install.packages("renv")
renv::install("stan-dev/cmdstanr")
install.packages("rstan", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
install.packages("bayesplot")
install.packages("here")
install.packages("ggpubr")
install.packages("Formula")
install.packages("htmltools")
install.packages("digest")
install.packages("fastmap")
install.packages("miniUI")
install.packages("httpuv")
install.packages("shiny")
install.packages("later")
install.packages("xtable")
install.packages("promises")
install.packages("ggExtra")
install.packages("mime")
install.packages("patchwork")


# First, need to load in all of the model runs and all of the packages.
source("analysis/analysis/00-load-model-runs.R")

##### Run diagnostic summaries #####

# A first summary: how many parameters are we estimating?
# we can drop duplicated parameters that have the exact same estimates, because
# we often have parameters stored in multiple places
# we should alo plot these parameters
UCW_fit_summary %>% 
  filter(!(duplicated(median) & duplicated(mean))) %>% 
  filter(!(variable == "lp__")) -> UCW_unique_parameters_fit_summary

UCH_fit_summary %>% 
  filter(!(duplicated(median) & duplicated(mean))) %>% 
  filter(!(variable == "lp__")) -> UCH_unique_parameters_fit_summary

MCW_fit_summary %>% 
  filter(!(duplicated(median) & duplicated(mean))) %>% 
  filter(!(variable == "lp__")) -> MCW_unique_parameters_fit_summary

MCH_fit_summary %>% 
  filter(!(duplicated(median) & duplicated(mean))) %>% 
  filter(!(variable == "lp__")) -> MCH_unique_parameters_fit_summary

SRW_fit_summary %>% 
  filter(!(duplicated(median) & duplicated(mean))) %>% 
  filter(!(variable == "lp__")) -> SRW_unique_parameters_fit_summary

SRH_fit_summary %>% 
  filter(!(duplicated(median) & duplicated(mean))) %>% 
  filter(!(variable == "lp__")) -> SRH_unique_parameters_fit_summary

dim(UCW_unique_parameters_fit_summary)
dim(UCH_unique_parameters_fit_summary)
dim(MCW_unique_parameters_fit_summary)
dim(MCH_unique_parameters_fit_summary)
dim(SRW_unique_parameters_fit_summary)
dim(SRH_unique_parameters_fit_summary)



# Create six-panel diagnostic plots for rhat, ess_bulk, and ess_tail

create_rhat_hist <- function(fit_summary, population){
  rhat_plot <- ggplot(fit_summary, aes(x = rhat)) +
    geom_histogram() + 
    # annotate("text", label = "Range:", x = max(fit_summary$rhat, na.rm = T), y = 600, hjust = 1) +
    # annotate("text", label = round(range(fit_summary$rhat, na.rm = T),3)[1], x = max(fit_summary$rhat, na.rm = T), y = 550,hjust = 1) +
    # annotate("text", label = paste0(" - ", round(range(fit_summary$rhat, na.rm = T),3)[2]), x = max(fit_summary$rhat, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population)
  
  return(rhat_plot)
}

create_ess_bulk_hist <- function(fit_summary, population){
  ess_bulk_plot <- ggplot(fit_summary, aes(x = ess_bulk)) +
    geom_histogram() + 
    # annotate("text", label = "Range:", x = max(fit_summary$ess_bulk, na.rm = T), y = 600, hjust = 1) +
    # annotate("text", label = round(range(fit_summary$ess_bulk, na.rm = T),3)[1], x = max(fit_summary$ess_bulk, na.rm = T), y = 550,hjust = 1) +
    # annotate("text", label = paste0(" - ", round(range(fit_summary$ess_bulk, na.rm = T),3)[2]), x = max(fit_summary$ess_bulk, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population) +
    scale_x_continuous(lim = c(0, 4200))
  
  return(ess_bulk_plot)
}

create_ess_tail_hist <- function(fit_summary, population){
  ess_tail_plot <- ggplot(fit_summary, aes(x = ess_tail)) +
    geom_histogram() + 
    # annotate("text", label = "Range:", x = max(fit_summary$ess_tail, na.rm = T), y = 600, hjust = 1) +
    # annotate("text", label = round(range(fit_summary$ess_tail, na.rm = T),3)[1], x = max(fit_summary$ess_tail, na.rm = T), y = 550,hjust = 1) +
    # annotate("text", label = paste0(" - ", round(range(fit_summary$ess_tail, na.rm = T),3)[2]), x = max(fit_summary$ess_tail, na.rm = T), y = 500, hjust = 1) +
    ylab("N parameters") +
    ggtitle(population) +
    scale_x_continuous(lim = c(0, 4200))
  
  return(ess_tail_plot)
}

# Run functions for all populations
UCW_rhat_plot <- create_rhat_hist(fit_summary = UCW_unique_parameters_fit_summary, population = "UCW")
UCH_rhat_plot <- create_rhat_hist(fit_summary = UCH_unique_parameters_fit_summary, population = "UCH")
MCW_rhat_plot <- create_rhat_hist(fit_summary = MCW_unique_parameters_fit_summary, population = "MCW")
MCH_rhat_plot <- create_rhat_hist(fit_summary = MCH_unique_parameters_fit_summary, population = "MCH")
SRW_rhat_plot <- create_rhat_hist(fit_summary = SRW_unique_parameters_fit_summary, population = "SRW")
SRH_rhat_plot <- create_rhat_hist(fit_summary = SRH_unique_parameters_fit_summary, population = "SRH")

UCW_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = UCW_unique_parameters_fit_summary, population = "UCW")
UCH_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = UCH_unique_parameters_fit_summary, population = "UCH")
MCW_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = MCW_unique_parameters_fit_summary, population = "MCW")
MCH_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = MCH_unique_parameters_fit_summary, population = "MCH")
SRW_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = SRW_unique_parameters_fit_summary, population = "SRW")
SRH_ess_bulk_plot <- create_ess_bulk_hist(fit_summary = SRH_unique_parameters_fit_summary, population = "SRH")

UCW_ess_tail_plot <- create_ess_tail_hist(fit_summary = UCW_unique_parameters_fit_summary, population = "UCW")
UCH_ess_tail_plot <- create_ess_tail_hist(fit_summary = UCH_unique_parameters_fit_summary, population = "UCH")
MCW_ess_tail_plot <- create_ess_tail_hist(fit_summary = MCW_unique_parameters_fit_summary, population = "MCW")
MCH_ess_tail_plot <- create_ess_tail_hist(fit_summary = MCH_unique_parameters_fit_summary, population = "MCH")
SRW_ess_tail_plot <- create_ess_tail_hist(fit_summary = SRW_unique_parameters_fit_summary, population = "SRW")
SRH_ess_tail_plot <- create_ess_tail_hist(fit_summary = SRH_unique_parameters_fit_summary, population = "SRH")

# arrange them
rhat_comb_plots <- ggarrange(plotlist = list(UCW_rhat_plot, UCH_rhat_plot, MCW_rhat_plot, MCH_rhat_plot, SRW_rhat_plot, SRH_rhat_plot),
                             ncol = 3, nrow = 2)

ess_bulk_comb_plots <- ggarrange(plotlist = list(UCW_ess_bulk_plot, UCH_ess_bulk_plot, MCW_ess_bulk_plot, MCH_ess_bulk_plot, SRW_ess_bulk_plot, SRH_ess_bulk_plot),
                                 ncol = 3, nrow = 2)

ess_tail_comb_plots <- ggarrange(plotlist = list(UCW_ess_tail_plot, UCH_ess_tail_plot, MCW_ess_tail_plot, MCH_ess_tail_plot, SRW_ess_tail_plot, SRH_ess_tail_plot),
                                 ncol = 3, nrow = 2)

ggsave(file = "stan_actual/output/diagnostics/rhat_comb_plots.png", plot = rhat_comb_plots,
       height = 10, width = 10)
ggsave(file = "stan_actual/output/diagnostics/ess_tail_comb_plots.png", plot = ess_tail_comb_plots,
       height = 10, width = 10)
ggsave(file = "stan_actual/output/diagnostics/ess_bulk_comb_plots.png", plot = ess_bulk_comb_plots,
       height = 10, width = 10)

# Ok, some SRH and SRW parameters look a bit iffy. 
arrange(SRH_fit_summary, ess_tail)[1:100,] -> bad_ess_tail_SRH
bad_ess_tail_SRH$variable
# They're mostly yearxorigin parameters, but not exclusively...

# Draw some traceplots
subset_draws(SRH_fit, variable = c("btemp1xorigin4_matrix_9_39_DE")) -> SRH_fit_draws_to_plot
mcmc_trace(SRH_fit_draws_to_plot, facet_args = list(ncol = 1, strip.position = "left")) -> SRH_traceplots
# chain 1 jsut gets stuck in a weird parameter space before returning to reality




#### Diagnostics: Traceplots for selected parameters ####
subset_draws(UCW_fit, variable = c("b0_matrix_2_3","bspillwindow_matrix_2_1",
                                               "btemp0_matrix_2_3")) -> UCW_fit_draws_to_plot
mcmc_trace(UCW_fit_draws_to_plot, facet_args = list(ncol = 1, strip.position = "left")) -> UCW_traceplots

subset_draws(UCH_fit, variable = c("b0_matrix_2_3","bspillwindow_matrix_2_1",
                                   "btemp0_matrix_2_3")) -> UCH_fit_draws_to_plot
mcmc_trace(UCH_fit_draws_to_plot, facet_args = list(ncol = 1, strip.position = "left")) -> UCH_traceplots

subset_draws(MCW_fit, variable = c("b0_matrix_4_5","bspillwindow_matrix_2_1",
                                   "btemp1xorigin1_matrix_2_3")) -> MCW_fit_draws_to_plot
mcmc_trace(MCW_fit_draws_to_plot, facet_args = list(ncol = 1, strip.position = "left")) -> MCW_traceplots

subset_draws(MCH_fit, variable = c("b0_matrix_4_5","bspillwindow_matrix_2_1",
                                   "btemp1xorigin1_matrix_2_3")) -> MCH_fit_draws_to_plot
mcmc_trace(MCH_fit_draws_to_plot, facet_args = list(ncol = 1, strip.position = "left")) -> MCH_traceplots

subset_draws(SRW_fit, variable = c("b0_matrix_2_3","bspillwindow_matrix_2_1",
                                   "btemp0_matrix_2_3")) -> SRW_fit_draws_to_plot
mcmc_trace(SRW_fit_draws_to_plot, facet_args = list(ncol = 1, strip.position = "left")) -> SRW_traceplots

subset_draws(SRH_fit, variable = c("b0_matrix_2_3","bspillwindow_matrix_2_1",
                                   "btemp0_matrix_2_3")) -> SRH_fit_draws_to_plot
mcmc_trace(SRH_fit_draws_to_plot, facet_args = list(ncol = 1, strip.position = "left")) -> SRH_traceplots

ggsave(file = "stan_actual/output/traceplots/UCW_traceplots.png", plot = UCW_traceplots,
       height = 8, width = 8)
ggsave(file = "stan_actual/output/traceplots/UCH_traceplots.png", plot = UCH_traceplots,
       height = 8, width = 8)
ggsave(file = "stan_actual/output/traceplots/MCW_traceplots.png", plot = MCW_traceplots,
       height = 8, width = 8)
ggsave(file = "stan_actual/output/traceplots/MCH_traceplots.png", plot = MCH_traceplots,
       height = 8, width = 8)
ggsave(file = "stan_actual/output/traceplots/SRW_traceplots.png", plot = SRW_traceplots,
       height = 8, width = 8)
ggsave(file = "stan_actual/output/traceplots/SRH_traceplots.png", plot = SRH_traceplots,
       height = 8, width = 8)


#### Run diagnostics: Pairs plots ####

# For one specific movement
one_movement_param_pairs <- function(fit, fit_summary, movement, all_movements){
  parameters <- unique(fit_summary$variable)
  
  if (length(grep(paste0(movement), parameters)) > 1) {
    one_movement_params <- parameters[grep(paste0(movement), parameters)]
    # remove NDE parameters
    one_movement_params <- one_movement_params[!(grepl("NDE", one_movement_params))]
    # remove year parameters
    # one_movement_params <- one_movement_params[!(grepl("\\[", one_movement_params))]
    one_movement_params <- one_movement_params[!(grepl("year", one_movement_params))]
    # remove any parameters from movements that are not supposed to be in there
    
    # if/else statement ensures that if we're going to be dropping all parameters, that's a mistake so don't do it
    if (sum(grepl(paste(all_movements[!(all_movements == movement)], collapse = "|"), one_movement_params)) == length(one_movement_params)) {
      one_movement_params <- one_movement_params
    } else {
      one_movement_params <- one_movement_params[!(grepl(paste(all_movements[!(all_movements == movement)], collapse = "|"), one_movement_params))]
    }
    
    # generate the plot
    pairs_plot <- mcmc_pairs(fit, pars = one_movement_params,
                             off_diag_args = list(size = 0.5))
  } else {
    print(paste0("One or no parameters for movement ", movement))
  }
  
  return(pairs_plot)
}


# run for all states, all populations
populations <- c("UCW", "UCH", "MCW", "MCH", "SRW", "SRH")

## generate all plots
UCW_movements_plotlist <- vector(mode = "list", length = length(UCW_movements))
UCH_movements_plotlist <- vector(mode = "list", length = length(UCH_movements))
MCW_movements_plotlist <- vector(mode = "list", length = length(MCW_movements))
MCH_movements_plotlist <- vector(mode = "list", length = length(MCH_movements))
SRW_movements_plotlist <- vector(mode = "list", length = length(SRW_movements))
SRH_movements_plotlist <- vector(mode = "list", length = length(SRH_movements))


for (i in 1:length(populations)){
  for (j in 1:length(eval(parse(text = paste0(populations[i], "_movements"))))){
    pairs_plot <- one_movement_param_pairs(fit = eval(parse(text = paste0(populations[i], "_fit"))),
                                           fit_summary = eval(parse(text = paste0(populations[i], "_fit_summary"))),
                                           movement = eval(parse(text = paste0(populations[i], "_movements[j]"))),
                                           all_movements = eval(parse(text = paste0(populations[i], "_movements"))))
    
    eval(parse(text = paste0(populations[i], "_movements_plotlist[[j]] <- pairs_plot")))
    print(paste0("i = ", i))
    print(paste0("j = ", j))
  }
}


## save all plots
for (i in 1:length(populations)){
  for (j in 1:length(eval(parse(text = paste0(populations[i], "_movements_plotlist"))))){
    ggsave(file = paste0("stan_actual/output/pairs_plots/", populations[i], "/", eval(parse(text = paste0(populations[i], "_movements[j]"))), 
                         "_pairs_plot.png"), plot = eval(parse(text = paste0(populations[i], "_movements_plotlist[[j]]"))),
           height = 15, width = 15)
  }
}


## Pairs plot, for all intercept parameters movements out of one state
one_state_intercept_param_pairs <- function(fit, fit_summary, state){
  parameters <- unique(fit_summary$variable)
  
    one_state_params <- parameters[grep(paste0("b0_matrix_", state, "_"), parameters)]
    # remove NDE parameters
    one_state_params <- one_state_params[!(grepl("NDE", one_state_params))]
    # remove year parameters
    # one_state_params <- one_state_params[!(grepl("\\[", one_state_params))]
    one_state_params <- one_state_params[!(grepl("year", one_state_params))]
    
    if (length(one_state_params) > 1) {
      # generate the plot
      pairs_plot <- mcmc_pairs(fit, pars = one_state_params,
                               off_diag_args = list(size = 0.5))
      
    } else {
      print(paste0("only one parameter for state ", state))
    }
    
  
  return(pairs_plot)
}


# run for all states, all populations
populations <- c("UCW", "UCH", "MCW", "MCH", "SRW", "SRH")
states <- 2:9

## generate all plots
UCW_states_plotlist <- vector(mode = "list", length = length(states))
UCH_states_plotlist <- vector(mode = "list", length = length(states))
MCW_states_plotlist <- vector(mode = "list", length = length(states))
MCH_states_plotlist <- vector(mode = "list", length = length(states))
SRW_states_plotlist <- vector(mode = "list", length = length(states))
SRH_states_plotlist <- vector(mode = "list", length = length(states))


for (i in 1:length(populations)){
  for (j in 1:length(states)){
    pairs_plot <- one_state_intercept_param_pairs(fit = eval(parse(text = paste0(populations[i], "_fit"))),
                                           fit_summary = eval(parse(text = paste0(populations[i], "_fit_summary"))),
                                           state = states[j])
    
    eval(parse(text = paste0(populations[i], "_states_plotlist[[j]] <- pairs_plot")))
    print(paste0("i = ", i))
    print(paste0("j = ", j))
  }
}


## save all plots
for (i in 1:length(populations)){
  for (j in 1:length(states)){
    ggsave(file = paste0("stan_actual/output/pairs_plots_intercept/", populations[i], "/", populations[i],"_", states[j], 
                         "_pairs_plot.png"), plot = eval(parse(text = paste0(populations[i], "_states_plotlist[[j]]"))),
           height = 15, width = 15)
  }
}




#### Pairs plot, comparing temperature and spill parameters ####
plot_spillwindow_v_temp1_by_state_origin <- function(fit, fit_summary, from_state, origin_numeric){
  parameters <- unique(fit_summary$variable)
  
  spillwindow_param <- parameters[grep(paste0("bspillwindow_matrix_", from_state), parameters)]
  temp1_params <- parameters[grep(paste0("btemp1xorigin", origin_numeric, "_matrix_", from_state), parameters)]
  pars_to_compare <- c(spillwindow_param, temp1_params) 
  
  # remove NDE parameters
  pars_to_compare <- pars_to_compare[!(grepl("NDE", pars_to_compare))]
  
  if (length(pars_to_compare) > 1) {
    # generate the plot
    pairs_plot <- mcmc_pairs(fit, pars = pars_to_compare,
                             off_diag_args = list(size = 0.5))
    
  } else {
    print(paste0("Temp1 and spillwindow do not both affect movements out of this state:", from_state))
  }
  
  
  return(pairs_plot)
}


# run this for all mainstem states and all origins
# run for UCW
for (i in 2:9){ # for all mainstem states
  for (j in 1:max(origin_param_map[8:11,"wild"], na.rm = T)) { # for all origins
    plot <- plot_spillwindow_v_temp1_by_state_origin(fit = UCW_fit, fit_summary = UCW_fit_summary, from_state = i, origin_numeric = j)
    
    ggsave(file = paste0("stan_actual/output/pairs_plots/temp_v_spill/", "UCW_from_", i, "_origin_", j, 
                         "temp_v_spill_pairs_plot.png"), plot,
           height = 15, width = 15)
  } 
}

# run for UCH
for (i in 2:9){ # for all mainstem states
  for (j in 1:max(origin_param_map[8:11,"hatchery"], na.rm = T)) { # for all origins
    plot <- plot_spillwindow_v_temp1_by_state_origin(fit = UCH_fit, fit_summary = UCH_fit_summary, from_state = i, origin_numeric = j)
    
    ggsave(file = paste0("stan_actual/output/pairs_plots/temp_v_spill/", "UCH_from_", i, "_origin_", j, 
                         "temp_v_spill_pairs_plot.png"), plot,
           height = 15, width = 15)
  } 
}

# run for MCW
for (i in 2:9){ # for all mainstem states
  for (j in 1:max(origin_param_map[1:7,"wild"], na.rm = T)) { # for all origins
    plot <- plot_spillwindow_v_temp1_by_state_origin(fit = MCW_fit, fit_summary = MCW_fit_summary, from_state = i, origin_numeric = j)
    
    ggsave(file = paste0("stan_actual/output/pairs_plots/temp_v_spill/", "MCW_from_", i, "_origin_", j, 
                         "temp_v_spill_pairs_plot.png"), plot,
           height = 15, width = 15)
  } 
}

# run for MCH
for (i in 2:9){ # for all mainstem states
  for (j in 1:max(origin_param_map[1:7,"hatchery"], na.rm = T)) { # for all origins
    plot <- plot_spillwindow_v_temp1_by_state_origin(fit = MCH_fit, fit_summary = MCH_fit_summary, from_state = i, origin_numeric = j)
    
    ggsave(file = paste0("stan_actual/output/pairs_plots/temp_v_spill/", "MCH_from_", i, "_origin_", j, 
                         "temp_v_spill_pairs_plot.png"), plot,
           height = 15, width = 15)
  } 
}

# run for SRW
for (i in 2:9){ # for all mainstem states
  for (j in 1:max(origin_param_map[12:17,"wild"], na.rm = T)) { # for all origins
    plot <- plot_spillwindow_v_temp1_by_state_origin(fit = SRW_fit, fit_summary = SRW_fit_summary, from_state = i, origin_numeric = j)
    
    ggsave(file = paste0("stan_actual/output/pairs_plots/temp_v_spill/", "SRW_from_", i, "_origin_", j, 
                         "temp_v_spill_pairs_plot.png"), plot,
           height = 15, width = 15)
  } 
}

# run for SRH
for (i in 2:9){ # for all mainstem states
  for (j in 1:max(origin_param_map[12:17,"hatchery"], na.rm = T)) { # for all origins
    plot <- plot_spillwindow_v_temp1_by_state_origin(fit = SRH_fit, fit_summary = SRH_fit_summary, from_state = i, origin_numeric = j)
    
    ggsave(file = paste0("stan_actual/output/pairs_plots/temp_v_spill/", "SRH_from_", i, "_origin_", j, 
                         "temp_v_spill_pairs_plot.png"), plot,
           height = 15, width = 15)
  } 
}




