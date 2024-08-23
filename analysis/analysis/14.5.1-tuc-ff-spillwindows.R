# 14.5.1-tuc-ff-spillwindows.R

# This script runs final fates + covariates simulations.

#### Load libraries, state information ####
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(rlang)
library(tibble)
library(forcats)
library(lubridate)

# set the working directory so that here() will cooperate
setwd("/gscratch/scrubbed/mmin/")
library(here)

source(here::here("analysis", "analysis", "14.0-ff-cov-functions.R"))
print("sourcing scripts complete")

#### Spill volume comparison ####

### Tucannon River ###

# Here, we will run multiple spill scenarios at different projects at the same time.
# We will examine the effect of increased summer/fall spill at BON and spring spill at MCN

# Tucannon River Steelhead:
# compare homing based on changing spill during July 1 - September 15 period at BON
# 0, 50, 100, 150 kcfs
# compare homing based on changing spill during March 1 - March 31 at MCN
# 0, 50, 100, 150 kcfs
# compare homing based on changing spill during March 1 - March 31 at ICH
# 0, 25, 50, 75 kcfs
BON_spillvolume_values <- c(0, 0.5, 1, 1.5)
MCN_spillvolume_values <- c(0, 0.5, 1, 1.5)
ICH_spillvolume_values <- c(0, 0.25, 0.5, 0.75)

fix_run_years <- rep_years$fix_run_year[2:4]

TUC_BON_MCN_ICH_spillvolume_homing <- data.frame()

for (i in 1:length(BON_spillvolume_values)){
  for (j in 1:length(fix_run_years)){
    TUC_FF <- compare_final_fate_fixcov_times_rear_type_SR(niter = ff_iter, nsim = ff_nsim,
                                                                       origin_select = "Tucannon River",
                                                                       fix_state = NA, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                                       fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                                                       fix_spillwindow_states = c(F, T, T, rep(F, 4), T, F),
                                                                       fix_spillwindow_values = c(NA, BON_spillvolume_values[i], MCN_spillvolume_values[i], rep(NA, 4), ICH_spillvolume_values[i], NA),
                                                                       fix_spillwindow_start_days = c(NA, yday("2024-07-01"), yday("2024-03-01"), rep(NA, 4), yday("2024-03-01"), 0),
                                                                       fix_spillwindow_end_days = c(NA, yday("2024-09-15"), yday("2024-03-31"), rep(NA, 4), yday("2024-03-31"), 0),
                                                                       fix_winterspill_value = NA,
                                                                       fix_run_year = fix_run_years[j])
    TUC_FF$BON_spillvolume <-  BON_spillvolume_values[i]
    TUC_FF$MCN_spillvolume <-  MCN_spillvolume_values[i]
    TUC_FF$ICH_spillvolume <-  ICH_spillvolume_values[i]
    TUC_FF$fix_run_year <-  fix_run_years[j]
    
    TUC_homing <- subset(TUC_FF, state == "Tucannon River")
    
    TUC_BON_MCN_ICH_spillvolume_homing %>% 
      bind_rows(., TUC_homing) -> TUC_BON_MCN_ICH_spillvolume_homing
    
  }
}

TUC_BON_MCN_ICH_spillvolume_homing %>% 
  mutate(BON_spillvolume_actual = BON_spillvolume*100) %>%
  mutate(MCN_spillvolume_actual = MCN_spillvolume*100) %>% 
  mutate(ICH_spillvolume_actual = ICH_spillvolume*100) -> TUC_BON_MCN_ICH_spillvolume_homing

TUC_BON_MCN_ICH_spillvolume_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> TUC_BON_MCN_ICH_spillvolume_homing

# save data
save(TUC_BON_MCN_ICH_spillvolume_homing, file = here::here("stan_actual", "output", "final_fates_covariates", "TUC_BON_MCN_ICH_spillvolume_homing.rda"))

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
TUC_BON_MCN_ICH_spillvolume_homing_plot <- ggplot(TUC_BON_MCN_ICH_spillvolume_homing, aes(x = BON_spillvolume_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = rear_type,
                                                                            shape = conditions)) +
  geom_point(size = 3.5, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3)) +
  scale_shape_manual(values = c(15,16,17)) +
  ylab("Homing Probability") +
  xlab("Spill volume at Bonneville Dam from 7/1-9/15 (kcfs)") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = rear_colors) +
  theme(plot.title = elemTUC_text(size = 12),
        # axis.text.y = elemTUC_text(color = rev(state_significance_colors)),
        axis.title = elemTUC_text(size = 14),
        axis.text.x = elemTUC_text(size = 12)) +
  ggtitle("Homing by Tucannon River Steelhead under variable spill conditions at BON, MCN, and ICH")

ggsave(here::here("stan_actual", "output", "final_fates_covariates", "TUC_BON_MCN_ICH_spillvolume_homing_plot_temps.png"), TUC_BON_MCN_ICH_spillvolume_homing_plot, height = 8, width = 8)