# 14.1-JDR-ff-cov.R

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

### John Day River ###

# John Day River Steelhead:
# compare homing based on changing spill during April 3 - June 30 period at BON
# 50, 100, 150 kcfs
BON_spillvolume_values <- c(0.5, 1, 1.5)

fix_run_years <- rep_years$fix_run_year[2:4]

JDR_BON_spillvolume_homing <- data.frame()

for (i in 1:length(BON_spillvolume_values)){
  for (j in 1:length(fix_run_years)){
    JDR_FF <- compare_final_fate_fixcov_times_rear_type_MC(niter = ff_iter, nsim = ff_nsim,
                                                                       origin_select = "John Day River",
                                                                       fix_state = NA, fix_temp_season0_value = NA, fix_temp_season1_value = NA,
                                                                       fix_spillwindow_season0_value = NA, fix_spillwindow_season1_value = NA,
                                                                       fix_spillwindow_states = c(F, T, rep(F, 7)),
                                                                       fix_spillwindow_values = c(NA, BON_spillvolume_values[i], rep(NA, 7)),
                                                                       fix_spillwindow_start_days = c(0, yday("2024-04-03"), rep(0, 7)),
                                                                       fix_spillwindow_end_days = c(0, yday("2024-06-30"), rep(0, 7)),
                                                                       fix_winterspill_value = NA,
                                                                       fix_run_year = fix_run_years[j])
    JDR_FF$BON_spillvolume <-  BON_spillvolume_values[i]
    JDR_FF$fix_run_year <-  fix_run_years[j]
    
    JDR_homing <- subset(JDR_FF, state == "John Day River")
    
    JDR_BON_spillvolume_homing %>% 
      bind_rows(., JDR_homing) -> JDR_BON_spillvolume_homing
    
  }
}

JDR_BON_spillvolume_homing %>% 
  mutate(BON_spillvolume_actual = BON_spillvolume*100) -> JDR_BON_spillvolume_homing

JDR_BON_spillvolume_homing %>% 
  left_join(rep_years, by = "fix_run_year") -> JDR_BON_spillvolume_homing

# save data
save(JDR_BON_spillvolume_homing, file = here::here("stan_actual", "output", "final_fates_covariates", "JDR_BON_spillvolume_homing.rda"))

rear_colors <- c(hatchery = "#ff7f00", wild = "#33a02c")
JDR_BON_spillvolume_homing_plot <- ggplot(JDR_BON_spillvolume_homing, aes(x = BON_spillvolume_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`, color = rear_type,
                                                                            shape = conditions)) +
  geom_point(size = 3.5, position=position_dodge(width=3)) +
  geom_linerange(position=position_dodge(width=3)) +
  scale_shape_manual(values = c(15,16,17)) +
  ylab("Homing Probability") +
  xlab("Spill volume at Bonneville Dam from April 3 - June 30 (kcfs)") +
  # ggtitle(" ") +
  # Create a scale common to all
  scale_y_continuous(lim = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = rear_colors) +
  theme(plot.title = element_text(size = 12),
        # axis.text.y = element_text(color = rev(state_significance_colors)),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12)) +
  ggtitle("Homing by John Day River Steelhead under different spill conditions at Bonneville Dam during April 3 - June 30")

ggsave(here::here("stan_actual", "output", "final_fates_covariates", "JDR_BON_spillvolume_homing_plot_temps.png"), JDR_BON_spillvolume_homing_plot, height = 8, width = 8)