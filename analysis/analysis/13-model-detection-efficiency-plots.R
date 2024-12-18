# 13-model-detection-efficiency-plots

# This script takes the output from the stan model runs in 05-stan-runs and
# plots how detection efficiency varies through time

# First, need to load in all of the model runs and all of the packages.
source("analysis/analysis/00-load-model-runs.R")

#### Extract all parameter values from the model fit objects ####

# Here, you need to make sure to check the origin params, because those change by DPS

# function to take a parameter type (fixed effect) and store all of them in an array
extract_DE_parameters <- function(fit, fit_summary){
  # extract b0 as an array
  parameters <- fit_summary$variable
  parameters[grepl("_alpha|_beta" , parameters)] -> param_subset
  
  # add in fifteenmile_beta and imnaha_beta
  DE_params <- c(param_subset[1:24], "fifteenmile_beta", param_subset[25], "imnaha_beta", param_subset[26:length(param_subset)])
  
  
  # arrange all parameter values together in a df
  DE_param_matrix <- matrix(data = 0, nrow = length(DE_params),
                                         length(as.matrix(fit[,,1])))
  
  rownames(DE_param_matrix) <- DE_params
  
  
  
  # extract the draws for each and store in a matrix
  for(i in 1:nrow(DE_param_matrix)){
    if (DE_params[i] %in% c("fifteenmile_beta", "imnaha_beta")){
      DE_param_matrix[i, ] <- 0
    } else {
      DE_param_matrix[i, ] <- as.matrix(fit[,,DE_params[i]])
    }
    
  }
  
  return(DE_param_matrix)
}

# function to reformat param matrix and calculate DE using those params + discharge
calculate_DE <- function(DE_param_matrix, tributary, tributary_design_matrices_array){
  
  # get the index for the tributary state (needs to be mouth since that's where we're correcting for DE)
  tributary_state <- intersect(grep(tributary, model_states, ignore.case = TRUE), grep("Mouth", model_states, ignore.case = TRUE))
  
  # use the state indexing to get the right design matrix
  trib_design_matrix <- tributary_design_matrices_array[,,tributary_state]
  
  # create a matrix to store DE per year, per iter
  niter <- 4000
  DE_matrix <- matrix(nrow = nrow(trib_design_matrix), ncol = niter)
  
  # for each run year, get a confidence interval around detection efficiency by using the different draws
  for (i in 1:nrow(trib_design_matrix)){
    # if there is no intercept term, that means there is no DE correction for that year - so skip and leave as NA
    if(sum(trib_design_matrix[i,1:21]) == 0){
      
      
    } else {
      for (j in 1:niter){
        eta <- sum(trib_design_matrix[i,] * DE_param_matrix[,j])
        DE_matrix[i,j] <- exp(eta)/(1 + exp(eta))
      }
      
    }
    
    
  }
  
  return(DE_matrix)
  
}

# function to plot DE plus credible intervals
plot_DE_rear <- function(DE_by_year_wild = NULL, DE_by_year_hatchery = NULL, plot_title){
  rear_colors <- c(hatchery = "#ff7f00", natural = "#33a02c")
  
  if(is.null(DE_by_year_hatchery)){
    niter <- 4000
    colnames(DE_by_year_wild) <- paste0("iter", 1:niter)
    
    DE_by_year_wild %>% 
      as.data.frame() %>% 
      rownames_to_column("year") %>% 
      mutate(year = as.numeric(year)) %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
      group_by(year) %>% 
      summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "natural") -> DE_by_year_wild_long
    
    DE_by_year_wild_long$year_actual <- 2004:2021
    
    DE_by_year_wild_long -> DE_by_year_rear_long

  } else if (is.null(DE_by_year_wild)){
    niter <- 4000
    colnames(DE_by_year_hatchery) <- paste0("iter", 1:niter)
    
    DE_by_year_hatchery %>% 
      as.data.frame() %>% 
      rownames_to_column("year") %>% 
      mutate(year = as.numeric(year)) %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
      group_by(year) %>% 
      summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> DE_by_year_hatchery_long
    
    DE_by_year_hatchery_long$year_actual <- 2004:2021
    
    DE_by_year_hatchery_long -> DE_by_year_rear_long
    
  } else {
    niter <- 4000
    colnames(DE_by_year_wild) <- paste0("iter", 1:niter)
    colnames(DE_by_year_hatchery) <- paste0("iter", 1:niter)
    
    DE_by_year_wild %>% 
      as.data.frame() %>% 
      rownames_to_column("year") %>% 
      mutate(year = as.numeric(year)) %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
      group_by(year) %>% 
      summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "natural") -> DE_by_year_wild_long
    
    DE_by_year_wild_long$year_actual <- 2004:2021
    
    DE_by_year_hatchery %>% 
      as.data.frame() %>% 
      rownames_to_column("year") %>% 
      mutate(year = as.numeric(year)) %>% 
      pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
      group_by(year) %>% 
      summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
      pivot_wider(names_from = q, values_from = prob) %>% 
      mutate(rear = "hatchery") -> DE_by_year_hatchery_long
    
    DE_by_year_hatchery_long$year_actual <- 2004:2021
    
    DE_by_year_wild_long %>% 
      bind_rows(., DE_by_year_hatchery_long) -> DE_by_year_rear_long
    
  }
  
  
  DE_by_year_plot <- ggplot(DE_by_year_rear_long, aes(x = year_actual, y = `0.5`, ymin = `0.025`, ymax = `0.975`,
                                                      color = rear, fill = rear)) +
    geom_line() +
    scale_color_manual(values = rear_colors) +
    scale_fill_manual(values = rear_colors) +
    geom_ribbon(alpha = 0.2, color = NA) +
    scale_y_continuous(lim = c(0,1), expand = c(0,0)) +
    scale_x_continuous(lim = c(2003,2022), expand = c(0,0)) +
    xlab("Year") +
    ylab("Detection probability") +
    ggtitle(plot_title) +
    guides(color=guide_legend(title="Rearing Type"), fill=guide_legend(title="Rearing Type")) +
    theme(panel.grid.major = element_line(color = "gray90"),
          panel.background = element_rect(fill = "white", color = NA),
          panel.border = element_rect(color = NA, fill=NA, linewidth=0.4),
          # these plot margins are to leave space for the population name on the big figure
          plot.margin = unit(c(0.2, 0.2, 0.2, 0.2),"cm"))
  
  return(DE_by_year_plot)
}



### extract all params ###
UCW_DE_param_matrix <- extract_DE_parameters(fit = UCW_fit, fit_summary = UCW_fit_summary)
UCH_DE_param_matrix <- extract_DE_parameters(fit = UCH_fit, fit_summary = UCH_fit_summary)
MCW_DE_param_matrix <- extract_DE_parameters(fit = MCW_fit, fit_summary = MCW_fit_summary)
MCH_DE_param_matrix <- extract_DE_parameters(fit = MCH_fit, fit_summary = MCH_fit_summary)
SRW_DE_param_matrix <- extract_DE_parameters(fit = SRW_fit, fit_summary = SRW_fit_summary)
SRH_DE_param_matrix <- extract_DE_parameters(fit = SRH_fit, fit_summary = SRH_fit_summary)


#### calculate and plot annual detection efficiency by tributary ####

# first some quick sanity checks to make sure that our data is formatted properly
dim(MCW_envir$data$tributary_design_matrices_array)
# run years, parameters, states
MCW_envir$data$tributary_design_matrices_array[,,34] # this should be discharge per year for Asotin Creek


### Middle Columbia

deschutes_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
deschutes_DE_plot <- plot_DE_rear(DE_by_year_wild = deschutes_DE_by_year_wild, plot_title = "Deschutes River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "deschutes_DE_plot.png"), deschutes_DE_plot, height = 8, width = 8)

john_day_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "john", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
john_day_DE_plot <- plot_DE_rear(DE_by_year_wild = john_day_DE_by_year_wild, plot_title = "John Day River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "john_day_DE_plot.png"), john_day_DE_plot, height = 8, width = 8)

fifteenmile_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "fifteenmile", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
fifteenmile_DE_plot <- plot_DE_rear(DE_by_year_wild = fifteenmile_DE_by_year_wild, plot_title = "Fifteenmile Creek Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "fifteenmile_DE_plot.png"), fifteenmile_DE_plot, height = 8, width = 8)

umatilla_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "umatilla", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
umatilla_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "umatilla", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array)
umatilla_DE_plot <- plot_DE_rear(DE_by_year_wild = umatilla_DE_by_year_wild, 
                            DE_by_year_hatchery = umatilla_DE_by_year_hatchery, 
                            plot_title = "Umatilla River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "umatilla_DE_plot.png"), umatilla_DE_plot, height = 8, width = 8)

yakima_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "yakima", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
yakima_DE_plot <- plot_DE_rear(DE_by_year_wild = yakima_DE_by_year_wild,
                          plot_title = "Yakima River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "yakima_DE_plot.png"), yakima_DE_plot, height = 8, width = 8)

walla_walla_DE_by_year_wild <- calculate_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "walla", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array)
walla_walla_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "walla", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array)
walla_walla_DE_plot <- plot_DE_rear(DE_by_year_wild = walla_walla_DE_by_year_wild,
                                    DE_by_year_hatchery = walla_walla_DE_by_year_hatchery,
                                    plot_title = "Walla Walla River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "walla_walla_DE_plot.png"), walla_walla_DE_plot, height = 8, width = 8)

### Upper Columbia
wenatchee_DE_by_year_wild <- calculate_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "wenatchee", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array)
wenatchee_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "wenatchee", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array)
wenatchee_DE_plot <- plot_DE_rear(DE_by_year_wild = wenatchee_DE_by_year_wild, 
                                  DE_by_year_hatchery = wenatchee_DE_by_year_hatchery,
                                  plot_title = "Wenatchee River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "wenatchee_DE_plot.png"), wenatchee_DE_plot, height = 8, width = 8)

entiat_DE_by_year_wild <- calculate_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "entiat", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array)
entiat_DE_plot <- plot_DE_rear(DE_by_year_wild = entiat_DE_by_year_wild, plot_title = "Entiat River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "entiat_DE_plot.png"), entiat_DE_plot, height = 8, width = 8)

okanogan_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "okanogan", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array)
okanogan_DE_plot <- plot_DE_rear(DE_by_year_hatchery = okanogan_DE_by_year_hatchery, plot_title = "Okanogan River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "okanogan_DE_plot.png"), okanogan_DE_plot, height = 8, width = 8)

methow_DE_by_year_wild <- calculate_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "methow", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array)
methow_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "methow", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array)
methow_DE_plot <- plot_DE_rear(DE_by_year_wild = methow_DE_by_year_wild,
                               DE_by_year_hatchery = methow_DE_by_year_hatchery,
                               plot_title = "Methow River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "methow_DE_plot.png"), methow_DE_plot, height = 8, width = 8)

### Snake River
asotin_DE_by_year_wild <- calculate_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "asotin", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array)
asotin_DE_plot <- plot_DE_rear(DE_by_year_wild = asotin_DE_by_year_wild, plot_title = "Asotin Creek Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "asotin_DE_plot.png"), asotin_DE_plot, height = 8, width = 8)

imnaha_DE_by_year_wild <- calculate_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "imnaha", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array)
imnaha_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = SRH_DE_param_matrix, tributary = "imnaha", tributary_design_matrices_array = SRH_envir$data$tributary_design_matrices_array)
imnaha_DE_plot <- plot_DE_rear(DE_by_year_wild = imnaha_DE_by_year_wild,
                               DE_by_year_hatchery = imnaha_DE_by_year_hatchery,
                               plot_title = "Imnaha River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "imnaha_DE_plot.png"), imnaha_DE_plot, height = 8, width = 8)

tucannon_DE_by_year_wild <- calculate_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "tucannon", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array)
tucannon_DE_by_year_hatchery <- calculate_DE(DE_param_matrix = SRH_DE_param_matrix, tributary = "tucannon", tributary_design_matrices_array = SRH_envir$data$tributary_design_matrices_array)
tucannon_DE_plot <- plot_DE_rear(DE_by_year_wild = tucannon_DE_by_year_wild,
                                 DE_by_year_hatchery = tucannon_DE_by_year_hatchery,
                                 plot_title = "Tucannon River Estimated Detection Probability")
ggsave(here::here("stan_actual", "output", "detection_efficiency", "tucannon_DE_plot.png"), tucannon_DE_plot, height = 8, width = 8)









#### Calculate average detection efficiency by tributary ####

# This is the same function from 09-model-temperature-plots.R

extract_covariate_experiences <- function(envir, rear, origin_select){
  origin_vector <- vector(length = nrow(envir$data$cat_X_mat))
  for(i in 1:nrow(envir$data$cat_X_mat)){
    origin_vector[i] <- which(envir$data$cat_X_mat[i,]==1)
  }
  
  
  
  # for spill days - include the winter post-overshoot vector, which contains
  # info on whether they could have experienced winter spill conditions or not
  # add a new fish_ID column, which is not the same as tag code but will allow us to differentiate between fish
  pop_states_dates <- data.frame(fish_ID = rep(1:length(origin_vector), each = ncol(envir$data$y)),
                                 state = as.vector(t(envir$data$y)),
                                 date = as.vector(t(envir$data$transition_dates)),
                                 year = ceiling(as.vector(t(envir$data$transition_dates))/365.25)+1,
                                 origin = rep(origin_vector, each = ncol(envir$data$y)))
  
  
  # add mainstem dam for each state
  dam_index <- data.frame(dam = c("BON", "MCN", "PRA", "RIS", "RRE", "WEL", "ICH", "LGR"),
                          state = seq(2,9))
  pop_states_dates %>% 
    left_join(., dam_index, by = "state") -> pop_states_dates
  
  
  # reformat covariates so that they can be joined
  as.data.frame(envir$data$spill_window_data) %>% 
    rownames_to_column("date") %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(cols = -c(date), names_to = "dam", values_to = "spill_window") -> spill_window_long
  
  as.data.frame(envir$data$temperature_data) %>% 
    rownames_to_column("date") %>% 
    mutate(date = as.numeric(date)) %>% 
    pivot_longer(cols = -c(date), names_to = "dam", values_to = "temperature") -> temp_long
  
  as.data.frame(envir$data$winter_spill_days_data) %>% 
    rownames_to_column("year") %>% 
    mutate(year = as.numeric(year)) %>% 
    pivot_longer(cols = -c(year), names_to = "dam", values_to = "winter_spill") -> spill_days_long
  
  
  # add temperature, spill window, winter spill days
  pop_states_dates %>% 
    left_join(., spill_window_long, by = c("date", "dam")) %>% 
    left_join(., temp_long, by = c("date", "dam")) %>% 
    left_join(., spill_days_long, by = c("year", "dam")) -> pop_states_dates_covariates
  
  # drop observations in the loss state and with index 0
  pop_states_dates_covariates %>% 
    filter(!(state %in% c(0,43))) -> pop_states_dates_covariates
  
  # now add winter post-overshoot vector
  pop_states_dates_covariates$winter_post_overshoot_vector <- as.vector(envir$data$winter_post_overshoot_vector)
  
  
  # Now, keep only the origin selected
  if(rear == "wild"){
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$wild
  } else {
    origin_numeric <- subset(origin_param_map, natal_origin == origin_select)$hatchery
  }
  
  subset(pop_states_dates_covariates, origin == origin_numeric) -> origin_states_dates_covariates
  
  return(origin_states_dates_covariates)
  
}



calculate_weighted_DE <- function(DE_param_matrix, tributary, tributary_design_matrices_array,
                                 covariate_experiences){
  
  ## calculate estimated detection efficiency by year
  # get the index for the tributary state (needs to be mouth since that's where we're correcting for DE)
  tributary_state <- intersect(grep(tributary, model_states, ignore.case = TRUE), grep("Mouth", model_states, ignore.case = TRUE))
  
  tributary_from_state_df <- data.frame(from = c(rep(2,5),3,3,5,6,7,7,8,9,9),
                                        to = grep("Mouth", model_states, ignore.case = FALSE))
  
  tributary_mainstem_state <- subset(tributary_from_state_df, to == tributary_state)$from
  
  # use the state indexing to get the right design matrix
  trib_design_matrix <- tributary_design_matrices_array[,,tributary_state]
  
  # create a matrix to store DE per year, per iter
  niter <- 4000
  DE_matrix <- matrix(nrow = nrow(trib_design_matrix), ncol = niter)
  
  # for each run year, get a confidence interval around detection efficiency by using the different draws
  for (i in 1:nrow(trib_design_matrix)){
    # if there is no intercept term, that means there is no DE correction for that year - so skip and leave as NA
    if(sum(trib_design_matrix[i,1:21]) == 0){
      
      
    } else {
      for (j in 1:niter){
        eta <- sum(trib_design_matrix[i,] * DE_param_matrix[,j])
        DE_matrix[i,j] <- exp(eta)/(1 + exp(eta))
      }
      
    }
    
    
  }
  
  colnames(DE_matrix) <- paste0("iter", 1:niter)
  
  DE_matrix %>% 
    as.data.frame() %>% 
    rownames_to_column("year") %>% 
    mutate(year = as.numeric(year)) %>% 
    pivot_longer(cols = starts_with("iter"), names_to = "iter", values_to = "detection_probability") %>% 
    group_by(year) %>% 
    summarise(prob = quantile(detection_probability, c(0.025, 0.5, 0.975), na.rm = T), q = c(0.025, 0.5, 0.975)) %>% 
    pivot_wider(names_from = q, values_from = prob) %>% 
    mutate(rear = "wild") -> DE_matrix_long
  
  DE_matrix_long$year_actual <- 2004:2021
  
  ## calculate number of transitions into that tributary by year
  covariate_experiences %>% 
    mutate(next_state = ifelse(lead(fish_ID) == fish_ID, lead(state), 43)) -> covariate_experiences
  
  # get a table of counts by run year, into that tributary
  as.data.frame(table(subset(covariate_experiences, state == tributary_mainstem_state & next_state == tributary_state)$year)) %>% 
    dplyr::rename(index = Var1, count = Freq) %>% 
    mutate(index = as.numeric(as.character(index))) -> trib_entries_by_year
  
  # add in the actual year (not index year). 2004 = year 1
  year_indices <- data.frame(index = 1:19, year_actual = 2004:2022)
  trib_entries_by_year %>% 
    left_join(year_indices, by = "index") -> trib_entries_by_year
  
  DE_matrix_long %>% 
    left_join(dplyr::select(trib_entries_by_year, count, year_actual), by = "year_actual") -> DE_by_year
  
  # create a weighted average
  # make sure to drop any years where DE isn't estimated in the model
  DE_by_year %>% 
    filter(!(is.na(`0.5`))) %>% 
    group_by(rear) %>% 
    mutate(weight = count/sum(count, na.rm = T)) %>% 
    summarise(weighted_DE = sum(`0.5`*weight, na.rm = T)) -> weighted_DE
  
  return(list(DE_by_year = DE_by_year, weighted_average = weighted_DE$weighted_DE))
}


# first, get all covariate experiences
DES_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Deschutes River")
JDR_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "John Day River")
FIF_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Fifteenmile Creek")
UMA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Umatilla River")
UMA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Umatilla River")
YAK_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Yakima River")
WAWA_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Walla Walla River")
WAWA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Walla Walla River")
WEN_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Wenatchee River")
WEN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Wenatchee River")
ENT_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Entiat River")
OKA_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Okanogan River")
MET_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Methow River")
MET_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Methow River")
ASO_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Asotin Creek")
IMN_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Imnaha River")
IMN_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Imnaha River")
TUC_hatchery_covariate_experiences <- extract_covariate_experiences(envir = MCH_envir, rear = "hatchery", origin_select = "Tucannon River")
TUC_wild_covariate_experiences <- extract_covariate_experiences(envir = MCW_envir, rear = "wild", origin_select = "Tucannon River")



### Middle Columbia
des_des_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                            covariate_experiences = DES_wild_covariate_experiences)
uma_hatchery_des_average_DE <- calculate_weighted_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array,
                                                     covariate_experiences = UMA_hatchery_covariate_experiences)
uma_wild_des_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "deschutes", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                                 covariate_experiences = UMA_wild_covariate_experiences)


john_day_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "john", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                             covariate_experiences = JDR_wild_covariate_experiences)

fifteenmile_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "fifteenmile", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                                covariate_experiences = FIF_wild_covariate_experiences)



umatilla_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "umatilla", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                             covariate_experiences = UMA_wild_covariate_experiences)
umatilla_DE_by_year_hatchery <- calculate_weighted_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "umatilla", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array,
                                                      covariate_experiences = UMA_wild_covariate_experiences)


yakima_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "yakima", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                           covariate_experiences = YAK_wild_covariate_experiences)


walla_walla_average_DE <- calculate_weighted_DE(DE_param_matrix = MCW_DE_param_matrix, tributary = "walla", tributary_design_matrices_array = MCW_envir$data$tributary_design_matrices_array,
                                                covariate_experiences = WAWA_wild_covariate_experiences)

walla_walla_DE_by_year_hatchery <- calculate_weighted_DE(DE_param_matrix = MCH_DE_param_matrix, tributary = "walla", tributary_design_matrices_array = MCH_envir$data$tributary_design_matrices_array,
                                                         covariate_experiences = WAWA_wild_covariate_experiences)

### Upper Columbia
wenatchee_average_DE <- calculate_weighted_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "wenatchee", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array,
                                              covariate_experiences = WEN_wild_covariate_experiences)

wenatchee_DE_by_year_hatchery <- calculate_weighted_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "wenatchee", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array,
                                                       covariate_experiences = WEN_wild_covariate_experiences)

entiat_average_DE <- calculate_weighted_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "entiat", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array,
                                           covariate_experiences = ENT_wild_covariate_experiences)

okanogan_DE_by_year_hatchery <- calculate_weighted_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "okanogan", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array,
                                                      covariate_experiences = OKA_wild_covariate_experiences)

methow_average_DE <- calculate_weighted_DE(DE_param_matrix = UCW_DE_param_matrix, tributary = "methow", tributary_design_matrices_array = UCW_envir$data$tributary_design_matrices_array,
                                           covariate_experiences = MET_wild_covariate_experiences)

methow_DE_by_year_hatchery <- calculate_weighted_DE(DE_param_matrix = UCH_DE_param_matrix, tributary = "methow", tributary_design_matrices_array = UCH_envir$data$tributary_design_matrices_array,
                                                    covariate_experiences = MET_wild_covariate_experiences)

### Snake River
asotin_average_DE <- calculate_weighted_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "asotin", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array,
                                           covariate_experiences = ASO_wild_covariate_experiences)

imnaha_average_DE <- calculate_weighted_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "imnaha", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array,
                                           covariate_experiences = IMN_wild_covariate_experiences)

imnaha_DE_by_year_hatchery <- calculate_weighted_DE(DE_param_matrix = SRH_DE_param_matrix, tributary = "imnaha", tributary_design_matrices_array = SRH_envir$data$tributary_design_matrices_array,
                                                    covariate_experiences = IMN_wild_covariate_experiences)


tucannon_average_DE <- calculate_weighted_DE(DE_param_matrix = SRW_DE_param_matrix, tributary = "tucannon", tributary_design_matrices_array = SRW_envir$data$tributary_design_matrices_array,
                                             covariate_experiences = TUC_wild_covariate_experiences)


tucannon_DE_by_year_hatchery <- calculate_weighted_DE(DE_param_matrix = SRH_DE_param_matrix, tributary = "tucannon", tributary_design_matrices_array = SRH_envir$data$tributary_design_matrices_array,
                                                      covariate_experiences = TUC_wild_covariate_experiences)


