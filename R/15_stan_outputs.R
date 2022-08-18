# 15 stan_actual_outputs


# This script will load and investigate the outputs from our stan models.
library(cmdstanr)
library(rstan)
library(bayesplot)

# Plot the traceplots
# They look okay honestly, not horribly autocorrelated
mcmc_trace(fit$draws(), pars = c("borigin1_matrix_21_8"))

# Inspect parameter estimates
snake_fit_summary <- fit$summary()



##### Extract our parameters of interest #####

# We would like to reformat this, so that we have a dataframe where rows = transitions and columns = different parameters
# The values - we want a 
snake_fit_summary %>% 
  filter(.,!(grepl(",", variable))) %>% 
  subset(., variable != "lp__") -> snake_params

# Create a transition column
snake_params %>% 
  # mutate(transition = str_extract(variable, "[^_]*_[^_]*"))
  mutate(transition = sub('.*matrix_', "", variable)) %>% 
  mutate(parameter = sub("_.*", "", variable)) -> snake_params

# Reformat this df
snake_params %>% 
  dplyr::select(transition, parameter, mean) %>% 
  pivot_wider(id_cols = "transition", names_from = "parameter", values_from = "mean") %>% 
  mutate(from = sub("_.*", "", transition)) %>% 
  mutate(to = sub(".*_", "", transition)) -> snake_mean

# reorder
snake_mean %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) -> snake_mean

# replace those NAs with zeros for calculations later
snake_mean[is.na(snake_mean)] <- 0


##### Calculate movement probabilities ####

# write a function to calculate mlogit probabilities
unique_from <- unique(snake_mean$from)

to_states_list <- vector(mode = "list", length = length(unique_from))
for (i in 1:length(unique_from)){
  j <- unique_from[i]
  subset(snake_mean, from == j) %>% 
    distinct(to) -> to_states_df
  
  to_states_list[[i]] <- as.vector(to_states_df$to)
}

# Create an empty df to store these
snake_movement_probabilities <- data.frame(from = snake_mean$from,
                                           to = snake_mean$to,
                                           Tucannon_River = NA, 
                                           Asotin_Creek = NA, 
                                           Clearwater_River = NA, 
                                           Salmon_River = NA, 
                                           Grande_Ronde_River = NA, 
                                           Imnaha_River = NA)

  for (i in 1:length(to_states_list)){
    from_parameter_estimates <- subset(parameter_estimates, from == unique_from[i]) 
    
    for (j in 1:length(to_states_list[[i]])){
      # origin 1 - "Tucannon_River"
      # numerator
      numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin1[j])
      denominator <- 1
      for (k in 1:length(to_states_list[[i]])){
        denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin1[k]) 
      }
      # populate with numerator/denominator
      snake_movement_probabilities[snake_movement_probabilities$from == unique_from[i],"Tucannon_River"][j] <- numerator/denominator
      
      # origin 2 - "Asotin_Creek"
      # numerator
      numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin2[j])
      denominator <- 1
      for (k in 1:length(to_states_list[[i]])){
        denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin2[k]) 
      }
      # populate with numerator/denominator
      snake_movement_probabilities[snake_movement_probabilities$from == unique_from[i],"Asotin_Creek"][j] <- numerator/denominator
      
      # origin 3 -  "Clearwater_River"
      # numerator
      numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin3[j])
      denominator <- 1
      for (k in 1:length(to_states_list[[i]])){
        denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin3[k]) 
      }
      # populate with numerator/denominator
      snake_movement_probabilities[snake_movement_probabilities$from == unique_from[i],"Clearwater_River"][j] <- numerator/denominator
      
      # origin 4 - "Salmon_River"
      # numerator
      numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin4[j])
      denominator <- 1
      for (k in 1:length(to_states_list[[i]])){
        denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin4[k]) 
      }
      # populate with numerator/denominator
      snake_movement_probabilities[snake_movement_probabilities$from == unique_from[i],"Salmon_River"][j] <- numerator/denominator
      
      # origin 5 - "Grande_Ronde_River"
      # numerator
      numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin5[j])
      denominator <- 1
      for (k in 1:length(to_states_list[[i]])){
        denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin5[k]) 
      }
      # populate with numerator/denominator
      snake_movement_probabilities[snake_movement_probabilities$from == unique_from[i],"Grande_Ronde_River"][j] <- numerator/denominator
      
      # origin 6 - "Imnaha_River"
      # numerator
      numerator <- exp(from_parameter_estimates$b0[j] - from_parameter_estimates$borigin1[j] -
                        from_parameter_estimates$borigin2[j] -
                         from_parameter_estimates$borigin3[j] -
                         from_parameter_estimates$borigin4[j] -
                         from_parameter_estimates$borigin5[j]
                       )
      denominator <- 1
      for (k in 1:length(to_states_list[[i]])){
        denominator <- denominator + exp(from_parameter_estimates$b0[k] - from_parameter_estimates$borigin1[k] -
                                           from_parameter_estimates$borigin2[k] -
                                           from_parameter_estimates$borigin3[k] -
                                           from_parameter_estimates$borigin4[k] -
                                           from_parameter_estimates$borigin5[k]) 
      }
      # populate with numerator/denominator
      snake_movement_probabilities[snake_movement_probabilities$from == unique_from[i],"Imnaha_River"][j] <- numerator/denominator
    }
  }

# check on these numbers
snake_movement_probabilities %>% 
  group_by(from) %>% 
  summarise(sum(Tucannon_River)) %>% 
  dplyr::rename(non_loss = `sum(Tucannon_River)`) %>% 
  mutate(loss = 1 - non_loss) -> snake_tuc_loss_probs

# Convert into percentage probabilities for each natal origin

# Tucannon_River
snake_movement_probabilities %>% 
  dplyr::select(from, to, Tucannon_River) %>% 
  mutate(Tucannon_River = round(Tucannon_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Tucannon_River") %>% 
  column_to_rownames("from") -> Tucannon_River_probability_matrix
# add loss probability
Tucannon_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Tucannon_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Tucannon_River_probability_matrix

# Asotin_Creek
snake_movement_probabilities %>% 
  dplyr::select(from, to, Asotin_Creek) %>% 
  mutate(Asotin_Creek = round(Asotin_Creek*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Asotin_Creek") %>% 
  column_to_rownames("from") -> Asotin_Creek_probability_matrix
# add loss probability
Asotin_Creek_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Asotin_Creek_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Asotin_Creek_probability_matrix

# Clearwater_River
snake_movement_probabilities %>% 
  dplyr::select(from, to, Clearwater_River) %>% 
  mutate(Clearwater_River = round(Clearwater_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Clearwater_River") %>% 
  column_to_rownames("from") -> Clearwater_River_probability_matrix
# add loss probability
Clearwater_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Clearwater_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Clearwater_River_probability_matrix

# Salmon_River
snake_movement_probabilities %>% 
  dplyr::select(from, to, Salmon_River) %>% 
  mutate(Salmon_River = round(Salmon_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Salmon_River") %>% 
  column_to_rownames("from") -> Salmon_River_probability_matrix
# add loss probability
Salmon_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Salmon_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Salmon_River_probability_matrix

# Grande_Ronde_River
snake_movement_probabilities %>% 
  dplyr::select(from, to, Grande_Ronde_River) %>% 
  mutate(Grande_Ronde_River = round(Grande_Ronde_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Grande_Ronde_River") %>% 
  column_to_rownames("from") -> Grande_Ronde_River_probability_matrix
# add loss probability
Grande_Ronde_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Grande_Ronde_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Grande_Ronde_River_probability_matrix

# Imnaha_River
snake_movement_probabilities %>% 
  dplyr::select(from, to, Imnaha_River) %>% 
  mutate(Imnaha_River = round(Imnaha_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Imnaha_River") %>% 
  column_to_rownames("from") -> Imnaha_River_probability_matrix
# add loss probability
Imnaha_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Imnaha_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Imnaha_River_probability_matrix



# States 1 and 2, all the same for this ESU
Tucannon_River_probability_matrix[1,]
Tucannon_River_probability_matrix[2,]

# Let's check transitions out of the MCN to ICH to PRA reach
Tucannon_River_probability_matrix[3,]
Asotin_Creek_probability_matrix[3,]
Clearwater_River_probability_matrix[3,]
Salmon_River_probability_matrix[3,]
Grande_Ronde_River_probability_matrix[3,]
Imnaha_River_probability_matrix[3,]

# Let's check transitions out of the ICH to LGR reach
Tucannon_River_probability_matrix[8,]
Asotin_Creek_probability_matrix[8,]
Clearwater_River_probability_matrix[8,]
Salmon_River_probability_matrix[8,]
Grande_Ronde_River_probability_matrix[8,]
Imnaha_River_probability_matrix[8,]

# Let's check transitions out of the upstream LGR reach
Tucannon_River_probability_matrix[9,]
Asotin_Creek_probability_matrix[9,]
Clearwater_River_probability_matrix[9,]
Salmon_River_probability_matrix[9,]
Grande_Ronde_River_probability_matrix[9,]
Imnaha_River_probability_matrix[9,]


