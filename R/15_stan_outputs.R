# 15 stan_actual_outputs


# This script will load and investigate the outputs from our stan models.
library(cmdstanr)
library(rstan)
library(bayesplot)
library(posterior)
# library(rstan)

##### Snake River #####

# Plot the traceplots
# They look okay honestly, not horribly autocorrelated
# fit <- readRDS(here::here("stan_actual", "ESU_models", "snake", "100iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))
snake_fit <- readRDS(here::here("stan_actual", "ESU_models", "snake", "200iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))
# mcmc_trace(snake_fit$draws(), pars = c("borigin1_matrix_21_8"))

# Inspect parameter estimates
snake_fit_summary <- snake_fit$summary()

##### Extract our parameters of interest #####

##### Derived parameters/probabilities #####
# Look at the probabilities
snake_fit_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat)) -> snake_derived_probabilities

# see if these make sense
snake_derived_probabilities %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable))) -> snake_derived_probabilities

snake_derived_probabilities %>% 
  group_by(from) %>% 
  summarise(sum(mean))
# for six origins, this makes sense that it sums to six

snake_derived_probabilities %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) -> snake_derived_probabilities

snake_derived_probabilities %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Tucannon_River",
                         ifelse(origin == "origin2", "Asotin_Creek",
                                ifelse(origin == "origin3",  "Clearwater_River",
                                       ifelse(origin == "origin4", "Salmon_River",
                                              ifelse(origin == "origin5", "Grande_Ronde",
                                                     ifelse(origin == "origin6", "Imnaha_River", "error"))))))) -> snake_derived_probabilities 

subset(snake_derived_probabilities, origin == "Tucannon_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> tuc_mean_probabilities

subset(snake_derived_probabilities, origin == "Asotin_Creek",) %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> aso_mean_probabilities

subset(snake_derived_probabilities, origin == "Clearwater_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> clr_mean_probabilities

subset(snake_derived_probabilities, origin == "Salmon_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> sal_mean_probabilities

subset(snake_derived_probabilities, origin == "Grande_Ronde") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> grro_mean_probabilities

subset(snake_derived_probabilities, origin == "Imnaha_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> imn_mean_probabilities

# Inspect these probabilities
# Out of BON to MCN
dplyr::select(tuc_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(aso_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(clr_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(sal_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(grro_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(imn_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]

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
    # from_parameter_estimates <- subset(parameter_estimates, from == unique_from[i]) 
    from_parameter_estimates <- subset(snake_params, from == unique_from[i])
    
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


##### Upper Columbia #####

# Plot the traceplots
# They look okay honestly, not horribly autocorrelated
# fit <- readRDS(here::here("stan_actual", "ESU_models", "upper_columbia", "100_iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))
# Load the three chains
upper_columbia_seed101 <- readRDS(here::here("stan_actual", "ESU_models", "upper_columbia", "seed101_200_iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))
upper_columbia_seed102 <- readRDS(here::here("stan_actual", "ESU_models", "upper_columbia", "seed102_200_iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))
upper_columbia_seed103 <- readRDS(here::here("stan_actual", "ESU_models", "upper_columbia", "seed103_200_iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))
upper_columbia_3chains <- cmdstanr::read_cmdstan_csv(c(upper_columbia_seed101$output_files(), upper_columbia_seed102$output_files(), upper_columbia_seed103$output_files()))
# upper_columbia_fit <- readRDS(here::here("stan_actual", "ESU_models", "upper_columbia", "200_iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))
# mcmc_trace(fit$draws(), pars = c("borigin1_matrix_3_2"))
upper_columbia_seed101$draws()
merge_chains(upper_columbia_seed101$draws(), upper_columbia_seed102$draws(), upper_columbia_seed103$draws()) -> upcol_3chains_draws

# Inspect parameter estimates
upper_columbia_fit_summary <- upper_columbia_seed101$summary()

##### Derived parameters/probabilities #####
# Look at the probabilities
upper_columbia_fit_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat)) -> upper_columbia_derived_probabilities

# see if these make sense
upper_columbia_derived_probabilities %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable))) -> upper_columbia_derived_probabilities

upper_columbia_derived_probabilities %>% 
  group_by(from) %>% 
  summarise(sum(mean))
# for six origins, this makes sense that it sums to six

upper_columbia_derived_probabilities %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) -> upper_columbia_derived_probabilities

upper_columbia_derived_probabilities %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Wenatchee_River",
                         ifelse(origin == "origin2", "Entiat_River",
                                ifelse(origin == "origin3",  "Okanogan_River",
                                       ifelse(origin == "origin4", "Methow_River","error"))))) -> upper_columbia_derived_probabilities 

subset(upper_columbia_derived_probabilities, origin == "Wenatchee_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> wen_mean_probabilities

subset(upper_columbia_derived_probabilities, origin == "Entiat_River",) %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> ent_mean_probabilities

subset(upper_columbia_derived_probabilities, origin == "Okanogan_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> oka_mean_probabilities

subset(upper_columbia_derived_probabilities, origin == "Methow_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> met_mean_probabilities

# Inspect these probabilities
# Out of BON to MCN
dplyr::select(wen_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(ent_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(oka_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(met_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]


##### Extract our parameters of interest #####

# We would like to reformat this, so that we have a dataframe where rows = transitions and columns = different parameters
# The values - we want a 
upper_columbia_fit_summary %>% 
  filter(.,!(grepl(",", variable))) %>% 
  subset(., variable != "lp__") -> upper_columbia_params

# Create a transition column
upper_columbia_params %>% 
  # mutate(transition = str_extract(variable, "[^_]*_[^_]*"))
  mutate(transition = sub('.*matrix_', "", variable)) %>% 
  mutate(parameter = sub("_.*", "", variable)) -> upper_columbia_params

# Reformat this df
upper_columbia_params %>% 
  dplyr::select(transition, parameter, mean) %>% 
  pivot_wider(id_cols = "transition", names_from = "parameter", values_from = "mean") %>% 
  mutate(from = sub("_.*", "", transition)) %>% 
  mutate(to = sub(".*_", "", transition)) -> upper_columbia_mean

# reorder
upper_columbia_mean %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) -> upper_columbia_mean

# replace those NAs with zeros for calculations later
upper_columbia_mean[is.na(upper_columbia_mean)] <- 0


##### Calculate movement probabilities ####

# write a function to calculate mlogit probabilities
unique_from <- unique(upper_columbia_mean$from)

to_states_list <- vector(mode = "list", length = length(unique_from))
for (i in 1:length(unique_from)){
  j <- unique_from[i]
  subset(upper_columbia_mean, from == j) %>% 
    distinct(to) -> to_states_df
  
  to_states_list[[i]] <- as.vector(to_states_df$to)
}

# Create an empty df to store these
upper_columbia_movement_probabilities <- data.frame(from = upper_columbia_mean$from,
                                                    to = upper_columbia_mean$to,
                                                    Wenatchee_River = NA, 
                                                    Entiat_River = NA, 
                                                    Okanogan_River = NA, 
                                                    Methow_River = NA)

for (i in 1:length(to_states_list)){
  from_parameter_estimates <- subset(upper_columbia_mean, from == unique_from[i]) 
  
  for (j in 1:length(to_states_list[[i]])){
    # origin 1 - "Wenatchee_River"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin1[j])
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin1[k]) 
    }
    # populate with numerator/denominator
    upper_columbia_movement_probabilities[upper_columbia_movement_probabilities$from == unique_from[i],"Wenatchee_River"][j] <- numerator/denominator
    
    # origin 2 - "Entiat_River"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin2[j])
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin2[k]) 
    }
    # populate with numerator/denominator
    upper_columbia_movement_probabilities[upper_columbia_movement_probabilities$from == unique_from[i],"Entiat_River"][j] <- numerator/denominator
    
    # origin 3 -  "Okanogan_River"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin3[j])
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin3[k]) 
    }
    # populate with numerator/denominator
    upper_columbia_movement_probabilities[upper_columbia_movement_probabilities$from == unique_from[i],"Okanogan_River"][j] <- numerator/denominator
    
    # origin 4 - "Methow_River"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] - from_parameter_estimates$borigin1[j] -
                       from_parameter_estimates$borigin2[j] -
                       from_parameter_estimates$borigin3[j]
    )
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] - from_parameter_estimates$borigin1[k] -
                                         from_parameter_estimates$borigin2[k] -
                                         from_parameter_estimates$borigin3[k]) 
    }
    # populate with numerator/denominator
    upper_columbia_movement_probabilities[upper_columbia_movement_probabilities$from == unique_from[i],"Methow_River"][j] <- numerator/denominator
  }
}

# Look into ESU-level averages (based on intercept term, only) for certain movement probabilities
# 3 (MCN to ICH or PRA) to 4 (PRA to RIS)
exp(subset(upper_columbia_mean, transition == "3_4")$b0)/(1 + sum(exp(subset(upper_columbia_mean, from == 3)$b0)))

# 3 (MCN to ICH or PRA) to 8 (ICH to LGR)
exp(subset(upper_columbia_mean, transition == "3_8")$b0)/(1 + sum(exp(subset(upper_columbia_mean, from == 3)$b0)))

# 3 (MCN to ICH or PRA) to 2 (BON to MCN)
exp(subset(upper_columbia_mean, transition == "3_2")$b0)/(1 + sum(exp(subset(upper_columbia_mean, from == 3)$b0)))



# check on these numbers
upper_columbia_movement_probabilities %>% 
  group_by(from) %>% 
  summarise(sum(Wenatchee_River)) %>% 
  dplyr::rename(non_loss = `sum(Wenatchee_River)`) %>% 
  mutate(loss = 1 - non_loss) -> upper_columbia_tuc_loss_probs

# Convert into percentage probabilities for each natal origin

# Wenatchee_River
upper_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Wenatchee_River) %>% 
  mutate(Wenatchee_River = round(Wenatchee_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Wenatchee_River") %>% 
  column_to_rownames("from") -> Wenatchee_River_probability_matrix
# add loss probability
Wenatchee_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Wenatchee_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Wenatchee_River_probability_matrix

# Entiat_River
upper_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Entiat_River) %>% 
  mutate(Entiat_River = round(Entiat_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Entiat_River") %>% 
  column_to_rownames("from") -> Entiat_River_probability_matrix
# add loss probability
Entiat_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Entiat_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Entiat_River_probability_matrix

# Okanogan_River
upper_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Okanogan_River) %>% 
  mutate(Okanogan_River = round(Okanogan_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Okanogan_River") %>% 
  column_to_rownames("from") -> Okanogan_River_probability_matrix
# add loss probability
Okanogan_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Okanogan_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Okanogan_River_probability_matrix

# Methow_River
upper_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Methow_River) %>% 
  mutate(Methow_River = round(Methow_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Methow_River") %>% 
  column_to_rownames("from") -> Methow_River_probability_matrix
# add loss probability
Methow_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Methow_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Methow_River_probability_matrix


# States 1 and 2, all the same for this ESU
Wenatchee_River_probability_matrix[1,]
Wenatchee_River_probability_matrix[2,]

# Let's check transitions out of the MCN to ICH to PRA reach
Wenatchee_River_probability_matrix[3,]
Entiat_River_probability_matrix[3,]
Okanogan_River_probability_matrix[3,]
Methow_River_probability_matrix[3,]

# Let's check transitions out of the ICH to LGR reach
Wenatchee_River_probability_matrix[8,]
Entiat_River_probability_matrix[8,]
Okanogan_River_probability_matrix[8,]
Methow_River_probability_matrix[8,]

# Let's check transitions out of the upstream LGR reach
Wenatchee_River_probability_matrix[9,]
Entiat_River_probability_matrix[9,]
Okanogan_River_probability_matrix[9,]
Methow_River_probability_matrix[9,]




##### Middle Columbia #####

# Plot the traceplots
# They look okay honestly, not horribly autocorrelated
# fit <- readRDS(here::here("stan_actual", "ESU_models", "middle_columbia", "100_iter_parallel_middle_columbia_stan_actual_int_origin_stan_fit.rds"))
middle_columbia_fit <- readRDS(here::here("stan_actual", "ESU_models", "middle_columbia", "200_iter_parallel_middle_columbia_stan_actual_int_origin_stan_fit.rds"))
# mcmc_trace(fit$draws(), pars = c("borigin1_matrix_3_2"))

# Inspect parameter estimates
middle_columbia_fit_summary <- middle_columbia_fit$summary()



##### Extract our parameters of interest #####

##### Derived parameters/probabilities #####
# Look at the probabilities
middle_columbia_fit_summary %>% 
  filter(.,grepl("probs", variable))  %>% 
  filter(!is.na(rhat)) -> middle_columbia_derived_probabilities

# see if these make sense
middle_columbia_derived_probabilities %>% 
  # mutate(from = sub("\\[", "", variable))
  mutate(from = sub(",.*", "", sub(".*\\[", "", variable))) %>% 
  mutate(to = sub(".*,", "", sub("\\]", "", variable))) -> middle_columbia_derived_probabilities

middle_columbia_derived_probabilities %>% 
  group_by(from) %>% 
  summarise(sum(mean))
# for six origins, this makes sense that it sums to six

middle_columbia_derived_probabilities %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) -> middle_columbia_derived_probabilities

middle_columbia_derived_probabilities %>% 
  mutate(origin = sub("_.*", "", variable)) %>% 
  mutate(origin = ifelse(origin == "origin1", "Deschutes_River",
                         ifelse(origin == "origin2", "Fifteenmile_Creek",
                                ifelse(origin == "origin3",  "John_Day_River",
                                       ifelse(origin == "origin4", "Umatilla_River",
                                              ifelse(origin == "origin5", "Yakima_River",
                                                     ifelse(origin == "origin6", "Walla_Walla_River", "error"))))))) -> middle_columbia_derived_probabilities 

subset(middle_columbia_derived_probabilities, origin == "Deschutes_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> des_mean_probabilities

subset(middle_columbia_derived_probabilities, origin == "Fifteenmile_Creek",) %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> fif_mean_probabilities

subset(middle_columbia_derived_probabilities, origin == "John_Day_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> jdr_mean_probabilities

subset(middle_columbia_derived_probabilities, origin == "Umatilla_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> uma_mean_probabilities

subset(middle_columbia_derived_probabilities, origin == "Yakima_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> yak_mean_probabilities

subset(middle_columbia_derived_probabilities, origin == "Walla_Walla_River") %>% 
  dplyr::select(mean, to, from) %>% 
  arrange(to) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "mean") %>% 
  arrange(from) -> wawa_mean_probabilities

# Inspect these probabilities
# Out of BON to MCN
dplyr::select(des_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(fif_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(jdr_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(uma_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(yak_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]
dplyr::select(wawa_mean_probabilities, "1", "3", "10", "11", "12", "13", "14", "27", "29")[2,]

# Out of MCN to ICH or PRA
dplyr::select(des_mean_probabilities, "2", "4", "8", "15", "16", "29")[3,]
dplyr::select(fif_mean_probabilities, "2", "4", "8", "15", "16", "29")[3,]
dplyr::select(jdr_mean_probabilities, "2", "4", "8", "15", "16", "29")[3,]
dplyr::select(uma_mean_probabilities, "2", "4", "8", "15", "16", "29")[3,]
dplyr::select(yak_mean_probabilities, "2", "4", "8", "15", "16", "29")[3,]
dplyr::select(wawa_mean_probabilities, "2", "4", "8", "15", "16", "29")[3,]


##### Model parameters #####

# We would like to reformat this, so that we have a dataframe where rows = transitions and columns = different parameters
# The values - we want a 
middle_columbia_fit_summary %>% 
  filter(.,!(grepl(",", variable))) %>% 
  subset(., variable != "lp__") -> middle_columbia_params

# Create a transition column
middle_columbia_params %>% 
  # mutate(transition = str_extract(variable, "[^_]*_[^_]*"))
  mutate(transition = sub('.*matrix_', "", variable)) %>% 
  mutate(parameter = sub("_.*", "", variable)) -> middle_columbia_params

# Reformat this df
middle_columbia_params %>% 
  dplyr::select(transition, parameter, mean) %>% 
  pivot_wider(id_cols = "transition", names_from = "parameter", values_from = "mean") %>% 
  mutate(from = sub("_.*", "", transition)) %>% 
  mutate(to = sub(".*_", "", transition)) -> middle_columbia_mean

# reorder
middle_columbia_mean %>% 
  mutate(from = as.numeric(from), to = as.numeric(to)) %>% 
  arrange(from, to) -> middle_columbia_mean

# replace those NAs with zeros for calculations later
middle_columbia_mean[is.na(middle_columbia_mean)] <- 0


##### Calculate movement probabilities ####

# write a function to calculate mlogit probabilities
unique_from <- unique(middle_columbia_mean$from)

to_states_list <- vector(mode = "list", length = length(unique_from))
for (i in 1:length(unique_from)){
  j <- unique_from[i]
  subset(middle_columbia_mean, from == j) %>% 
    distinct(to) -> to_states_df
  
  to_states_list[[i]] <- as.vector(to_states_df$to)
}

# Create an empty df to store these
middle_columbia_movement_probabilities <- data.frame(from = middle_columbia_mean$from,
                                           to = middle_columbia_mean$to,
                                           Deschutes_River = NA, 
                                           Fifteenmile_Creek = NA, 
                                           John_Day_River = NA, 
                                           Umatilla_River = NA, 
                                           Yakima_River = NA, 
                                           Walla_Walla_River = NA)

for (i in 1:length(to_states_list)){
  from_parameter_estimates <- subset(middle_columbia_mean, from == unique_from[i]) 
  
  for (j in 1:length(to_states_list[[i]])){
    # origin 1 - "Deschutes_River"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin1[j])
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin1[k]) 
    }
    # populate with numerator/denominator
    middle_columbia_movement_probabilities[middle_columbia_movement_probabilities$from == unique_from[i],"Deschutes_River"][j] <- numerator/denominator
    
    # origin 2 - "Fifteenmile_Creek"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin2[j])
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin2[k]) 
    }
    # populate with numerator/denominator
    middle_columbia_movement_probabilities[middle_columbia_movement_probabilities$from == unique_from[i],"Fifteenmile_Creek"][j] <- numerator/denominator
    
    # origin 3 -  "John_Day_River"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin3[j])
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin3[k]) 
    }
    # populate with numerator/denominator
    middle_columbia_movement_probabilities[middle_columbia_movement_probabilities$from == unique_from[i],"John_Day_River"][j] <- numerator/denominator
    
    # origin 4 - "Umatilla_River"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin4[j])
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin4[k]) 
    }
    # populate with numerator/denominator
    middle_columbia_movement_probabilities[middle_columbia_movement_probabilities$from == unique_from[i],"Umatilla_River"][j] <- numerator/denominator
    
    # origin 5 - "Yakima_River"
    # numerator
    numerator <- exp(from_parameter_estimates$b0[j] + from_parameter_estimates$borigin5[j])
    denominator <- 1
    for (k in 1:length(to_states_list[[i]])){
      denominator <- denominator + exp(from_parameter_estimates$b0[k] + from_parameter_estimates$borigin5[k]) 
    }
    # populate with numerator/denominator
    middle_columbia_movement_probabilities[middle_columbia_movement_probabilities$from == unique_from[i],"Yakima_River"][j] <- numerator/denominator
    
    # origin 6 - "Walla_Walla_River"
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
    middle_columbia_movement_probabilities[middle_columbia_movement_probabilities$from == unique_from[i],"Walla_Walla_River"][j] <- numerator/denominator
  }
}

# Check on specific ones - which way are middle columbia steelhead overshooting?
# 3 (MCN to ICH or PRA) to 4 (PRA to RIS)
exp(subset(middle_columbia_mean, transition == "3_4")$b0)/(1 + sum(exp(subset(middle_columbia_mean, from == 3)$b0)))

# 3 (MCN to ICH or PRA) to 8 (ICH to LGR)
exp(subset(middle_columbia_mean, transition == "3_8")$b0)/(1 + sum(exp(subset(middle_columbia_mean, from == 3)$b0)))

# 3 (MCN to ICH or PRA) to 2 (BON to MCN)
exp(subset(middle_columbia_mean, transition == "3_2")$b0)/(1 + sum(exp(subset(middle_columbia_mean, from == 3)$b0)))

# check on these numbers
middle_columbia_movement_probabilities %>% 
  group_by(from) %>% 
  summarise(sum(Deschutes_River)) %>% 
  dplyr::rename(non_loss = `sum(Deschutes_River)`) %>% 
  mutate(loss = 1 - non_loss) -> middle_columbia_des_loss_probs

# Convert into percentage probabilities for each natal origin

# Deschutes_River
middle_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Deschutes_River) %>% 
  mutate(Deschutes_River = round(Deschutes_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Deschutes_River") %>% 
  column_to_rownames("from") -> Deschutes_River_probability_matrix
# add loss probability
Deschutes_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Deschutes_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Deschutes_River_probability_matrix

# Fifteenmile_Creek
middle_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Fifteenmile_Creek) %>% 
  mutate(Fifteenmile_Creek = round(Fifteenmile_Creek*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Fifteenmile_Creek") %>% 
  column_to_rownames("from") -> Fifteenmile_Creek_probability_matrix
# add loss probability
Fifteenmile_Creek_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Fifteenmile_Creek_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Fifteenmile_Creek_probability_matrix

# John_Day_River
middle_columbia_movement_probabilities %>% 
  dplyr::select(from, to, John_Day_River) %>% 
  mutate(John_Day_River = round(John_Day_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "John_Day_River") %>% 
  column_to_rownames("from") -> John_Day_River_probability_matrix
# add loss probability
John_Day_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(John_Day_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> John_Day_River_probability_matrix

# Umatilla_River
middle_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Umatilla_River) %>% 
  mutate(Umatilla_River = round(Umatilla_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Umatilla_River") %>% 
  column_to_rownames("from") -> Umatilla_River_probability_matrix
# add loss probability
Umatilla_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Umatilla_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Umatilla_River_probability_matrix

# Yakima_River
middle_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Yakima_River) %>% 
  mutate(Yakima_River = round(Yakima_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Yakima_River") %>% 
  column_to_rownames("from") -> Yakima_River_probability_matrix
# add loss probability
Yakima_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Yakima_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Yakima_River_probability_matrix

# Walla_Walla_River
middle_columbia_movement_probabilities %>% 
  dplyr::select(from, to, Walla_Walla_River) %>% 
  mutate(Walla_Walla_River = round(Walla_Walla_River*100,1)) %>% 
  pivot_wider(id_cols = "from", names_from = "to", values_from = "Walla_Walla_River") %>% 
  column_to_rownames("from") -> Walla_Walla_River_probability_matrix
# add loss probability
Walla_Walla_River_probability_matrix %>% 
  bind_cols(., non_loss = rowSums(Walla_Walla_River_probability_matrix, na.rm = TRUE)) %>% 
  mutate(loss = 100 - non_loss) %>% 
  dplyr::select(-non_loss) -> Walla_Walla_River_probability_matrix



# States 1 and 2, all the same for this ESU
Deschutes_River_probability_matrix[1,]
Deschutes_River_probability_matrix[2,]

# Let's check transitions out of the MCN to ICH to PRA reach
Deschutes_River_probability_matrix[3,]
Fifteenmile_Creek_probability_matrix[3,]
John_Day_River_probability_matrix[3,]
Umatilla_River_probability_matrix[3,]
Yakima_River_probability_matrix[3,]
Walla_Walla_River_probability_matrix[3,]

# Let's check transitions out of the ICH to LGR reach
Deschutes_River_probability_matrix[8,]
Fifteenmile_Creek_probability_matrix[8,]
John_Day_River_probability_matrix[8,]
Umatilla_River_probability_matrix[8,]
Yakima_River_probability_matrix[8,]
Walla_Walla_River_probability_matrix[8,]

# Let's check transitions out of the upstream LGR reach
Deschutes_River_probability_matrix[9,]
Fifteenmile_Creek_probability_matrix[9,]
John_Day_River_probability_matrix[9,]
Umatilla_River_probability_matrix[9,]
Yakima_River_probability_matrix[9,]
Walla_Walla_River_probability_matrix[9,]


##### Compare ESUs at different locations #####

# state 3 - the branch state

# 3 (MCN to ICH or PRA) to 4 (PRA to RIS)
exp(subset(middle_columbia_mean, transition == "3_4")$b0)/(1 + sum(exp(subset(middle_columbia_mean, from == 3)$b0)))
exp(subset(upper_columbia_mean, transition == "3_4")$b0)/(1 + sum(exp(subset(upper_columbia_mean, from == 3)$b0)))
exp(subset(snake_mean, transition == "3_4")$b0)/(1 + sum(exp(subset(snake_mean, from == 3)$b0)))


# 3 (MCN to ICH or PRA) to 8 (ICH to LGR)
exp(subset(middle_columbia_mean, transition == "3_8")$b0)/(1 + sum(exp(subset(middle_columbia_mean, from == 3)$b0)))
exp(subset(upper_columbia_mean, transition == "3_8")$b0)/(1 + sum(exp(subset(upper_columbia_mean, from == 3)$b0)))
exp(subset(snake_mean, transition == "3_8")$b0)/(1 + sum(exp(subset(snake_mean, from == 3)$b0)))


# 3 (MCN to ICH or PRA) to 2 (BON to MCN)
exp(subset(middle_columbia_mean, transition == "3_2")$b0)/(1 + sum(exp(subset(middle_columbia_mean, from == 3)$b0)))
exp(subset(upper_columbia_mean, transition == "3_2")$b0)/(1 + sum(exp(subset(upper_columbia_mean, from == 3)$b0)))
exp(subset(snake_mean, transition == "3_2")$b0)/(1 + sum(exp(subset(snake_mean, from == 3)$b0)))

