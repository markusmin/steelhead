# stan DE model diagnostics

library(cmdstanr)
library(posterior)
library(tidyverse)
library(lubridate)
library(bayesplot)
library(here)
library(shinystan)

# Load the model run
# snake_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "diagnostics", "snake6", "50iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))

# 2022-11-30 - load the run with only 40iter, but that doesn't have divergent transitions
# snake_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "diagnostics", "snake22", "100iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))

# 2022-12-03 - load the vectorized 1stan run
snake_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "diagnostics", "snake_vectorized_1", "20iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))

# 2022-12-04 - load 200 iter runs from each ESU
snake_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "snake", "200iter_parallel_snake_stan_actual_int_origin_stan_fit.rds"))
middle_columbia_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "middle_columbia", "200iter_parallel_middle_columbia_stan_actual_int_origin_stan_fit.rds"))
upper_columbia_fit <- readRDS(here::here("stan_actual", "deteff_ESU_models", "results", "200iter", "upper_columbia", "200iter_parallel_upper_columbia_stan_actual_int_origin_stan_fit.rds"))

# plot some traces
# acf(fit$draws()[,2,3])
# mcmc_trace(fit$draws(), pars = c("b0_matrix_1_2"))
# snake_summary$variable

# Note all of the parameters that we want to trace:
snake_model_params = c("b0_matrix_1_2", "b0_matrix_2_1", "b0_matrix_2_3", "b0_matrix_2_10_DE", "b0_matrix_2_10_NDE"      ,
                 "b0_matrix_2_12_DE", "b0_matrix_2_12_NDE", "b0_matrix_2_14_DE", "b0_matrix_2_14_NDE", "b0_matrix_2_16_DE", "b0_matrix_2_16_NDE"   ,
                 "b0_matrix_2_18_DE", "b0_matrix_2_18_NDE", "b0_matrix_2_41", "b0_matrix_3_2", "b0_matrix_3_4", "b0_matrix_3_8"           ,
                 "b0_matrix_3_20_DE", "b0_matrix_3_20_NDE", "b0_matrix_3_22_DE", "b0_matrix_3_22_NDE", "b0_matrix_4_3", "b0_matrix_4_5",
                 "b0_matrix_5_4", "b0_matrix_5_6", "b0_matrix_5_24_DE", "b0_matrix_5_24_NDE", "b0_matrix_6_5", "b0_matrix_6_7"           ,
                 "b0_matrix_6_26_DE", "b0_matrix_6_26_NDE", "b0_matrix_7_6", "b0_matrix_7_28_DE", "b0_matrix_7_28_NDE", "b0_matrix_7_30_DE"       ,
                 "b0_matrix_7_30_NDE", "b0_matrix_7_42", "b0_matrix_8_3", "b0_matrix_8_9", "b0_matrix_8_32_DE", "b0_matrix_8_32_NDE"      ,
                 "b0_matrix_9_8", "b0_matrix_9_34_DE", "b0_matrix_9_34_NDE", "b0_matrix_9_36", "b0_matrix_9_37", "b0_matrix_9_38",
                 "b0_matrix_9_39_DE", "b0_matrix_9_39_NDE", "b0_matrix_10_2", "b0_matrix_12_2", "b0_matrix_14_2", "b0_matrix_16_2"          ,
                 "b0_matrix_18_2", "b0_matrix_20_3", "b0_matrix_22_3", "b0_matrix_24_5", "b0_matrix_26_6", "b0_matrix_28_7"          ,
                 "b0_matrix_30_7", "b0_matrix_32_8", "b0_matrix_34_9", "b0_matrix_36_9", "b0_matrix_37_9", "b0_matrix_38_9"          ,
                 "b0_matrix_39_9", "b0_matrix_41_2", "b0_matrix_42_7", "borigin1_matrix_3_8", "borigin1_matrix_8_3", "borigin1_matrix_8_9"     ,
                 "borigin1_matrix_8_32_DE", "borigin1_matrix_8_32_NDE", "borigin1_matrix_9_8", "borigin1_matrix_9_34_DE", "borigin1_matrix_9_34_NDE", "borigin1_matrix_9_36"    ,
                 "borigin1_matrix_9_37", "borigin1_matrix_9_38", "borigin1_matrix_9_39_DE", "borigin1_matrix_9_39_NDE", "borigin1_matrix_32_8", "borigin1_matrix_34_9"    ,
                 "borigin1_matrix_36_9", "borigin1_matrix_37_9", "borigin1_matrix_38_9", "borigin1_matrix_39_9", "borigin2_matrix_3_8", "borigin2_matrix_8_3"     ,
                 "borigin2_matrix_8_9", "borigin2_matrix_8_32_DE", "borigin2_matrix_8_32_NDE", "borigin2_matrix_9_8", "borigin2_matrix_9_34_DE", "borigin2_matrix_9_34_NDE",
                 "borigin2_matrix_9_36", "borigin2_matrix_9_37", "borigin2_matrix_9_38", "borigin2_matrix_9_39_DE", "borigin2_matrix_9_39_NDE", "borigin2_matrix_32_8"    ,
                 "borigin2_matrix_34_9", "borigin2_matrix_36_9", "borigin2_matrix_37_9", "borigin2_matrix_38_9", "borigin2_matrix_39_9", "borigin3_matrix_3_8"     ,
                 "borigin3_matrix_8_3", "borigin3_matrix_8_9", "borigin3_matrix_8_32_DE", "borigin3_matrix_8_32_NDE", "borigin3_matrix_9_8", "borigin3_matrix_9_34_DE" ,
                 "borigin3_matrix_9_34_NDE", "borigin3_matrix_9_36", "borigin3_matrix_9_37", "borigin3_matrix_9_38", "borigin3_matrix_9_39_DE", "borigin3_matrix_9_39_NDE",
                 "borigin3_matrix_32_8", "borigin3_matrix_34_9", "borigin3_matrix_36_9", "borigin3_matrix_37_9", "borigin3_matrix_38_9", "borigin3_matrix_39_9"    ,
                 "borigin4_matrix_3_8", "borigin4_matrix_8_3", "borigin4_matrix_8_9", "borigin4_matrix_8_32_DE", "borigin4_matrix_8_32_NDE", "borigin4_matrix_9_8"   ,
                 "borigin4_matrix_9_34_DE", "borigin4_matrix_9_34_NDE", "borigin4_matrix_9_36", "borigin4_matrix_9_37", "borigin4_matrix_9_38", "borigin4_matrix_9_39_DE" ,
                 "borigin4_matrix_9_39_NDE", "borigin4_matrix_32_8", "borigin4_matrix_34_9", "borigin4_matrix_36_9", "borigin4_matrix_37_9", "borigin4_matrix_38_9"    ,
                 "borigin4_matrix_39_9", "borigin5_matrix_3_8", "borigin5_matrix_8_3", "borigin5_matrix_8_9", "borigin5_matrix_8_32_DE", "borigin5_matrix_8_32_NDE",
                 "borigin5_matrix_9_8", "borigin5_matrix_9_34_DE", "borigin5_matrix_9_34_NDE", "borigin5_matrix_9_36", "borigin5_matrix_9_37", "borigin5_matrix_9_38" ,
                 "borigin5_matrix_9_39_DE", "borigin5_matrix_9_39_NDE", "borigin5_matrix_32_8", "borigin5_matrix_34_9", "borigin5_matrix_36_9", "borigin5_matrix_37_9",
                 "borigin5_matrix_38_9", "borigin5_matrix_39_9", "asotin_alpha1", "asotin_alpha2", "deschutes_alpha1", "entiat_alpha1"           ,
                 "fifteenmile_alpha1", "hood_alpha1", "imnaha_alpha1", "john_day_alpha1", "methow_alpha1", "methow_alpha2"           ,
                 "okanogan_alpha1", "tucannon_alpha1", "tucannon_alpha2", "umatilla_alpha1", "umatilla_alpha2", "walla_walla_alpha1",
                 "walla_walla_alpha2", "walla_walla_alpha3", "wenatchee_alpha1", "yakima_alpha1", "asotin_beta", "deschutes_beta"         ,
                 "entiat_beta", "hood_beta", "john_day_beta", "methow_beta"  ,
                 "okanogan_beta", "tucannon_beta", "umatilla_beta", "walla_walla_beta", "wenatchee_beta", "yakima_beta")

middle_columbia_model_params = c("b0_matrix_1_2", 
                                 "b0_matrix_2_1", 
                                 "b0_matrix_2_3", 
                                 "b0_matrix_2_10_DE", 
                                 "b0_matrix_2_10_NDE", 
                                 "b0_matrix_2_12_DE", 
                                 "b0_matrix_2_12_NDE", 
                                 "b0_matrix_2_14_DE", 
                                 "b0_matrix_2_14_NDE", 
                                 "b0_matrix_2_16_DE", 
                                 "b0_matrix_2_16_NDE", 
                                 "b0_matrix_2_18_DE", 
                                 "b0_matrix_2_18_NDE", 
                                 "b0_matrix_2_41", 
                                 "b0_matrix_3_2", 
                                 "b0_matrix_3_4", 
                                 "b0_matrix_3_8", 
                                 "b0_matrix_3_20_DE", 
                                 "b0_matrix_3_20_NDE", 
                                 "b0_matrix_3_22_DE", 
                                 "b0_matrix_3_22_NDE", 
                                 "b0_matrix_4_3", 
                                 "b0_matrix_4_5", 
                                 "b0_matrix_5_4", 
                                 "b0_matrix_5_6", 
                                 "b0_matrix_5_24_DE", 
                                 "b0_matrix_5_24_NDE", 
                                 "b0_matrix_6_5", 
                                 "b0_matrix_6_7", 
                                 "b0_matrix_6_26_DE", 
                                 "b0_matrix_6_26_NDE", 
                                 "b0_matrix_7_6", 
                                 "b0_matrix_7_28_DE", 
                                 "b0_matrix_7_28_NDE", 
                                 "b0_matrix_7_30_DE", 
                                 "b0_matrix_7_30_NDE", 
                                 "b0_matrix_7_42", 
                                 "b0_matrix_8_3", 
                                 "b0_matrix_8_9", 
                                 "b0_matrix_8_32_DE", 
                                 "b0_matrix_8_32_NDE", 
                                 "b0_matrix_9_8", 
                                 "b0_matrix_9_34_DE", 
                                 "b0_matrix_9_34_NDE", 
                                 "b0_matrix_9_36", 
                                 "b0_matrix_9_37", 
                                 "b0_matrix_9_38", 
                                 "b0_matrix_9_39_DE", 
                                 "b0_matrix_9_39_NDE", 
                                 "b0_matrix_10_2", 
                                 "b0_matrix_12_2", 
                                 "b0_matrix_14_2", 
                                 "b0_matrix_16_2", 
                                 "b0_matrix_18_2", 
                                 "b0_matrix_20_3", 
                                 "b0_matrix_22_3", 
                                 "b0_matrix_24_5", 
                                 "b0_matrix_26_6", 
                                 "b0_matrix_28_7", 
                                 "b0_matrix_30_7", 
                                 "b0_matrix_32_8", 
                                 "b0_matrix_34_9", 
                                 "b0_matrix_36_9", 
                                 "b0_matrix_37_9", 
                                 "b0_matrix_38_9", 
                                 "b0_matrix_39_9", 
                                 "b0_matrix_41_2", 
                                 "b0_matrix_42_7", 
                                 
                                 "borigin1_matrix_1_2", 
                                 "borigin1_matrix_2_1", 
                                 "borigin1_matrix_2_3", 
                                 "borigin1_matrix_2_10_DE", 
                                 "borigin1_matrix_2_10_NDE", 
                                 "borigin1_matrix_2_12_DE", 
                                 "borigin1_matrix_2_12_NDE", 
                                 "borigin1_matrix_2_14_DE", 
                                 "borigin1_matrix_2_14_NDE", 
                                 "borigin1_matrix_2_16_DE", 
                                 "borigin1_matrix_2_16_NDE", 
                                 "borigin1_matrix_2_18_DE", 
                                 "borigin1_matrix_2_18_NDE", 
                                 "borigin1_matrix_2_41", 
                                 "borigin1_matrix_3_2", 
                                 "borigin1_matrix_3_4", 
                                 "borigin1_matrix_3_8", 
                                 "borigin1_matrix_3_20_DE", 
                                 "borigin1_matrix_3_20_NDE", 
                                 "borigin1_matrix_3_22_DE", 
                                 "borigin1_matrix_3_22_NDE", 
                                 "borigin1_matrix_4_3", 
                                 "borigin1_matrix_8_3", 
                                 "borigin1_matrix_10_2", 
                                 "borigin1_matrix_12_2", 
                                 "borigin1_matrix_14_2", 
                                 "borigin1_matrix_16_2", 
                                 "borigin1_matrix_18_2", 
                                 "borigin1_matrix_20_3", 
                                 "borigin1_matrix_22_3", 
                                 "borigin1_matrix_41_2", 
                                 
                                 "borigin2_matrix_1_2", 
                                 "borigin2_matrix_2_1", 
                                 "borigin2_matrix_2_3", 
                                 "borigin2_matrix_2_10_DE", 
                                 "borigin2_matrix_2_10_NDE", 
                                 "borigin2_matrix_2_12_DE", 
                                 "borigin2_matrix_2_12_NDE", 
                                 "borigin2_matrix_2_14_DE", 
                                 "borigin2_matrix_2_14_NDE", 
                                 "borigin2_matrix_2_16_DE", 
                                 "borigin2_matrix_2_16_NDE", 
                                 "borigin2_matrix_2_18_DE", 
                                 "borigin2_matrix_2_18_NDE", 
                                 "borigin2_matrix_2_41", 
                                 "borigin2_matrix_3_2", 
                                 "borigin2_matrix_3_4", 
                                 "borigin2_matrix_3_8", 
                                 "borigin2_matrix_3_20_DE", 
                                 "borigin2_matrix_3_20_NDE", 
                                 "borigin2_matrix_3_22_DE", 
                                 "borigin2_matrix_3_22_NDE", 
                                 "borigin2_matrix_4_3", 
                                 "borigin2_matrix_8_3", 
                                 "borigin2_matrix_10_2", 
                                 "borigin2_matrix_12_2", 
                                 "borigin2_matrix_14_2", 
                                 "borigin2_matrix_16_2", 
                                 "borigin2_matrix_18_2", 
                                 "borigin2_matrix_20_3", 
                                 "borigin2_matrix_22_3", 
                                 "borigin2_matrix_41_2", 
                                 
                                 "borigin3_matrix_1_2", 
                                 "borigin3_matrix_2_1", 
                                 "borigin3_matrix_2_3", 
                                 "borigin3_matrix_2_10_DE", 
                                 "borigin3_matrix_2_10_NDE", 
                                 "borigin3_matrix_2_12_DE", 
                                 "borigin3_matrix_2_12_NDE", 
                                 "borigin3_matrix_2_14_DE", 
                                 "borigin3_matrix_2_14_NDE", 
                                 "borigin3_matrix_2_16_DE", 
                                 "borigin3_matrix_2_16_NDE", 
                                 "borigin3_matrix_2_18_DE", 
                                 "borigin3_matrix_2_18_NDE", 
                                 "borigin3_matrix_2_41", 
                                 "borigin3_matrix_3_2", 
                                 "borigin3_matrix_3_4", 
                                 "borigin3_matrix_3_8", 
                                 "borigin3_matrix_3_20_DE", 
                                 "borigin3_matrix_3_20_NDE", 
                                 "borigin3_matrix_3_22_DE", 
                                 "borigin3_matrix_3_22_NDE", 
                                 "borigin3_matrix_4_3", 
                                 "borigin3_matrix_8_3", 
                                 "borigin3_matrix_10_2", 
                                 "borigin3_matrix_12_2", 
                                 "borigin3_matrix_14_2", 
                                 "borigin3_matrix_16_2", 
                                 "borigin3_matrix_18_2", 
                                 "borigin3_matrix_20_3", 
                                 "borigin3_matrix_22_3", 
                                 "borigin3_matrix_41_2", 
                                 
                                 "borigin4_matrix_1_2", 
                                 "borigin4_matrix_2_1", 
                                 "borigin4_matrix_2_3", 
                                 "borigin4_matrix_2_10_DE", 
                                 "borigin4_matrix_2_10_NDE", 
                                 "borigin4_matrix_2_12_DE", 
                                 "borigin4_matrix_2_12_NDE", 
                                 "borigin4_matrix_2_14_DE", 
                                 "borigin4_matrix_2_14_NDE", 
                                 "borigin4_matrix_2_16_DE", 
                                 "borigin4_matrix_2_16_NDE", 
                                 "borigin4_matrix_2_18_DE", 
                                 "borigin4_matrix_2_18_NDE", 
                                 "borigin4_matrix_2_41", 
                                 "borigin4_matrix_3_2", 
                                 "borigin4_matrix_3_4", 
                                 "borigin4_matrix_3_8", 
                                 "borigin4_matrix_3_20_DE", 
                                 "borigin4_matrix_3_20_NDE", 
                                 "borigin4_matrix_3_22_DE", 
                                 "borigin4_matrix_3_22_NDE", 
                                 "borigin4_matrix_4_3", 
                                 "borigin4_matrix_8_3", 
                                 "borigin4_matrix_10_2", 
                                 "borigin4_matrix_12_2", 
                                 "borigin4_matrix_14_2", 
                                 "borigin4_matrix_16_2", 
                                 "borigin4_matrix_18_2", 
                                 "borigin4_matrix_20_3", 
                                 "borigin4_matrix_22_3", 
                                 "borigin4_matrix_41_2", 
                                 "borigin5_matrix_1_2", 
                                 "borigin5_matrix_2_1", 
                                 "borigin5_matrix_2_3", 
                                 "borigin5_matrix_2_10_DE", 
                                 "borigin5_matrix_2_10_NDE", 
                                 "borigin5_matrix_2_12_DE", 
                                 "borigin5_matrix_2_12_NDE", 
                                 "borigin5_matrix_2_14_DE", 
                                 "borigin5_matrix_2_14_NDE", 
                                 "borigin5_matrix_2_16_DE", 
                                 "borigin5_matrix_2_16_NDE", 
                                 "borigin5_matrix_2_18_DE", 
                                 "borigin5_matrix_2_18_NDE", 
                                 "borigin5_matrix_2_41", 
                                 "borigin5_matrix_3_2", 
                                 "borigin5_matrix_3_4", 
                                 "borigin5_matrix_3_8", 
                                 "borigin5_matrix_3_20_DE", 
                                 "borigin5_matrix_3_20_NDE", 
                                 "borigin5_matrix_3_22_DE", 
                                 "borigin5_matrix_3_22_NDE", 
                                 "borigin5_matrix_4_3", 
                                 "borigin5_matrix_8_3", 
                                 "borigin5_matrix_10_2", 
                                 "borigin5_matrix_12_2", 
                                 "borigin5_matrix_14_2", 
                                 "borigin5_matrix_16_2", 
                                 "borigin5_matrix_18_2", 
                                 "borigin5_matrix_20_3", 
                                 "borigin5_matrix_22_3", 
                                 "borigin5_matrix_41_2", 
                       "asotin_alpha1", "asotin_alpha2", "deschutes_alpha1", "entiat_alpha1"           ,
                       "fifteenmile_alpha1", "hood_alpha1", "imnaha_alpha1", "john_day_alpha1", "methow_alpha1", "methow_alpha2"           ,
                       "okanogan_alpha1", "tucannon_alpha1", "tucannon_alpha2", "umatilla_alpha1", "umatilla_alpha2", "walla_walla_alpha1",
                       "walla_walla_alpha2", "walla_walla_alpha3", "wenatchee_alpha1", "yakima_alpha1", "asotin_beta", "deschutes_beta"         ,
                       "entiat_beta", "hood_beta", "john_day_beta", "methow_beta"  ,
                       "okanogan_beta", "tucannon_beta", "umatilla_beta", "walla_walla_beta", "wenatchee_beta", "yakima_beta")

upper_columbia_model_params = c("b0_matrix_1_2",
                                "b0_matrix_2_1",
                                "b0_matrix_2_3",
                                "b0_matrix_2_10_DE",
                                "b0_matrix_2_10_NDE",
                                "b0_matrix_2_12_DE",
                                "b0_matrix_2_12_NDE",
                                "b0_matrix_2_14_DE",
                                "b0_matrix_2_14_NDE",
                                "b0_matrix_2_16_DE",
                                "b0_matrix_2_16_NDE",
                                "b0_matrix_2_18_DE",
                                "b0_matrix_2_18_NDE",
                                "b0_matrix_2_41",
                                "b0_matrix_3_2",
                                "b0_matrix_3_4",
                                "b0_matrix_3_8",
                                "b0_matrix_3_20_DE",
                                "b0_matrix_3_20_NDE",
                                "b0_matrix_3_22_DE",
                                "b0_matrix_3_22_NDE",
                                "b0_matrix_4_3",
                                "b0_matrix_4_5",
                                "b0_matrix_5_4",
                                "b0_matrix_5_6",
                                "b0_matrix_5_24_DE",
                                "b0_matrix_5_24_NDE",
                                "b0_matrix_6_5",
                                "b0_matrix_6_7",
                                "b0_matrix_6_26_DE",
                                "b0_matrix_6_26_NDE",
                                "b0_matrix_7_6",
                                "b0_matrix_7_28_DE",
                                "b0_matrix_7_28_NDE",
                                "b0_matrix_7_30_DE",
                                "b0_matrix_7_30_NDE",
                                "b0_matrix_7_42",
                                "b0_matrix_8_3",
                                "b0_matrix_8_9",
                                "b0_matrix_8_32_DE",
                                "b0_matrix_8_32_NDE",
                                "b0_matrix_9_8",
                                "b0_matrix_9_34_DE",
                                "b0_matrix_9_34_NDE",
                                "b0_matrix_9_36",
                                "b0_matrix_9_37",
                                "b0_matrix_9_38",
                                "b0_matrix_9_39_DE",
                                "b0_matrix_9_39_NDE",
                                "b0_matrix_10_2",
                                "b0_matrix_12_2",
                                "b0_matrix_14_2",
                                "b0_matrix_16_2",
                                "b0_matrix_18_2",
                                "b0_matrix_20_3",
                                "b0_matrix_22_3",
                                "b0_matrix_24_5",
                                "b0_matrix_26_6",
                                "b0_matrix_28_7",
                                "b0_matrix_30_7",
                                "b0_matrix_32_8",
                                "b0_matrix_34_9",
                                "b0_matrix_36_9",
                                "b0_matrix_37_9",
                                "b0_matrix_38_9",
                                "b0_matrix_39_9",
                                "b0_matrix_41_2",
                                "b0_matrix_42_7",
                                
                                "borigin1_matrix_3_4",
                                "borigin1_matrix_4_3",
                                "borigin1_matrix_4_5",
                                "borigin1_matrix_5_4",
                                "borigin1_matrix_5_6",
                                "borigin1_matrix_5_24_DE",
                                "borigin1_matrix_5_24_NDE",
                                "borigin1_matrix_6_5",
                                "borigin1_matrix_6_7",
                                "borigin1_matrix_6_26_DE",
                                "borigin1_matrix_6_26_NDE",
                                "borigin1_matrix_7_6",
                                "borigin1_matrix_7_28_DE",
                                "borigin1_matrix_7_28_NDE",
                                "borigin1_matrix_7_30_DE",
                                "borigin1_matrix_7_30_NDE",
                                "borigin1_matrix_7_42",
                                "borigin1_matrix_24_5",
                                "borigin1_matrix_26_6",
                                "borigin1_matrix_28_7",
                                "borigin1_matrix_30_7",
                                "borigin1_matrix_42_7",
                                
                                "borigin2_matrix_3_4",
                                "borigin2_matrix_4_3",
                                "borigin2_matrix_4_5",
                                "borigin2_matrix_5_4",
                                "borigin2_matrix_5_6",
                                "borigin2_matrix_5_24_DE",
                                "borigin2_matrix_5_24_NDE",
                                "borigin2_matrix_6_5",
                                "borigin2_matrix_6_7",
                                "borigin2_matrix_6_26_DE",
                                "borigin2_matrix_6_26_NDE",
                                "borigin2_matrix_7_6",
                                "borigin2_matrix_7_28_DE",
                                "borigin2_matrix_7_28_NDE",
                                "borigin2_matrix_7_30_DE",
                                "borigin2_matrix_7_30_NDE",
                                "borigin2_matrix_7_42",
                                "borigin2_matrix_24_5",
                                "borigin2_matrix_26_6",
                                "borigin2_matrix_28_7",
                                "borigin2_matrix_30_7",
                                "borigin2_matrix_42_7",
                                
                                "borigin3_matrix_3_4",
                                "borigin3_matrix_4_3",
                                "borigin3_matrix_4_5",
                                "borigin3_matrix_5_4",
                                "borigin3_matrix_5_6",
                                "borigin3_matrix_5_24_DE",
                                "borigin3_matrix_5_24_NDE",
                                "borigin3_matrix_6_5",
                                "borigin3_matrix_6_7",
                                "borigin3_matrix_6_26_DE",
                                "borigin3_matrix_6_26_NDE",
                                "borigin3_matrix_7_6",
                                "borigin3_matrix_7_28_DE",
                                "borigin3_matrix_7_28_NDE",
                                "borigin3_matrix_7_30_DE",
                                "borigin3_matrix_7_30_NDE",
                                "borigin3_matrix_7_42",
                                "borigin3_matrix_24_5",
                                "borigin3_matrix_26_6",
                                "borigin3_matrix_28_7",
                                "borigin3_matrix_30_7",
                                "borigin3_matrix_42_7", 
                                 "asotin_alpha1", "asotin_alpha2", "deschutes_alpha1", "entiat_alpha1"           ,
                                 "fifteenmile_alpha1", "hood_alpha1", "imnaha_alpha1", "john_day_alpha1", "methow_alpha1", "methow_alpha2"           ,
                                 "okanogan_alpha1", "tucannon_alpha1", "tucannon_alpha2", "umatilla_alpha1", "umatilla_alpha2", "walla_walla_alpha1",
                                 "walla_walla_alpha2", "walla_walla_alpha3", "wenatchee_alpha1", "yakima_alpha1", "asotin_beta", "deschutes_beta"         ,
                                 "entiat_beta", "hood_beta", "john_day_beta", "methow_beta"  ,
                                 "okanogan_beta", "tucannon_beta", "umatilla_beta", "walla_walla_beta", "wenatchee_beta", "yakima_beta")




# reduce number of parameters for each of the three ESU models
snake_fit_shinystan <- as.shinystan(snake_fit, pars = snake_model_params)
middle_columbia_fit_shinystan <- as.shinystan(middle_columbia_fit, pars = middle_columbia_model_params)
upper_columbia_fit_shinystan <- as.shinystan(upper_columbia_fit, pars = upper_columbia_model_params)

launch_shinystan(snake_fit_shinystan)
launch_shinystan(middle_columbia_fit_shinystan)
launch_shinystan(upper_columbia_fit_shinystan)
