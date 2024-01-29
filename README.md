# steelhead
Steelhead overshoot and fallback in the Columbia and Snake River Basins

*last updated: 2024-01-29*

Note: The data files are very large and therefore are not uploaded to GitHub. If you would like access to the data files, please contact me directly.

To reproduce this analysis, you must run the following scripts **sequentially**:

**Note**: Once the model has been finalized, these scripts will be re-organized such that all of the final versions of each script (i.e., the ones needed to reproduce the results) are together in a single folder. 

## Part 1: Prepare the fish movement data for inclusion in the model

#### Step 1: Query from PTAGIS

The code to setup/explain the queries is can be found in `R/01_data_querying.Rmd`.

All of the PTAGIS files that were queried can be found in `PTAGIS_queries/complete_tag_histories/antenna_info/`

#### Step 2: Convert raw detections into detections at specific sites
This also notes whether detections at fish ladders are in fact, non-ascents. Diagrams of the various site configurations can be found in the `site_configuration` folder.

For computational reasons, the full dataset was split into four, with the same script (except for the files imported and exported) used on each of the four subsets.

* `from_hyak_transfer/2022-07-27_det_hist/CTH1-4/CTH1-4_02_hyak_detection_histories_v2.R`
* `from_hyak_transfer/2022-07-27_det_hist/CTH5-8/CTH5-8_02_hyak_detection_histories_v2.R`
* `from_hyak_transfer/2022-07-27_det_hist/CTH9-11/CTH9-11_02_hyak_detection_histories_v2.R`
* `from_hyak_transfer/2022-07-27_det_hist/CTH12-14/CTH12-14_02_hyak_detection_histories_v2.R`

Once the above four scripts were run, this script was run to join them and export a single file that could then be fed into the next step:
`from_hyak_transfer/2022-07-27_det_hist/02.5_join_det_hist.R`. Note that this script also makes a few fixes to the files, and then exports them here: `from_hyak_transfer/2022-07-27_det_hist/complete_det_hist_postprocessed.csv`.

#### Step 3: Convert detections at specific sites into state occupancy/transitions
This script takes the `complete_det_hist_postprocessed.csv` file that was exported during step 2 and converts this series of detections at different sites into a history of transitions between different states in our model.

For computational reasons, the full dataset was once again split into four, with the same script (except for the files imported and exported) used on each of the four subsets.

* `from_hyak_transfer/2022-11-02-complete_det_hist_ckpt/03_hyak_complete_detection_histories_v4_part1.R`
* `from_hyak_transfer/2022-11-02-complete_det_hist_ckpt/03_hyak_complete_detection_histories_v4_part2.R`
* `from_hyak_transfer/2022-11-02-complete_det_hist_ckpt/03_hyak_complete_detection_histories_v4_part3.R`
* `from_hyak_transfer/2022-11-02-complete_det_hist_ckpt/03_hyak_complete_detection_histories_v4_part4.R`

Each of these files exports one CSV (`states_complete_part1.csv`, `states_complete_part2.csv`, `states_complete_part3.csv`, and `states_complete_part4.csv`) of the state occupancy history of each fish as defined by our model states.

#### Step 4: Clean up the state occupancy

The script `R/03.5_stepwise_states_cleaning.R` takes as inputs the four CSVs from step 3 and cleans them up to remove any juvenile movements, identify kelt movements, and identify repeat spawners. This script then exports the full state histories for all fish in `stan_actual/adults_states_complete.csv`, and also exports the state histories for the different DPSs (Snake, Upper Columbia, and Middle Columbia) to the `stan_actual/` folder in a few of the sub-folders for models (eventually, this will be moved to a different location).

## Part 2: Prepare the covariate data for inclusion in the model

### Temperature

#### Step 1: Filter out outlier temperatures
The script `R/22_temperature_for_modeling.Rmd` filters out outlier temperatures manually (based on a visual inspection of temperature data), based on interannual averages (data points >4 degrees from the average value for that day of year are removed), and based on a moving average (data points >2 degrees from the 7-day moving average are removed).

#### Step 2: Fit MARSS model to estimate temperature
The script `R/23_temperature_MAR_model.R` uses as inputs the forebay and tailrace temperatures from the previous step and uses them to fit a MARSS model to estimate temperature (this addresses the missing data). 

#### Step 3: Estimate temperatures within windows
The script `R/24_model_based_temperature_windows.Rmd` uses as inputs the MARSS-estimated temperatures from the previous steps, calculates the median residence time for in each state, and then uses that median residence time to calculate the mean temperature using a window for each date (defined as the window of time from the date plus the median residence time in that state). These window temperatures are exported in `/covariate_data/for_model/window_temps_for_stan.csv`.

### Spill

The script `R/25_spill_data_prep.Rmd` takes the daily average spill at each dam and processes it in two ways: by calculating a window of the volume of spill around a date (the same approach as is used for temperature) and by calculating the total number of days of spill per year in the winter months (January, February, and March). The spill window data is exported in `/covariate_data/for_model/window_spill_for_stan.csv`. and the spill days data is exported in `/covariate_data/for_model/january_spill_df_for_stan.csv`, `/covariate_data/for_model/february_spill_df_for_stan.csv`, and `/covariate_data/for_model/march_spill_df_for_stan.csv`. Note that although the spill days data is exported separately for each of the winter months, it is later combined into a single value for the total number of days with any amount of spill across all months.

## Part 3: Fit the model 


