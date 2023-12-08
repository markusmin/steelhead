# steelhead
Steelhead overshoot and fallback in the Columbia and Snake River Basins

Note: The data files are very large and therefore are not uploaded to GitHub. If you would like access to the data files, please contact me directly.


## Overview

To reproduce this analysis, you must run the following scripts **sequentially**:

**Note**: Once the model has been finalized, these scripts will be re-organized such that all of the final versions of each script (i.e., the ones needed to reproduce the results) are together in a single folder. 

#### Step 1: Query from PTAGIS

The code to setup/explain the queries is can be found in `R/01_data_querying.Rmd`.

All of the PTAGIS files that were queried can be found in `PTAGIS_queries/complete_tag_histories/antenna_info/`

#### Step 2: Convert raw detections into detections at specific sites
This also notes whether detections at fish ladders are in fact, non-ascents. Diagrams of the various site configurations can be found in the `site_configuration` folder.

For computational reasons, the full dataset was split into four, with the same script (except for the files imported and exported) used on each of the four subsets.

`from_hyak_transfer/2022-07-27_det_hist/CTH1-4/CTH1-4_02_hyak_detection_histories_v2.R`
`from_hyak_transfer/2022-07-27_det_hist/CTH5-8/CTH5-8_02_hyak_detection_histories_v2.R`
`from_hyak_transfer/2022-07-27_det_hist/CTH9-11/CTH9-11_02_hyak_detection_histories_v2.R`
`from_hyak_transfer/2022-07-27_det_hist/CTH12-14/CTH12-14_02_hyak_detection_histories_v2.R`

Once the above four scripts were run, this script was run to join them and export a single file that could then be fed into the next step:
`from_hyak_transfer/2022-07-27_det_hist/02.5_join_det_hist.R`. Note that this script also makes a few fixes to the files, and then exports them here: `from_hyak_transfer/2022-07-27_det_hist/complete_det_hist_postprocessed.csv`.

#### Step 3: Convert detections at specific sites into state occupancy/transitions
This script takes the `complete_det_hist_postprocessed.csv` file that was exported during step 2 and converts this series of detections at different sites into a history of transitions between different states in our model.

For computational reasons, the full dataset was once again split into four, with the same script (except for the files imported and exported) used on each of the four subsets.

`from_hyak_transfer/2022-11-02-complete_det_hist_ckpt/03_hyak_complete_detection_histories_v4_part1.R`
`from_hyak_transfer/2022-11-02-complete_det_hist_ckpt/03_hyak_complete_detection_histories_v4_part2.R`
`from_hyak_transfer/2022-11-02-complete_det_hist_ckpt/03_hyak_complete_detection_histories_v4_part3.R`
`from_hyak_transfer/2022-11-02-complete_det_hist_ckpt/03_hyak_complete_detection_histories_v4_part4.R`

Each of these files exports one CSV (`states_complete_part1.csv`, `states_complete_part2.csv`, `states_complete_part3.csv`, and `states_complete_part4.csv`) of the state occupancy history of each fish as defined by our model states.

#### Step 4: Clean up the state occupancy

The script `R/03.5_stepwise_states_cleaning.R` takes as inputs the four CSVs from step 3 and cleans them up to remove any juvenile movements, identify kelt movements, and identify repeat spawners. This script then exports the full state histories for all fish in `stan_actual/adults_states_complete.csv`, and also exports the state histories for the different DPSs (Snake, Upper Columbia, and Middle Columbia) to the `stan_actual/` folder in a few of the sub-folders for models (eventually, this will be moved to a different location).


