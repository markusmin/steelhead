# steelhead
Steelhead overshoot and fallback in the Columbia and Snake River Basins

Note: The data files are very large and therefore are not uploaded to GitHub. If you would like access to the data files, please contact me directly.


## Overview

To reproduce this analysis, you must run the following scripts **sequentially**:

**Note**: Once the model has been finalized, these scripts will be re-organized such that all of the final versions of each script (i.e., the ones needed to reproduce the results) are together in a single folder. 

#### Step 1: Query from PTAGIS


#### Step 2: Convert raw detections into detections at specific sites
This also notes whether detections at fish ladders are in fact, non-ascents.

For computational reasons, the full dataset was split into four, with the same script (except for the files imported and exported) used on each of the four subsets.



Once the above four scripts were run, this script was run to join them and export a single file that could then be fed into the next step:
`from_hyak_transfer/2022-07-27_det_hist/02.5_join_det_hist.R`.





