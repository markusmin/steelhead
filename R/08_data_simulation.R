# 08: Data simulation

# Here we are generating a simulated dataset, in order to test our model

# Load libraries
library(here)
library(janitor)
library(tidyverse)
library(lubridate)

# We will be using a reduced number of states:

# 5 mainstem sites:
# 1: Mouth to BON
# 2: BON to MCN
# 3: MCN to ICH or PRA
# 4: PRA to RIS
# 5: ICH to LGR

# 4 tributary sites:
# 1: John Day River
# 2: Deschutes River
# 3: Yakima River
# 4: Tucannon River

# We will have 3 natal origins:
# John Day River, Yakima River, Tucannon River


# We will have two continuous covariates:
# Temperature, flow

# We will have two categorical covariates:
# Natal origin, rear type


# Our movement probabilities will follow the following functional relationships:

# For the mouth to BON state
# ascend BON: 


# All other relationships will not be a product of covariates














