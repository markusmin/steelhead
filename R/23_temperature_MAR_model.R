### Temperature MAR model

# This script will fit a MAR model to temperature data from different dams
# to estimate temperature for every date in our time series.

# The MAR model that we will fit is a state-space model, with forebay and tailrace
# temperatures as observations of a "Columbia River Basin" temperature, with
# each dam estimated as an offset of it
# We will have sixteen observations of the same process (forebay and tailrace temperatures for eight dams)

# Load libraries
library(tidyverse)
library(MARSS)

# Load data
# note that these data files were produced by script 22, which filtered out bad data
# in addition to some other steps. The only columns that we need from these
# files are date, tailrace_temp, and forebay_temp
BON_temp <- read.csv(here::here("covariate_data", "temp_processed", "BON_temp_processed.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
MCN_temp <- read.csv(here::here("covariate_data", "temp_processed", "MCN_temp_processed.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
PRA_temp <- read.csv(here::here("covariate_data", "temp_processed", "PRA_temp_processed.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
RIS_temp <- read.csv(here::here("covariate_data", "temp_processed", "RIS_temp_processed.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
RRE_temp <- read.csv(here::here("covariate_data", "temp_processed", "RRE_temp_processed.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
ICH_temp <- read.csv(here::here("covariate_data", "temp_processed", "ICH_temp_processed.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
WEL_temp <- read.csv(here::here("covariate_data", "temp_processed", "WEL_temp_processed.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)
LGR_temp <- read.csv(here::here("covariate_data", "temp_processed", "LGR_temp_processed.csv")) %>% dplyr::select(date, tailrace_temp, forebay_temp)

# make the column data names more informative by including dam
colnames(BON_temp) <- paste0("BON_", colnames(BON_temp))
colnames(MCN_temp) <- paste0("MCN_", colnames(MCN_temp))
colnames(PRA_temp) <- paste0("PRA_", colnames(PRA_temp))
colnames(RIS_temp) <- paste0("RIS_", colnames(RIS_temp))
colnames(RRE_temp) <- paste0("RRE_", colnames(RRE_temp))
colnames(ICH_temp) <- paste0("ICH_", colnames(ICH_temp))
colnames(WEL_temp) <- paste0("WEL_", colnames(WEL_temp))
colnames(LGR_temp) <- paste0("LGR_", colnames(LGR_temp))


# reformat data for model - drop date and transpose, bind all dams together
t(BON_temp[,2:3]) %>% 
  rbind(., t(MCN_temp[,2:3])) %>% 
  rbind(., t(PRA_temp[,2:3])) %>% 
  rbind(., t(RIS_temp[,2:3])) %>% 
  rbind(., t(RRE_temp[,2:3])) %>% 
  rbind(., t(ICH_temp[,2:3])) %>% 
  rbind(., t(WEL_temp[,2:3])) %>% 
  rbind(., t(LGR_temp[,2:3])) -> temp_for_MAR

dat = temp_for_MAR
# fit the model


# A are the offsets/biases for each dam - eight elements, each repeated twice
# BON is the reference
# a <- matrix(data = c("a_BON", "a_BON", "a_MCN", "a_MCN", "a_PRA", "a_PRA", "a_RIS", "a_RIS",
#        "a_RRE", "a_RRE", "a_ICH", "a_ICH", "a_WEL", "a_WEL", "a_LGR", "a_LGR"),
#        nrow = 16, ncol = 1)
a <- matrix(data = as.list(rep(0, 16)),
            nrow = 16, ncol = 1)
a[3:16] <- c("a_MCN", "a_MCN", "a_PRA", "a_PRA", "a_RIS", "a_RIS",
       "a_RRE", "a_RRE", "a_ICH", "a_ICH", "a_WEL", "a_WEL", "a_LGR", "a_LGR")
A <- a

# no bias term for Columbia river temp
U <- "zero"

# just one process and no interactions
B <- matrix(1)

# Z is sixteen observations of the same process
Z <- matrix(1, 16, 1)

# we need eight elements to correspond to the 
# eight dams, and each element will be represented twice (once for forebay and once for tailrace)
r <- matrix(data = as.list(rep(0, 256)), nrow = 16, ncol = 16)
diag(r) <- c("r_BON", "r_BON", "r_MCN", "r_MCN", "r_PRA", "r_PRA", "r_RIS", "r_RIS",
             "r_RRE", "r_RRE", "r_ICH", "r_ICH", "r_WEL", "r_WEL", "r_LGR", "r_LGR")
# R <- r
# variance among dams is probably the same, given all the same equipment
R <- "diagonal and equal"
Q <- matrix("q") # this matrix is a 1x1 matrix, because there's only one trend we're estimating
model.list <- list(B = B, U = U, Q = Q, Z = Z, A = A, R = R)
kem <- MARSS(dat, model = model.list)

# plot MARSS model fit
fit_df <- data.frame(date = seq(ymd('2005-06-01'),ymd('2022-12-31'), by = 'days'),
                     fit = kem$states)


# some have considerably larger r (variance) values. Why?
# is it because of the difference between forebay and tailrace temperatures?
BON_temp[,2] - BON_temp[,3] -> BON_fore_tail_diff
summary(BON_fore_tail_diff)
MCN_temp[,2] - MCN_temp[,3] -> MCN_fore_tail_diff
summary(MCN_fore_tail_diff)
PRA_temp[,2] - PRA_temp[,3] -> PRA_fore_tail_diff
summary(PRA_fore_tail_diff)
RIS_temp[,2] - RIS_temp[,3] -> RIS_fore_tail_diff
summary(RIS_fore_tail_diff)
RRE_temp[,2] - RRE_temp[,3] -> RRE_fore_tail_diff
summary(RRE_fore_tail_diff)
ICH_temp[,2] - ICH_temp[,3] -> ICH_fore_tail_diff
summary(ICH_fore_tail_diff)
WEL_temp[,2] - WEL_temp[,3] -> WEL_fore_tail_diff
summary(WEL_fore_tail_diff)
LGR_temp[,2] - LGR_temp[,3] -> LGR_fore_tail_diff
summary(LGR_fore_tail_diff)

# plot Columbia River temperature
dates <- seq(ymd('2005-06-01'),ymd('2022-12-31'), by = 'days')

par(mfrow = c(1,1))


plot(dates, kem$states, ylab = "Columbia River Temperature", 
     xlab = "", type = "l")
lines(dates, kem$states - 1.96 * kem$states.se, type = "l", lwd = 1, lty = 2, col = "red")
lines(dates, kem$states + 1.96 * kem$states.se, type = "l", lwd = 1, lty = 2, col = "red")
title("Columbia River Temperature")


# extract the y hat values
hat_yt <- MARSShatyt(kem)

# plot the yhat for each dam
# add points for actual data

par(mfrow = c(2,2))
for (i in 1:16) {
  plot(dates, kem$ytT[i, ], ylab = "Temperature", 
       xlab = "", type = "l")
  lines(dates, kem$ytT[i, ] - 1.96 * kem$ytT.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  lines(dates, kem$ytT[i, ] + 1.96 * kem$ytT.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  # points(dates, dat[i, ], cex = 0.2)
  title(rownames(dat)[i])
}

# plot just the data
par(mfrow = c(1,1))
for (i in 1:16) {
  plot(dates, dat[i, ], cex = 0.4)
  title(rownames(dat)[i])
}

# export the yhat values for tailrace temperatures for modeling
temp_yhat <- data.frame(date = dates,
                        BON = kem$ytT[1,],
                        MCN = kem$ytT[3,],
                        PRA = kem$ytT[5,],
                        RIS = kem$ytT[7,],
                        RRE = kem$ytT[9,],
                        ICH = kem$ytT[11,],
                        WEL = kem$ytT[13,],
                        LGR = kem$ytT[15,])


# extract A values and add to the state
a_est <- as.data.frame(coef(kem)$A)
a_MCN <- a_est["a_MCN",]
a_PRA <- a_est["a_PRA",]
a_RIS <- a_est["a_RIS",]
a_RRE <- a_est["a_RRE",]
a_ICH <- a_est["a_ICH",]
a_WEL <- a_est["a_WEL",]
a_LGR <- a_est["a_LGR",]


# BON estimate is just the state
BON_temp_model_est <- kem$states

# every other is just plus the bias term
MCN_temp_model_est <- kem$states + a_MCN
PRA_temp_model_est <- kem$states + a_PRA
RIS_temp_model_est <- kem$states + a_RIS
RRE_temp_model_est <- kem$states + a_RRE
ICH_temp_model_est <- kem$states + a_ICH
WEL_temp_model_est <- kem$states + a_WEL
LGR_temp_model_est <- kem$states + a_LGR

temp_model_est <- data.frame(date = dates,
                             MOUTH = rep(NA, 6423),
                             BON = BON_temp_model_est[1,],
                             MCN = MCN_temp_model_est[1,],
                             PRA = PRA_temp_model_est[1,],
                             RIS = RIS_temp_model_est[1,],
                             RRE = RRE_temp_model_est[1,],
                             ICH = ICH_temp_model_est[1,],
                             WEL = WEL_temp_model_est[1,],
                             LGR = LGR_temp_model_est[1,])


write.csv(temp_model_est, here::here("covariate_data", "for_model", "temp", "temp_mod_est.csv"), row.names = FALSE)






