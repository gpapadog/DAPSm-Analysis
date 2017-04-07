# ImputeTS for our data

config_path <- '~/Github/DAPSm-Analysis/config.R'

# Loading libraries.
library(data.table)
library(corrplot)
library(MatchIt)
library(Matching)
library(fields)
library(ggplot2)
library(DAPSm)
library(stringr)
library(arepa)
library(optmatch)
library(imputeTS)

source(config_path)
# Setting the working directory.
setwd(wd)

# Sourcing the functions that I will be using.
source('CaliperEst_function.R')
source('CreateNOxControlsFunction.R')
source('Data_analysis_functions.R')
source('Data_analysis_models_functions.R')
source('dist.caliper_function.R')
source('expit.R')
source('GBMPropScores_function.R')
source('make_data_functions.R')
source('OptPSmatch_function.R')
source('predict_variable_functions.R')
source('PredictHeatInput.R')
source('PSmatchEst_function.R')
source('StandDiff_function.R')
source('CleanPPunits_function.R')
source('AverageSulfurContent_function.R')
source('ImputeHeatInput_function.R')

# ------------------- PART 1------------------- #
# ------ CREATING THE ANALYSIS DATA SET ------ #


# ---- STEP 1. Loading unit level data and predicting heat input.
full_data <- LoadUnitLevelData(data_dir)
full_data[, V1 := NULL]
dat_unit <- CleanPPunits(full_data, year = year, month = month)
dat_unit <- AverageSulfurContent(dat_unit)
setkeyv(dat_unit, c('uID', 'Year', 'Month'))

if (impute_with_ts) {
  subdta <- ImputeHeatInput(dat_unit, year, month, method = 'kalman')
  subdta_ym <- subset(subdta, Year == year & Month %in% month)
} else {
  lm_pred <- NewPredictHeatInput(dat_unit, year, month, time_use)
  print(lm_pred$total)
  print(lm_pred$missing)
  print(lm_pred$rsquared)
  print(lm_pred$num_pred)
  subdta_ym <- lm_pred$data
  wh <- which(subdta_ym$Heat.Input..MMBtu. < 0)
  print(paste('Setting', length(wh), 'negative heat input entries to 0.'))
  subdta_ym$Heat.Input..MMBtu.[wh] <- 0
}

# ---- STEP 2: Aggregating to the facility level.
dat_facility <- UnitToFacility(dat_unit = subdta_ym)

print(paste('Dropping', sum(dat_facility$totHeatInput == 0, na.rm = TRUE),
            'facilities for heat input = 0'))
dat_facility <- subset(dat_facility, totHeatInput > 0 | is.na(totHeatInput))

# ---- STEP 3: Linking the aggregated data to ozone monitors.
dat <- LinkPPtoMonitors(dat_facility, within_km, year, month, OzTempCen = data_dir)

# Dropping facilities with missing data for at least one month.
wh <- which(dat$nmonths != length(month))
print(paste('Dropping', length(wh), 'out of', length(unique(dat$FacID)),
            'facilities due to missing information on at least one month'))
if (length(wh) > 0) {
  dat <- dat[- wh, ]
}

analysis_dat <- CleanData(dat, plotcor = FALSE)
analysis_dat <- ReformData(analysis_dat)


