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



