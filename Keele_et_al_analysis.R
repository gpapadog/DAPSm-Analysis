# ImputeTS for our data

config_path <- '~/Github/DAPSm-Analysis/config.R'

# Loading libraries.
library(data.table)
library(corrplot)
library(fields)
library(ggplot2)
library(DAPSm)
library(stringr)
library(arepa)
library(optmatch)
library(Rcplex)


source(config_path)
# Setting the working directory.
setwd(wd)

# Sourcing the functions that I will be using.
source('CreateNOxControlsFunction.R')
source('Data_analysis_functions.R')
source('Data_analysis_models_functions.R')
source('expit.R')
source('make_data_functions.R')
source('StandDiff_function.R')
source('CleanPPunits_function.R')
source('AverageSulfurContent_function.R')
source('PredictHeatInput_function.R')
source('keele_match_function.R')
source('Keele_et_al_functions/01_subsetmatch2.R')
source('Keele_et_al_functions/02_errorhandling.R')
source('Keele_et_al_functions/03_problemparameters2.R')
source('Keele_et_al_functions/04_constraintmatrix2.R')


# ------------------- PART 1------------------- #
# ------ CREATING THE ANALYSIS DATA SET ------ #


# ---- STEP 1. Loading unit level data and predicting heat input.
full_data <- LoadUnitLevelData(data_dir)
full_data[, V1 := NULL]
dat_unit <- CleanPPunits(full_data, year = year, month = month)
dat_unit <- AverageSulfurContent(dat_unit)
setkeyv(dat_unit, c('uID', 'Year', 'Month'))

# Imputing heat input at the unit level.
lm_pred <- NewPredictHeatInput(dat_unit, year, month, time_use)
print(lm_pred$total)
print(lm_pred$missing)
print(lm_pred$rsquared)
print(lm_pred$num_pred)
subdta_ym <- lm_pred$data
# Negative predicted values set to 0.
wh <- which(subdta_ym$Heat.Input..MMBtu. < 0)
print(paste('Setting', length(wh), 'negative heat input entries to 0.'))
subdta_ym$Heat.Input..MMBtu.[wh] <- 0


# ---- STEP 2: Aggregating to the facility level.
dat_facility <- UnitToFacility(dat_unit = subdta_ym)

print(paste('Dropping', sum(dat_facility$totHeatInput == 0, na.rm = TRUE),
            'facilities for heat input = 0'))
print(paste('Dropping', sum(is.na(dat_facility$totHeatInput)),
            'facilities for missing heat input'))
dat_facility <- subset(dat_facility, totHeatInput > 0)

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


# ------------------- PART 2------------------- #
# --------- SETTING UP THE ANALYSIS --------- #
outcome_names <- c('mean4maxOzone', 'totNOxemissions')
out.col <- which(names(analysis_dat) == outcome_analysis)
trt.col <- which(names(analysis_dat) == 'SnCR')
coord.cols <- which(names(analysis_dat) %in% c('Fac.Longitude', 'Fac.Latitude'))
drop_col <- which(names(analysis_dat) %in% setdiff(outcome_names, outcome_analysis))
conf <- setdiff(1:ncol(analysis_dat), c(out.col, trt.col, coord.cols, drop_col))
subdta <- analysis_dat[, unique(c(trt.col, out.col, coord.cols, conf)), with = FALSE]


n_matches = 1  # Number of matched per treated.
use_controls = NULL  # Whether specific controls need to be used
enforce_constraints <- FALSE
trt.col <- 1
out.col <- 2
coord.cols <- c(4, 3)

naive_ps <- glm(as.formula(paste('SnCR ~ . - Fac.Latitude - Fac.Longitude -',
                                 outcome_analysis)), data = subdta, family = 'binomial')
subdta[, prop.scores := fitted(naive_ps)]
ignore.cols <- which(names(subdta) == 'prop.scores')
ignore.cols.coords <- c(ignore.cols, coord.cols)
cols.balance <- setdiff(1:dim(subdta)[2],
                        c(trt.col, out.col, ignore.cols.coords,
                          which(names(subdta) == 'prop.scores')))

# ------------------- PART 3 ------------------- #
# --------- ANALYSIS WITH MATCHING ON PS --------- #

mom_covs <- NULL
mom_tols <- NULL
exact_breaks <- seq(min(subdta$prop.scores), max(subdta$prop.scores), length.out = 6)
exact_covs <- matrix(as.numeric(cut(subdta$prop.scores, exact_breaks)), ncol = 1)


# Fitting Keele.
subset_weight <- seq(10, 1500, by = 20)
num_match <- numeric(length(subset_weight))
bal <- matrix(NA, nrow = length(subset_weight), ncol = length(cols.balance))
distance <- numeric(length(subset_weight))

for (ss in 1:length(subset_weight)) {
  keele <- keele_match(subdta, trt_col = trt.col, out_col = out.col,
                       coords.columns = coord.cols, subset_weight = subset_weight[ss],
                       n_matches = n_matches, use_controls = use_controls,
                       enforce_constraints = enforce_constraints, pairsRet = TRUE,
                       cols.balance = cols.balance, mom_tols = mom_tols,
                       mom_covs = mom_covs, exact_covs = exact_covs)
  num_match[ss] <- keele$num_match
  distance[ss] <- mean(keele$distance)
  bal[ss, ] <- keele$balance[2, ]
}

plot(subset_weight, num_match)
plot(subset_weight, distance)
plot(subset_weight, apply(bal, 1, function(x) sum(abs(x) > cutoff)))
plot(subset_weight, apply(bal, 1, function(x) mean(abs(x))))
plot(subset_weight, apply(bal, 1, function(x) max(abs(x))))



# ------------------- PART 3 ------------------- #
# --------- ANALYSIS WITH MATCHING ON COVARIATES --------- #

mom_covs_ind <- c(2, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
exact_covs_ind <- c(6, 20, 21, 22)
mom_covs <- as.matrix(subdta)[, mom_covs_ind]
mom_tols <- caliper * apply(mom_covs, 2, sd)
exact_covs <- as.matrix(subdta)[, exact_covs_ind]


# Fitting Keele.
subset_weight <- seq(900, 1200, by = 10)
num_match <- numeric(length(subset_weight))
bal <- matrix(NA, nrow = length(subset_weight), ncol = length(cols.balance))
distance <- numeric(length(subset_weight))

for (ss in 1:length(subset_weight)) {
  if (ss %% 20 == 0) {
    print(ss)
  }
  keele <- keele_match(subdta, trt_col = trt.col, out_col = out.col,
                       coords.columns = coord.cols, subset_weight = subset_weight[ss],
                       n_matches = n_matches, use_controls = use_controls,
                       enforce_constraints = enforce_constraints, pairsRet = TRUE,
                       cols.balance = cols.balance, mom_tols = mom_tols,
                       mom_covs = mom_covs, exact_covs = exact_covs)
  num_match[ss] <- keele$num_match
  distance[ss] <- mean(keele$distance)
  bal[ss, ] <- keele$balance[2, ]
}

plot(subset_weight, num_match)
plot(subset_weight, distance)
plot(subset_weight, apply(bal, 1, function(x) sum(abs(x) > cutoff)))
plot(subset_weight, apply(bal, 1, function(x) mean(abs(x))))
plot(subset_weight, apply(bal, 1, function(x) max(abs(x))))




