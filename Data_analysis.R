# Author: Georgia Papadogeorgou
# Date: 5/15/2016
# Desc: Using the facility level data, we fit the models and perform the analysis
#       using the naive, and the spatial propensity score matching methods.

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


# ------------------- PART 1------------------- #
# ------ CREATING THE ANALYSIS DATA SET ------ #


# ---- STEP 1. Loading unit level data and predicting heat input.
dat_unit <- LoadUnitLevelData(data_dir)
subdta <- PredictHeatInput(dat_unit, year, month, time_use)

# ---- STEP 2: Aggregate to the facility level.
dat_facility <- UnitToFacility(dat_unit = subdta)

# Dropping facilities that did not operate.
print(paste('Dropping', sum(dat_facility$totHeatInput == 0, na.rm = TRUE),
            'facilities for heat input = 0'))
dat_facility <- subset(dat_facility, totHeatInput > 0 | is.na(totHeatInput))

# ---- STEP 3: Linking the aggregated data to ozone monitors.
# To run the next command you need to run the script that links ozone, temperature
# and Census information.
dat <- LinkPPtoMonitors(dat_facility, within_km, year, month, OzTempCen = data_dir)

# Dropping facilities with missing data for at least one month.
wh <- which(dat$nmonths != length(month))
print(paste('Dropping', length(wh), 'out of', length(unique(dat$FacID)),
            'facilities due to missing information on at least one month'))
dat <- dat[- wh, ]

analysis_dat <- CleanData(dat, plotcor = FALSE)
analysis_dat <- ReformData(analysis_dat)
# analysis_dat includes data on the observations we will use.
analysis_dat[, meanOzone := NULL]
analysis_dat[, meanmaxOzone := NULL]


# ------------------- PART 2------------------- #
# --------- SETTING UP THE ANALYSIS --------- #
outcome_names <- c('mean4maxOzone', 'totNOxemissions')
out.col <- which(names(analysis_dat) == outcome_analysis)
trt.col <- which(names(analysis_dat) == 'SnCR')
coord.cols <- which(names(analysis_dat) %in% c('Fac.Longitude', 'Fac.Latitude'))
drop_col <- which(names(analysis_dat) %in% setdiff(outcome_names, outcome_analysis))
conf <- setdiff(1:ncol(analysis_dat), c(out.col, trt.col, coord.cols, drop_col))
subdta <- analysis_dat[, unique(c(trt.col, out.col, coord.cols, conf)), with = FALSE]


# ------------------- PART 3 ------------------- #
# --------- PERFORMING THE ANALYSIS --------- #

# ---  3a. Initializing the analysis.

set.seed(1234)

trt.col <- 1
out.col <- 2
coord.cols <- c(4, 3)

# Dropping missing data
wh <- is.na(subdta[, names(subdta)[2], with = FALSE])
if (any(wh)) {
  print(paste('Dropping', sum(wh), 'observations due to missing outcome.'))
  subdta <- subdta[- which(wh), ]
}

# Defining the results matrix.
methods <- c('Naive', 'GBM', 'Distance Caliper', 'DAPSm', 'Keele et al')
result <- array(NA, dim = c(length(methods), 3))
dimnames(result) <- list(method = methods,
                         statistic = c('LB', 'Estimate', 'UB'))
num_match <- numeric(length(methods))
names(num_match) <- methods
distance <- num_match

naive_ps <- glm(as.formula(paste('SnCR ~ . - Fac.Latitude - Fac.Longitude -',
                                 outcome_analysis)), data = subdta, family = 'binomial')
subdta[, prop.scores := fitted(naive_ps)]
ignore.cols <- which(names(subdta) == 'prop.scores')
ignore.cols.coords <- c(ignore.cols, coord.cols)
cols.balance <- setdiff(1:dim(subdta)[2],
                        c(trt.col, out.col, ignore.cols.coords,
                          which(names(subdta) == 'prop.scores')))
bal <- array(NA, dim = c(length(methods) + 1, length(cols.balance)))
dimnames(bal) <- list(method = c(methods, 'Full-data'),
                      variable = names(subdta)[cols.balance])


# ----- 3b. Starting the analysis.

# Fitting the naive.
naive.match <- NaiveModel(subdta, trt.col = trt.col, out.col = out.col,
                          caliper = caliper, coord.cols = coord.cols,
                          cols.balance = cols.balance,
                          matching_algorithm = matching_algorithm,
                          remove_unmatchables = remove_unmatchables)
result[1, ] <- naive.match$result
num_match[1] <- naive.match$num_match
distance[1] <- naive.match$distance
bal[c(1, 5), ] <- naive.match$balance[c(2, 1), ]

# Fitting GBM
GBM.match <- GBMmodel(subdta, trt.col, out.col, caliper, coord.cols,
                      cols.balance = cols.balance, seed = 1234,
                      matching_algorithm = matching_algorithm,
                      remove_unmatchables = remove_unmatchables)
result[2, ] <- GBM.match$result
num_match[2] <- GBM.match$num_match
distance[2] <- GBM.match$distance
bal[2, ] <- GBM.match$balance[2, ]

# Fitting Distance Caliper
cal.match <- DistCalModel(subdta, caliper, dist.caliper = 0.2, coord.cols,
                          ignore.cols = ignore.cols.coords, trt.col, out.col,
                          cols.balance = cols.balance, coord_dist = TRUE,
                          matching_algorithm = matching_algorithm,
                          remove_unmatchables = remove_unmatchables)
result[3, ] <- cal.match$result
num_match[3] <- cal.match$num_match
distance[3] <- cal.match$distance
bal[3, ] <- cal.match$balance[2, ]

# Fitting DAPSm.
w_bal <- CalcDAPSWeightBalance(subdta, weights, cols.balance, trt.col, out.col,
                               coords.columns = coord.cols, caliper,
                               coord_dist = TRUE,
                               matching_algorithm = matching_algorithm,
                               remove.unmatchables = remove_unmatchables)

dapsm <- DAPSchoiceModel(balance = w_bal$balance, cutoff = cutoff,
                         dataset = subdta, pairs = w_bal$pairs,
                         full_pairs = w_bal$full_pairs,
                         distance_DAPS = w_bal$distance_DAPS,
                         out.col = out.col, weights = weights,
                         trt.col = trt.col)
result[4, ] <- dapsm$est + c(- 1, 0, 1) * 1.96 * dapsm$se
num_match[4] <- dapsm$num_match
distance[4] <- dapsm$distance
bal[4, ] <- dapsm$balance[2, ]

# Plotting the standardized difference of means as a function of weight.
PlotWeightBalance(w_bal$balance, full_data = -5, weights, cutoff, inset = -0.5)
# Power plant and area-level characteristics separately.
PlotWeightBalance(abs(w_bal$balance[, , c(1, 2, 14, 16, 17, 18)]),
                  full_data = -5, weights, cutoff, axis_cex = 0.6,
                  mar = c(4, 4, 2, 4))
PlotWeightBalance(abs(w_bal$balance[, , - c(1, 2, 14, 16, 17, 18)]),
                  full_data = -5, weights, cutoff, axis_cex = 0.6,
                  mar = c(4, 4, 2, 7), inset = -0.35)



# Fitting Keele.
source('Keele_et_al_functions/01_subsetmatch2.R')
source('Keele_et_al_functions/02_errorhandling.R')
source('Keele_et_al_functions/03_problemparameters2.R')
source('Keele_et_al_functions/04_constraintmatrix2.R')
n_matches = 1  # Number of matched per treated.
subset_weight = 100  # Whether all the treated units need to be used
use_controls = NULL  # Whether specific controls need to be used
enforce_constraints <- FALSE
library(Rcplex)


keele_dta <- subdta[order(subdta$SnCR, decreasing = TRUE), ]
t_ind <- as.numeric(keele_dta$SnCR)
coords <- cbind(keele_dta$Fac.Longitude, keele_dta$Fac.Latitude)
dist_mat <- rdist.earth(coords[t_ind == 1, ], coords[t_ind == 0, ])

out <- subsetmatch(dist_mat = dist_mat, t_ind = t_ind, n_matches = n_matches, 
                   mom_covs = NULL, mom_weights = NULL, mom_tols = NULL,
                   exact_covs = NULL, near_exact_covs = NULL, near_exact_devs = NULL, 
                   fine_covs = NULL, near_fine_covs = NULL, near_fine_devs = NULL, 
                   subset_weight, use_controls, enforce_constraints)
keele_dta <- keele_dta[c(out$t_id, out$c_id), ]
keele_dta <- as.data.frame(keele_dta)
lmod <- lm(keele_dta[, out.col] ~ keele_dta[, trt.col])
result[5, ] <- lmod$coef[2] + 1.96 * summary(lmod)$coef[2, 2] * c(- 1, 0, 1)
num_match[5] <- length(out$t_id)




# Plotting the results.
PlotResults(result, title = paste(outcome_analysis, paste(month, collapse = ','),
                                  '/', year))

apply(bal, 1, function(x) c(sum = sum(abs(x) > cutoff),
                            mean = mean(abs(x)),
                            max = max(abs(x))))

num_match
distance
sum(subdta$SnCR)
dapsm$weight


# Plotting maps of the matched pairs
MatchedDataMap(naive.match$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'Naive pairs', point_data = FALSE)
MatchedDataMap(cal.match$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'Distance Caliper pairs', point_data = FALSE)
MatchedDataMap(dapsm$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'DAPSm pairs', point_data = FALSE)



