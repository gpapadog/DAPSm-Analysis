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
library(Rcplex)


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
source('PSmatchEst_function.R')
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
analysis_dat[, hasNOxControl_byHI := NULL]

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

trt.col <- 1
out.col <- 2
coord.cols <- c(4, 3)

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
set.seed(1241)

# Fitting the naive.
naive.match <- NaiveModel(subdta, trt.col = trt.col, out.col = out.col,
                          caliper = caliper, coord.cols = coord.cols,
                          cols.balance = cols.balance,
                          matching_algorithm = matching_algorithm,
                          remove_unmatchables = remove_unmatchables)
result[1, ] <- naive.match$result
num_match[1] <- naive.match$num_match
distance[1] <- naive.match$distance
bal[c(1, 6), ] <- naive.match$balance[c(2, 1), ]
sum(abs(bal[1, ]) > cutoff)

# Fitting GBM
GBM.match <- GBMmodel(subdta, trt.col, out.col, gbm.caliper, coord.cols,
                      cols.balance = cols.balance, seed = 1234,
                      matching_algorithm = matching_algorithm,
                      remove_unmatchables = remove_unmatchables)
result[2, ] <- GBM.match$result
num_match[2] <- GBM.match$num_match
distance[2] <- GBM.match$distance
bal[2, ] <- GBM.match$balance[2, ]

# Fitting Distance Caliper
cal.match <- DistCalModel(subdta, caliper, dist.caliper = dist_cal, coord.cols,
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


# Fitting Keele.
subdta <- subdta[order(subdta$SnCR, decreasing = TRUE), ]
n_trt <- sum(subdta$SnCR)

mom_covs_ind <- c(5, 7 : 19)
exact_covs_ind <- c(6, 20, 21, 22)
mom_covs <- as.matrix(subdta)[, mom_covs_ind]
mom_tols <- keele_caliper * apply(mom_covs[1 : n_trt, ], 2, sd)
exact_covs <- as.matrix(subdta)[, exact_covs_ind]

keele <- keele_match(subdta, trt_col = trt.col, out_col = out.col,
                     coords.columns = coord.cols, subset_weight = subset_weight,
                     n_matches = n_matches, use_controls = use_controls,
                     enforce_constraints = enforce_constraints, pairsRet = TRUE,
                     cols.balance = cols.balance, mom_tols = mom_tols,
                     mom_covs = mom_covs, exact_covs = exact_covs)
result[5, ] <- keele$CI
num_match[5] <- keele$num_match
distance[5] <- mean(keele$distance)
bal[5, ] <- keele$balance[2, ]





# ----- PART 4. Looking at the results. ------ #


# Plotting the standardized difference of means as a function of weight.
PlotWeightBalance(w_bal$balance, full_data = -5, weights, cutoff, inset = -0.5)
# Power plant and area-level characteristics separately.
cov_names <- c('% Operating Capacity', 'ARP Phase 2', '4th Max Temp', '% Urban',
               '% White', '% Black', '% Hispanic', '% High School',
               'Household Income', '% Poor', '% Occupied', '% 5-year residents',
               'House Value', 'Heat Input', 'Population / square mile',
               'Gas facility', 'Small sized facility', 'Medium sized facility')
dimnames(w_bal$balance)[[3]] <- cov_names

cols <- paste0('gray', c(80, 60, 30, 5))
PlotWeightBalance(abs(w_bal$balance[, , c(1, 2, 14, 16, 17, 18)]),
                  full_data = -5, weights, cutoff, axis_cex = 0.8,
                  mar = c(4, 4, 2, 3), leg_cex = 0.7, inset = - 0,
                  cols = cols)
title(main = 'Power plant characteristics')
PlotWeightBalance(abs(w_bal$balance[, , - c(1, 2, 14, 16, 17, 18)]),
                  full_data = -5, weights, cutoff, axis_cex = 0.8,
                  mar = c(4, 4, 2, 3), inset = - 0, leg_cex = 0.7,
                  cols = cols)
title(main = 'Area level characteristics')


apply(bal, 1, function(x) c(sum = sum(abs(x) > cutoff),
                            mean = mean(abs(x)),
                            max = max(abs(x))))

num_match
distance
sum(subdta$SnCR)
dapsm$weight


# Plotting the results.
PlotResults(result[c(1, 3, 5, 4), ], color = c(rep('grey65', 3), 'grey35'),
            title = paste(outcome_analysis, paste(month, collapse = ','), '/', year))

# Plotting maps of the matched pairs
MatchedDataMap(naive.match$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'Naive pairs', point_data = FALSE)
MatchedDataMap(GBM.match$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'GBM pairs', point_data = FALSE)
MatchedDataMap(cal.match$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'Distance Caliper pairs', point_data = FALSE)
MatchedDataMap(dapsm$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'DAPSm pairs', point_data = FALSE)
MatchedDataMap(keele$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'Keele et al pairs', point_data = FALSE)


DAPSWeightCE(dataset = subdta, out.col = out.col, trt.col = trt.col,
             weights = weights, pairs = w_bal$pairs, chosen_w = dapsm$weight)