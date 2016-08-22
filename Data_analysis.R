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

source('~/Github/DAPSm-Analysis/config.R')

# Setting the working directory.
setwd(wd)

# Sourcing the functions that I will be using.
source(paste0(source.path, 'PSmatchEst_function.R'))
source(paste0(source.path, 'DistCal_functions.R'))
source(paste0(source.path, 'GBMPropScores_function.R'))
source(paste0(source.path, 'expit.R'))
source(paste0(source.path, 'StandDiff_function.R'))
source(paste0(source.path, 'Data_analysis_functions.R'))
source(paste0(source.path, 'Data_analysis_models_functions.R'))
source(paste0(source.path, 'make_data_functions.R'))
source(paste0(source.path, 'predict_variable_functions.R'))
source(paste0(source.path, 'CreateNOxControlsFunction.R'))
source(paste0(source.path, 'PredictHeatInput.R'))


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
dat <- LinkPPtoMonitors(dat_facility, within_km, year, month, OzTempCensus)

# Dropping facilities with missing data for at least one month.
wh <- which(dat$nmonths != length(month))
print(paste('Dropping', length(wh), 'out of', length(unique(dat$FacID)),
            'facilities due to missing information on at least one month'))
dat <- dat[- wh, ]

analysis_dat <- CleanData(dat, plotcor = FALSE)
analysis_dat <- ReformData(analysis_dat)
# analysis_dat includes data on the observations we will use.



# ------------------- PART 2------------------- #
# --------- SETTING UP THE ANALYSIS --------- #
outcome_names <- c('meanOzone', 'meanmaxOzone', 'mean4maxOzone', 'totNOxemissions')
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
result <- array(NA, dim = c(4, 3))
dimnames(result) <- list(methods = c('Naive', 'GBM', 'Distance Caliper', 'DAPSm'),
                         statistic = c('LB', 'Estimate', 'UB'))

num_match <- numeric(4)
names(num_match) <- dimnames(result)[[1]]
distance <- num_match

naive_ps <- glm(as.formula(paste('SnCR ~ . - Fac.Latitude - Fac.Longitude -',
                                 outcome_analysis)), data = subdta, family = 'binomial')
subdta[, prop.scores := fitted(naive_ps)]
ignore.cols <- which(names(subdta) == 'prop.scores')
ignore.cols.coords <- c(ignore.cols, coord.cols)
cols.balance <- setdiff(1:dim(subdta)[2],
                        c(trt.col, out.col, ignore.cols.coords,
                          which(names(subdta) == 'prop.scores')))
bal <- array(NA, dim = c(5, length(cols.balance)))
dimnames(bal) <- list(method = c(dimnames(result)[[1]], 'Full-data'),
                      variable = names(subdta)[cols.balance])


# ----- 3b. Starting the analysis.

# Fitting the naive.
naive.match <- NaiveModel(subdta, trt.col, out.col, caliper, coord.cols,
                          cols.balance = cols.balance)
result[1, ] <- naive.match$result
num_match[1] <- naive.match$num_match
distance[1] <- naive.match$distance
bal[c(1, 5), ] <- naive.match$balance[c(2, 1), ]

# Fitting GBM
GBM.match <- GBMmodel(subdta, trt.col, out.col, caliper, coord.cols,
                      cols.balance = cols.balance, seed = 1234)
result[2, ] <- GBM.match$result
num_match[2] <- GBM.match$num_match
distance[2] <- GBM.match$distance
bal[2, ] <- GBM.match$balance[2, ]

# Fitting Distance Caliper
cal.match <- DistCalModel(subdta, caliper, dist.caliper = 0.2, coord.cols,
                          ignore.cols = ignore.cols.coords, trt.col, out.col,
                          cols.balance = cols.balance, coord_dist = TRUE)
result[3, ] <- cal.match$result
num_match[3] <- cal.match$num_match
distance[3] <- cal.match$distance
bal[3, ] <- cal.match$balance[2, ]

# Fitting DAPSm.
w_bal <- CalcDAPSWeightBalance(subdta, weights, cols.balance, trt.col, out.col,
                               coords.columns = coord.cols, caliper,
                               coord_dist = TRUE)

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


# Plotting the results.
PlotResults(result, title = paste(outcome_analysis, paste(month, collapse = ','),
                                  '/', year))

apply(bal, 1, function(x) c(sum = sum(abs(x) > cutoff),
                            mean = mean(abs(x)),
                            max = max(abs(x))))

num_match
distance
sum(subdta$SnCR)

DAPS.match.opt$weight
DAPS.match.choice$weight


# Plotting maps of the matched pairs
MatchedDataMap(naive.match$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'Naive pairs')
MatchedDataMap(cal.match$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'Distance Caliper pairs')
MatchedDataMap(DAPS.match.choice$pairs, trt_coords = c(3, 4), con_coords = c(7, 8),
               plot.title = 'DAPSm pairs')

# Getting causal effect estimate for different weight.
CEweight <- DAPSWeightCE(subdta, trt.col = trt.col, weights = weights,
                         pairs = w_bal$pairs, out.col = out.col,
                         chosen_w = dapsm$weight)
CEweight$plot


