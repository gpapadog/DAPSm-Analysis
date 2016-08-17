# Author: Georgia Papadogeorgou
# Date: 5/15/2016
# Desc: Using the facility level data, we fit the models and perform the analysis
#       using the naive, and the spatial propensity score matching methods. I also
#       perform an analysis for the choice of w for DAPS.

# Loading libraries.
library(data.table)
library(corrplot)
library(MatchIt)
library(Matching)
library(fields)
library(ggplot2)
library(DAPSm)

setwd("/Users/georgiapapadogeorgou//Documents/ARP/Application")

# Sourcing in functions.
# Sourcing the functions that I will be using
source.path <- paste0('/Users/georgiapapadogeorgou/Documents/ARP/',
                      'DAPS_Simulations/functions/')
source(paste0(source.path, 'estimating_functions.R'))
source(paste0(source.path, 'DistCal_functions.R'))
source(paste0(source.path, 'GBM_functions.R'))
source(paste0(source.path, 'expit.R'))
source(paste0(source.path, 'balance_functions.R'))
source('Facility_level_Analysis/Data_analysis_functions.R')
source('Facility_level_Analysis/Data_analysis_models_functions.R')

# ---- 1. Loading and forming the data.


# If we want to do this for different time period, we need to change time_use in the
# predict_aggregate_link.R file, for time periods used for prediction of heat input.
#year <- 2005
month <- 6:8
within_km <- 100
#time_use <- list('2004', '2003', '2006', '2007')

# year <- 2014
# time_use <- list('2013', '2012', '2011', '2010')

year <- 2004
time_use <- list('2003', '2002', '2005', '2006')

# Sourcing the code the loads unit level data, predicts unit heat input for units with
# missing heat input information, aggregates unit level data to facility level, and links
# with monitoring data.
source('Make Data Code/predict_aggregate_link.R')

analysis_dat <- CleanData(dat, plotcor = FALSE)
analysis_dat <- ReformData(analysis_dat)

# Which temperature do we want to keep:
temp_keep <- 'mean4MaxTemp'
temp_names <- c('avgTemp', 'mean4MaxTemp', 'meanMaxTemp')
analysis_dat[, setdiff(temp_names, temp_keep) := NULL, with = FALSE]

# for (ii in 1:ncol(subdta)) {
#   plot(as.data.frame(subdta)[, ii], subdta$meanOzone, cex = 0.8,
#        col = ifelse(subdta[, SnCR] == 1, 1, 12), main = colnames(subdta)[ii])
# }


# ------ 2. Setting up the analysis.

outcome_analysis <- 'totNOxemissions'

outcome_names <- c('meanOzone', 'meanmaxOzone', 'mean4maxOzone', 'totNOxemissions')
out.col <- which(names(analysis_dat) == outcome_analysis)
trt.col <- which(names(analysis_dat) == 'SnCR')
coord.cols <- which(names(analysis_dat) %in% c('Fac.Longitude', 'Fac.Latitude'))
drop_col <- which(names(analysis_dat) %in% setdiff(outcome_names, outcome_analysis))
conf <- setdiff(1:ncol(analysis_dat), c(out.col, trt.col, coord.cols, drop_col))
subdta <- analysis_dat[, unique(c(trt.col, out.col, coord.cols, conf)), with = FALSE]


# ----- 3. Analysis.


# ---  3a. Initializing the analysis.

set.seed(1234)

trt.col <- 1
out.col <- 2
caliper <- 1
cutoff <- 0.15
subsampling <- FALSE
subsamples <- 200
weights <- seq(0, 1, length.out = 40)
coord.cols <- c(4, 3)

# Dropping missing data
wh <- is.na(subdta[, names(subdta)[2], with = FALSE])
if (any(wh)) {
  print(paste('Dropping', sum(wh), 'observations due to missing outcome.'))
  subdta <- subdta[- which(wh), ]
}

# Defining the results matrix.
result <- array(NA, dim = c(5, 3))
dimnames(result) <- list(methods = c('Naive', 'GBM', 'Distance Caliper',
                                     'DAPS optimal','DAPS choice'),
                         statistic = c('LB', 'Estimate', 'UB'))
result_corSE <- result # Saving the results with SEs calculated by subsampling.
rownames(result_corSE) <- paste0(rownames(result_corSE), '-CorSE')

num_match <- numeric(5)
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
bal <- array(NA, dim = c(6, length(cols.balance)),
             dimnames = list(method = c(dimnames(result)[[1]], 'Full-data'),
                             variable = names(subdta)[cols.balance]))


# ----- 3b. Starting the analysis.

# Fitting the naive.
naive.match <- NaiveModel(subdta, trt.col, out.col, caliper, coord.cols,
                          subsampling = subsampling, subsamples = subsamples,
                          cols.balance = cols.balance)
result[1, ] <- naive.match$result
num_match[1] <- naive.match$num_match
distance[1] <- naive.match$distance
bal[c(1, 6), ] <- naive.match$balance[c(2, 1), ]

# Fitting GBM
GBM.match <- GBMmodel(subdta, trt.col, out.col, caliper, coord.cols,
                      subsampling = subsampling, subsamples = subsamples,
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


# Fitting DAPS optimal
DAPS.match.opt <- DAPSoptModel(subdta, trt.col, out.col, coord.cols, caliper,
                               cols.balance, cutoff, coord_dist = TRUE)
result[4, ] <- DAPS.match.opt$result
num_match[4] <- DAPS.match.opt$num_match
distance[4] <- DAPS.match.opt$distance
bal[4, ] <- DAPS.match.opt$balance


# Fitting DAPS choice.
w_bal <- CalcDAPSWeightBalance(subdta, weights, cols.balance, trt.col, out.col,
                               coords.columns = coord.cols, caliper,
                               coord_dist = TRUE)

DAPS.match.choice <- DAPSchoiceModel(balance = w_bal$balance, cutoff = cutoff,
                                     dataset = subdta, pairs = w_bal$pairs,
                                     full_pairs = w_bal$full_pairs,
                                     distance_DAPS = w_bal$distance_DAPS,
                                     out.col = out.col, weights = weights)
result[5, ] <- DAPS.match.choice$est + c(- 1, 0, 1) * 1.96 * DAPS.match.choice$se
num_match[5] <- DAPS.match.choice$num_match
distance[5] <- DAPS.match.choice$distance
bal[5, ] <- DAPS.match.choice$balance[2, ]

# Plotting the standardized difference of means as a function of weight.
PlotWeightBalance(w_bal$balance, full_data = -5, weights, cutoff, inset = -0.5)
# Power plant and area-level characteristics separately.
PlotWeightBalance(abs(w_bal$balance[, , c(1, 2, 14, 16, 17, 18)]),
                  full_data = -5, weights, cutoff, axis_cex = 0.6,
                  mar = c(4, 4, 2, 4))
PlotWeightBalance(abs(w_bal$balance[, , - c(1, 2, 14, 16, 17, 18)]),
                  full_data = -5, weights, cutoff, axis_cex = 0.6,
                  mar = c(4, 4, 2, 7), inset = -0.35)

if (subsampling) {
  result_corSE[1, ] <- naive.match$result_corSE
  result_corSE[2, ] <- GBM.match$result_corSE
  result_corSE[3, ] <- cal.match$result_corSE
  result_corSE[4, ] <- DAPS.match.opt$result_corSE
  result_corSE[5, ] <- DAPS.match.choice$result_corSE
}

# Plotting the results.
PlotResults(result, title = paste(outcome_analysis,
                                  paste(month, collapse = ','), '/', year))

apply(bal, 1, function(x) c(sum = sum(abs(x) > cutoff),
                            mean = mean(abs(x)),
                            max = max(abs(x))))

num_match
sum(subdta$SnCR)

DAPS.match.opt$weight
DAPS.match.choice$weight


# Getting causal effect estimate for different weight.
CEweight <- DAPSWeightCE(subdta, weights, w_bal$pairs, out.col = out.col,
                         chosen_w = DAPS.match.choice$weight)
CEweight$plot


