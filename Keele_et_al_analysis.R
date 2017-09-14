# Loading libraries.
library(mipmatch)
library(foreign)
library(fields)
library(data.table)

source_path <- paste0('~/Documents/ARP/DAPS_Simulations/functions/Keele/')
source(paste0(source_path, '01_subsetmatch2.R'))
source(paste0(source_path, '02_errorhandling.R'))
source(paste0(source_path, '03_problemparameters2.R'))
source(paste0(source_path, '04_constraintmatrix2.R'))
source('~/Github/DAPSm-Analysis/keele_match_function.R')
source('~/Github/DAPSm/R/CalculateBalance_function.R')

load('~/Dropbox/DATAverse/subdta.dat')
subdta[, SnCR := as.numeric(SnCR)]
naive_ps <- glm(as.formula(paste('SnCR ~ . - Fac.Latitude - Fac.Longitude -',
                                 'mean4maxOzone')), data = subdta,
                family = 'binomial')
subdta[, prop.scores := fitted(naive_ps)]

subdta <- as.data.frame(subdta)
subdta <- subdta[order(subdta$SnCR, decreasing = TRUE), ]
subdta$MedianHValue <- subdta$MedianHValue / 1000

t_ind <- subdta$SnCR
n_trt <- sum(t_ind)
n_matches <- 1  # Number of matched per treated.
use_controls <- NULL  # Whether specific controls need to be used
enforce_constraints <- FALSE
coord.cols <- c(4, 3)

trt.col <- 1
out.col <- 2
ignore.cols <- which(names(subdta) == 'prop.scores')
ignore.cols.coords <- c(ignore.cols, coord.cols)
cols.balance <- setdiff(1:dim(subdta)[2],
                        c(trt.col, out.col, ignore.cols.coords,
                          which(names(subdta) == 'prop.scores')))
caliper <- 1

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

mom_covs_ind <- setdiff(which(names(subdta) %in% names(naive_ps$coefficients)[- 1]),
                        c(6, 20, 21, 22))
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




