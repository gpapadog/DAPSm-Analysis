#' Function that fits matching in distance calipers and returns the causal effect
#' estimate under a specification of caliper, and distance caliper.
#' 
#' @param dataset
#' Data frame including treatment, outcome, coordinates, and observed confounders.
#' @out.col If outcome column name is not 'Y', out.col should be the index of the
#' outcome column.
#' @trt.col If treatment is not named 'X', trt.col should be set to the index of the
#' treatment column.
#' @param ps.caliper
#' A caliper for the PS difference of matched pairs. Defaults to 0.1. Scalar. Caliper
#' is defined as the number of sd of the ps distribution.
#' @param dist.caliper
#' The quantile of all pairwise treated-control distances that should be used as
#' distance caliper.
#' @param coords.columns
#' If the columns of coordinates are not named 'Longitude' and Latitude',
#' coords.columns are the column indeces corresponding to longitude and latitude
#' accordingly.
#' @param ignore.cols All column indeces that should not be included in the linear
#' model. Often, this should be set to all column indeces corresponding to columns
#' other than outcome, treatment, and observed confounders.
#' @param SEreturn
#' Whether a standard error estimate is returned. Defaults to FALSE.
#' @param pairsRet
#' Whether we want to return the information on the matched pairs. Logical. Defaults to
#' FALSE.
#' @param coord_dist
#' Set to true when we want to use a distance function that calculates the spherical
#' distance of points instead of euclidean. Defaults to FALSE.
#' @param true_value
#' Numeric. If provided, an indicator of whether the CI covers the true value is
#' returned.
#' @param matching_algorithm Argument with options 'optimal', or 'greedy'. The optimal
#' choice uses the optmatch R package to acquire the matches based on propensity score
#' difference and a caliper on distance. The greedy option matches treated and control
#' units sequentially, starting from the ones with the smallest propensity score
#' difference. Defaults to 'optimal'.
#' @param remove.unmatchables Logical. Argument of the optmatch function. Defaults to
#' FALSE. If set to FALSE, the matching fails unless all treated units are matched. If
#' set to TRUE, matching might return matches only for some of the treated units.
#' 
#' @return A vector of length 2. The first element is the causal effect estimated from
#' a linear model on the matched pairs, adjusting for no observed confounding. The
#' second element is the causal effect estimate from a linear model on the matched
#' pairs that adjusts for all covariates that are not in ignore.cols.
#' 
CaliperEst <- function(dataset, out.col = NULL, trt.col = NULL, ps.caliper = 0.1,
                       dist.caliper = 0.25, coords.columns = NULL, ignore.cols = NULL,
                       SEreturn = FALSE, pairsRet = FALSE, coord_dist = FALSE,
                       true_value = NULL,
                       matching_algorithm = c('optimal', 'greedy'),
                       remove.unmatchables = FALSE) {
  
  matching_algorithm <- match.arg(matching_algorithm)
  r <- NULL
  
  dataset <- as.data.frame(dataset)
  # Naming outcome and treatment as 'Y', 'X'.
  dataset <- FormDataset(dataset, ignore.cols = NULL,
                         out.col = out.col, trt.col = trt.col)
  
  cal.daps <- dist.caliper(treated = dataset[dataset$X == 1, ],
                           control = dataset[dataset$X == 0, ],
                           ps.caliper = ps.caliper, dist.quan = dist.caliper,
                           coords.columns = coords.columns, coord_dist = coord_dist,
                           matching_algorithm = matching_algorithm, 
                           remove.unmatchables = remove.unmatchables)
  
  # If no matches were acheived, return missing values.
  if (nrow(cal.daps) == 0) {
    warning(paste0('No matches were acheived for distance caliper = ', dist.caliper,
                   ', and propensity score caliper = ', ps.caliper))
    r$est <- NA
    r$se <- NA
    r$pairs <- matrix(NA, nrow = 0, ncol = 12)
    r$cover <- NA
    return(r)
  }
  
  pairs.out <- cal.daps$match
  names(pairs.out) <- rownames(cal.daps)
  pairs.out <- na.omit(pairs.out)
  pairs.daps <- dataset[c(as.numeric(names(pairs.out)),
                          as.numeric(pairs.out)), ]
  # Dropping the columns that will not be included.
  # Outcome and treatment columns have already been renamed.
  pairs_small <- FormDataset(pairs.daps, ignore.cols = ignore.cols)
  est <- numeric(2)
  lmod1 <- lm(Y~X, data = pairs_small)
  lmod2 <- lm(Y~X + ., data = pairs_small)
  est[1] <- lmod1$coef[2] 
  est[2] <- lmod2$coef[2]
  
  r$est <- est
  
  se_est <- c(summary(lmod1)$coef[2, 2], summary(lmod2)$coef[2, 2])
  if (!is.null(true_value)) {
    r$cover <- (abs(true_value - est) < qnorm(0.975) * se_est)
  }
  
  if (SEreturn) {
    r$SE <- se_est
  }
  if (pairsRet) {
    which_cols <- c(which(names(dataset) %in% c('X', 'Y', 'prop.scores')))
    which_cols <- c(which_cols, coords.columns)
    pairs <- pairs.daps[1:(nrow(pairs.daps) / 2), which_cols]
    pairs <- cbind(pairs, pairs.daps[(nrow(pairs.daps) / 2 + 1):
                                       nrow(pairs.daps), which_cols])
    names(pairs) <- paste0(rep(c('Trt.', 'Con.'),
                               each = length(which_cols)), names(pairs))
    pairs$IDtrt <- as.numeric(names(pairs.out))
    pairs$IDcon <- as.numeric(pairs.out)
    r$pairs <- as.matrix(pairs)[, c(1, 3:6, 8:12, 2, 7)]
  }
  
  return(r)
}
