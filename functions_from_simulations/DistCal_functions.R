# Created:      3/13/2016
# Description:  NEW distance caliper functions that use the matrix matching
#               technique.
# Author:       Georgia Papadogeorgou
#
# Date: 4/9/2016
# Update: I got rid of the option to apply a caliper on the absolute difference
# of PS.
#
# Date: 7/1/2016
# Update: I return the outcome information on the pairs also.
# Update: (7/29/2016) Added coord_dist option.


dist.caliper <- function(treated, control, ps.caliper = 0.1,
                         dist.quan = 0.25, coords.columns = NULL,
                         coord_dist = FALSE) {
  # Function that takes in a two data frames one with treatment
  # and one with control units and returns a matrix of the matched
  # pairs, that are matched using a distance caliper and a PS caliper.
  #
  # Args:
  #  treated:        A data frame include the treated units and the variables:
  #                  longitude, latitude and propensity scores (must be named
  #                  'prop.scores'). The rownames of treated should be the unit
  #                  ids.
  #  control:        Control units. Same variables as in treated. Rownames
  #                  should be the unit ids of the controls
  #  ps.caliper:     A caliper of propensity score difference for matching.
  #                  Caliper is set as number of sd of the ps distribution.
  #  dist.quan:      The distance of all treated-control pairs is calculated.
  #                  dist.quan is a scalar between 0, and 1 describing the
  #                  quantile of all the pairwise distances that should be used
  #                  as a distance caliper.
  #  coords.columns: If the columns of coordinates are not named 'Longitude',
  #                  'Latitude', coords.cols should be the column indeces
  #                  corresponding to longitude and latitude accordingly.
  #  coord_dist:     Set to true when we want to use a distance function that
  #                  calculates the spherical distance of points instead of
  #                  euclidean. Defaults to FALSE.
  #
  # Returns:
  #  A dataframe, where each row corresponds to each treated unit, and includes
  #  the control unit to which it was matched, their propensity score difference,
  #  their DAPS difference, their distance, their standardized distance.
  
  require(fields)  # For rdist().
  
  # Setting the caliper.
  caliper <- caliper * sd(c(treated$prop.scores, control$prop.scores))
  
  if (!is.null(coords.columns)) {
    names(treated)[coords.columns] <- c('Longitude', 'Latitude')
    names(control)[coords.columns] <- c('Longitude', 'Latitude')
  }
  
  if (coord_dist) {
    dist.mat <- rdist.earth(cbind(treated$Longitude, treated$Latitude),
                            cbind(control$Longitude, control$Latitude))
  } else {
    dist.mat <- rdist(cbind(treated$Longitude, treated$Latitude),
                      cbind(control$Longitude, control$Latitude))
  }
  colnames(dist.mat) <- rownames(control)
  rownames(dist.mat) <- rownames(treated)
  
  cut.off <- quantile(dist.mat, probs = dist.quan)
  indic_dist <- ifelse(dist.mat > cut.off, Inf, 1)

  mat <- data.frame(match = rep(NA, dim(treated)[1]),
                    distance = rep(NA, dim(treated)[1]),
                    # Matching on PS_difference.
                    prop.diff = rep(NA, dim(treated)[1]))
  rownames(mat) <- rownames(treated)

  prop_diff <- t(sapply(treated$prop.scores,
                        function(x) rep(x, nrow(control))))
  prop_diff <- prop_diff - sapply(control$prop.scores,
                                  function(x) rep(x, nrow(treated)))
  prop_diff <- abs(prop_diff)
  
  # Putting distance = Inf for large distances.
  prop_diff <- prop_diff * indic_dist
  
  pairs <- MinDistMatch(prop_diff, caliper)
  matched_trt <- pairs[, 1]
  matched_con <- pairs[, 2]
  
  mat$match[matched_trt] <- rownames(control)[matched_con]
  for (ii in 1:length(matched_trt)) {
    wh_trt <- matched_trt[ii]
    wh_con <- matched_con[ii]
    mat$prop.diff[wh_trt] <- treated$prop.scores[wh_trt] -
      control$prop.scores[wh_con]
    mat$distance[wh_trt] <- dist.mat[wh_trt, wh_con]
  }
  return(mat)
}



CaliperEst <- function(dataset, out.col = NULL, trt.col = NULL,
                       ps.caliper = 0.1, dist.caliper = 0.25,
                       coords.columns = NULL, ignore.cols = NULL,
                       SEreturn = FALSE, pairsRet = FALSE,
                       coord_dist = FALSE, true_value = NULL) {
  # Function that fits matching in distance calipers and returns the
  # causal effect estimate under a specification of caliper, distance
  # caliper, and number of permutations of matching.
  #
  # Args:
  #  dataset:       Data frame including treatment, outcome, coordinates, and
  #                 observed confounders.
  #  out.col:       If outcome column name is not 'Y', out.col should be the
  #                 index of the outcome column.
  #  trt.col:       If treatment is not named 'X', trt.col should be set to the
  #                 index of the treatment column.
  #  ps.caliper:    A caliper for the PS difference of matched pairs. Defaults
  #                 to 0.1. Scalar. Caliper is defined as the number of sd of
  #                 the ps distribution.
  #  dist.caliper:  The quantile of all pairwise treated-control distances that
  #                 should be used as distance caliper.
  #  coords.columns:If the columns of coordinates are not named 'Longitude' and
  #                 Latitude', coords.columns are the column indeces corresponding
  #                 to longitude and latitude accordingly.
  #  ignore.cols:   All column indeces that should not be included in the linear
  #                 model. Often, this should be set to all column indeces
  #                 corresponding to columns other than outcome, treatment, and
  #                 and observed confounders.
  #  SEreturn:      Whether a standard error estimate is returned. Defaults to
  #                 FALSE.
  #  pairsRet:      Whether we want to return the information on the matched
  #                 pairs. Logical. Defaults to FALSE.
  #  coord_dist:    Set to true when we want to use a distance function that
  #                 calculates the spherical distance of points instead of
  #                 euclidean. Defaults to FALSE.
  #  true_value:    Numeric. If provided, an indicator of whether the CI covers
  #                 the true value is returned.
  #
  # Returns:
  #  A vector of length 2. The first element is the causal effect
  #  estimated from a linear model on the matched pairs, adjusting
  #  for no observed confounding. The second element is the causal
  #  effect estimate from a linear model on the matched pairs that
  #  adjusts for all covariates that are not in ignore.cols.
  
  dataset <- as.data.frame(dataset)
  # Naming outcome and treatment as 'Y', 'X'.
  dataset <- FormDataset(dataset, ignore.cols = NULL,
                         out.col = out.col, trt.col = trt.col)
  
  cal.daps <- dist.caliper(treated = dataset[dataset$X == 1, ],
                           control = dataset[dataset$X == 0, ],
                           ps.caliper = ps.caliper,
                           dist.quan = dist.caliper,
                           coords.columns = coords.columns,
                           coord_dist = coord_dist)
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
  
  r <- NULL
  r$est <- est
  
  if (!is.null(true_value)) {
    se_est <- c(summary(lmod1)$coef[2, 2], summary(lmod2)$coef[2, 2])
    r$cover <- (abs(true_value) - est < qnorm(0.975) * se_est)
  }
  
  if (SEreturn) {
    r$SE <- c(summary(lmod1)$coef[2, 2], summary(lmod2)$coef[2, 2])
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
