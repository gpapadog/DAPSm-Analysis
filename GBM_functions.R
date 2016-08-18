# Created: 10/1/2015
# Last Updated: 10/20/2015
# Description: Functions that perform estimation of the propensity
#              score. Based on this propensity score estimates, we
#              estimate the causal effect using 1-1 nearest neighbor
#              matching.
# Author: Georgia Papadogeorgou

GBMPropScores <- function(dta, ignore.cols = NULL, interaction.depth = 3,
                          trt.col = NULL) {
  # Function that receives a dataset, and calculates the propensity score
  # of the treatment X using generalized boosted models.
  #
  # Args:
  #  dta:        The dataset on which the PS and linear models are fit.
  #  ignore.cols:The indeces of the columns that should not be included
  #              in the PS model.
  #  trt.col:    Column index of treatment. Required if treatment is
  #              not named 'X'.
  #  interaction.depth: GBM interaction depth.
  #
  # Returns:
  #  Vector of propensity score estimates.
  require(gbm)
  
  dta <- as.data.frame(dta)
  dta <- FormDataset(dta, trt.col = trt.col)
  dta <- FormDataset(dta, ignore.cols = ignore.cols)
  dta <- dta[, which(names(dta) != 'Y')]

  # If we have observed confounding.
  if (dim(dta)[2] > 1) {
    gbmod <- gbm(X ~ ., data = dta, distribution = 'bernoulli',
                 interaction.depth = interaction.depth)
    prop.scores <- gbmod$fit
  } else {
    prop.scores <- rep(0, dim(dta)[1])
  }
  prop.scores <- expit(prop.scores)
  return(prop.scores)
}


FitMatchedGBM <- function(dta, ignore.cols = NULL, caliper = 0.1,
                          true_value = NULL) {
  # Function that receives a dataset, fits a propensity score model of
  # treatment on all variables that are not set in ignore.cols using
  # Generalized Boosted Model, and using the Match() function of the
  # Matching package estimates the effect of treated from a linear
  # model on the matched dataset.
  #
  # Args:
  #  dta:         The dataset on which the PS and linear models are fit.
  #               The outcome of interest should be named Y, and the
  #               treatment X.
  #  ignore.cols: The indeces of the columns that should not be included
  #               in the PS model.
  #  caliper:     Caliper for propensity score differnce of matched pairs
  #               Set to 0.1.
  #  true_value:  Numeric. If provided, an indicator of whether the CI
  #               covers the true value is returned.
  #
  # Returns:
  #  Effect estimate of treatment X through GBM-estimates
  #  propensity score matching.
  require(Matching)
  
  prop.scores <- GBMPropScores(dta, ignore.cols)
  match.mod <- Match(Y = dta$Y, Tr = dta$X,
                     X = prop.scores, estimand = 'ATT',
                     caliper = caliper)
  if (!is.null(true_value)) {
    r$cover <- (abs(true_value - r$est) <
                  qnorm(0.975) * as.numeric(match.mod$est))
  }
  return(match.mod$est)
}
