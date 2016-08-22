# Created: 10/1/2015
# Description: Functions that perform estimation of the propensity
#              score using GBM.
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

