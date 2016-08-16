# Author: Georgia Papadogeorgou
# Date: 1/6/2015
# Description: Functions that calculate balance of variables.
# Last Updated: 2/27/2016
# Update: Added Mean distance.

StandDiff <- function(dataset, trt, con, col) {
  # Function that calculates the standardized difference in means
  # for two subsets of a dataset.
  #
  # Args:
  #  dataset: Includes the information on the covariates whose
  #           standardized difference we want to calculate.
  #  trt:     Indeces of 1st subset of the dataset. This subset
  #           will be used for calculating the standard deviation.
  #  con:     Indeces of the 2nd subset.
  #  col:     1 or more column indeces for the columns whose
  #           standardized difference we want to calculate.
  #
  # Returns:
  #  A vector where each element is the standardized difference
  #  of each column in col, for the two data subsets.
  
  sub1 <- as.matrix(dataset[trt, col], nrow = length(trt),
                    ncol = length(col))
  sub2 <- as.matrix(dataset[con, col], nrow = length(con),
                    ncol = length(col))
  
  stand_diff <- apply(sub1, 2, mean) - apply(sub2, 2, mean)
  stand_diff <- stand_diff / apply(sub1, 2, sd)
  return(stand_diff)
}

MeanDistance <- function(dataset, trt, con, coords.cols){
  # Function that calculates the mean Eucleadean distance of
  # the coordinate columns coords.cols, for pairs in trt, con.
  #
  # Args:
  #  dataset:     Includes the coordinates of the points.
  #  trt:         Indeces of 1st subset of dataset.
  #  con:         Indeces of 2nd subset. Equal length to trt.
  #  coords.cols: Columns of coordinates.
  
  sub1 <- as.matrix(dataset[trt, coords.cols], nrow = length(trt),
                    ncol = length(coords.cols))
  sub2 <- as.matrix(dataset[con, coords.cols], nrow = length(con),
                    ncol = length(coords.cols))
  
  distances <- (sub1 - sub2) ^ 2
  distances <- mean(apply(distances, 1, function(x) sqrt(sum(x))))
  return(distances)
}


