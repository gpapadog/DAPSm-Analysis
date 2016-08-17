# Author: Georgia Papadogeorgou
# Date Commented: 1/26/2016


expit <- function(x) {
  # Function that returns the expit of a given number
  # or vector.
  # Args:
  #  x: Numeric of vector
  # Returns:
  #  The expit of the argument.
  return(exp(x) / (1 + exp(x)))
}


ColDrop <- function(dataset, ignore.cols) {
  # Function that takes a dataset and drop columns of the data
  # if ignore.col is not null.
  #
  # Args:
  #  dataset:    A matrix or dataframe.
  #  ignore.col: Indeces of columns to be dropped. If set to NULL
  #              dataset is returned with no change.
  
  if (!is.null(dataset)) {
    dataset <- dataset[, - ignore.cols]
  }
  return(dataset)
}
