# Author: Georgia Papadogeorgou
# Date Commented: 1/26/2016
# Last Updated: 3/20/2016
# Update: Moved FormDataset() in this file.


FormDataset <- function(dataset, ignore.cols = NULL,
                        out.col = NULL, trt.col = NULL) {
  # Function that takes a dataset, and reforms it if any of the
  # ignore.cols, out.col, or trt.col is not NULL. Specificially,
  # it drops columns in ignore.cols, renames the outcome column
  # to 'Y', and the treatment column to 'X'.
  #
  # Args:
  #  dataset:     Data frame that we want to reform.
  #  ignore.cols: Indeces of columns that should be dropped from
  #               the data frame. If not specified, no columns
  #               are dropped.
  #  out.col:     If out.col is not NULL, the column of index
  #               out.col will be renamed to 'Y'.
  #  trt.col:     If trt.col is not NULL, the column of index
  #               trt.col will be renamed to 'X'.
  #
  # Returns:
  #  A data frame of same number of rows, after dropping columns
  #  in ignore.cols, and renaming columns out.col, trt.col.
  
  if (!is.null(out.col)) {
    names(dataset)[out.col] <- 'Y'
  }
  if (!is.null(trt.col)) {
    names(dataset)[trt.col] <- 'X'
  }
  if (!is.null(ignore.cols)) {
    dataset <- dataset[, - ignore.cols]
  }
  return(dataset)
}


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