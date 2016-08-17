
expit <- function(x) {
  # Function that returns the expit of a given number
  # or vector.
  # Args:
  #  x: Numeric of vector
  # Returns:
  #  The expit of the argument.
  return(exp(x) / (1 + exp(x)))
}
