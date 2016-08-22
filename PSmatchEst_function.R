# Date Created: 1/18/2016
# Author: Georgia Papadogeorgou

PSmatchEst <- function(dataset, out.col = NULL, trt.col = NULL,
                       pscores = NULL, estimand = 'ATT',
                       caliper = 0.1, SEreturn = FALSE,
                       pairsRet = FALSE, replace = FALSE,
                       coord.cols = NULL, true_value = NULL) {
  # Function uses the given propesity scores to estimate the
  # causal effect of a treatment.
  # 
  # Args:
  #  dataset: Data frame including all the variables that we want to include
  #           in the linear model. dataset should include at least an
  #           outcome and treatment column.
  #  out.col: If outcome column name is not 'Y', out.col should be the index
  #           of the outcome column.
  #  trt.col: If treatment is not named 'X', trt.col should be set to the
  #           index of the treatment column.
  #  pscores: If there is a column in dataset names prop.scores that includes
  #           the propensity score values to be used in matching pscores can
  #           be left NULL. If not, pscores can be numeric representing the
  #           index of dataset column that is the propensity score column.
  #           Otherwise, it can be a vector of propensity scores, and of
  #           length equal to the number of observations in dataset.
  # estimand: A character string for the estimand. The default estimand is
  #           "ATT", the sample average treatment effect for the treated.
  #           "ATE" is the sample average treatment effect, and 'ATC' is the
  #           sample average treatment effect for the controls.
  # caliper:  A caliper for the PS difference of matched pairs. Defaults to
  #           0.1.
  # SEreturn: Logical. Whether we want to return the standard error. Defaults
  #           to FALSE.
  # pairsRet: Whether we want to return the information on the matched pairs.
  #           Logical. Defaults to FALSE.
  # replace:  Whether we want matching with replacement. Logical. Defaults to
  #           FALSE.
  # coord.cols: Only necessary when pairsRet is set to TRUE. It is the
  #           indeces of the columns that include the coordinates.
  # true_value: Numeric. If provided, an indicator of whether the CI covers
  #           the true value is returned.
  #
  # Returns:
  #  The estimate from the Match() function of the Matching R package, using
  #  the given propensity scores. It also returns the SE estimated from the
  #  Match() function, as well as the matched pairs.
  
  # Reforming the dataset so that the names of outcome and treatment columns
  # are 'Y', 'X' accordingly.
  dataset <- as.data.frame(dataset)
  dataset <- FormDataset(dataset = dataset, ignore.cols = NULL,
                         out.col = out.col, trt.col = trt.col)
  
  # Saving the index of the outcome column if it has not been specified.
  if (is.null(out.col)) {
    out.col <- which(names(dataset) == 'Y')
  }
  
  if (!is.null(pscores)) {
    if (length(pscores) == 1) {
      names(dataset)[pscores] <- 'prop.scores'
    } else {
      dataset$prop.scores <- pscores
    }
  }
  
  match.mod <- Match(Y = dataset$Y, Tr = dataset$X,
                     X = dataset$prop.scores, replace = replace,
                     estimand = estimand, caliper = caliper)
  
  r <- NULL
  r$est <- as.numeric(match.mod$est)
  if (SEreturn) {
    r$SE <- c(SE = match.mod$se, SE.stand = match.mod$se.standard)
  }
  
  if (!is.null(true_value)) {
    r$cover <- (abs(true_value - r$est) <
                  qnorm(0.975) * as.numeric(match.mod$est))
  }
  
  ### Add pairs here from match.mod in the form that we want.
  if (pairsRet) {
    which_cols <- c(which(names(dataset) == 'X'), out.col,
                    which(names(dataset) == 'prop.scores'), coord.cols)
    pairs <- cbind(dataset[match.mod$index.treated, which_cols],
                   dataset[match.mod$index.control, which_cols])
    names(pairs) <- paste0(rep(c('Trt.', 'Con.'),
                               each = length(which_cols)), names(pairs))
    pairs$IDtrt <- match.mod$index.treated
    pairs$IDcon <- match.mod$index.control
    # Rearranging to ensure that columns 9, 10 are still the IDs.
    r$pairs <- as.matrix(pairs)[, c(1, 3:6, 8:12, 2, 7)]
  }
  return(r)
}



