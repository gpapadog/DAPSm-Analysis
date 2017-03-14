# Date Created: 2/9/2017
# Description: Using the optimal matching methods for acquiring matches.
# Author: Georgia Papadogeorgou
#

#' Optimal matching based on the propensity score difference and effect estimation.
#' 
#' @param dataset Data frame including at least an outcome and treatment column.
#' @param out.col If outcome column name is not 'Y', out.col should be the index of the
#' outcome column.
#' @param trt.col If treatment is not named 'X', trt.col should be set to the index of
#' the treatment column.
#' @param pscores If there is a column in dataset names prop.scores that includes the
#' propensity score values to be used in matching pscores can be left NULL. If not,
#' pscores can be numeric representing the index of dataset column that is the
#' propensity score column. Otherwise, it can be a vector of propensity scores, and of
#' length equal to the number of observations in dataset.
#' @param caliper A caliper for the PS difference of matched pairs. Defaults to 0.1.
#' @param SEreturn Logical. Whether we want to return the standard error. Defaults to
#' FALSE.
#' @param pairsRet Whether we want to return the information on the matched pairs.
#' Logical. Defaults to FALSE.
#' @param coord.cols Only necessary when pairsRet is set to TRUE. It is the indeces of
#' the columns that include the coordinates.
#' @param true_value Numeric. If provided, an indicator of whether the CI covers the
#' the true value is returned.
#' @param remove.unmatchables Logical. Argument of the optmatch function. Defaults to
#' FALSE. If set to FALSE, the matching fails unless all treated units are matched. If
#' set to TRUE, matching might return matches only for some of the treated units.
OptPSmatch <- function(dataset, out.col = NULL, trt.col = NULL, pscores = NULL,
                       caliper = 0.1, SEreturn = FALSE, pairsRet = FALSE,
                       coord.cols = NULL, true_value = NULL,
                       remove.unmatchables = FALSE) {
  
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
  
  # Creating distance matrix based on the propensity score.
  wh_trt <- which(dataset$X == 1)
  wh_con <- which(dataset$X == 0)
  
  D <- match_on(dataset$X ~ dataset$prop.scores, method = 'euclidean')
  D <- D + caliper(D, caliper * sd(dataset$prop.scores))
  
  opt_match <- pairmatch(as.matrix(D), data = dataset,
                         remove.unmatchables = remove.unmatchables)
  subdta <- dataset[!is.na(opt_match), ]
  
  # If no pairs are returned, return NA for all arguments.
  if (nrow(subdta) == 0) {
    warning('No matches found.')
    r <- NULL
    r$est <- NA
    r$SE <- NA
    r$cover <- NA
    r$pairs <- NA
  } else {
    lmod <- lm(Y ~ X, data = subdta)
    
    r <- NULL
    r$est <- lmod$coefficients[2]
    if (SEreturn) {
      r$SE <- summary(lmod)$coef[2, 2]
    }
    
    if (!is.null(true_value)) {
      se <- as.numeric(summary(lmod)$coef[2, 2])
      r$cover <- (abs(true_value - r$est) < qnorm(0.975) * se)
    }
    
    if (pairsRet) {
      which_cols <- c(which(names(dataset) == 'X'), out.col,
                      which(names(dataset) == 'prop.scores'), coord.cols)
      pairs_ids <- sort(as.character(unique(opt_match[!is.na(opt_match)])))
      pairs_ids <- data.frame(group = pairs_ids)
      match_trt <- cbind(wh_trt, group = as.character(opt_match[wh_trt]))
      match_con <- cbind(wh_con, group = as.character(opt_match[wh_con]))
      pairs_ids <- merge(pairs_ids, match_trt, by = 'group')
      pairs_ids <- merge(pairs_ids, match_con, by = 'group')
      pairs_ids <- pairs_ids[, - 1]
      pairs_ids[, 1] <- as.numeric(as.character(pairs_ids[, 1]))
      pairs_ids[, 2] <- as.numeric(as.character(pairs_ids[, 2]))
      pairs_ids <- pairs_ids[order(pairs_ids[, 1]), ]
      
      pairs <- cbind(dataset[pairs_ids$wh_trt, which_cols],
                     dataset[pairs_ids$wh_con, which_cols])
      names(pairs) <- paste0(rep(c('Trt.', 'Con.'),
                                 each = length(which_cols)), names(pairs))
      pairs$IDtrt <- pairs_ids$wh_trt
      pairs$IDcon <- pairs_ids$wh_con
      
      # Rearranging to ensure that columns 9, 10 are still the IDs.
      r$pairs <- as.matrix(pairs)[, c(1, 3:6, 8:12, 2, 7)]
    }
  }
  
  return(r)
}



