# Description: Code that fits the models of the data analysis.
# Author: Georgia Papadogeorgou
# Date: 5/31/2016

GetBalancePlots <- function(subdta, matched_ind) {
  x <- as.data.frame(subdta)
  x <- x[matched_ind, ]
  plots <- NULL
  ind <- 1
  for (ii in 1:ncol(x)) {
    if (class(x[, ii]) %in% c('integer', 'numeric')) {
      if (names(x)[ii] != 'SnCR') {
        plots[[ind]] <- ggplot(x, aes_string(names(x)[ii], fill = 'SnCR')) +
          geom_histogram(bins = 40, alpha = 0.8, position = 'identity') +
          ggtitle(names(x)[ii])
        ind <- ind + 1
      }
    }
  }
  return(plots)
}

NaiveModel <- function(subdta, trt.col, out.col, caliper, coord.cols,
                       cols.balance, matching_algorithm = c('optimal', 'greedy'),
                       remove_unmatchables = TRUE) {

  matching_algorithm <- match.arg(matching_algorithm)
  r <- NULL

  if (matching_algorithm == 'greedy') {
    naive.match <- PSmatchEst(subdta, trt.col = trt.col, out.col = out.col,
                              SEreturn = TRUE, caliper = caliper,
                              replace = FALSE, estimand = 'ATT',
                              pairsRet = TRUE, coord.cols = coord.cols)
  } else {
    naive.match <- OptPSmatch(subdta, out.col = out.col, trt.col = trt.col,
                              caliper = caliper, SEreturn = TRUE, pairsRet = TRUE,
                              coord.cols = coord.cols,
                              remove.unmatchables = remove_unmatchables)
  }
  r$result <- naive.match$est + naive.match$SE * 1.96 * c(- 1, 0, 1)
  if (is.null(dim(naive.match$pairs))) {
    r$num_match <- 0
    r$distance <- Inf
    r$balance <- rep(NA, length(cols.balance))
  } else {
    r$num_match <- dim(naive.match$pairs)[1]
    r$distance <- rdist.earth(naive.match$pairs[, c(3, 4)], naive.match$pairs[, c(7, 8)])
    r$distance <- diag(r$distance)
    r$distance <- mean(r$distance)
    dtaAfter <- subdta[c(naive.match$pairs[, 9], naive.match$pairs[, 10]), ]
    r$balance <- CalculateBalance(dtaBef = subdta, dtaAfter = dtaAfter, trt = trt.col,
                                  cols = cols.balance)
    r$balance_plots <- GetBalancePlots(subdta, c(naive.match$pairs[, 9],
                                                 naive.match$pairs[, 10]))
  }
  r$pairs <- naive.match$pairs
  return(r)
}


GBMmodel <- function(subdta, trt.col, out.col, caliper, coord.cols,
                     cols.balance, interaction.depth = 3, seed = 1234,
                     matching_algorithm = c('optimal', 'greedy'),
                     remove_unmatchables = TRUE) {

  r <- NULL
  set.seed(seed)
  gbm.ps <- GBMPropScores(subdta, trt.col = trt.col, ignore.cols =
                            c(out.col, which(names(subdta) == 'prop.scores')),
                          interaction.depth = interaction.depth)
  
  if (matching_algorithm == 'greedy') {
    GBM.match <- PSmatchEst(subdta, pscores = gbm.ps, trt.col = trt.col,
                            out.col = out.col, SEreturn = TRUE, caliper = caliper,
                            estimand = 'ATT', pairsRet = TRUE, coord.cols = coord.cols)
  } else {
    GBM.match <- OptPSmatch(subdta, out.col = out.col, trt.col = trt.col,
                            pscores = gbm.ps, caliper = caliper, SEreturn = TRUE,
                            pairsRet = TRUE, coord.cols = coord.cols,
                            remove.unmatchables = remove_unmatchables)
  }
  
  r$result <- GBM.match$est + GBM.match$SE * 1.96 * c(- 1, 0, 1)
  if (is.null(dim(GBM.match$pairs))) {
    r$num_match <- 0
    r$distance <- Inf
    r$balance <- rep(NA, length(cols.balance))
  } else {
    r$num_match <- dim(GBM.match$pairs)[1]
    r$distance <- rdist.earth(GBM.match$pairs[, c(3, 4)], GBM.match$pairs[, c(7, 8)])
    r$distance <- diag(r$distance)
    r$distance <- mean(r$distance)
    dtaAfter <- subdta[c(GBM.match$pairs[, 9], GBM.match$pairs[, 10]), ]
    r$balance <- CalculateBalance(dtaBef = subdta, dtaAfter = dtaAfter, trt = trt.col,
                                  cols = cols.balance)
    r$balance_plots <- GetBalancePlots(subdta, c(GBM.match$pairs[, 9],
                                                 GBM.match$pairs[, 10]))
  }
  r$pairs <- GBM.match$pairs
  return(r)
}


DistCalModel <- function(subdta, caliper, dist.caliper, coord.cols, ignore.cols,
                         trt.col, out.col, cols.balance, coord_dist = TRUE,
                         matching_algorithm = c('optimal', 'greedy'),
                         remove_unmatchables = TRUE) {

  r <- NULL
  cal.match <- CaliperEst(dataset = subdta, ps.caliper = caliper,
                          dist.caliper = dist.caliper,
                          coords.columns = coord.cols,
                          ignore.cols = ignore.cols.coords,
                          SEreturn = TRUE, pairsRet = TRUE, trt.col = trt.col,
                          out.col = out.col, coord_dist = coord_dist,
                          matching_algorithm = matching_algorithm,
                          remove.unmatchables = remove_unmatchables)
  r$result <- cal.match$est[1] + 1.96 * c(-1, 0, 1) * cal.match$SE[1]
  r$num_match <- dim(cal.match$pairs)[1]
  D <- rdist.earth(cal.match$pairs[, c(3, 4)], cal.match$pairs[, c(7, 8)])
  r$distance <- mean(diag(D))
  r$balance <- CalculateBalance(dtaBef = subdta,
                                dtaAfter = subdta[c(cal.match$pairs[, 9],
                                                    cal.match$pairs[, 10])],
                                trt = trt.col, cols = cols.balance)
  r$balance_plots <- GetBalancePlots(subdta, c(cal.match$pairs[, 9],
                                               cal.match$pairs[, 10]))
  r$pairs <- cal.match$pairs
  
  return(r)
}



DAPSoptModel <- function(subdta, trt.col, out.col, coord.cols,
                         caliper, cols.balance, cutoff, coord_dist = TRUE) {

  r <- NULL
  out_name <- names(subdta)[out.col]
  
  # For optimal weight:
  DAPS.match <- DAPSest(subdta, out.col = out.col, trt.col = trt.col,
                        coords.columns = coord.cols, caliper = caliper,
                        weight = 'optimal', cov.cols = cols.balance,
                        w_tol = 0.001, caliper_type = 'DAPS',
                        cutoff = cutoff, pairsRet = TRUE,
                        coord_dist = coord_dist)
  lmod <- lm(as.formula(paste(out_name, '~ SnCR')),
             data = subdta[c(DAPS.match$ind_trt, DAPS.match$ind_cnt), ])
  r$result <- summary(lmod)$coef[2, 1] + summary(lmod)$coef[2, 2] * 1.96 * c(-1, 0, 1)
  r$num_match <- dim(DAPS.match$pairs)[1]
  r$distance <- mean(rdist(DAPS.match$pairs[, c(3, 4)], DAPS.match$pairs[, c(7, 8)]))
  r$balance <- DAPS.match$stand_diff
  r$weight <- DAPS.match$weight
  
  r$balance_plots <- GetBalancePlots(subdta, c(DAPS.match$pairs[, 9],
                                               DAPS.match$pairs[, 10]))
  r$pairs <- DAPS.match$pairs
  
  return(r)
}








