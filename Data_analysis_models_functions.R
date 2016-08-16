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
                       subsampling = FALSE, subsamples = NULL, cols.balance) {
  
  r <- NULL

  naive.match <- PSmatchEst(subdta, trt.col = trt.col, out.col = out.col,
                            SEreturn = TRUE, caliper = caliper,
                            replace = FALSE, estimand = 'ATT',
                            pairsRet = TRUE, coord.cols = coord.cols)
  r$result <- naive.match$est + naive.match$SE * 1.96 * c(-1, 0, 1)
  r$num_match <- dim(naive.match$pairs)[1]
  r$distance <- mean(rdist(naive.match$pairs[, c(3, 4)],
                         naive.match$pairs[, c(7, 8)]))
  r$balance <- CalculateBalance(dtaBef = subdta,
                                dtaAfter = subdta[c(naive.match$pairs[, 9],
                                                    naive.match$pairs[, 10])],
                                trt = trt.col, cols = cols.balance)
  r$balance_plots <- GetBalancePlots(subdta, c(naive.match$pairs[, 9],
                                               naive.match$pairs[, 10]))
  r$pairs <- naive.match$pairs
  
  if (subsampling) {
    b <- nrow(subdta) * 0.5
    # Subsampling for standard errors of the naive.
    match_sub <- numeric(subsamples)
    est_sub <- numeric(subsamples)
    for (ii in 1:subsamples) {
      D <- subdta[sample(1:nrow(subdta), b, replace = FALSE), ]
      naive_ps <- glm(as.formula(paste('SnCR ~ . - Fac.Latitude - Fac.Longitude - ',
                                       out_name)), data = D, family = 'binomial')
      D[, prop.scores := fitted(naive_ps)]
      naive.match <- PSmatchEst(D, trt.col = trt.col, out.col = out.col,
                                SEreturn = TRUE, caliper = caliper,
                                replace = FALSE, estimand = 'ATT',
                                pairsRet = TRUE, coord.cols = coord.cols)
      est_sub[ii] <- naive.match$est
      match_sub[ii] <- nrow(naive.match$pairs) * 2
    }
    se_naive <- sqrt(mean(match_sub) / r$num_match) * sd(est_sub)
    r$result_corSE <- r$result[2] + se_naive * 1.96 * c(-1, 0, 1)
  }
  return(r)
}



GBMmodel <- function(subdta, trt.col, out.col, caliper, coord.cols,
                     subsampling = FALSE, subsamples = NULL, cols.balance,
                     interaction.depth = 3, seed = 1234) {

  r <- NULL
  set.seed(seed)
  gbm.ps <- GBMPropScores(subdta, trt.col = trt.col, ignore.cols =
                            c(out.col, which(names(subdta) == 'prop.scores')),
                          interaction.depth = interaction.depth)
  GBM.match <- PSmatchEst(subdta, pscores = gbm.ps, trt.col = trt.col,
                          out.col = out.col, SEreturn = TRUE, caliper = caliper,
                          estimand = 'ATT', pairsRet = TRUE, coord.cols = coord.cols)
  r$result <- GBM.match$est + GBM.match$SE * 1.96 * c(-1, 0, 1)
  r$num_match <- dim(GBM.match$pairs)[1]
  r$distance <- mean(rdist(GBM.match$pairs[, c(3, 4)], GBM.match$pairs[, c(7, 8)]))
  r$balance <- CalculateBalance(dtaBef = subdta,
                                dtaAfter = subdta[c(GBM.match$pairs[, 9],
                                                    GBM.match$pairs[, 10])],
                                trt = trt.col, cols = cols.balance)
  r$balance_plots <- GetBalancePlots(subdta, c(GBM.match$pairs[, 9],
                                               GBM.match$pairs[, 10]))
  r$pairs <- GBM.match$pairs
  
  if (subsampling) {
    b <- nrow(subdta) * 0.5
    # Subsampling for GBM:
    # Subsampling for standard errors of the naive.
    match_sub <- numeric(subsamples)
    est_sub <- numeric(subsamples)
    for (ii in 1:subsamples) {
      D <- subdta[sample(1:nrow(subdta), b, replace = FALSE), ]
      gbm.ps <- GBMPropScores(D, trt.col = trt.col, ignore.cols =
                                c(out.col, which(names(subdta) == 'prop.scores')))
      GBM.match <- PSmatchEst(D, pscores = gbm.ps, trt.col = trt.col,
                              out.col = out.col, SEreturn = TRUE, caliper = caliper,
                              estimand = 'ATT', pairsRet = TRUE, coord.cols = coord.cols)
      est_sub[ii] <- GBM.match$est
      match_sub[ii] <- nrow(GBM.match$pairs) * 2
    }
    se_naive <- sqrt(mean(match_sub) / r$num_match) * sd(est_sub)
    r$result_corSE <- r$result[2] + se_naive * 1.96 * c(-1, 0, 1)
  }
  return(r)
}


DistCalModel <- function(subdta, caliper, dist.caliper, coord.cols, ignore.cols,
                         trt.col, out.col, subsampling = FALSE, subsamples = NULL,
                         cols.balance, coord_dist = TRUE) {

  r <- NULL
  cal.match <- CaliperEst(dataset = subdta, ps.caliper = caliper, dist.caliper = dist.caliper,
                          coords.columns = coord.cols, ignore.cols = ignore.cols.coords,
                          SEreturn = TRUE, pairsRet = TRUE, trt.col = trt.col,
                          out.col = out.col, coord_dist = coord_dist)
  r$result <- cal.match$est[1] + 1.96 * c(-1, 0, 1) * cal.match$SE[1]
  r$num_match <- dim(cal.match$pairs)[1]
  r$distance <- mean(rdist(cal.match$pairs[, c(3, 4)], cal.match$pairs[, c(7, 8)]))
  r$balance <- CalculateBalance(dtaBef = subdta,
                                dtaAfter = subdta[c(cal.match$pairs[, 9],
                                                    cal.match$pairs[, 10])],
                                trt = trt.col, cols = cols.balance)
  r$balance_plots <- GetBalancePlots(subdta, c(cal.match$pairs[, 9],
                                               cal.match$pairs[, 10]))
  r$pairs <- cal.match$pairs
  
  if (subsampling) {
    b <- nrow(subdta) * 0.5
    # Subsampling for matching in distance caliper.
    match_sub <- numeric(subsamples)
    est_sub <- numeric(subsamples)
    for (ii in 1:subsamples) {
      D <- subdta[sample(1:nrow(subdta), b, replace = FALSE), ]
      cal.match <- CaliperEst(dataset = D, ps.caliper = caliper, dist.caliper = dist.caliper,
                              coords.columns = coord.cols, ignore.cols = ignore.cols.coords,
                              SEreturn = TRUE, pairsRet = TRUE, trt.col = trt.col,
                              out.col = out.col, coord_dist = coord_dist)
      est_sub[ii] <- cal.match$est[1]
      match_sub[ii] <- nrow(cal.match$pairs) * 2
    }
    se_naive <- sqrt(mean(match_sub) / r$num_match) * sd(est_sub)
    r$result_corSE <- r$result[2] + se_naive * 1.96 * c(-1, 0, 1)
  }
  return(r)
}



DAPSoptModel <- function(subdta, trt.col, out.col, coord.cols,
                         caliper, cols.balance, cutoff, subsampling = FALSE,
                         subsamples = NULL, coord_dist = TRUE) {

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
  
  if (subsampling) {
    b <- nrow(subdta) * 0.5
    # Subsampling for DAPS with optimal weight.
    match_sub <- numeric(subsamples)
    est_sub <- numeric(subsamples)
    for (ii in 1:subsamples) {
      if (ii %% 50 == 0) print(ii)
      D <- subdta[sample(1:nrow(subdta), b, replace = FALSE), ]
      DAPS.match <- DAPSest(D, out.col = out.col, trt.col = trt.col,
                            weight = 'optimal',
                            coords.columns = coord.cols, caliper = caliper,
                            cov.cols = cols.balance, w_tol = 0.01, caliper_type = 'DAPS',
                            cutoff = cutoff, pairsRet = TRUE, quiet = TRUE,
                            coord_dist = coord_dist)
      est_sub[ii] <- DAPS.match$est
      match_sub[ii] <- nrow(DAPS.match$pairs) * 2
    }
    se_naive <- sqrt(mean(match_sub) / r$num_match) * sd(est_sub)
    r$result_corSE <- r$result[2] + se_naive * 1.96 * c(-1, 0, 1)
  }
  return(r)
}




