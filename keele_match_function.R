keele_match <- function(dta, trt_col, out_col, coords.columns, exact_covs = NULL,
                        mom_covs = NULL, mom_tols = NULL, near_fine_covs = NULL,
                        near_fine_devs = NULL, subset_weight, true_value = NULL,
                        pairsRet = FALSE, use_controls = NULL, cols.balance = NULL,
                        enforce_constraints = FALSE, n_matches = 1, miles = TRUE,
                        distance_measure = c('geo', 'euclidean')) {
  
  require(Rcplex)
  require(fields)

  distance_measure <- match.arg(distance_measure)
  
  dta <- as.data.frame(dta)
  dataset <- dta[order(dta[, trt_col], decreasing = TRUE), ]
  
  t_ind <- as.numeric(dataset[, trt_col])
  coords <- cbind(dataset[, coords.columns])
  if (distance_measure == 'geo') {
    dist_mat <- rdist.earth(coords[t_ind == 1, ], coords[t_ind == 0, ],
                            miles = miles)
  } else {
    dist_mat <- rdist(coords[t_ind == 1, ], coords[t_ind == 0, ])
  }

  matched <- subsetmatch(dist_mat = dist_mat, t_ind = t_ind,
                         exact_covs = exact_covs, mom_covs = mom_covs,
                         mom_tols = mom_tols, near_fine_covs = near_fine_covs,
                         near_fine_devs = near_fine_devs,
                         subset_weight = subset_weight,
                         use_controls = use_controls, n_matches = n_matches,
                         enforce_constraints = enforce_constraints)
  
  t_id <- matched$t_id
  c_id <- matched$c_id
  
  matched_data <- dataset[c(t_id, c_id), ]
  names(matched_data)[c(out_col, trt_col)] <- c('Y', 'X')
  lmod <- lm(Y ~ X, data = matched_data)
  
  r <- NULL
  r$est <- lmod$coef[2]
  r$SE <- summary(lmod)$coef[2, 2]
  r$CI <- r$est + c(- 1, 0, 1) * r$SE * qnorm(0.975)
  
  if (!is.null(true_value)) {
    r$cover <- (abs(true_value - r$est) < qnorm(0.975) * r$SE)
  }
  
  if (pairsRet) {
    which_cols <- c(which(names(matched_data) %in% c('X', 'Y', 'prop.scores')))
    which_cols <- c(which_cols, coords.columns)
    
    pairs <- matched_data[, which_cols]
    pairs <- cbind(pairs[1:length(t_id), ], pairs[- c(1:length(t_id)), ])
    names(pairs) <- paste0(rep(c('Trt.', 'Con.'),
                               each = length(which_cols)), names(pairs))
    pairs$IDtrt <- t_id
    pairs$IDcnt <- c_id
    r$pairs <- as.matrix(pairs[, c(1, 3:6, 8:12, 2, 7)])
    
    if (distance_measure == 'geo') {
      r$distance <- diag(rdist.earth(r$pairs[, c(3, 4)], r$pairs[, c(7, 8)],
                                     miles = miles))
    } else {
      r$distance <- diag(rdist(r$pairs[, c(3, 4)], r$pairs[, c(7, 8)]))
    }
    
    r$num_match <- length(r$distance)
    r$balance <- CalculateBalance(dtaBef = dataset, dtaAfter = matched_data,
                                  trt = trt_col, cols = cols.balance)
    r$diff_means <- CalculateBalance(dtaBef = dataset, dtaAfter = matched_data,
                                     trt = trt_col, cols = cols.balance,
                                     diff_means = TRUE)
    
  }
  
  return(r)
}

