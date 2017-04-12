keele_match <- function(dataset, dist_mat, t_ind, exact_covs = NULL, mom_covs = NULL,
                        mom_tols = NULL, near_fine_covs = NULL, near_fine_devs = NULL,
                        subset_weight, true_value = NULL, pairsRet = FALSE,
                        coords.cols = NULL) {
  
  keele_match <- subsetmatch(dist_mat = dist_mat, t_ind = t_ind,
                             exact_covs = exact_covs,
                             mom_covs = mom_covs, mom_tols = mom_tols,
                             near_fine_covs = near_fine_covs,
                             near_fine_devs = near_fine_devs,
                             subset_weight = subset_weight)
  
  t_id <- keele_match$t_id
  c_id = keele_match$c_id
  
  matched_data <- dataset[c(t_id, c_id), ]
  lmod <- lm(Y ~ X, data = matched_data)
  
  r <- NULL
  r$est <- lmod$coef[2]
  r$SE <- summary(lmod)$coef[2, 2]
  if (!is.null(true_value)) {
    r$cover <- (abs(true_value - r$est) < qnorm(0.975) * r$SE)
  }
  
  if (pairsRet) {
    which_cols <- c(which(names(dataset) %in% c('X', 'Y', 'prop.scores')))
    which_cols <- c(which_cols, coords.cols)
    
    pairs <- matched_data[, which_cols]
    pairs <- cbind(pairs[1:length(t_id), ], pairs[- c(1:length(t_id)), ])
    names(pairs) <- paste0(rep(c('Trt.', 'Con.'),
                               each = length(which_cols)), names(pairs))
    pairs$IDtrt <- t_id
    pairs$IDcnt <- c_id
    r$pairs <- as.matrix(pairs[, c(1, 3:6, 8:12, 2, 7)])
  }
  return(r)
}