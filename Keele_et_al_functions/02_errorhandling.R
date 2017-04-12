#! Rank based Mahalanobis distance, from Paul Rosenbaum's Design of Observational Studies, p. 251
errorhandling = function(dist_mat, t_ind, n_matches,
                         mom_covs, mom_weights, mom_tols,
                         ks_covs, ks_n_grid, ks_weights, ks_tols,
                         exact_covs, 
                         near_exact_covs, near_exact_devs,
                         fine_covs,
                         near_fine_covs, near_fine_devs,
                         use_controls,
                         enforce_constraints) { 
  
  t_ind_sort = t_ind[order(t_ind, decreasing = TRUE)]
  if (!identical(t_ind, t_ind_sort)) {
    stop("The data needs to be sorted in decreasing order by the treatment indicator")
  }
  
  if (!is.null(mom_covs)) {
    if ((is.null(mom_weights) & is.null(mom_tols)) | (!is.null(mom_weights) & !is.null(mom_tols))) {
      stop("With mom_covs either mom_weights or mom_tols needs to be specified")
    }
  }
  if (!is.null(mom_covs)) {
    if (!is.null(mom_weights)) {
      if (ncol(mom_covs)!=length(mom_weights)) {
        stop("The number of columns in mom_covs needs to be equal to the length of mom_weights")
      }
    }	
    if (!is.null(mom_tols)) {
      if (ncol(mom_covs)!=length(mom_tols)) {
        stop("The number of columns in mom_covs needs to be equal to the length of mom_tols")
      }
    }
  }
  
  if (!is.null(ks_covs)) {
    if ((is.null(ks_weights) & is.null(ks_tols)) | (!is.null(ks_weights) & !is.null(ks_tols))) {
      stop("With ks_covs either ks_weights or ks_tols needs to be specified")
    }
  }
  if (!is.null(ks_covs)) {
    if (!is.null(ks_weights)) {
      if (ncol(ks_covs)!=length(ks_weights)) {
        stop("The number of columns in ks_covs needs to be equal to the length of ks_weights")
      }
    }
    if (!is.null(ks_tols)) {
      if (ncol(ks_covs)!=length(ks_tols)) {
        stop("The number of columns in ks_covs needs to be equal to the length of ks_tols")
      }
    }
  }
  
  if (!is.null(near_exact_covs)) {
    if (!is.null(near_exact_devs)) {
      if (length(near_exact_devs)!=1) {
        if (ncol(near_exact_covs)!=length(near_exact_devs)) {
          stop("If different to 1, the length of near_exact_devs has to be equal to the number of columns in near_exact_covs")
        }	
      }
    }
  }
  
  if (!is.null(near_fine_covs)) {
    if (!is.null(near_fine_devs)) {
      if (length(near_fine_devs)!=1) {
        if (ncol(near_fine_covs)!=length(near_fine_devs)) {
          stop("If different to 1, the length of near_fine_devs has to be equal to the number of columns in near_fine_covs")
        }	
      }
    }
  }
  
}
