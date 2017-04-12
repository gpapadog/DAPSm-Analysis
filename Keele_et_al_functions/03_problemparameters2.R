problemparameters = function(dist_mat, t_ind, n_matches,
                             mom_covs, mom_weights, mom_tols,
                             ks_covs, ks_n_grid, ks_weights, ks_tols,
                             exact_covs, 
                             near_exact_covs, near_exact_dev,
                             fine_covs,
                             near_fine_covs, near_fine_dev,
                             dir_covs, dir_deltas,
                             subset_weight,
                             use_controls) {
  
  n_t = sum(t_ind)
  n_c = length(t_ind)-n_t
  
  n_mom_covs = 0
  if(!is.null(mom_covs)) {
    n_mom_covs = ncol(mom_covs)
  }
  
  n_ks_covs = 0
  if(!is.null(ks_covs)) {
    n_ks_covs = ncol(ks_covs)
  }
  
  n_dec_vars = (n_t*n_c)+n_mom_covs+n_ks_covs
  
  ks_covs_aux = NULL
  if (!is.null(ks_covs)) {	
    ks_grid = matrix(0, nrow = ks_n_grid, ncol = n_ks_covs)
    for (i in 1:n_ks_covs) {
      ks_covs_t_aux = ks_covs[, i][t_ind==1]
      ks_grid[, i] = quantile(ks_covs_t_aux, probs = seq(1/ks_n_grid, 1, 1/ks_n_grid))
    }	
    ks_covs_aux = matrix(0, nrow = length(t_ind), ncol = ks_n_grid*n_ks_covs)
    for (i in 1:n_ks_covs) {
      k = (i-1)*ks_n_grid
      for (j in 1:ks_n_grid) {
        ks_covs_aux[, j+k][ks_covs[, i]<ks_grid[j, i]] = 1
      }
    }
  }
  
  if (is.null(subset_weight)) {
    cvec = as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE))
  }
  if (!is.null(subset_weight)) {
    cvec = as.vector(matrix(t(dist_mat), nrow = 1, byrow = TRUE)) -
      (subset_weight * rep(1, n_t*n_c))
  }
  
  constraintmat_out = constraintmatrix(t_ind, n_matches,
                                       mom_covs, mom_weights, mom_tols,
                                       ks_covs, ks_covs_aux, ks_n_grid, ks_weights, ks_tols,
                                       exact_covs,
                                       near_exact_covs, near_exact_dev,
                                       fine_covs,
                                       near_fine_covs, near_fine_dev,
                                       dir_covs, dir_deltas,
                                       subset_weight,
                                       use_controls)
  
  cnstrn_mat = constraintmat_out$cnstrn_mat
  bvec_7 = constraintmat_out$bvec_7
  bvec_8 = constraintmat_out$bvec_8
  
  bvec = c(rep(n_matches, n_t), rep(1, n_c))
  if (!is.null(mom_covs)) {
    bvec = c(bvec, rep(0, 2*n_mom_covs))	
  }
  if (!is.null(mom_covs)) {
    bvec = c(bvec, rep(0, 2*n_ks_covs*ks_n_grid))
  }
  
  if (!is.null(exact_covs)) {
    bvec = c(bvec, rep(0, ncol(exact_covs))) 
  }	
  if (!is.null(near_exact_covs)) {
    if (length(near_exact_devs)==1) {
      bvec = c(bvec, near_exact_dev) 
    }
    if (length(near_exact_devs)>1) {
      bvec = c(bvec, near_exact_devs+1)
    } 
  }
  if (!is.null(fine_covs)) {
    bvec = c(bvec, bvec_7) 
  }
  if (!is.null(near_fine_covs)) {
    bvec = c(bvec, bvec_8) 
  }
  if (!is.null(dir_covs)) {
    bvec = c(bvec, rep(0, ncol(dir_covs))) 
  }	
  if (!is.null(use_controls)) {
    bvec = c(bvec, sum(use_controls)) 
  }
  
  ub = c(rep(1, n_t*n_c), rep(Inf, n_mom_covs), rep(Inf, n_ks_covs))
  if (!is.null(mom_tols) | !is.null(ks_tols)) { 
    ub = rep(1, n_t*n_c)
  }					 						 
  sense = c(rep("L", n_t), rep("L", n_c), rep("L", 2*n_mom_covs), rep("L", 2*n_ks_covs*ks_n_grid))
  if (!is.null(exact_covs)) {
    sense = c(sense, rep("E", ncol(exact_covs))) 
  }
  if (!is.null(near_exact_covs)) {
    sense = c(sense, rep("L", ncol(near_exact_covs))) 
  }
  if (!is.null(fine_covs)) {
    sense = c(sense, rep("E", length(bvec_7))) 
  }
  if (!is.null(near_fine_covs)) {
    sense = c(sense, rep(c("L"), length(bvec_8))) 
  }
  if (!is.null(dir_covs)) {
    sense = c(sense, rep("E", ncol(dir_covs))) 
  }
  if (!is.null(use_controls)) {
    sense = c(sense, "E") 
  }
  
  vtype = c(rep("B", n_t*n_c), rep("C", n_mom_covs), rep("C", n_ks_covs))
  if (!is.null(mom_tols) | !is.null(ks_tols)) {					 
    vtype = rep("B", n_t*n_c)
  }						 
  
  c_index = rep(1:n_c, n_t)	
  
  return(list(n_t = n_t, n_c = n_c, 
              n_matches = n_matches, 
              n_mom_covs = n_mom_covs, mom_weights = mom_weights, 
              n_ks_covs = n_ks_covs, ks_n_grid = ks_n_grid, ks_weights = ks_weights, 
              cvec = cvec, 
              Amat = cnstrn_mat, 
              bvec = bvec, 
              ub = ub, 
              sense = sense,
              vtype = vtype, 
              c_index = c_index))	
  
}
