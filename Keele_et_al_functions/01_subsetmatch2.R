subsetmatch = function(dist_mat, t_ind, n_matches = 1,
                       mom_covs = NULL, mom_weights = NULL, mom_tols = NULL,
                       exact_covs = NULL, 
                       near_exact_covs = NULL, near_exact_devs = NULL,
                       fine_covs = NULL,
                       near_fine_covs = NULL, near_fine_devs = NULL,
                       subset_weight = NULL,
                       use_controls = NULL,
                       enforce_constraints = FALSE,
                       quiet = TRUE) {
  
  ks_covs = NULL
  ks_n_grid = 10
  ks_weights = NULL
  ks_tols = NULL
  
  dir_covs = NULL
  dir_deltas = NULL
  
  errorhandling(dist_mat, t_ind, n_matches,
                mom_covs, mom_weights, mom_tols,
                ks_covs, ks_n_grid, ks_weights, ks_tols,
                exact_covs, 
                near_exact_covs, near_exact_devs,
                fine_covs,
                near_fine_covs, near_fine_devs,
                use_controls,
                enforce_constraints)
  
  if (!quiet) {
    cat(format("  Generating the parameters..."), "\n")
  }
  prmtrs = problemparameters(dist_mat, t_ind, n_matches,
                             mom_covs, mom_weights, mom_tols,
                             ks_covs, ks_n_grid, ks_weights, ks_tols,
                             exact_covs, 
                             near_exact_covs, near_exact_devs,
                             fine_covs,
                             near_fine_covs, near_fine_devs,
                             dir_covs, dir_deltas,
                             subset_weight,
                             use_controls)
  n_t = prmtrs$n_t
  n_c = prmtrs$n_c
  n_matches = prmtrs$n_matches
  
  n_covs_tot = 0
  weights = rep(0, n_covs_tot)	
  
  if (!is.null(mom_weights) | !is.null(ks_weights)) {
    n_covs_tot = prmtrs$n_mom_covs+prmtrs$n_ks_covs
    if (prmtrs$n_mom_covs != 0 & prmtrs$n_ks_covs == 0) {
      weights = prmtrs$mom_weights
    }
    if (prmtrs$n_mom_covs == 0 & prmtrs$n_ks_covs != 0) {
      weights = prmtrs$ks_weights
    }					
    if (prmtrs$n_mom_covs != 0 & prmtrs$n_ks_covs != 0) {
      weights = c(prmtrs$mom_weights, prmtrs$ks_weights)
    }
  }
  
  cvec = prmtrs$cvec
  Amat = prmtrs$Amat
  bvec = prmtrs$bvec
  ub = prmtrs$ub 
  sense = prmtrs$sense
  vtype = prmtrs$vtype
  c_index = prmtrs$c_index
  
  if (!quiet) {
    cat(format("  Finding the optimal matches..."), "\n")
  }
  ptm = proc.time()
  out = Rcplex(cvec, Amat, bvec, ub = ub, sense = sense, vtype = vtype, n = 1,
               control = list(trace = 0, round = enforce_constraints*1))
  time = (proc.time()-ptm)[3]
  
  if (is.na(out$obj)) {
    cat(format("  Error: problem infeasible!"), "\n")
    obj_total = NA
    obj_dist_mat = NA
    obj_covs = NA
    c_id = NA
    group_id = NA
    time = NA
  }
  
  if (!is.na(out$obj)) {
    if (!quiet) {
      cat(format("  Optimal matches found"), "\n")
    }
    if (n_covs_tot==0) {			
      t_id = sort(rep(1:n_t, n_c))[out$xopt==1]
      c_id = (c_index+n_t)[out$xopt==1]	
    }
    if (n_covs_tot!=0) {			
      aux = (length(out$xopt)-(n_covs_tot-1)):length(out$xopt)
      t_id = sort(rep(1:n_t, n_c))[out$xopt[-aux]==1]
      c_id = (c_index+n_t)[out$xopt[-aux]==1]
    }
    
    group_id_t = 1:(length(t_id))
    group_id_c = sort(rep(1:length(t_id), n_matches))
    group_id = c(group_id_t, group_id_c)
    
    obj_total = out$obj
    obj_covs = 0
    
    if (!is.null(mom_weights) | !is.null(ks_weights)) {
      obj_covs = sum(weights*out$xopt[-(1:(length(out$xopt)-n_covs_tot))])
    }
    obj_dist_mat = obj_total-obj_covs
  }
  
  return(list(obj_total = obj_total, obj_dist_mat = obj_dist_mat, obj_covs = obj_covs, t_id = t_id, c_id = c_id, group_id = group_id, time = time))
  
}
