constraintmatrix = function(t_ind, n_matches,
                            mom_covs, mom_weights, mom_tols,
                            ks_covs, ks_covs_aux, ks_n_grid, ks_weights, ks_tols,
                            exact_covs, 
                            near_exact_covs, near_exact_dev,
                            fine_covs,
                            near_fine_covs, near_fine_dev,
                            dir_covs, dir_deltas,
                            subset_weight,
                            use_controls) {
  
  n_t = sum(t_ind)
  n_c = length(t_ind)-n_t	
  
  n_tot = n_t*n_c
  
  row_ind_1 = sort(rep(1:n_t, n_c))
  col_ind_1 = 1:n_tot
  ones_1 = rep(1, n_tot)
  
  row_ind_2 = sort(rep(1:n_c, n_t))+n_t
  col_ind_2 = rep(seq(1, n_t*n_c, n_c), n_c)+(sort(rep(1:n_c, n_t))-1)
  ones_2 = rep(1, n_tot)
  row_ind_cur	= max(row_ind_2)
  
  mom_ks_covs = NULL
  if (!is.null(mom_covs) | !is.null(ks_covs)) {
    row_ind_3.4 = 0
    n_mom_covs = 0
    if(!is.null(mom_covs)) {
      n_mom_covs = ncol(mom_covs)
      
    }
    n_ks_covs = 0
    if(!is.null(ks_covs)) {
      n_ks_covs = ncol(ks_covs)
    }
    if(!is.null(mom_covs) & is.null(ks_covs_aux)) {
      mom_ks_covs = mom_covs
      mom_ks_tols = mom_tols
    }
    if(is.null(mom_covs) & !is.null(ks_covs_aux)) {
      mom_ks_covs = ks_covs_aux
      mom_ks_tols = NA
      for (i in 1:ncol(ks_covs)) {
        mom_ks_tols = c(mom_ks_tols, rep(ks_tols[i], ks_n_grid))
      }
      mom_ks_tols = mom_ks_tols[-1]
    }
    if(!is.null(mom_covs) & !is.null(ks_covs_aux)) {
      mom_ks_covs = cbind(mom_covs, ks_covs_aux)
      mom_ks_tols = mom_tols
      for (i in 1:ncol(ks_covs)) {
        mom_ks_tols = c(mom_ks_tols, rep(ks_tols[i], ks_n_grid))
      }
    }	
  }					 
  if (!is.null(mom_ks_covs)) {
    n_mom_ks_covs = ncol(mom_ks_covs)
    if (!is.null(mom_weights) | !is.null(ks_weights)) {
      row_ind_3.4 = sort(rep(1:(2*n_mom_ks_covs)+n_t+n_c, n_tot+1))
    }	
    if (!is.null(mom_tols) | !is.null(ks_tols)) {
      row_ind_3.4 = sort(rep(1:(2*n_mom_ks_covs)+n_t+n_c, n_tot))
    }	
    col_ind_3.4 = NA
    mom_ks_vals_3.4 = NA
    j = 1
    k = 0
    for (i in 1:n_mom_ks_covs) {
      if (n_mom_covs != 0 & i <= n_mom_covs) {	
        if (!is.null(mom_weights) | !is.null(ks_weights)) {
          col_ind_3.4 = c(col_ind_3.4, rep(c(1:n_tot, n_tot+i), 2))
        }	
        if (!is.null(mom_tols) | !is.null(ks_tols)) {
          col_ind_3.4 = c(col_ind_3.4, rep(1:n_tot, 2))
        }
      }
      if (n_ks_covs != 0 & i > n_mom_covs) {	
        if (!is.null(mom_weights) | !is.null(ks_weights)) {
          col_ind_3.4 = c(col_ind_3.4, rep(c(1:n_tot, n_tot+n_mom_covs+j), 2))
          k = k+1
          if (k >= ks_n_grid) {
            j = j+1
            k = 0	
          }
        }	
        if (!is.null(mom_tols) | !is.null(ks_tols)) {
          col_ind_3.4 = c(col_ind_3.4, rep(1:n_tot, 2))
          k = k+1
          if (k >= ks_n_grid) {
            j = j+1
            k = 0	
          }
        }
      }
      temp_mean_1 = rep(mom_ks_covs[t_ind==0, i], n_t)-(mom_ks_covs[t_ind==1, i])[sort(rep(1:n_t, n_c))]
      if (!is.null(mom_weights) | !is.null(ks_weights)) {
        temp_mean_2 = c(rep(temp_mean_1, n_t), -1)
        temp_mean_3 = c(-rep(temp_mean_1, n_t), -1)
      }	
      if (!is.null(mom_tols) | !is.null(ks_tols)) {
        temp_mean_2 = temp_mean_1-(mom_ks_tols[i]*rep(1, n_t*n_c))
        temp_mean_3 = -temp_mean_1-(mom_ks_tols[i]*rep(1, n_t*n_c))
      }	
      mom_ks_vals_3.4 = c(mom_ks_vals_3.4, temp_mean_2, temp_mean_3)
      if (i == 1) {
        col_ind_3.4 = col_ind_3.4[-1]
        mom_ks_vals_3.4 = mom_ks_vals_3.4[-1]
      }
    }
    row_ind_cur	= max(row_ind_3.4)
  }	
  
  rows_exact = NULL
  cols_exact = NULL
  vals_exact = NULL
  if (!is.null(exact_covs)) {
    n_exact_cats = ncol(exact_covs)
    j = 1
    for (i in 1:n_exact_cats) {
      rows_exact = c(rows_exact, rep(row_ind_cur+j, n_t*n_c))
      cols_exact = c(cols_exact, 1:(n_t*n_c))
      dist_exact_cov = abs(outer(exact_covs[t_ind==1, i], exact_covs[t_ind==0, i], "-"))
      dist_exact_cov = t(dist_exact_cov)
      vals_exact = c(vals_exact, as.vector(dist_exact_cov))
      if (j == 1) {
        rows_exact = rows_exact[-1]
        cols_exact = cols_exact[-1]
        vals_exact = vals_exact[-1]
      }
      j = j+1
    }	
    row_ind_5 = rows_exact
    col_ind_5 = cols_exact
    exact_vals_5 = vals_exact
    row_ind_cur	= max(row_ind_5)
  }	
  
  rows_near_exact = NULL
  cols_near_exact = NULL
  vals_near_exact = NULL
  if (!is.null(near_exact_covs)) {
    n_near_exact_cats = ncol(near_exact_covs)
    j = 1
    for (i in 1:n_near_exact_cats) {
      rows_near_exact = c(rows_near_exact, rep(row_ind_cur+j, n_t*n_c))
      cols_near_exact = c(cols_near_exact, 1:(n_t*n_c))
      dist_near_exact_cov = abs(outer(near_exact_covs[t_ind==1, i], near_exact_covs[t_ind==0, i], "-"))
      dist_near_exact_cov = t(dist_near_exact_cov)
      vals_near_exact = c(vals_near_exact, as.vector(dist_near_exact_cov))
      j = j+1
    }	
    row_ind_6 = rows_near_exact
    col_ind_6 = cols_near_exact
    near_exact_vals_6 = vals_near_exact
    row_ind_cur	= max(row_ind_6)
  }
  
  bvec_7 = NA
  rows_fine = NULL
  cols_fine = NULL
  vals_fine = NULL
  if (!is.null(fine_covs)) {
    fine_covs_2 = rep(NA, nrow(fine_covs))
    n_fine_covs = ncol(fine_covs)
    j = 1
    for (i in 1:n_fine_covs) {	
      aux = factor(fine_covs[, i])
      fine_covs_2 = cbind(fine_covs_2, diag(nlevels(aux))[aux,])
      if (j == 1) {
        fine_covs_2 = fine_covs_2[, -1]
      }
      j = j+1
    }
    n_fine_cats = ncol(fine_covs_2)
    j = 1
    for (i in 1:n_fine_cats) {
      rows_fine = c(rows_fine, rep(row_ind_cur+j, n_t*n_c))
      cols_fine = c(cols_fine, 1:(n_t*n_c))
      dist_fine_cov = outer(fine_covs_2[t_ind==1, i], fine_covs_2[t_ind==0, i], "-")
      dist_fine_cov = t(dist_fine_cov)
      vals_fine = c(vals_fine, as.vector(dist_fine_cov))
      if (j == 1) {
        rows_fine = rows_fine[-1]
        cols_fine = cols_fine[-1]
        vals_fine = vals_fine[-1]
      }
      j = j+1
    }	
    row_ind_7 = rows_fine
    col_ind_7 = cols_fine
    fine_vals_7 = vals_fine
    bvec_7 = rep(0, n_fine_cats)
    row_ind_cur	= max(row_ind_7)
  }
  
  bvec_8 = NA
  rows_near_fine = NULL
  cols_near_fine = NULL
  vals_near_fine = NULL
  vec_n_cats_near_fine_covs = NA
  if (!is.null(near_fine_covs)) {	
    near_fine_covs_2 = rep(NA, nrow(near_fine_covs))
    n_near_fine_covs = ncol(near_fine_covs)
    j = 1
    for (i in 1:n_near_fine_covs) {	
      
      cats_near_fine_cov = as.numeric(names(table(near_fine_covs[, i])))
      vec_n_cats_near_fine_covs = c(vec_n_cats_near_fine_covs, length(cats_near_fine_cov))
      
      aux = factor(near_fine_covs[, i])
      near_fine_covs_2 = cbind(near_fine_covs_2, diag(nlevels(aux))[aux,])
      if (j == 1) {
        near_fine_covs_2 = near_fine_covs_2[, -1]
      }
      j = j+1
    }
    n_near_fine_cats = ncol(near_fine_covs_2)
    j = 1
    for (i in 1:n_near_fine_cats) {
      rows_near_fine = c(rows_near_fine, rep(row_ind_cur+j, n_t*n_c))
      cols_near_fine = c(cols_near_fine, 1:(n_t*n_c))
      dist_near_fine_cov = outer(near_fine_covs_2[t_ind==1, i], near_fine_covs_2[t_ind==0, i], "-")
      dist_near_fine_cov = t(dist_near_fine_cov)
      vals_near_fine = c(vals_near_fine, as.vector(dist_near_fine_cov))
      if (j == 1) {
        rows_near_fine = rows_near_fine[-1]
        cols_near_fine = cols_near_fine[-1]
        vals_near_fine = vals_near_fine[-1]
      }
      j = j+1
    }	
    row_ind_8 = rows_near_fine
    col_ind_8 = cols_near_fine
    near_fine_vals_8 = vals_near_fine
    
    vec_n_cats_near_fine_covs = vec_n_cats_near_fine_covs[-1]
    
    if (length(near_fine_devs) == 1) {
      bvec_8 = rep(near_fine_devs, n_near_fine_cats)
    }
    if (length(near_fine_devs) > 1) {
      bvec_8 = rep(near_fine_devs, vec_n_cats_near_fine_covs)
    }
    row_ind_cur	= max(row_ind_8)
  }
  
  if (!is.null(dir_covs)) {
    row_ind_a = NA
    col_ind_a = NA
    dir_covs_vals_a = NA
    for (i in 1:ncol(dir_covs)) {
      aux = outer(dir_covs[t_ind==1, i], dir_covs[t_ind==0, i], "-")
      temp = aux
      temp[aux>(dir_deltas[i])] = 1
      temp[aux<=(dir_deltas[i]) & aux>=(-dir_deltas[i])] = 0
      temp[aux<(-dir_deltas[i])] = -1
      temp = as.vector(t(temp))
      dir_covs_vals_a = c(dir_covs_vals_a, temp)
      row_ind_a = c(row_ind_a, rep(row_ind_cur+i, n_tot))
      col_ind_a = c(col_ind_a, 1:n_tot)
    }
    dir_covs_vals_a = dir_covs_vals_a[-1]
    row_ind_a = row_ind_a[-1]
    col_ind_a = col_ind_a[-1]
    row_ind_cur	= max(row_ind_a)
  }	
  
  if (!is.null(use_controls)) {
    use_controls = use_controls[(n_t+1):(n_t+n_c)]
    use_controls_aux = rep(use_controls, n_t)
    col_ind_9 = (1:n_tot)[use_controls_aux==1]
    row_ind_9 = rep(row_ind_cur+1, length(col_ind_9))
    use_controls_vals_9 = rep(1, length(col_ind_9))
  }
  
  row_ind = c(row_ind_1, row_ind_2)
  col_ind = c(col_ind_1, col_ind_2)
  vals = c(ones_1, ones_2)
  if (!is.null(mom_ks_covs)) {
    row_ind = c(row_ind, row_ind_3.4)
    col_ind = c(col_ind, col_ind_3.4)
    vals = c(vals, mom_ks_vals_3.4)
  }
  if (!is.null(exact_covs)) {
    row_ind = c(row_ind, row_ind_5)
    col_ind = c(col_ind, col_ind_5)
    vals = c(vals, exact_vals_5)
  }
  if (!is.null(near_exact_covs)) {
    row_ind = c(row_ind, row_ind_6)
    col_ind = c(col_ind, col_ind_6)
    vals = c(vals, near_exact_vals_6)
  }
  if (!is.null(fine_covs)) {
    row_ind = c(row_ind, row_ind_7)
    col_ind = c(col_ind, col_ind_7)
    vals = c(vals, fine_vals_7)
  }
  if (!is.null(near_fine_covs)) {
    row_ind = c(row_ind, row_ind_8)
    col_ind = c(col_ind, col_ind_8)
    vals = c(vals, near_fine_vals_8)
  }
  if (!is.null(dir_covs)) {
    row_ind = c(row_ind, row_ind_a)
    col_ind = c(col_ind, col_ind_a)
    vals = c(vals, dir_covs_vals_a)	
  }	
  if (!is.null(use_controls)) {
    row_ind = c(row_ind, row_ind_9)
    col_ind = c(col_ind, col_ind_9)
    vals = c(vals, use_controls_vals_9)	
  }
  
  aux = cbind(row_ind, col_ind, vals)[order(col_ind), ]
  cnstrn_mat = simple_triplet_matrix(i = aux[, 1], j = aux[, 2], v = aux[, 3])
  
  return(list(cnstrn_mat = cnstrn_mat, bvec_7 = bvec_7, bvec_8 = bvec_8))
  
}


