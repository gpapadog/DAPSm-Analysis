#' @param mean_na_cutoff A unit needs to have up to mean_na_cutoff percentage of
#' missingness in order to impute its values.
ImputeHeatInput <- function(dat_unit, year, month, method = c('kalman', 'seadec'),
                            mean_na_cutoff = 0.5, plot_imputation = FALSE) {
  
  method <- match.arg(method)
  impute_function <- na.kalman
  if (method == 'seadec') {
    impute_function <- na.seadec
  }
  
  missing_start_which <- which(is.na(dat_unit$Heat.Input..MMBtu.))
  missing_start <- length(missing_start_which)
  
  # Finding the unit ids with missing heat input.
  unique_uID <- with(dat_unit, unique(uID))
  missing_info <- matrix(0, nrow = length(unique_uID), ncol = 2)
  colnames(missing_info) <- c('n_obs', 'n_miss')
  rownames(missing_info) <- unique_uID
  dat_unit_year_month <- subset(dat_unit, Year == year & Month %in% month)
  
  for (ii in 1:length(unique_uID)) {
    D <- subset(dat_unit_year_month, uID == unique_uID[ii])
    missing_info[ii, 1] <- dim(D)[1]
    missing_info[ii, 2] <- sum(with(D, is.na(Heat.Input..MMBtu.)))
  }
  missing_info[, 2] <- missing_info[, 2] + length(month) - missing_info[, 1]
  
  print('Number of observations per unit.')
  print(table(missing_info[, 1]))
  print('Number of missing observations')
  print(table(missing_info[, 2]))
  
  uIDs_with_missing <- rownames(missing_info)[missing_info[, 2] > 0]
  
  time_points <- paste0(year, '_', month[c(1, length(month))])
  
  uID_imputed <- NULL
  
  for (ii in 1:length(uIDs_with_missing)) {
    D <- subset(dat_unit, uID == uIDs_with_missing[ii])
    D <- D[, c('Year', 'Month', 'year_month', 'Heat.Input..MMBtu.')]
    setorderv(D, c('Year', 'Month'))

    # We will not impute if we have too many missing data.
    impute <- (mean(is.na(D$Heat.Input..MMBtu.)) <= mean_na_cutoff)
    
    # We need to have equally distanced data.
    time_data <- table(D$Month)
    all_time_data <- all(as.numeric(names(time_data)) == 1:12)
    if (all_time_data) {
      all_time_data <- length(unique(time_data)) == 1
    }
    if (!all_time_data) {  # If data are not equally distanced, make them.
      D <- TSequallyDistanced(D)
    }
    wh_points <- which(D$year_month %in% time_points)
    
    if (length(time_data) > 6) {  # Operation was NOT only during some months.
      # Recalculate the percentage of missing.
      impute <- (mean(is.na(D$Heat.Input..MMBtu.)) <= mean_na_cutoff)
    }
      
    # We will only impute if there are observed data on either side.
    if (impute) {
      impute <- (any(!is.na(D$Heat.Input..MMBtu.[1 : wh_points[1]])) &
                   any(!is.na(D$Heat.Input..MMBtu.[- c(1 : wh_points[2])])))
    }

    if (impute) {
      uID_imputed <- c(uID_imputed, uIDs_with_missing[ii])
      imp <- impute_function(D$Heat.Input..MMBtu.)
      if (plot_imputation) {
        plotNA.distribution(D$Heat.Input..MMBtu.)
        plotNA.imputations(D$Heat.Input..MMBtu., imp)
        abline(v = wh_points, col = 'purple')
      }
      dat_unit[dat_unit$uID == uIDs_with_missing[ii]]$Heat.Input..MMBtu. <- imp
    }
  }
  
  missing_end_which <- which(is.na(dat_unit$Heat.Input..MMBtu.))
  missing_end <- length(missing_end_which)
  
  predicted <- setdiff(missing_start_which, missing_end_which)
  
  print(paste('Initial number of missing entries:', missing_start))
  print(paste('Number of imputed entries:', missing_start - missing_end))
  print(paste('Number of imputed less than 0, that are set to 0:',
              sum(dat_unit$Heat.Input..MMBtu.[predicted] < 0)))
  dat_unit$Heat.Input..MMBtu.[predicted] <- 0
  
  return(dat_unit)
}






TSequallyDistanced <- function(D) {
  min_year <- D$Year[1] # Since ordered.
  min_month <- D$Month[1]
  max_year <- D$Year[nrow(D)]
  max_month <- D$Month[nrow(D)]
  
  all_months <- paste0(rep(seq(min_year, max_year), each = 12), '_',
                       rep(1:12, max_year - min_year + 1))
  if (max_month < 12) {
    all_months <- all_months[- c((length(all_months) + max_month - 12 + 1) :
                                   length(all_months))]
  }
  if (min_month > 1) {
    all_months <- all_months[- c(1:(min_month - 1))]
  }
  all_months <- data.table(year_month = all_months)
  all_months[, Year := as.numeric(substr(year_month, 1, 4))]
  all_months[, Month := as.numeric(substr(year_month, 6, length(year_month)))]
  D <- merge(all_months, D, all.x = TRUE, by = 'year_month')
  D[, Year := Year.x]
  D[, Month := Month.x]
  D <- D[, c('Year', 'Month', 'year_month', 'Heat.Input..MMBtu.')]
  setorderv(D, c('Year', 'Month'))
  return(D)
}
