# Description: Using data from a different time period to predict a variable on
#              the month of interest.
# Date:        5/13/2016


library(data.table)


GetPredictData <- function(time_pred, time_use, variable, dat_unit) {
  # Get the data in a form that is used for predicting. Each row corresponds to
  # a unique unit, with variable information across time in the columns.
  #
  # Args:
  #  time_pred: The time (month of what year) we are interested in prediciting.
  #             should be of the form 'YYYY_M' or 'YYYY_MM'.
  #  time_use:  A vector of times in the above form to be used for prediction.
  #  variable:  The variable name we want to predict.
  #  dat_unit:  The data frame including the variables: 'uID' as the unique
  #             unit id, the variable in the argument variable, and 'year_month'
  #             describing what month the data refer to for each observation.
  #
  # Returns a data frame.
  
  dat_unit <- as.data.table(dat_unit)
  all_variables <- c('uID', variable, 'year_month')
  
  
  subdta <- subset(dat_unit, year_month %in% c(time_use, time_pred))
  subdta <- subdta[, all_variables, with = FALSE]
  setnames(subdta, variable, "Var")
  
  # Dropping units without an entry on time_pred.
  subdta[, entry_at_timepred := any(year_month == time_pred), by = uID]
  subdta <- subset(subdta, entry_at_timepred == 1)
  subdta[, entry_at_timepred := NULL]

  
  subdta <- reshape(subdta, idvar = "uID", timevar = "year_month", direction = "wide")
  X <- as.data.frame(subdta)[, - 1]
  rownames(X) <- subdta[, uID]
  
  wh <- which(colnames(X) == paste0('Var.', time_pred))
  X <- X[, c(wh, setdiff(1:ncol(X), wh))]
  
  return(X)
}


PredictivePower <- function(time_pred, time_use, variable, dat_unit) {
  # Function that checks whether we can predict a variable in our data using
  # the same variable from different time point.
  #
  # Args:
  #  time_pred: The time (month of what year) we are interested in prediciting.
  #             should be of the form 'YYYY_M' or 'YYYY_MM'.
  #  time_use:  A vector of times in the above form to be used for prediction.
  #  variable:  The variable name we want to predict.
  #  dat_unit:  The data frame including the variables: 'uID' as the unique
  #             unit id, the variable in the argument variable, and 'year_month'
  #             describing what month the data refer to for each observation.
  #
  # Returns:
  #  A list, including lm, linear model of the predicting variable on the rest,
  #  missing as the number of observations with at least one of the entries
  #  missing, num_obs as the total number of observations (including missing).
  #  Also, we plot the variable of time_pred against all time_use time points.
  
  X <- GetPredictData(time_pred, time_use, variable, dat_unit)
  form <- as.formula(paste(paste0('Var.', time_pred), ' ~ .'))
  lmod <- lm(form, data = X)
  
  r <- NULL
  r$lm <- lmod
  r$adj.rsquared <- summary(lmod)$adj.r
  r$missing <- sum(apply(X, 1, function(x) any(is.na(x))))
  missing <- which(is.na(X[, 1]))
  r$missing_pred <- length(missing)
  r$num_obs <- nrow(X)
  remain <- X[, - 1, drop = FALSE]
  r$cant_predict <- sum(apply(remain, 1, function(x) any(is.na(x)))[missing])
  for (ii in 1:length(time_use)) {
    plot(X[, ii + 1], X[, 1], xlab = colnames(X)[ii + 1],
         ylab = colnames(X)[1])
  }
  return(r)
}





PredictVariable <- function(time_pred, time_use, variable, dat_unit) {
  # Function that fits linear model using the time_use variable inputs, and
  # predicts the value of the variable at the missing locations of time_pred.
  # Multiple models can be used.
  #
  # Args:
  #  time_pred: The time (month of what year) we are interested in prediciting.
  #             should be of the form 'YYYY_M' or 'YYYY_MM'.
  #  time_use:  A list of models we want to fit. Each list element must be a
  #             vector of months in the above format.
  #  variable:  The variable name we want to predict.
  #  dat_unit:  The data frame including the variables: 'uID' as the unique
  #             unit id, the variable in the argument variable, and 'year_month'
  #             describing what month the data refer to for each observation.
  #
  # Returns:
  #  A list with 1) data frame where the unit ids are the rownames, and some
  #  values are imputed, 2) the number of imputed values for each model, 3)
  #  a list of imputed indeces per model.
  
  r <- NULL
  
  if (class(time_use) == 'character') {
    time_use <- list(time_use)
  }
  
  # How many each model predicted.
  predicted_num <- rep(NA, length(time_use))
  # List of predicted indeces.
  predicted <- NULL
  # Model r-squared
  R2 <- rep(NA, length(time_use))
  
  # What are the data that we used and on which we want to predict:
  X <- GetPredictData(time_pred, time_use[[1]], variable, dat_unit)
  row_names <- rownames(X)
  dta <- X[, 1]
  names(dta) <- row_names
  missing <- which(is.na(dta))
  
  # For each row of time_use we will predict the maximum amount of observations.
  for (ii in 1:length(time_use)) {
    
    # What are the data that we used and on which we want to predict:
    X <- GetPredictData(time_pred, time_use[[ii]], variable, dat_unit)
    
    # Using the function to get the fitted linear model.
    p <- PredictivePower(time_pred, time_use[[ii]], variable, dat_unit)
    
    R2[ii] <- p$adj.rsquared
    dta[missing] <- predict(p$lm, newdata = X[missing, - 1, drop = FALSE])
    predicted_num[ii] <- length(missing) - sum(is.na(dta))
    missing_new <- which(is.na(dta))
    predicted[[ii]] <- setdiff(missing, missing_new)
    missing <- missing_new
  }
  r$data <- dta
  r$predicted <- predicted
  r$predicted_num <- predicted_num
  r$r_squared <- R2
  return(r)
}



PredictVariableMonths <- function(year, month, time_use, variable, dat_unit) {
  # Function that fits linear model using the time_use variable inputs, and
  # predicts the value of the variable at the missing locations of time_pred.
  # Multiple models can be used. time_inter can include multiple time points.
  #
  # Args:
  #  year:      The year we want to predict the variable for. (Numeric)
  #  month:     The months within the year we want to predict. (Vector)
  #  time_use:  A list of years we want to use for every month we want to predict.
  #  variable:  The variable name we want to predict.
  #  dat_unit:  The data frame including the variables: 'uID' as the unique
  #             unit id, the variable in the argument variable, and 'year_month'
  #             describing what month the data refer to for each observation.
  #
  # Returns:
  #  A list with 1) data frame where the unit ids are the rownames, and some
  #  values are imputed, 2) the number of imputed values for each model, 3)
  #  a list of imputed indeces per model.
  
  time_pred <- paste0(year, '_', month)
  
  r <- NULL
  r$predicted_num <- array(NA, dim = c(length(time_pred), length(time_use)))
  dimnames(r$predicted_num) <- list(time_pred = time_pred, time_use = time_use)
  r$r_squared <- r$predicted_num
  col_ind <- which(names(subdta) == variable)
  subdta <- as.data.frame(subset(dat_unit,  year_month %in% time_pred))
  r$data <- subdta[, col_ind]
  
  # For each month we predict from all previous ones.
  for (ii in 1:length(time_pred)) {
    time_use_ii <- time_use
    wh <- which(subdta$Month == month[ii])
    
    for (jj in 1:length(time_use)) {
      time_use_ii[[jj]] <- paste0(time_use[[jj]], '_', month[ii])
    }
    
    # For every month we predict the variable
    p <- PredictVariable(time_pred[ii], time_use_ii, variable, dat_unit)
    r$predicted_num[ii, ] <- p$predicted_num
    r$r_squared[ii, ] <- p$r_squared
    
    # Sanity check using dimensions.
    length(p$data)
    dim(subdta[wh, ])
    # The next two numbers add to length(p$data) which is equal to nrow(subdta).
    mean(p$data == subdta$Heat.Input..MMBtu.[wh], na.rm = TRUE)
    sum(is.na(subdta$Heat.Input..MMBtu.[wh]))
    
    # Assinging the new values (which include the predicted Heat input values).
    r$data[wh] <- p$data
  }
  return(r)
}



