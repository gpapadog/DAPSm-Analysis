NewPredictHeatInput <- function(dat_unit, year, month, time_use) {
  
  time_use <- as.numeric(time_use)
  other_years <- subset(dat_unit, Year %in% time_use & Month %in% month)

  rsquared <- array(NA, dim = c(length(time_use), length(month)))
  dimnames(rsquared) <- list(use = time_use, pred = paste0(year, '_', month))
  num_pred <- rsquared
  
  r <- NULL
  r_data <- subset(dat, Year == year & Month %in% month)
  r$total <- r_data[, list(missing = length(Heat.Input..MMBtu.)), by = 'Month']
  r$missing <- r_data[, list(missing = sum(is.na(Heat.Input..MMBtu.))), by = 'Month']
  
  for (mm in 1:length(month)) {
    wh_month <- month[mm]
    for (uu in 1:length(time_use)) {
      
      dat <- rbind(r_data, other_years)
      
      wh_year <- time_use[uu]
      D <- subset(dat, Year %in% c(year, wh_year) & Month == wh_month)
  
      D <- D[, c('uID', 'Year', 'Heat.Input..MMBtu.')]
      D <- D[, list(var_use = Heat.Input..MMBtu.[Year == wh_year],
                    var_pred = Heat.Input..MMBtu.[Year == year]), by = uID]
      
      curr_missing <- sum(is.na(D$var_pred))
      cant_predict <- sum(is.na(D$var_use) & is.na(D$var_pred))
      
      D <- subset(D, !is.na(D$var_use))
      
      lmod <- lm(D$var_pred ~ D$var_use)
      rsquared[uu, mm] <- summary(lmod)$r.squared
      
      newdata <- data.frame(Int = 1, X = D$var_use)
      D[, prediction := predict(lmod, newdata = newdata)]
      D[, val := ifelse(!is.na(var_pred), var_pred, prediction)]
      D <- D[, c('uID', 'val')]
      D[, Month := wh_month]
    
      r_data <- merge(r_data, D, by = c('uID', 'Month'), all.x = TRUE)
      r_data[, Heat.Input..MMBtu. := ifelse(!is.na(Heat.Input..MMBtu.),
                                            Heat.Input..MMBtu., val)]
      r_data[, val := NULL]
      
      num_pred[uu, mm] <- curr_missing - cant_predict - sum(is.na(D$val))
    }
  }
  r$data <- r_data
  r$num_pred <- num_pred
  r$rsquared <- rsquared
  return(r)
}

