PredictHeatInput <- function(dat_unit, year, month, time_use) {
  
  variable <- 'Heat.Input..MMBtu.'
  
  time_pred <- paste0(year, '_', month)
  subdta <- subset(dat_unit, year_month %in% time_pred)
  p <- PredictVariableMonths(year, month, time_use, variable, dat_unit)
  
  # Assinging the new values (which include the predicted Heat input values).
  print(paste(sum(p$data < 0, na.rm = TRUE), 'out of', 
              sum(is.na(subdta$Heat.Input..MMBtu.)), 'predicted values are',
              'negative, and are forced to 0.'))
  subdta$Heat.Input..MMBtu. <- p$data
  subdta$Heat.Input..MMBtu.[subdta$Heat.Input..MMBtu. < 0] <- 0
  
  print('The R-squared of the models used was:')
  print(p$r_squared)
  return(subdta)
}