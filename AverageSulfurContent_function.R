AverageSulfurContent <- function(dat_unit) {
  
  dat_unit[, Sulfur.Content := as.numeric(Sulfur.Content)]
  
  A = dat_unit[, list(FacID = FacID[1],
                      Initial.Year.of.Operation = Initial.Year.of.Operation[1],
                      AverSulfContent = mean(Sulfur.Content, na.rm = TRUE),
                      Facility.Latitude.x = Facility.Latitude.x[1],
                      Facility.Longitude.x = Facility.Longitude.x[1],
                      SO2..tons. = SO2..tons.[1],
                      NOx..tons. = NOx..tons.[1],
                      Heat.Input..MMBtu. = Heat.Input..MMBtu.[1],
                      Operating.Time = Operating.Time[1],
                      Status = Status[1],
                      Fuel1.IsCoal = Fuel1.IsCoal[1],
                      NOx.Scrub1 = NOx.Scrub1[1],
                      NOx.Scrub2 = NOx.Scrub2[1],
                      NOx.Scrub3 = NOx.Scrub3[1],
                      NOx.Scrub4 = NOx.Scrub4[1],
                      Is.Phase1 = Is.Phase1[1],
                      year_month = year_month[1]),
               by = c("uID", 'Year', 'Month')]
  print('Aggregation completed.')
  
  return(A)
}