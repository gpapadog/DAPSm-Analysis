AverageSulfurContent <- function(dat_unit) {
  
  dat_unit[, Sulfur.Content := as.numeric(Sulfur.Content)]
  
  A = dat_unit[, list(FacID = FacID[1],
                      Unit.ID = Unit.ID[1],
                      Initial.Year.of.Operation = Initial.Year.of.Operation[1],
                      Sulfur.Content = mean(Sulfur.Content, na.rm = TRUE),
                      Facility.Latitude.x = Facility.Latitude.x[1],
                      Facility.Longitude.x = Facility.Longitude.x[1],
                      SO2..tons. = SO2..tons.[1],
                      NOx..tons. = NOx..tons.[1],
                      CO2..short.tons. = CO2..short.tons.[1],
                      Gross.Load..MW.h. = Gross.Load..MW.h.[1],
                      FIPS = FIPS[1],
                      Heat.Input..MMBtu. = Heat.Input..MMBtu.[1],
                      Operating.Time = Operating.Time[1],
                      Status = Status[1],
                      Max.Hourly.HI.Rate..MMBtu.hr. = Max.Hourly.HI.Rate..MMBtu.hr.[1],
                      Fuel1.IsCoal = Fuel1.IsCoal[1],
                      Fuel.Type..Primary..x = Fuel.Type..Primary..x[1],
                      NOx.Scrub1 = NOx.Scrub1[1],
                      NOx.Scrub2 = NOx.Scrub2[1],
                      NOx.Scrub3 = NOx.Scrub3[1],
                      NOx.Scrub4 = NOx.Scrub4[1],
                      Is.Phase1 = Is.Phase1[1],
                      Is.Phase2 = Is.Phase2[1],
                      year_month = year_month[1]),
               by = c("uID", 'Year', 'Month')]
  print('Aggregation completed.')
  
  return(A)
}