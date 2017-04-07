#' Function that takes the whole data set with power plant information and drops units
#' and facilities that were either retired, not operating, or not not have any data
#' during the study period.
CleanPPunits <- function(dat_unit, year, month) {
  
  # Dropping units that were retired during the time period.
  retired_in_year <- subset(dat_unit, Year == year & Month %in% month & Is.Retired)
  drop_uID <- unique(retired_in_year$uID)
  print(paste(length(drop_uID), 'units in', length(unique(retired_in_year$FacID)),
              'facility were retired.'))
  dat_unit <- subset(dat_unit, ! (uID %in% drop_uID))
  rm(retired_in_year)
  rm(drop_uID)
  
  # Dropping units that had not started operating during the period.
  not_operating <- subset(dat_unit, Initial.Year.of.Operation > year)
  drop_uID <- unique(not_operating$uID)
  print(paste(length(drop_uID), 'units in', length(unique(not_operating$FacID)),
              'facility were not yet operating'))
  dat_unit <- subset(dat_unit, ! (uID %in% drop_uID))
  rm(not_operating)
  rm(drop_uID)
  
  # Dropping facilities that do not have data during the time period.
  facilities_in_year <- subset(dat_unit, Year == year & Month %in% month)
  facilities_in_year <- with(facilities_in_year, unique(FacID))
  drop_fac <- length(unique(dat_unit$FacID)) - length(facilities_in_year)
  print(paste(drop_fac, 'facilities do not have data during this period.'))
  dat_unit <- subset(dat_unit, FacID %in% facilities_in_year)
  
  return(dat_unit)
}