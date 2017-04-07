#' Function that takes the whole data set with power plant information and drops units
#' and facilities that were either retired, not operating, or not not have any data
#' during the study period.
CleanPPunits <- function(dat_unit, year, month) {
  
  # Dropping facilities that do not have data during the time period.
  facilities_in_year <- subset(dat_unit, Year == year & Month %in% month)
  facilities_in_year <- with(facilities_in_year, unique(FacID))
  drop_fac <- length(unique(dat_unit$FacID)) - length(facilities_in_year)
  print(paste(drop_fac, 'facilities do not have data during this period.'))
  dat_unit <- subset(dat_unit, FacID %in% facilities_in_year)
  
  # Dropping units that were retired during the time period.
  retired_in_year <- subset(dat_unit, Year == year & Month %in% month & Is.Retired)
  drop_uID <- unique(retired_in_year$uID)
  print(paste(length(drop_uID), 'units in', length(unique(retired_in_year$FacID)),
              'facilities were retired.'))
  dat_unit <- subset(dat_unit, ! (uID %in% drop_uID))
  rm(retired_in_year)
  rm(drop_uID)
  
  # Dropping units that started operating after the period.
  not_operating <- subset(dat_unit, Initial.Year.of.Operation > year)
  drop_uID <- unique(not_operating$uID)
  print(paste(length(drop_uID), 'units in', length(unique(not_operating$FacID)),
              'facilities were not yet operating'))
  dat_unit <- subset(dat_unit, ! (uID %in% drop_uID))
  rm(not_operating)
  rm(drop_uID)
  
  # Dropping units that were not operating during the period.
  not_operating <- subset(dat_unit, Status != 'Operating')
  not_operating <- subset(not_operating, Year == year & Month %in% month)
  drop_uID <- unique(not_operating$uID)
  print(paste(length(drop_uID), 'units in', length(unique(not_operating$FacID)),
              'facilities were not operating during the whole period.'))
  dat_unit <- subset(dat_unit, ! (uID %in% drop_uID))
  rm(not_operating)
  rm(drop_uID)
  
  # Dropping units with no data during or after the time period.
  no_data_uID <- dat_unit[, list(maxYear = max(Year),
                                 maxMonth = max(Month[Year == max(Year)]),
                                 FacID = FacID[1]), by = uID]
  no_data_uID <- subset(no_data_uID, maxYear < year |
                          (maxYear == year & maxMonth < month[length(month)]))
  drop_uID <- unique(no_data_uID$uID)
  print(paste(length(drop_uID), 'units in', length(unique(no_data_uID$FacID)),
              'facilities stopped having data before the study period.'))
  dat_unit <- subset(dat_unit, ! (uID %in% drop_uID))
  rm(no_data_uID)
  rm(drop_uID)
  
  # Dropping units with no data during or before the time period.
  no_data_uID <- dat_unit[, list(minYear = min(Year),
                                 minMonth = min(Month[Year == min(Year)]),
                                 FacID = FacID[1]), by = uID]
  no_data_uID <- subset(no_data_uID, minYear > year |
                          (minYear == year & minMonth > month[1]))
  drop_uID <- unique(no_data_uID$uID)
  print(paste(length(drop_uID), 'units in', length(unique(no_data_uID$FacID)),
              'facilities started having data after the study period.'))
  dat_unit <- subset(dat_unit, ! (uID %in% drop_uID))
  rm(no_data_uID)
  rm(drop_uID)
  
  return(dat_unit)
}