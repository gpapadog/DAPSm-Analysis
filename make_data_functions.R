# Author: Georgia Papadogeorgou
# Date: 5/14/2016
# Desc: Basic functions for load and manipulating the power plant data in order
#       to predict heat input, aggregate to the facility level, and link with 
#       ozone, temperature and Census.
#
# Update: 7/1/2016
# I updated the functions to work for a period longer that a specific month.
# Update: 7/29/2016
# I no longer drop the facilities with missing or 0 operating time, we chose to
# ignore operating time because it was well predicted by heat input and gross load.

LoadUnitLevelData <- function(data_dir = '~/Dropbox/') {
  # Function that loads the unit level emissions for coal or gas units.
  #
  # Returns:
  #  The data frame (data.table format) with emissions and many other variables.
  
  # Read in Data that is monthly emissions data merged with other stuff
  dat <- fread(paste0(data_dir, "AMPD_Unit_with_Sulfur_Content_and_Regulations_",
                      "with_Facility_Attributes.csv"))
  dat$uID = paste(dat$Facility.ID..ORISPL., dat$Unit.ID, sep = "_")
  dat$FacID = dat$Facility.ID..ORISPL.
  dat$Year = as.numeric(dat$Year)
  dat$Month = as.numeric(dat$Month)
  dat[, year_month := paste(Year, Month, sep="_")]
  setkeyv(dat, c("uID", "Year", "Month"))
  setorderv(dat, c("uID", "Year", "Month"))
  dat <- unique(dat)
  
  # Restrict to Coal-fired or Natural Gas-fired units only
  dat_unit = subset(dat, grepl('Coal', Fuel.Type..Primary..x) |
                      grepl('Natural Gas', Fuel.Type..Primary..x))
  setkeyv(dat_unit, c("uID", "Year", "Month"))

  print('Unit Level Data Loading Complete.')
  return(dat_unit)
}


UnitToFacility <- function(dat_unit) {
  # Aggregating the power plant data from unit to facility, defining SO2, NOx controls.
  #
  # Args:
  #  dat_unit: The power plant data at the unit level.
  #
  # Returns:
  #  Data frame of the facility level data.
  
  ## -- Define NOx Control Strategies at the Unit Level
  print('Creating NOx control technologies variables')
  dat_unit = NOxcontroltechnologies(dat_unit)
  dat_unit$NumNOxControls = 0
  
  dat_unit[, idx := seq(nrow(dat_unit))]
  
  dat_unit[, SCR := as.numeric(SCR)]
  dat_unit[, SNCR := as.numeric(SNCR)]
  dat_unit[, LowNOxBurner := as.numeric(LowNOxBurner)]
  dat_unit[, OverFire := as.numeric(OverFire)]
  dat_unit[, Ammonia := as.numeric(Ammonia)]
  dat_unit[, CombustMod := as.numeric(CombustMod)]
  dat_unit[, WaterInj := as.numeric(WaterInj)]
  dat_unit[, OtherNOx := as.numeric(OtherNOx)]
  
  dat_unit[, NumNOxControls := sum(SCR, SNCR, LowNOxBurner, OverFire,
                                   Ammonia, CombustMod, WaterInj, OtherNOx),
           by = idx]
  dat_unit[, idx := NULL]

  # Group SCR and SnCR as one category
  dat_unit[, S_n_CR := as.numeric(SCR + SNCR >= 1)]
  dat_unit[, hasNOxControl := as.numeric(NumNOxControls > 0)]

  ## -- Define Operating Capacity as the maximum number of hours per month times the
  print('Creating Capacity related variables.')
  # Max Hourly Heat Input Rate
  ## maxopp is the max total number of hours operating time for each of 12 months
  maxopp <- with(dat_unit, tapply(Operating.Time, year_month, max, na.rm = TRUE))
  ### Max.Hourly.HI.Rate..MMBtu.hr. is the heat input capacity (MMBtu/hour)
  for (i in 1:length(maxopp))
    dat_unit[year_month == names(maxopp)[i],
             Capacity:= maxopp[i] * Max.Hourly.HI.Rate..MMBtu.hr.] # MMBtu/month.
  
  # Percent Capacity is Heat Input / Capacity
  dat_unit$PctCapacity <- dat_unit$Heat.Input..MMBtu / dat_unit$Capacity
  with(dat_unit, mean(PctCapacity > 1.5, na.rm = TRUE))
  dat_unit$PctCapacity[dat_unit$PctCapacity > 1.5] = NA ## > 150 % capacity
  
  
  ## -- Make Facility Level Data -- ##
  dat_unit[, Sulfur.Content := as.numeric(Sulfur.Content)]
  
  dat_facility = dat_unit[, list(nunits = length(unique(Unit.ID)),
                                 pctCoal = sum(grepl('Coal', Fuel.Type..Primary..x)) / length(unique(Unit.ID)),
                                 pctGas = sum(grepl('Natural Gas', Fuel.Type..Primary..x)) / length(unique(Unit.ID)),
                                 pctGas_byHI = sum(grepl('Natural Gas', Fuel.Type..Primary..x) * Heat.Input..MMBtu.) /
                                   sum(Heat.Input..MMBtu.),
                                 nunits_withsulfur = sum(!is.na(Sulfur.Content)),
                                 pctunits_withsulfur = sum(!is.na(Sulfur.Content)) / length(unique(Unit.ID)),
                                 meanSulfur_narm = sum(Sulfur.Content * Heat.Input..MMBtu., na.rm = TRUE) /
                                   sum(Heat.Input..MMBtu.[!is.na(Heat.Input..MMBtu.) & !is.na(Sulfur.Content)]),
                                 initialYear = Initial.Year.of.Operation[1],
                                 pctS_n_CR = sum(S_n_CR) / length(unique(Unit.ID)),
                                 S_n_CR_byHI_narm = sum(S_n_CR * Heat.Input..MMBtu., na.rm = TRUE) /
                                   sum(Heat.Input..MMBtu., na.rm = TRUE),
                                 S_n_CR_byHI = sum(S_n_CR * Heat.Input..MMBtu.) /
                                   sum(Heat.Input..MMBtu., na.rm = TRUE),
                                 # Percentage of heat input to units with NOx control.
                                 hasNOxControl_byHI = sum(hasNOxControl * Heat.Input..MMBtu.) /
                                   sum(Heat.Input..MMBtu., na.rm = TRUE),
                                 totOpTime_narm = sum(Operating.Time, na.rm = TRUE),
                                 totOpTime = sum(Operating.Time),
                                 totSO2emissions = sum(SO2..tons., na.rm = TRUE),
                                 totNOxemissions = sum(NOx..tons., na.rm = TRUE),
                                 totCO2emissions = sum(CO2..short.tons., na.rm = TRUE),
                                 totLoad = sum(Gross.Load..MW.h., na.rm = TRUE),
                                 totHeatInput_narm = sum(Heat.Input..MMBtu., na.rm = TRUE),
                                 totHeatInput = sum(Heat.Input..MMBtu.),
                                 pctCapacity = sum(Heat.Input..MMBtu.) / sum(Capacity),
                                 pctCapacity_byHI = sum(Heat.Input..MMBtu. * PctCapacity) /
                                   sum(Heat.Input..MMBtu., na.rm = TRUE),
                                 Phase2 = Is.Phase2[1],
                                 Fac.Latitude = Facility.Latitude.x[1],
                                 Fac.Longitude = Facility.Longitude.x[1],
                                 Fac.FIPS = FIPS[1],
                                 nmonths = length(unique(Month))),
                          by = "FacID"]
  print('Aggregation completed.')
  
  print('Percentage of facilities with pctCapacity > 1.5, set to NA.')
  print(with(dat_facility, mean(pctCapacity > 1.5, na.rm = TRUE)))
  dat_facility$pctCapacity[dat_facility$pctCapacity > 1.5] = NA ## > 150% capacity
  
  setkeyv(dat_facility, "FacID")
  setorderv(dat_facility, "FacID")

  return(dat_facility)
}




LinkPPtoMonitors <- function(dat, within_km, year, month, OzTempCen) {
  # Function that links the power plant data of a specific month to monitoring
  # data, after dropping facilities that were not operating.
  #
  # Args:
  #  dat:          Power plant data with coordinate information.
  #  within_km:    How many kilometers we want to perform the linkage at.
  #  year, month:  Year and month of the time period we are interested in.
  #  OzTempCen:    The file path where the file OzTempCen678_04.dat is saved.
  #
  # Returns:
  #  Data frame of the initial power plant data (after dropping non-operating
  #  ones), with additional ozone, temperature and Census information.
  

  load(paste0(OzTempCen, 'OzTempCen678_04.dat')) # OzTempCensus
  
  # Link the power plant data to ozone monitor data.
  ozpp_link <- spatial_link_index(dat, "Fac.Latitude", "Fac.Longitude", "FacID",
                                  OzTempCensus, "Latitude", "Longitude", "Monitor",
                                  within = within_km, closest = TRUE)

  print(paste('Monitors linked:', length(unique(ozpp_link$Monitor))))
  print(paste('Power plants linked:', length(unique(ozpp_link$FacID))))
  print(paste('Power plants dropped:', nrow(dat) - length(unique(ozpp_link$FacID))))
  
  
  # Merge linkage key with Ozone data, then aggragate to the facility level.
  
  # Merge monthly data with the linked file  
  OZmerge = merge(ozpp_link, OzTempCensus, by = "Monitor", all.x = TRUE)
  
  # Aggragate to the Facility Level
  oz_facility = OZmerge[, list(meanOzone = mean(meanOzone, na.rm = TRUE),
                               meanmaxOzone = mean(maxOzone, na.rm = TRUE),
                               mean4maxOzone = mean(fourthmaxOzone, na.rm = TRUE),
                               avgTemp = mean(avgTemp, na.rm = TRUE),
                               mean4MaxTemp = mean(mean4MaxTemp, na.rm = TRUE),
                               meanMaxTemp = mean(meanMaxTemp, na.rm = TRUE),
                               TotPop = sum(TotPop, na.rm = TRUE),
                               PctUrban = sum(PctUrban * TotPop, na.rm = TRUE) / sum(TotPop, na.rm = TRUE),
                               PctInUAs = sum(PctInUAs * TotPop, na.rm = TRUE) / sum(TotPop, na.rm = TRUE),
                               PctInUCs = sum(PctInUCs * TotPop, na.rm = TRUE) / sum(TotPop, na.rm = TRUE),
                               PctRural = sum(PctRural * TotPop, na.rm = TRUE) / sum(TotPop, na.rm = TRUE),
                               PctWhite = sum(PctWhite * TotPop, na.rm = TRUE) / sum(TotPop, na.rm = TRUE),
                               PctBlack = sum(PctBlack * TotPop, na.rm = TRUE) / sum(TotPop, na.rm = TRUE),
                               PctHisp = sum(PctHisp * TotPop, na.rm = TRUE) / sum(TotPop, na.rm = TRUE),
                               PctHighSchool = mean(PctHighSchool, na.rm = TRUE),
                               MedianHHInc = mean(MedianHHInc, na.rm = TRUE),
                               PctPoor = mean(PctPoor, na.rm = TRUE),
                               PctFemale = sum(PctFemale * TotPop, na.rm = TRUE) / sum(TotPop, na.rm = TRUE),
                               TotHUs = sum(TotHUs, na.rm = TRUE),
                               PctUrbanHUs = sum(PctUrbanHUs * TotHUs, na.rm = TRUE) / sum(TotHUs, na.rm = TRUE),
                               PctRuralHUs = sum(PctRuralHUs * TotHUs, na.rm = TRUE) / sum(TotHUs, na.rm = TRUE),
                               PctOccupied = sum(PctOccupied * TotHUs, na.rm = TRUE) / sum(TotHUs, na.rm = TRUE),
                               PctMovedIn5 = mean(PctMovedIn5, na.rm = TRUE),
                               MedianHValue = mean(MedianHValue, na.rm = TRUE),
                               PopPerSQM = mean(PopPerSQM, na.rm = TRUE),
                               nmonitors = length(Monitor)),
                        by = "FacID"]
  #565 Facilities with Ozone data
  
  ## Retain only those that are in both data sets (ie, no missing Ozone)
  dat <- merge(dat, oz_facility, by = "FacID")
  
  return(dat)
}





