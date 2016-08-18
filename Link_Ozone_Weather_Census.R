# Date: 2/7/2016
# Description: Getting the weather, and Census variables at the Ozone
#              monitor level.
# Update: 7/13/2016
# I updated the function to accommodate situation where the number of months is a
# vector.

library(arepa)
library(data.table)

GetOzoneTempCensus <- function(year, month) {
  # Function that can be run to load and link ozone temperature and census data.
  # Temp-Ozone within 150 km, Census-Ozone within 6 miles
  
  setwd('/Users/georgiapapadogeorgou/Documents/ARP/Application/AQS')
  source(paste0('/Users/georgiapapadogeorgou/Documents/ARP/Application/',
                'Make Data Code/subset_monitors_daily.R'))
  
  # Loading temperature data.
  print('Loading temperature data.')
  years = (year - 1) : year
  year_inter <- year  # Year we are interested in keeping.
  observation_percent = 0 # as per EPA suggestion
  
  Temp = load_daily_data(parameter = 'TEMP', years)
  Temp[, Year := as.numeric(substr(Date.Local, 1, 4))]
  Temp[, Month:= as.numeric(substr(Date.Local, 6, 7))]
  Temp[, Day:= as.numeric(substr(Date.Local, 9, 10))]
  Temp <- subset(Temp, Year == year_inter & Month %in% month)
  
  # Drop Microscale and for each day keep the observation with smallest POC.
  temp_daily = subset_monitors_daily(MONITORS = Temp, 62101, 
                                     observation_percent = observation_percent, 
                                     monitor_info_file = "~/Google Drive/monitor_list.csv")
  
  # Aggregate to Monthly
  # Choose the max, 4th max and mean over the period of months.
  temp_monthly = temp_daily[, list(Latitude = Latitude[1], Longitude = Longitude[1], 
                                   maxmaxTemp = max(X1st.Max.Value),
                                   fourthmaxTemp = X1st.Max.Value[order(X1st.Max.Value)[4]],
                                   meanTemp = mean(Arithmetic.Mean),
                                   meanMaxTemp = mean(X1st.Max.Value)),
                            by = "Monitor"]
  setkeyv(temp_monthly, "Monitor")
  
  Temp <- copy(temp_monthly)
  rm(temp_daily, temp_monthly)
  names(Temp)[names(Temp) == 'Monitor'] <- 'station'
  
  # Dropping the ones that are too low.
  wh <- which(Temp$meanTemp < 30)
  Temp <- Temp[- wh, ]
  print('Temperature data loaded.')
  
  
  # Loading the Ozone data.
  print('Loading Ozone data.')
  parameter_code = 44201 # for Ozone.
  
  aqsdat = load_daily_data(parameter = parameter_code, years)
  aqsdat[, Year := as.numeric(substr(Date.Local, 1, 4))]
  aqsdat[, Month:= as.numeric(substr(Date.Local, 6, 7))]
  aqsdat[, Day:= as.numeric(substr(Date.Local, 9, 10))]
  aqsdat <- subset(aqsdat, Year == year_inter & Month %in% month)
  
  OZ_daily = subset_monitors_daily(MONITORS = aqsdat, parameter_code, 
                                   observation_percent = observation_percent, 
                                   monitor_info_file = "~/Google Drive/monitor_list.csv")
  
  ### -- Aggregate to Monthly
  OZ_monthly = OZ_daily[, list(Latitude = Latitude[1], Longitude = Longitude[1],
                               fourthmaxOzone = X1st.Max.Value[order(X1st.Max.Value)[4]],
                               meanOzone = mean(Arithmetic.Mean),
                               maxOzone = max(X1st.Max.Value)),
                        by = "Monitor"]
  length(unique(OZ_monthly$Monitor)) ## 1319 Monitors
  setkeyv(OZ_monthly, "Monitor")
  
  OZ <- copy(OZ_monthly)
  rm(OZ_daily, OZ_monthly)
  print('Ozone data loaded.')
  
  
  # Linking Ozone data and weather data from aqs monitors.
  
  ## For every temp monitoring location, find all ozone monitors that are within 150km
  ## Note: this is NOT a unique linkage
  OZtemp_link = spatial_link_index(OZ, "Latitude", "Longitude", "Monitor",
                                   Temp, "Latitude", "Longitude", "station",
                                   within = 150, closest = FALSE)
  length(unique(OZtemp_link$station)) #715 temperature stations
  length(unique(OZtemp_link$Monitor)) #1180 monitors - some did not have a match.
  
  
  temp_withmonitor = merge(Temp, OZtemp_link, by = "station", allow.cartesian = TRUE)
  length(unique(temp_withmonitor$station)) #715 temperature stations
  length(unique(temp_withmonitor$Monitor)) #1180 monitors
  
  ## -- Aggregate to the monitor level
  temp_Monitor = temp_withmonitor[, list(avgTemp = mean(meanTemp, na.rm = TRUE),
                                         mean4MaxTemp = mean(fourthmaxTemp, na.rm = TRUE),
                                         meanMaxTemp = mean(meanMaxTemp, na.rm = TRUE),
                                         nmonitors = length(station)),
                                  by = "Monitor"]
  length(unique(temp_Monitor$Monitor)) # 1180 monitors
  mean(temp_Monitor$Monitor %in% OZ$Monitor)
  
  OZtemp <- merge(OZ, temp_Monitor, by = 'Monitor')
  length(unique(OZtemp$Monitor)) # 1171 monitors
  # OZtemp is the data frame including both ozone and temperature information,
  # at the ozone monitoring sites.
  
  
  
  
  # Now we have the Ozone and weather data, and we want to link it to the Census
  # data.
  
  ### 2000 Census
  censusdat = "~/Google Drive/ARP/Census_Dexter/2000/"
  d_census = fread(paste(censusdat, "subset_census_2000_with_zip_and_percentage.csv", sep=""))[, V1 := NULL]
  ##-----  Variables selected are determined in census_2000.R script and already subsetted here
  subset_census = d_census
  
  
  
  # Merge with Another Data Source (here, AQS Ozone and weather data are used)
  # We have the OZtemp data we want to merge census data with.
  
  # -- Get Zip Codes
  ZIP = get_zip_codes()
  
  ## -- Link AQS Monitors to all zip codes within a specified radius so they can be linked to Census data
  ## -- Note: This is not a unique linkage: a zip code can be assigned to more than one monitor
  OZzip_link <- spatial_link_index(OZtemp, "Latitude", "Longitude", "Monitor",
                                   ZIP, "Latitude.zip", "Longitude.zip", "zip",
                                   within = 9.656, closest = FALSE)
  length(unique(OZzip_link$Monitor))  # 1079 Monitors with a zip code located within within_km 
  length(unique(OZzip_link$zip))  # 9068 Zip codes located within within_km of a monitor
  OZzip_link$ZIP = as.character(OZzip_link$zip)
  
  ## -- Merge AQS data with Census data by zip code
  OZzip_link_census <- merge(OZzip_link, subset_census, by = "ZIP", all.x = TRUE)
  length(unique(OZzip_link_census$Monitor)) # 1088 Monitors
  length(unique(OZzip_link_census$ZIP))  # 9068 Zip Codes
  
  
  # Aggregate Census Data (for example, to the level of an AQS Monitor)
  
  ## - Note: Be careful here with how these aggregated variables are defined.  The denominator of percentages changes for some variables
  CENSUS_Monitor = OZzip_link_census[, list(TotPop = sum(TotPop, na.rm=TRUE), 
                                            PctUrban = sum(Urban, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                            PctInUAs = sum(InUAs, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                            PctInUCs = sum(InUCs, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                            PctRural = sum(Rural, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                            PctWhite = sum(White1, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                            PctBlack = sum(Black1, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                            PctHisp = sum(HispPop, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                            PctHighSchool = sum(HighSchool, na.rm = TRUE)/sum(Over25, na.rm = TRUE),
                                            MedianHHInc = mean(MedianHHInc, na.rm = TRUE),
                                            PctPoor = sum(Poor, na.rm = TRUE)/sum(PovUniverse, na.rm = TRUE),
                                            PctFemale = sum(Female, na.rm = TRUE)/sum(TotPop, na.rm = TRUE),
                                            TotHUs = sum(TotHUs, na.rm = TRUE),
                                            PctUrbanHUs = sum(UrbanHUs, na.rm = TRUE)/sum(TotHUs, na.rm = TRUE),
                                            PctRuralHUs = sum(RuralHUs, na.rm = TRUE)/sum(TotHUs, na.rm = TRUE),
                                            PctOccupied = sum(Occupied, na.rm = TRUE)/sum(TotHUs, na.rm = TRUE),
                                            PctMovedIn5 = sum(MovedIn5, na.rm = TRUE)/sum(Occupied, na.rm = TRUE),
                                            MedianHValue = mean(MedianHValue, na.rm = TRUE),
                                            PopPerSQM = sum(TotPop, na.rm = TRUE) / sum(LandSQMI, na.rm = TRUE)),
                                     by = "Monitor"]
  #                                          ### Note: PopPerSQM is missing for 1990 Census
  
  ## -- Note: places with missing TotPop and TotHUs will be given values of 0 with the above code.  Set them back to missing.
  CENSUS_Monitor$TotPop[CENSUS_Monitor$TotPop==0] = NA
  CENSUS_Monitor$TotHUs[CENSUS_Monitor$TotHUs==0] = NA
  colMeans(is.na(CENSUS_Monitor)) ## approx 0.01 (year 2000) missing census data
  length(unique(OZzip_link_census$Monitor))  #1088 Monitors
  
  
  # Some census variables are highly correlated.  Check which ones are close to redundant
  censvars = c( "TotPop", "PctUrban", "PctInUAs", "PctInUCs",                   
                "PctRural", "PctWhite", "PctBlack",                   
                "PctHisp", "PctHighSchool", "MedianHHInc",                
                "PctPoor", "PctFemale", "TotHUs",                     
                "PctUrbanHUs", "PctRuralHUs", "PctOccupied",                
                "PctMovedIn5", "MedianHValue", "PopPerSQM")
  
  M = cor(CENSUS_Monitor[, censvars, with = FALSE], use = "pairwise.complete.obs")
  par(mfrow = c(1,1))
  # pdf(file = "census_corplot.pdf")
  corrplot(M)
  dev.off()
  
  ## Check correlations in detail
  with(CENSUS_Monitor, plot(TotPop, TotHUs)) ## corr ~ 0.99
  with(CENSUS_Monitor, cor(TotPop, TotHUs, use = "pairwise.complete.obs"))
  
  
  OzTempCensus <- merge(OZtemp, CENSUS_Monitor, by = 'Monitor')
  
  OzTempCensus$Year <- year_inter
  OzTempCensus$Month <- paste(month, collapse = ', ')
  
  save(OzTempCensus, file = paste0('/Users/georgiapapadogeorgou/Documents/ARP/',
                                   'Application/Data/OzTempCen',
                                   paste(month, collapse = ''),
                                   '_', substr(as.character(year_inter), 3, 4), '.dat'))
  return(OzTempCensus)
}


