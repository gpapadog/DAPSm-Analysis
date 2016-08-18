# Author: Georgia Papadogeorgou
# Date: 5/14/2016
# Description: The goal is to have one code file using which we can create the data
#              for the SnCR on ozone analysis using power plant facilities. The tasks
#              performed here are:
#              1) Use other time points with adjusted R squared greated large enough
#              to predict heat input for units with missing heat input.
#              2) Aggregate to the facilitiy level, use new definitions of treatment
#              weighted by heat input.
#              3) Link with the ozone, census and temperature data.

time_pred <- paste0(year, '_', month)


# ---- STEP 1: Predicting heat input for unit level data.

# Loading unit level data.
dat_unit <- LoadUnitLevelData()
# Picking the subset which we are interested in.
subdta <- subset(dat_unit, year_month %in% time_pred)

# To check predictive power of the models, look at predicting_heat_input.R.
variable <- 'Heat.Input..MMBtu.'
p <- PredictVariableMonths(year, month, time_use, variable, dat_unit)

# Sanity check using dimensions.
length(p$data)
dim(subdta)
# The next two numbers add to length(p$data) which is equal to nrow(subdta).
sum(p$data == subdta$Heat.Input..MMBtu., na.rm = TRUE)
sum(is.na(subdta$Heat.Input..MMBtu.))

# Assinging the new values (which include the predicted Heat input values).
print(paste(sum(p$data < 0, na.rm = TRUE), 'out of', 
            sum(is.na(subdta$Heat.Input..MMBtu.)), 'predicted values are',
            'negative, and are forced to 0.'))
subdta$Heat.Input..MMBtu. <- p$data
subdta$Heat.Input..MMBtu.[subdta$Heat.Input..MMBtu. < 0] <- 0


# ---- STEP 2: Aggregate to the facility level.

# Aggregating only the data from this time period.
dat_facility <- UnitToFacility(dat_unit = subdta)
dim(dat_facility)  # 1318 facilities.
sum(dat_facility$nunits)  # 3955 units.
print(paste('Dropping', sum(dat_facility$totHeatInput == 0, na.rm = TRUE),
            'facilities for heat input = 0'))
dat_facility <- subset(dat_facility, totHeatInput > 0 | is.na(totHeatInput))
table(is.na(dat_facility$totHeatInput), is.na(dat_facility$S_n_CR_byHI))


# ---- STEP 3: Linking the aggregated data to ozone monitors with ozone
# temperature and Census information.
dat <- LinkPPtoMonitors(dat_facility, within_km, year, month)


# Facility IDS with missing data for at least one month.
wh <- which(dat$nmonths != length(month))
print(paste('Dropping', length(wh), 'out of', length(unique(dat$FacID)),
            'facilities due to missing information on at least one month'))
dat <- dat[- wh, ]
dat[, nmonths := NULL]

