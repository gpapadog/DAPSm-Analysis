# This is the configuration file for setting up the data analysis.

wd <- "~/Github/DAPSm-Analysis/"

# The folder where the scripts for performing the analysis are.
source.path <- paste0('~/Github/DAPSm-Analysis/')

# Specify the directory where the data on power plants are saved. The file must be
# named: "AMPD_Unit_with_Sulfur_Content_and_Regulations_with_Facility_Attributes.csv".
data_dir <- '~/Dropbox/'

# The file path where the linked ozone-temperature-Census data are saved.
OzTempCensus <- paste0('/Users/georgiapapadogeorgou/Documents/ARP/Application/',
                       'Data/OzTempCen678_04.dat')

# What is the outcome of interest
outcome_analysis <- 'mean4maxOzone'
# possible outcome_analysis are:
# 1. 'mean4maxOzone'
# 2. 'totNOxemissions'

# Specify the analysis details
caliper <- 1
cutoff <- 0.15
weights <- seq(0, 1, length.out = 40)  # values of w for the optimal DAPSm scan.


# Time period that we want to analyze.
year <- 2004
month <- 6:8
# Years we use in order to predict heat input when missing.
time_use <- list('2003', '2002', '2005', '2006')
# Within how many kilometers we want to link the power plants to the ozone sites.
within_km <- 100
