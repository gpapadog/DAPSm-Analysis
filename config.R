# This is the configuration file for setting up the data analysis.

wd <- "~/Github/DAPSm-Analysis/"
# This must be the directory where all the function scripts are saved.

# Specify the directory where the data on power plants, and covariate information are
# are saved. The files must be named:
# "AMPD_Unit_with_Sulfur_Content_and_Regulations_with_Facility_Attributes.csv", for
# the power plant data, and
# 'OzTempCen678_04.dat', for the covariate information data.
data_dir <- '~/Dropbox/DATAverse/'

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


# Optimal or greedy matching.
matching_algorithm <- 'greedy'
remove_unmatchables <- TRUE

# For the method of Keele et al.
n_matches <- 1  # Number of matched per treated.
use_controls <- NULL  # Whether specific controls need to be used
enforce_constraints <- FALSE
subset_weight <- 920

