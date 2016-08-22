# DAPSm-Analysis
Code to perform data analysis for the DAPSm paper

In general, repeating the data analysis requires that one would need to 

- Save all scripts at the same directory
- Save data sets
- Specify the directory where the scripts and data sets are saved
- Run the analysis script ```Data_analysis.R```.

But first...


## Before running the analysis

### Have all the pieces!

In order to run the analysis we need two data sets:

1. Data set on power plant characteristics available at ......
2. Data set on ozone, temperature, and Census information measured at the ozone AQS monitoring sites, available at ......

We further need the following R scripts:

1. ```config.R```: The file including the configuration of the analysis.
1. ```CreateNOxControlsFunction.R```: Function that created the different NOx emission control indicators.
1. ```Data_analysis_functions.R```: Functions that drop columns that we do not want to include in the analysis and reform others.
1. ```Data_analysis_models_functions.R```: Functions that fit the various methods that are fit on the data set. Depend on other scripts.
1. **```Data_analysis.R```**: Contains the main data analysis code which loads all data sets and sources in the necessary functions.
1. ```DistCal_functions.R```: Functions to fit matching in distance calipers. Sourced.
1. ```expit.R```: The expit functions.
1. ```GBMPropScores_function.R```: The function that is used to estimate the propensity scores using GBM.
1. ```make_data_functions.R```: This file includes functions that load the unit level data, aggregate unit level to facility level, and link to ozone monitoring sites.
1. ```predict_variable_functions.R```: Functions that are used when we predict heat input.
1. ```PredictHeatInput.R```: Function that is called to predict heat input at the unit level data.
1. ```PSmatchEst_function.R```: Function that is used to match data on a set of estimated propensity scores and return estimates of the causal effect, matched pairs, and balance.
1. ```StandDiff_function.R```: Function that calculates the standardized difference of means.


### Save scripts and data sets

You will need to specify the directory where all the scripts are (save all scripts in the same directory with the names used here), as well as the directory where the power plant data are, and the file path of the linked ozone, temperature, Census data.

### Installing necessary R packages

- In order to download, and link AQS data, Christine Choirat has created a R package called *arepa* and is available through Github.

- For the analysis using DAPSm, you will need the DAPSm R package which is also available through Github.

In order to download Github packages, you should install ```devtools``` on R studio by running ```install.packages('devtools')```, and load the corresponding library ```library(devtools)```.

Download the arepa and DAPSm packages using

```
devtools::install_github("czigler/arepa")
devtools::install_github("gpapadog/DAPSm")
```

Other R packages are also used in the analysis, and are all available on CRAN:
data.table, corrplot, MatchIt, Matching, fields, ggplot2, stringr.


### Updating the configuration file

This is something that needs to be done, but should not take more than a couple minutes! In the ```config.R``` file, one needs to specify:

- ```wd```: Your working directory.
- ```source.path```: The path to the folder where all the R scripts are saved.
- ```data_dir```: The path to the folder where the power plant data are saved. File name must match the name of the file provided.
- ```OzTempCensus```: The file path where the covariate information is saved.
- ```outcome_analysis```: The name of the variable which we consider as the outcome of interest. Options are:
    - ```mean4maxOzone```: The average of the fourth maximum ozone concentration.
    - ```totNOxemissions```: Total emissions of power plants.

The remaining variables should not be changed if one wants to replicate the results of the paper.

- ```caliper```: Caliper used in matching.
- ```cutoff```: The cutoff of ASDM below which a variable is considered balanced.
- ```weights```: The weights over which we fit DAPSm in order to investigate the optimal w.

- ```year```: Can only be changed if one has linked ozone, temperature, Census data on a year other than 2004.
- ```month```: Can only be changed if one has linked ozone, temperature, Census data on months other than June, July and August 2004.
- ```time_use```: A list of the years which will be used in order to predict heat input.
- ```within_km```: How many kilometers is the maximum linkage of power plants to ozone monitoring sites.


## Running the analysis script
The main data analysis script is named ```Data_analysis.R```.

- Open ```Data_analysis.R``` in R/Rstudio
- Specify the file path where the configuration file is saved as ```config_path```
- Run the full script