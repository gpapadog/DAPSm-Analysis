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

You will need to specify the directory where all the scripts are (save all scripts in the same directory with the names used here), as well as the file paths of the two data sets.

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


## Running the analysis script
The main data analysis script is named ```Data_analysis.R```.