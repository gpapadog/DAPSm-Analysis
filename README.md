# DAPSm-Analysis
Code to perform data analysis for the DAPSm paper

## Before running the analysis

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


### Linking ozone, temperature and Census information


## Running the analysis script
The main data analysis script is named ```Data_analysis.R```.