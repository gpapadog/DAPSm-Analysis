# Author: Georgia Papadogeorgou
# Date: 5/16/2016
# Desc: Functions to perform the data analysis. Includes functions that clean the data,
#       and functions that fit the models.
#
# Update: (7/1/2016) Added function that plots matched pairs on the map.
# Update: (7/29/2016) Dropping operating time.

CleanData <- function(dataset, plotcor = FALSE) {
  
  dataset[, nmonths := NULL]
  
  # Dropping data with missing heat input
  wh <- sum(is.na(dataset$totHeatInput))
  print(paste('Dropping', wh, 'facilities due to missing heat input.'))
  dataset <- subset(dataset, !is.na(totHeatInput))
  # Dropping facilities with NEGATIVE heat input
  wh <- sum(dataset$totHeatInput < 0)
  print(paste('Dropping', wh, 'facilities due to negative heat input.'))
  dataset <- subset(dataset, (totHeatInput >= 0))
  
  # Dropping year and month.
  if ('Year' %in% names(dataset)) {
    dataset[, Year := NULL]
  }
  if ('Month' %in% names(dataset)) {
    dataset[, Month := NULL]
  }
  if ('sdOzone' %in% names(dataset)) {
    dataset[, sdOzone := NULL]
  }
  
  # Dropping columns with sd 0.
  wh <- (apply(dataset, 2, function(x) sd(x, na.rm = TRUE)) < 0.0001)
  wh[is.na(wh)] <- FALSE
  if (any(wh, na.rm = TRUE)) {
    dataset[, names(dataset)[wh] := NULL, with = FALSE]
  }

  # When we have more control units we will estimate ATT:
  if (with(dataset, mean(S_n_CR_byHI >= 0.5)) < 0.5) {
    dataset[, SnCR := S_n_CR_byHI >= 0.5]
    print('Treated defined as SnCR_byHI >= 0.5')
  } else {
    dataset[, SnCR := S_n_CR_byHI < 0.5]
    print('Treated defined as SnCR_byHI < 0.5')
  }

  # Dropping columns we do not want to include in the analysis.
  dataset[, FacID := NULL]
  dataset[, Fac.FIPS := NULL]
#  dataset[, totNOxemissions := NULL]
  
  # Looking at the correlation of variables:
  if (plotcor) {
    C <- cor(dataset, use = 'pairwise.complete.obs')
    corrplot(C, tl.cex = 0.6)
  }

  # Dropping CO2 emissions and tot Load, keeping heat input with na.
  dataset[, totCO2emissions := NULL]
  dataset[, totLoad := NULL]
  dataset[, totHeatInput_narm := NULL]

  dataset[, pctCoal := NULL]
  dataset[, pctGas := NULL]
  dataset[, PctRuralHUs := NULL]
  dataset[, PctRural := NULL]
  dataset[, PctUrbanHUs := NULL]
  dataset[, pctS_n_CR := NULL]
  dataset[, S_n_CR_byHI_narm := NULL]
  dataset[, S_n_CR_byHI := NULL]
  dataset[, initialYear := NULL]
  
  dataset[, avgTemp := NULL]
  dataset[, meanMaxTemp := NULL]
  
  dataset[, meanSulfur_narm := NULL]

  dataset[, totOpTime_narm := NULL]
  # Dropping operating time since we will not use it in the analysis after all.
  dataset[, totOpTime := NULL]
  dataset[, pctCapacity := NULL]
  
  dataset[, TotHUs := NULL]
  dataset[, pctunits_withsulfur := NULL]
  
  dataset[, PctInUAs := NULL]  # 82% correlation with PctUrban
  dataset[, nmonitors := NULL]
  
  if (plotcor) {
    C <- cor(dataset, use = 'pairwise.complete.obs')
    corrplot(C, tl.cex = 0.6)
  }
  
  # Dropping the two observations with missing Census covariates.
  print(paste('Dropping', sum(is.na(dataset$pctCapacity_byHI)),
              'observations because of missing percent capacity.'))
  dataset <- dataset[!is.na(dataset$pctCapacity_byHI), ]
  
  # Dropping the two observations with missing Census covariates.
  print(paste('Dropping', sum(is.na(dataset$PctPoor)), 'observations',
              'because of missing Census information.'))
  dataset <- dataset[!is.na(dataset$PctFemale), ]
  
  return(dataset)
}


ReformData <- function(dataset) {
  
  # Taking the logarithm of heat input and tot operating time.
  dataset[, logHeatInput := log(totHeatInput)]
  # dataset[, logOpTime := log(totOpTime)]
  dataset[, logPopPerSQM := log(PopPerSQM)]
  # Changing PctGas_byHI to categories coal and gas
  dataset[, mostlyGas := as.numeric(pctGas_byHI >= 0.5)]
  # Creating categories for nunits
  dataset[, small_nunits := as.numeric(nunits %in% c(1, 2))]
  dataset[, med_nunits := as.numeric(nunits %in% c(3, 4, 5))]
  
  dataset[, totHeatInput := NULL]
  # dataset[, totOpTime := NULL]
  dataset[, PopPerSQM := NULL]
  dataset[, pctGas_byHI := NULL]
  dataset[, nunits := NULL]

  # Dropping pct variables. Plotted them against mean ozone with color of SnCR.
  dataset[, PctFemale := NULL]
  dataset[, PctInUCs := NULL]  # Kept PctUrban instead
  dataset[, TotPop := NULL]  # log of this has cor 88% with log PopPerSQM.
  
  dataset[, totSO2emissions := NULL] # Predictor of coal - gas?
  dataset[, nunits_withsulfur := NULL] # PctGas and nunits might be enugh.
  
  return(dataset)
}



PlotResults <- function(result, title = NULL, title.cex = 1, center = FALSE) {
  # Function that plots the data analysis results.
  #
  # Args:
  #  result: A matrix or data frame with 3 columns named 'LB'-lower bound,
  #          'Estimate', and 'UB'-upper bound for the estimate and CIs. The
  #          rownames of the matrix must correspond to the method used for
  #          each estimate.
  
  if (is.null(title)) {
    title <- 'Causal Effect estimates with 95% confidence intervals'
  }
  
  result <- as.data.frame(result)
  result$method <- rownames(result)
  result$method <- factor(result$method, levels = result$method)
  
  g <- ggplot(result, aes(x=method, y=Estimate, group=1)) +
    geom_errorbar(width=.1, aes(ymin = LB, ymax = UB), color = 'grey65',
                  data = result, cex = 1.5) +
    geom_point(shape=21, size=3, fill="grey65") +
    theme(
      panel.background = element_rect(fill = "grey92",
                                      size = 0.5, linetype = "solid"),
      panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                      colour = "white"), 
      panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                      colour = "white")
    ) +
    xlab('') +
    theme(panel.border = element_blank()) +
    ggtitle(title) +
    theme(plot.title = element_text(size = rel(title.cex)),
          axis.title = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.4)),
          legend.position = 'none')
  if (center) {
    g <- g + ylim(c(- max(abs(result[, 1:3])), max(abs(result[, 1:3]))))
  }
  print(g)
}


