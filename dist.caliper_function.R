#' Function that takes in a two data frames one with treatment and one with control
#' units and returns a matrix of the matched pairs, that are matched using a distance
#' caliper and a PS caliper.
#' 
#' @param treated A data frame include the treated units and the variables: longitude,
#' latitude and propensity scores (must be named 'prop.scores'). The rownames of
#' treated should be the unit ids.
#' @param control Control units. Same variables as in treated. Rownames should be the
#' unit ids of the controls
#' @param ps.caliper A caliper of propensity score difference for matching. Caliper is
#' set as number of sd of the ps distribution.
#' @param dist.quan The distance of all treated-control pairs is calculated. dist.quan
#' is a scalar between 0, and 1 describing the quantile of all the pairwise distances
#' that should be used as a distance caliper.
#' @param coords.columns
#' If the columns of coordinates are not named 'Longitude', 'Latitude', coords.cols
#' should be the column indeces corresponding to longitude and latitude accordingly.
#' @param coord_dist
#' Set to true when we want to use a distance function that calculates the spherical
#' distance of points instead of euclidean. Defaults to FALSE.
#' @param matching_algorithm Argument with options 'optimal', or 'greedy'. The optimal
#' choice uses the optmatch R package to acquire the matches based on propensity score
#' difference and a caliper on distance. The greedy option matches treated and control
#' units sequentially, starting from the ones with the smallest propensity score
#' difference. Defaults to 'optimal'.
#' 
#' @return A dataframe, where each row corresponds to each treated unit, and includes
#' the control unit to which it was matched, their propensity score difference, their
#' DAPS difference, their distance, their standardized distance.
#' 
dist.caliper <- function(treated, control, ps.caliper = 0.1, dist.quan = 0.25,
                         coords.columns = NULL, coord_dist = FALSE,
                         matching_algorithm = c('optimal', 'greedy')) {
  
  matching_algorithm <- match.arg(matching_algorithm)
  require(fields)  # For rdist().
  require(optmatch)  # For caliper() and match_on() functions.
  
  # Setting the caliper.
  caliper <- caliper * sd(c(treated$prop.scores, control$prop.scores))
  
  if (!is.null(coords.columns)) {
    names(treated)[coords.columns] <- c('Longitude', 'Latitude')
    names(control)[coords.columns] <- c('Longitude', 'Latitude')
  }
  
  if (coord_dist) {
    dist.mat <- rdist.earth(cbind(treated$Longitude, treated$Latitude),
                            cbind(control$Longitude, control$Latitude))
  } else {
    dist.mat <- rdist(cbind(treated$Longitude, treated$Latitude),
                      cbind(control$Longitude, control$Latitude))
  }
  colnames(dist.mat) <- rownames(control)
  rownames(dist.mat) <- rownames(treated)
  
  cut.off <- quantile(dist.mat, probs = dist.quan)
  # Using caliper() function from the optmatch package that returns 0 or Inf.
  infinite_distance <- caliper(dist.mat, cut.off)
  
  # Using match_on() from optmatch to get the propensity score difference.
  treatment_indicator <- c(rep(1, nrow(treated)), rep(0, nrow(control)))
  propensity_scores <- c(treated$prop.scores, control$prop.scores)
  D <- match_on(treatment_indicator ~ propensity_scores, method = 'euclidean')
  # Adding the distance caliper.
  D <- D + infinite_distance
  # Adding the propensity score caliper.
  D <- D + caliper(D, caliper)
  
  # The matrix of matches is mat.
  mat <- data.frame(match = rep(NA, dim(treated)[1]),
                    distance = rep(NA, dim(treated)[1]),
                    prop.diff = rep(NA, dim(treated)[1]))
  rownames(mat) <- rownames(treated)
  
  if (matching_algorithm == 'greedy') {
    pairs <- MinDistMatch(as.matrix(D), caliper = NULL)
  } else {  # Optimal matching using the optmatch R package.
    opt_match <- pairmatch(D, data = data.frame(treatment_indicator))
    
    pairs_ids <- sort(as.character(unique(opt_match[!is.na(opt_match)])))
    wh_trt <- 1:nrow(treated)
    wh_con <- (nrow(treated) + 1) : (nrow(treated) + nrow(control))
    match_trt <- cbind(wh_trt, group = as.character(opt_match[wh_trt]))
    match_con <- cbind(wh_con, group = as.character(opt_match[wh_con]))
    pairs <- merge(match_trt, match_con, by = 'group')
    pairs <- pairs[, - which(names(pairs) == 'group')]
    pairs <- na.omit(pairs)
    pairs[, 1] <- as.numeric(as.character(pairs[, 1]))
    pairs[, 2] <- as.numeric(as.character(pairs[, 2])) - nrow(treated)
  }
  
  matched_trt <- pairs[, 1]
  matched_con <- pairs[, 2]
  
  mat$match[matched_trt] <- rownames(control)[matched_con]
  for (ii in 1:length(matched_trt)) {
    wh_trt <- matched_trt[ii]
    wh_con <- matched_con[ii]
    mat$prop.diff[wh_trt] <- treated$prop.scores[wh_trt] -
      control$prop.scores[wh_con]
    mat$distance[wh_trt] <- dist.mat[wh_trt, wh_con]
  }
  return(mat)
}
