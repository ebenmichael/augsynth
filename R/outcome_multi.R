################################################################################
## Fitting outcome models for multiple treatment groups
################################################################################


#' Get fixed effects from pre-treatment data for each level
#'
#' @param X Matrix of outcomes
#' @param trt Vector of treatment status for each unit
#' @param mask Matrix of treatment statuses
#' @param force Fixed effects: 1="unit", 2="time", 3="two-way"
#' @param time_cohort Boolean indicating whether to use time cohorts
#' @noRd
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_feff <- function(X, trt, mask, force, time_cohort) {

    ttot <- ncol(X)
    n <- nrow(X)
    # grps <- trt[is.finite(trt)]
    # iterate over treatment cohorts
    grps <- unique(trt[is.finite(trt)])
    J <- length(grps)
    which_t <- (1:n)[is.finite(trt)]

    if(force %in% c(2,3)) {
        ## compute time fixed effects from pure controls
        time_eff <- matrix(colMeans(X[!is.finite(trt),, drop = F],
                            na.rm = TRUE),
                            nrow=nrow(X),
                            ncol=ncol(X),
                            byrow=T)
    } else {
      time_eff <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    }
    residuals <- X - time_eff
    y0hat <- time_eff
    if(force %in% c(1,3)) {

        ## compute unit fixed effects from pre-intervention outcomes
        unit_eff <- lapply(grps, 
                            function(tj) matrix(
                                            rowMeans(residuals[, 1:tj, drop = F],
                                                     na.rm = TRUE),
                                            nrow=nrow(X), ncol=ncol(X)))
        residuals <- lapply(1:J, function(j) residuals -
                                                unit_eff[[j]])
        y0hat <- unit_eff
    }

    if(force == 3) {
        y0hat <- lapply(unit_eff, function(ufj) time_eff + ufj)
    }
    
    # go from treatment cohorts to individuals
    if(force %in% c(1,3) & !time_cohort) {
      names(residuals) <- as.character(grps)
      residuals <- residuals[as.character(trt[is.finite(trt)])]
      names(y0hat) <- as.character(grps)
      y0hat <- y0hat[as.character(trt[is.finite(trt)])]
    }
    
    return(list(y0hat = y0hat,
                residuals = residuals))
    
}

