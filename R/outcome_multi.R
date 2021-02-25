################################################################################
## Fitting outcome models for multiple treatment groups
################################################################################


#' Use gsynth to fit factor model with 
#' @importFrom utils capture.output
#' @param long_df A long dataframe with 4 columns in the order unit, time, trt, outcome
#' @param X Matrix of outcomes
#' @param trt Vector of treatment status for each unit
#' @param r Number of factors to use (or start with if CV==1)
#' @param r.end Max number of factors to consider if CV==1
#' @param force Fixed effects (0=none, 1=unit, 2=time, 3=two-way)
#' @param CV Whether to do CV (0=no CV, 1=yes CV)
#' @noRd
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_gsynth_multi <- function(long_df, X, trt, r=0, force=3, CV=1) {
    if(!requireNamespace("gsynth", quietly = TRUE)) {
        stop("In order to fit generalized synthetic controls, you must install the gsynth package.")
    }
    ttot <- ncol(X)
    n <- nrow(X)

    labels <- colnames(long_df)
    gsyn <- gsynth::gsynth(data = long_df, Y = labels[4], D = labels[3], index = c(labels[1], labels[2]), force = force, CV = CV, r=r)
    
    y0hat <- matrix(0, nrow=n, ncol=ttot)
    y0hat[!is.finite(trt),]  <- t(gsyn$Y.co - gsyn$est.co$residuals)
    
    y0hat[is.finite(trt),] <- t(gsyn$Y.ct)
    
    ## add treated prediction for whole pre-period
    gsyn$est.co$Y.ct <- gsyn$Y.ct
    return(list(y0hat=y0hat,
                params=gsyn$est.co))
}



#' Get fixed effects from pre-treatment data for each level
#'
#' @param X Matrix of outcomes
#' @param trt Vector of treatment status for each unit
#' @param mask Matrix of treatment statuses
#' @param force Fixed effects: 1="unit", 2="time", 3="two-way"
#' @noRd
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_feff <- function(X, trt, mask, force) {

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
    if(force %in% c(1,3)) {
      names(residuals) <- as.character(grps)
    residuals <- residuals[as.character(trt[is.finite(trt)])]
    names(y0hat) <- as.character(grps)
    y0hat <- y0hat[as.character(trt[is.finite(trt)])]
    }
    
    return(list(y0hat = y0hat,
                residuals = residuals))
    
}

