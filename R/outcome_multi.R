################################################################################
## Fitting outcome models for multiple treatment groups
################################################################################


#' Use gsynth to fit factor model with 
#' @importFrom utils capture.output
#' @param X Matrix of outcomes
#' @param trt Vector of treatment status for each unit
#' @param r Number of factors to use (or start with if CV==1)
#' @param r.end Max number of factors to consider if CV==1
#' @param force Fixed effects (0=none, 1=unit, 2=time, 3=two-way)
#' @param CV Whether to do CV (0=no CV, 1=yes CV)
#'
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_gsynth_multi <- function(X, trt, r=0, r.end=5, force=3, CV=1) {

    if(!requireNamespace("gsynth", quietly = TRUE)) {
        stop("In order to fit generalized synthetic controls, you must install the gsynth package.")
    }


    ttot <- ncol(X)
    n <- nrow(X)
    
    ## observed matrix (everything observed)
    I <- matrix(1, ttot, n)

    ## treatment matrix
    trt_mat <- matrix(0, nrow=n, ncol=ttot)
    trt_mat[is.finite(trt),] <-
        t(vapply(trt[is.finite(trt)],
               function(ti) c(rep(0, ti), rep(1, ttot-ti)),
               numeric(ttot)))

    ## use internal gsynth function
    capture.output(gsyn <- gsynth:::synth.core(t(X), NULL, t(trt_mat), I,
                                               r=r, r.end=r.end,
                                               force=force, CV=CV,
                                               tol=0.001))
    ## get predicted outcomes
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
#' 
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_feff <- function(X, trt, mask, force) {

    ttot <- ncol(X)
    n <- nrow(X)
    grps <- trt[is.finite(trt)]
    J <- length(grps)
    which_t <- (1:n)[is.finite(trt)]        

    residuals <- lapply(1:J, function(j) X)
    if(force %in% c(2,3)) {
        ## compute time fixed effects from pure controls
        time_eff <- lapply(1:J, function(j) matrix(colMeans(X[!is.finite(trt),]),
                                                            nrow=nrow(X), ncol=ncol(X),
                                                            byrow=T))
        residuals <- lapply(1:J, function(j) residuals[[j]] - time_eff[[j]])
        y0hat <- time_eff
    }

    if(force %in% c(1,3)) {

        ## compute unit fixed effects from pre-intervention outcomes
        unit_eff <- lapply(1:J, function(j) matrix(rowMeans(residuals[[j]][, mask[j,]==1]),
                                                   nrow=nrow(X), ncol=ncol(X)))

        residuals <- lapply(1:J, function(j) residuals[[j]] - unit_eff[[j]])
        y0hat <- lapply(1:J, function(j) unit_eff[[j]])
    }

    if(force == 3) {
        y0hat <- lapply(1:J, function(j) time_eff[[j]] + unit_eff[[j]])
    }
    
    return(list(y0hat=y0hat,
                residuals=residuals))
    
}

