################################################################################
## Main functions for augmented synthetic controls Method
################################################################################

#' augsynth: A package implementing the Augmented Synthetic Controls Method
#' @docType package
#' @name augsynth
#' @import magrittr
#' @import dplyr
#' @import LowRankQP
#' @import ggplot2
#' @import tidyr
#' @import LiblineaR
#' @import glmnet
NULL


#' Fit Augmented SCM
#' @param formula outcome ~ treatment
#' @param unit Name of unit column
#' @param time Name of time column
#' @param t_int Time of intervention
#' @param data Panel data as dataframe
#' @param progfunc Outcome model to us, default is ridge regression
#' @param weightfunc Weighting function to use, default is SCM
#' @param opts.prog Optional options for fitting outcome model
#' @param opts.weights Optional options for fitting synth weights
#'
#' @return augsynth object that contains:
#'         \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#'          \item{"data"}{Panel data as matrices}
#'         }
#' @export
augsynth <- function(formula, unit, time, t_int, data,
                     progfunc="Ridge", weightfunc="SCM",
                     opts.prog=NULL, opts.weights=NULL) {

    unit <- enquo(unit)
    time <- enquo(time)
    
    ## format data
    outcome <- terms(formula)[[2]]
    trt <- terms(formula)[[3]]
    wide <- format_data(outcome, trt, unit, time, t_int, data)
    synth_data <- do.call(format_synth, wide)
    
    ## fit augsynth
    if(progfunc == "Ridge") {
        asyn <- do.call(fit_ridgeaug_formatted,
                        c(list(wide_data=wide, synth_data=synth_data), opts.weights))
    } else {
        asyn <- fit_augsyn(wide, synth_data, progfunc, weightfunc, opts.prog, opts.weights)
    }
    asyn$data <- wide
    asyn$progfunc <- progfunc
    asyn$weightfunc <- weightfunc
    ##format output
    class(asyn) <- "augsynth"
    return(asyn)
}

#' Get prediction of average outcome under control
#' @param asyn augsynth object
#'
#' @return Vector of predicted post-treatment control averages
predict.augsynth <- function(asyn) {

    y <- asyn$data$y
    trt <- asyn$data$trt
    mhat <- asyn$mhat
    
    m1 <- colMeans(mhat[trt==1,,drop=F])

    resid <- (y[trt==0,,drop=F] - mhat[trt==0,drop=F])

    return(m1 + t(resid) %*% asyn$weights)
}
