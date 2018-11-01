################################################################################
## Main functions for augmented synthetic controls Method
################################################################################

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
augsynth <- function(formula, unit, time, t_int, data,
                     progfunc="Ridge", weightfunc="SCM",
                     opts.prog=list(), opts.weights=list()) {

    unit <- enquo(unit)
    time <- enquo(time)
    
    ## format data
    outcome <- terms(formula)[[2]]
    trt <- terms(formula)[[3]]
    wide <- format_data(outcome, trt, unit, time, t_int, data)
    synth_data <- do.call(format_synth, wide)
    
    ## fit augsynth
    asyn <- fit_ridgeaug_formatted(wide, synth_data)

    asyn$data <- wide
    
    ##format output
    class(asyn) <- "augsynth"
    return(asyn)
}
