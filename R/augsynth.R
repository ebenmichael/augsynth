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
#' @param form outcome ~ treatment | auxillary covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param t_int Time of intervention
#' @param data Panel data as dataframe
#' @param progfunc What function to use to impute control outcomes
#'                 Ridge=Ridge regression (allows for standard errors),
#'                 None=No outcome model,
#'                 EN=Elastic Net, RF=Random Forest, GSYN=gSynth,
#'                 Comp=softImpute, MCP=MCPanel, CITS=CITS
#'                 CausalImpact=Bayesian structural time series with CausalImpact
#'                 seq2seq=Sequence to sequence learning with feedforward nets
#' @param weightfunc Weighting function to use, default is SCM
#' @param opts_prog Optional options for fitting outcome model
#' @param opts_weights Optional options for fitting synth weights
#' @param cov_agg Covariate aggregation functions, if NULL then use mean with NAs omitted
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
augsynth <- function(form, unit, time, t_int, data,
                     progfunc=c("Ridge", "None", "EN", "RF", "GSYN", "MCP","CITS", "CausalImpact", "seq2seq"),
                     weightfunc=c("SCM", "None"),
                     opts_prog=NULL, opts_weights=NULL,
                     cov_agg=NULL) {

    
    form <- Formula::Formula(form)
    unit <- enquo(unit)
    
    time <- enquo(time)
    
    ## format data
    outcome <- terms(formula(form, rhs=1))[[2]]
    trt <- terms(formula(form, rhs=1))[[3]]
    wide <- format_data(outcome, trt, unit, time, t_int, data)
    synth_data <- do.call(format_synth, wide)

    ## add covariates
    if(length(form)[2] == 2) {

        ## if no aggregation functions, use the mean (omitting NAs)
        cov_agg <- c(function(x) mean(x, na.rm=T))
        
        cov_form <- update(formula(terms(form, rhs=2, data=data)),
                           ~. - 1) ## ensure that there is no intercept

        ## pull out relevant covariates and aggregate
        model.matrix(cov_form,
                     model.frame(cov_form, data, na.action=NULL) ) %>%
            data.frame() %>%
            mutate(unit=pull(data, !!unit)) %>%
            group_by(unit) %>%
            summarise_all(cov_agg) %>%
            select(-unit) %>%
            as.matrix() -> Z
    } else {
        Z <- NULL
    }
    ## fit augsynth
    if(progfunc == "Ridge") {
        if(weightfunc == "SCM") {
            ## Ridge ASCM
            asyn <- do.call(fit_ridgeaug_formatted,
                            c(list(wide_data=wide, synth_data=synth_data, Z=Z),
                              opts_weights))
        } else if(weightfunc == "None") {
            ## Just ridge regression
            asyn <- do.call(fit_ridgeaug_formatted,
                            c(list(wide_data=wide, synth_data=synth_data,
                                   Z=Z, ridge=T, scm=F), opts_weights))
        }
    } else if(progfunc == "None") {
        ## Just SCM
        asyn <- do.call(fit_ridgeaug_formatted,
                        c(list(wide_data=wide, synth_data=synth_data, Z=Z, ridge=F, scm=T), opts_weights))
    } else {
        ## Other outcome models
        asyn <- fit_augsyn(wide, synth_data, progfunc, weightfunc, opts_prog, opts_weights)
    }
    asyn$data <- wide
    asyn$data$time <- data %>% distinct(!!time) %>% pull(!!time)
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


summary.augsynth <- function(augsynth) {

    summ <- list()
    if(augsynth$progfunc == "Ridge") {
        ## get standard errors
        synth_data <- format_synth(augsynth$data$X, augsynth$data$trt,
                                   augsynth$data$y)

        att_se <- loo_se_ridgeaug(augsynth$data, synth_data, lambda=augsynth$lambda)
        att <- data.frame(augsynth$data$time,
                          att_se$att,
                          att_se$se)
        names(att) <- c("Time", "Estimate", "Std. Error")

        summ$att <- att
        summ$sigma <- att_se$sigma
    }

    class(summ) <- "summary.augsynth"
    return(summ)
}
