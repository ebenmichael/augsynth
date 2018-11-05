################################################################################
## Methods to use flixbile outcome models
################################################################################

#' Fit E[Y(0)|X] and for each post-period and balance these
#'
#' @param wide_data Output of `format_data`
#' @param synth_data Output of `format_synth`
#' @param fit_progscore Function to fit prognostic score
#' @param fit_weights Function to fit synth weights
#' @param opts_prog Optional options for fitting prognostic score
#' @param opts_weights Optional options for fitting synth weights
#' 
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate (zero in this case)}
#' }
fit_progsyn_formatted <- function(wide_data, synth_data,
                                  fit_progscore, fit_weights,
                                  opts_prog=NULL, opts_weights=NULL) {

    X <- wide_data$X
    y <- wide_data$y
    trt <- wide_data$trt
    
    ## fit prognostic scores
    if(is.null(opts_prog)) {
        fitout <- fit_progscore(X, y, trt)
    } else {
        fitout <- do.call(fit_progscore,
                          c(list(X=X, y=y, trt=trt), opts_prog))
    }

    y0hat <- fitout$y0hat

    ## replace outcomes with fitted prognostic scores
    synth_data$Z0 <- t(as.matrix(y0hat[wide_data$trt == 0,,drop=FALSE]))
    synth_data$Z1 <- as.matrix(colMeans(as.matrix(y0hat[wide_data$trt == 1,,drop=FALSE])))

    ## fit synth/maxent weights
    if(is.null(opts_weights)) {
        syn <- fit_weights(synth_data)
    }
    syn <- do.call(fit_weights, c(list(synth_data=synth_data), opts_weights))

    syn$params <- fitout$params

    ## no outcome model
    syn$mhat <- matrix(0, nrow=nrow(y), ncol=ncol(y))
    
    return(syn)
}


#' Fit E[Y(0)|X] and for each post-period and balance these
#'
#' @param wide_data Output of `format_data`
#' @param synth_data Output of `format_synth`
#' @param progfunc What function to use to impute control outcomes
#'                 EN=Elastic Net, RF=Random Forest, GSYN=gSynth
#'                 CITS=Comparative interrupted time series
#'                 CausalImpact=Bayesian structural time series with CausalImpact
#'                 seq2seq=Sequence to sequence learning with feedforward nets
#' @param weightfunc What function to use to fit weights
#'                   SCM=Vanilla Synthetic Controls
#' @param fit_weights Function to fit synth weights
#' @param opts_prog Optional options for fitting prognostic score
#' @param opts_weights Optional options for fitting synth weights
#' 
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate (zero in this case)}
#' }
fit_progsyn <- function(wide_data, synth_data,
                        progfunc=c("EN", "RF", "GSYN", "CITS", "CausalImpact", "seq2seq"),
                        weightfunc=c("SCM"),
                        opts_prog=NULL, opts_weights=NULL) {
    ## prognostic score and weight functions to use
    if(progfunc == "EN") {
        progf <- fit_prog_reg
    } else if(progfunc == "RF") {
        progf <- fit_prog_rf
    } else if(progfunc == "GSYN"){
        progf <- fit_prog_gsynth
    } else if(progfunc == "CITS"){
        progf <- fit_prog_cits
    } else if(progfunc == "CausalImpact"){
        progf <- fit_prog_causalimpact
    } else if(progfunc == "seq2seq"){
        progf <- fit_prog_seq2seq
    } else {
        stop("progfunc must be one of 'EN', 'RF', 'GSYN', 'CITS', 'CausalImpact', 'seq2seq'")
    }

    
    if(weightfunc == "SCM") {
        weightf <- fit_synth_formatted
    } else {
        stop("weightfunc must be 'SCM'")
    }

    return(fit_progsyn_formatted(wide_data, synth_data,
                                 fit_progscore, fit_weights))
}

##### Augmented SCM with general outcome models

#' Fit E[Y(0)|X] and for each post-period and balance pre-period
#'
#' @param wide_data Output of `format_ipw`
#' @param synth_data Output of `synth_data`
#' @param fit_progscore Function to fit prognostic score
#' @param fit_weights Function to fit synth weights
#' @param opts_prog Optional options for fitting prognostic score
#' @param opts_weights Optional options for fitting synth weights
#' 
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#' }
fit_augsyn_formatted <- function(wide_data, synth_data,
                                fit_progscore, fit_weights,
                                opts_prog=NULL, opts_weights=NULL) {


    X <- wide_data$X
    y <- wide_data$y
    trt <- wide_data$trt
    
    ## fit prognostic scores
    if(is.null(opts_prog)) {
        fitout <- fit_progscore(X, y, trt)
    } else {
        fitout <- do.call(fit_progscore,
                          c(list(X=X, y=y, trt=trt),
                            opts_prog))
    }
    
    ## fit synth
    if(is.null(opts_weights)) {        
        syn <- fit_weights(synth_data)
    } else {
        syn <- do.call(fit_weights,
                       c(list(synth_data=synth_data),
                         opts_weights))
    }

    syn$params <- fitout$params

    syn$mhat <- fitout$y0hat
    
    return(syn)
}


#' Fit outcome model and balance pre-period
#' @param wide_data Output of `format_ipw`
#' @param synth_data Output of `synth_data`
#' @param progfunc What function to use to impute control outcomes
#'                 EN=Elastic Net, RF=Random Forest, GSYN=gSynth,
#'                 Comp=softImpute, MCP=MCPanel, CITS=CITS
#'                 CausalImpact=Bayesian structural time series with CausalImpact
#'                 seq2seq=Sequence to sequence learning with feedforward nets
#' @param weightfunc What function to use to fit weights
#'                   SC=Vanilla Synthetic Controls, ENT=Maximum Entropy
#' @param opts_prog Optional options for fitting prognostic score
#' @param opts_weights Optional options for fitting synth weights
#' 
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#' }
fit_augsyn <- function(wide_data, synth_data,
                       progfunc=c("EN", "RF", "GSYN", "MCP","CITS", "CausalImpact", "seq2seq"),
                       weightfunc=c("SCM"),
                       opts_prog = NULL,
                       opts_weights = NULL) {

    ## prognostic score and weight functions to use
    if(progfunc == "EN") {
        progf <- fit_prog_reg
    } else if(progfunc == "RF") {
        progf <- fit_prog_rf
    } else if(progfunc == "GSYN"){
        progf <- fit_prog_gsynth
    } else if(progfunc == "MCP"){
        progf <- fit_prog_mcpanel
    } else if(progfunc == "CITS") {
        progf <- fit_prog_cits
    } else if(progfunc == "CausalImpact") {
        progf <- fit_prog_causalimpact
    } else if(progfunc == "seq2seq"){
        progf <- fit_prog_seq2seq
    } else {
        stop("progfunc must be one of 'EN', 'RF', 'GSYN', 'MCP', 'CITS', 'CausalImpact', 'seq2seq'")
    }

    if(weightfunc == "SCM") {
        weightf <- fit_synth_formatted
    } else if(weightfunc == "NONE") {
        ## still fit synth even if none
        ## TODO: This is a dumb wasteful hack
        weightf <- fit_synth_formatted
    } else {
        stop("weightfunc must be one of `SCM`, `NONE`")
    }
    return(fit_augsyn_formatted(wide_data, synth_data,
                                progf, weightf,
                                opts_prog, opts_weights))
}



### Combine synth and gsynth by balancing pre-period residuals
#' Fit outcome model and balance residuals
#'
#' @param wide_data Output of `format_data`
#' @param synth_data Output of `format_synth`
#' @param fit_progscore Function to fit prognostic score
#' @param fit_weights Function to fit synth weights
#' @param opts.gsyn Optional options for gsynth
#' @param opts_weights Optional options for fitting synth weights
#' 
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#' }
fit_residaug_formatted <- function(wide_data, synth_data,
                                  fit_progscore, fit_weights,
                                  opts_prog=NULL, opts_weights=NULL) {


    X <- wide_data$X
    y <- wide_data$y
    trt <- wide_data$trt

    ## fit prognostic scores
    if(is.null(opts_prog)) {
        fitout <- fit_progscore(X, y, trt)
    } else {
        fitout <- do.call(fit_progscore,
                          c(list(X=X, y=y, trt=trt),
                            opts_prog))
    }

    
    y0hat <- fitout$y0hat

    ## get residuals
    ctrl_resids <- fitout$params$ctrl_resids
    trt_resids <- fitout$params$trt_resids
    
    ## replace outcomes with pre-period residuals
    t0 <- dim(X)[2]

    synth_data$Z0 <- ctrl_resids[1:t0, ]
    synth_data$Z1 <- as.matrix(trt_resids[1:t0])
    
    ## fit synth weights
    if(is.null(opts_weights)) {
        syn <- fit_weights(synth_data)
    } else {
        syn <- do.call(fit_weights,
                       c(list(synth_data=synth_data),
                         opts_weights))
    }

    syn$params <- fitout$params    

    ## return predicted values for treatment and control
    syn$mhat <- y0hat
    
    return(syn)
}
#' Fit outcome model and balance residuals
#'
#' @param wide_data Output of `format_data`
#' @param synth_data Output of `format_synth`
#' @param progfunc What function to use to impute control outcomes
#'                 GSYN=gSynth, MCP=MCPanel,
#'                 CITS=Comparative interrupted time series
#'                 CausalImpact=Bayesian structural time series with CausalImpact
#' @param weightfunc What function to use to fit weights
#'                   SCM=Vanilla Synthetic Controls
#'                   NONE=No reweighting, just outcome model
#' @param opts.gsyn Optional options for gsynth
#' @param opts_weights Optional options for fitting synth weights
#' 
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#' }
fit_residaug <- function(wide_data, synth_data,
                        progfunc=c("GSYN", "MCP", "CITS", "CausalImpact"),
                        weightfunc=c("SC","ENT", "SVD", "NONE"),
                        opts_prog = NULL,
                        opts_weights = NULL) {

    ## prognostic score and weight functions to use
    if(progfunc == "GSYN"){
        progf <- fit_prog_gsynth
    } else if(progfunc == "MCP"){
        progf <- fit_prog_mcpanel
    } else if(progfunc == "CITS") {
        progf <- fit_prog_cits
    } else if(progfunc == "CausalImpact") {
        progf <- fit_prog_causalimpact
    } else {
        stop("progfunc must be one of 'GSYN', 'MCP', 'CITS', 'CausalImpact'")
    }

    
    ## weight function to use
    if(weightfunc == "SCM") {
        weightf <- fit_synth_formatted
    } else if(weightfunc == "NONE") {
        ## still fit synth even if none
        ## TODO: This is a dumb wasteful hack
        weightf <- fit_synth_formatted
    } else {
        stop("weightfunc must be one of 'SCM', 'NONE'")
    }

    return(fit_residaug_formatted(wide_data, synth_data,
                                  progf, weightf,
                                  opts_prog, opts_weights))
}

