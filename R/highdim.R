################################################################################
## Methods to use flexible outcome models
################################################################################

##### Augmented SCM with general outcome models

#' Use zero weights, do nothing but output everything in the right way
#' @param synth_data Panel data in format of Synth::dataprep
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Synth weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#' }
fit_zero_weights <- function(X1, X0) {
    
    ## Imbalance is uniform weights imbalance
    uni_w <- matrix(1 / nrow(X0), nrow = nrow(X0), ncol = 1)
    unif_l2_imbalance <- sqrt(sum((t(X0) %*% uni_w - X1) ^ 2))
    scaled_l2_imbalance <- 1
    
    return(list(weights = rep(0, nrow(X0)),
                l2_imbalance = unif_l2_imbalance,
                scaled_l2_imbalance = scaled_l2_imbalance))
}



#' Fit E[Y(0)|X] and for each post-period and balance pre-period
#'
#' @param wide_data Output of `format_ipw`
#' @param synth_data Output of `synth_data`
#' @param fit_progscore Function to fit prognostic score
#' @param fit_weights Function to fit synth weights
#' @param ... optional arguments for outcome model
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#' }
fit_augsyn_formatted <- function(wide_data, X1, X0,
                                fit_progscore, fit_weights, ...) {


    X <- wide_data$X
    y <- wide_data$y
    trt <- wide_data$trt
    
    ## fit prognostic scores
    fitout <- do.call(fit_progscore,
                          list(X=X, y=y, trt=trt, ...))
    ## fit synth
    syn <- fit_weights(X1, X0)

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
#' @param scm Whether the SCM weighting function is used
#' @param ... optional arguments for outcome model
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#' }
fit_augsyn <- function(wide_data, X1, X0,
                       progfunc=c("EN", "RF", "GSYN", "MCP","CITS", "CausalImpact", "seq2seq"),
                       scm=T, ...) {
    ## prognostic score and weight functions to use
    progfunc = tolower(progfunc)
    if(progfunc == "en") {
        progf <- fit_prog_reg
    } else if(progfunc == "rf") {
        progf <- fit_prog_rf
    } else if(progfunc == "gsyn"){
        progf <- fit_prog_gsynth
    } else if(progfunc == "mcp"){
        progf <- fit_prog_mcpanel
    } else if(progfunc == "cits") {
        progf <- fit_prog_cits
    } else if(progfunc == "causalimpact") {
        progf <- fit_prog_causalimpact
    } else if(progfunc == "seq2seq"){
        progf <- fit_prog_seq2seq
    } else {
        stop("progfunc must be one of 'EN', 'RF', 'GSYN', 'MCP', 'CITS', 'CausalImpact', 'seq2seq'")
    }

    if(scm) {
        weightf <- fit_synth_formatted
    } else {
        ## still fit synth even if none
        ## TODO: This is a dumb wasteful hack
        weightf <- fit_zero_weights
    }
    return(fit_augsyn_formatted(wide_data, X1, X0,
                                progf, weightf, ...))
}