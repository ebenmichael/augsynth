#######################################################
# Helper scripts to fit synthetic controls to simulations
#######################################################

#' Fit synthetic controls on outcomes after formatting data
#' @param synth_data Panel data in format of Synth::dataprep
#'
#' @return \itemize{
#'          \item{"weights"}{Synth weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#' }
fit_synth_formatted <- function(synth_data, V = NULL) {

    if(!require("LowRankQP")) {
        stop("In order to use Synth, you must install LowRankQP")
    }
    ## if no  is supplied, set equal to 1
    if(is.null(V)) {
        custom.v <- rep(1, dim(synth_data$Z0)[1])
    } else {
        custom.v <- V
    }
    print(synth_data$Z1)
    print(custom.v)
    ## fit the weights    
    capture.output(synth_out <- Synth::synth(synth_data,
                                             custom.v=custom.v,
                                             quadopt="LowRankQP"))
    weights <- synth_out$solution.w
    loss <- synth_out$loss.w
    l2_imbalance <- sqrt(sum((synth_data$Z0 %*% weights - synth_data$Z1)^2))
    
    ## primal objective value scaled by least squares difference for mean
    uni_w <- matrix(1/ncol(synth_data$Z0), nrow=ncol(synth_data$Z0), ncol=1)
    unif_l2_imbalance <- sqrt(sum((synth_data$Z0 %*% uni_w - synth_data$Z1)^2))
    scaled_l2_imbalance <- l2_imbalance / unif_l2_imbalance
    
    return(list(weights=weights,
                l2_imbalance=l2_imbalance,
                scaled_l2_imbalance=scaled_l2_imbalance))
}
