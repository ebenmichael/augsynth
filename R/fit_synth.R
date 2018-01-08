#######################################################
# Helper scripts to fit synthetic controls to simulations
#######################################################

fit_synth_formatted <- function(data_out) {
    #' Fit synthetic controls on outcomes after formatting data
    #' @param data_out Panel data formatted by Synth::dataprep
    #' 
    is_treated <- data_out$is_treated
    data_out <- data_out$synth_data
    ## change the "predictors" to be the pre period outcomes
    data_out$X0 <- data_out$Z0
    data_out$X1 <- data_out$Z1
    ## set weights on predictors to be 0
    custom.v <- rep(1, dim(data_out$Z0)[1])

    ## fit the weights    
    capture.output(synth_out <- Synth::synth(data_out,
                                             custom.v=custom.v,
                                             quadopt="LowRankQP"))

    weights <- synth_out$solution.w
    loss <- synth_out$loss.w
    primal_obj <- sqrt(sum((data_out$Z0 %*% weights - data_out$Z1)^2))
    
    return(list(weights=weights,
                controls=data_out$Y0plot,
                is_treated=is_treated,
                primal_obj=primal_obj))
}


fit_synth <- function(outcomes, metadata, trt_unit=1) {
    #' Fit synthetic controls on outcomes, wrapper around fit_synth_formatted
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return Weights for synthetic controls, control outcomes as matrix,
    #'         and whether the unit is actually treated

    ## get the data into the right format
    data_out <- format_synth(outcomes, metadata, trt_unit)

    return(fit_synth_formatted(data_out))
}


impute_controls <- function(outcomes, fit, trt_unit) {
    #' Impute the controls after fitting synth
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param fit Output of fit_synth
    #' @param trt_unit Unit that is treated (target for regression)
    #'
    #' @return outcomes with additional synthetic control added and weights


    ## impute the controls
    syn_ctrl <- fit$controls %*% fit$weights
    ## replace true outcome with synthetic control
    syn_outcomes <- outcomes %>%
        filter(unit == trt_unit,
        (potential_outcome == "Y(1)" & treated == TRUE) |
        (potential_outcome == "Y(0)" & treated == FALSE)) %>%
        mutate(outcome = syn_ctrl,
               synthetic = "Y",
               potential_outcome = "Y(0)")
    return(list(outcomes=rbind(outcomes, syn_outcomes),
                weights=fit$weights))    
}


get_synth <- function(outcomes, metadata, trt_unit=1) {
    #' Fit synthetic controls on outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## get the synthetic controls weights
    out <- fit_synth(outcomes, metadata, trt_unit)

    ctrls <- impute_controls(outcomes, out, trt_unit)
    ctrls$primal_obj <- out$primal_obj
    return(ctrls)
}
