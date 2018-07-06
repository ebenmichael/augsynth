################################################################################
## Experimenting with other ideas for synth
################################################################################




get_svd_bal <- function(outcomes, metadata, trt_unit=1, r, hyperparam,
                    link=c("logit", "linear", "pos-linear"),
                    regularizer=c(NULL, "l1", "l2", "ridge", "linf"),
                    normalized=TRUE,
                    outcome_col=NULL,
                    cols=list(unit="unit", time="time",
                              outcome="outcome", treated="treated"),
                    opts=list()) {
    #' Find Balancing weights by solving the dual optimization problem
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param r Rank of matrix for SVD
    #' @param hyperparam Regularization hyperparameter
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param normalized Whether to normalize the weights
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #' @param opts Optimization options
    #'        \itemize{
    #'          \item{MAX_ITERS }{Maximum number of iterations to run}
    #'          \item{EPS }{Error tolerance}}
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## format data, reduce dim with SVD and get weights
    data_out <- format_ipw(outcomes, metadata, outcome_col, cols)

    data_out$X <- svd(data_out$X)$u[,1:r,drop=FALSE]
    out <- fit_balancer_formatted(data_out$X, data_out$trt, link, regularizer,
                                  hyperparam, normalized, opts)

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        data_out$outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        data_out$outcomes <- data_out$outcomes %>% dplyr::arrange_(outcome_col)
    }

    syndat <- format_data(outcomes, metadata, trt_unit, outcome_col, cols)
    out$controls <- syndat$synth_data$Y0plot
    ctrls <- impute_controls(syndat$outcomes, out, syndat$trt_unit)
    ctrls$dual <- out$dual
    ctrls$primal_obj <- out$primal_obj
    ctrls$l1_error <- out$l1_error
    ctrls$l2_error <- out$l2_error
    ctrls$linf_error <- out$linf_error
    ctrls$pscores <- out$pscores
    ctrls$eta <- out$eta
    ctrls$feasible <- out$feasible
    ctrls$scaled_primal_obj <- out$scaled_primal_obj
    ctrls$controls <- out$controls
    return(ctrls)
}
