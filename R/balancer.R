################################################################################
## Implement general problem with balancer package
################################################################################

fit_balancer_formatted <- function(X, trt,
                                   link=c("logit", "linear", "pos-linear"),
                                   regularizer=c(NULL, "l1", "l2", "ridge", "linf"),
                                   hyperparam, normalized=TRUE,
                                   opts=list()) {
    #' Find Balancing weights by solving the dual optimization problem
    #' @param X n x d matrix of covariates
    #' @param trt Vector of treatment status indicators
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param hyperparam Regularization hyperparameter
    #' @param normalized Whether to normalize the weights
    #' @param opts Optimization options
    #'        \itemize{
    #'          \item{MAX_ITERS }{Maximum number of iterations to run}
    #'          \item{EPS }{Error rolerance}}
    #'
    #' @return \itemize{
    #'          \item{theta }{Estimated dual propensity score parameters}
    #'          \item{weights }{Estimated primal weights}
    #'          \item{imbalance }{Imbalance in covariates}}

    ## fit weights with balancer
    out <- balancer::balancer(X, trt, type="att",
                              link=link, regularizer=regularizer,
                              hyperparam=hyperparam, normalized=normalized,
                              opts=opts)

    ## drop treated weights
    weights <- out$weights[trt == 0]
    x <- X[trt == 0,,drop=FALSE]
    y <- colMeans(X[trt == 1,,drop=FALSE])
    ## compute l1, l2, and linf error
    l2_error <- sqrt(sum((t(x) %*% weights - y)^2))
    l1_error <- sum(abs(t(x) %*% weights - y))
    linf_error <- max(abs(t(x) %*% weights - y))

    if(is.null(regularizer)) {
      primal_obj <- l2_error  
    } else if(regularizer=="l2") {
        primal_obj <- l2_error
    } else if(regularizer=="l1") {
        primal_obj <- linf_error
    } else if(regularizer == "linf") {
        primal_obj <- l1_error
    } else {
        primal_obj <- l2_error
    }

    ## primal objective value scaled by least squares difference for mean
    unif_primal_obj <- sqrt(sum((t(x) %*% rep(1/dim(x)[1], dim(x)[1]) - y)^2))
    scaled_primal_obj <- primal_obj / unif_primal_obj

    eta <- x %*% out$theta
    ## compute propensity scores
    pscores <- 1 / (1 + exp(-eta))

    ## get magnitude of vectors
    mag <- sqrt(sum(y^2))
    ## check for equality within 10^-3 * magnitude
    tol <- 1e-1
    if(hyperparam > 0) {
        equalfeasible <- abs(hyperparam - primal_obj) / hyperparam < tol
    } else {
        equalfeasible <- abs(hyperparam - primal_obj)  < tol^3
    }
    lessfeasible <- primal_obj < hyperparam
    feasible <- equalfeasible || lessfeasible
    
    return(list(weights=weights,
                dual=out$theta,
                controls=t(x),
                primal_obj=primal_obj,
                l1_error=l1_error,
                l2_error=l2_error,
                linf_error=linf_error,
                scaled_primal_obj=scaled_primal_obj,
                pscores=pscores,
                eta=eta,
                feasible=feasible))    
}


get_balancer <- function(outcomes, metadata, trt_unit=1, hyperparam,
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

    ## get the synthetic controls weights
    data_out <- format_ipw(outcomes, metadata, outcome_col, cols)
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
    ctrls <- impute_controls(syndat$outcomes, out, trt_unit)
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
