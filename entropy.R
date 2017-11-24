#############################################################
## Entropy regualrized synthetic controls (experimental)
#############################################################
library(alabama)


fit_entropy_formatted <- function(data_out, lam) {
    #' Fit entropy regularized synthetic controls
    #' @param data_out fomratted data from format_entropy
    #' @param lam regularization parameter, defaults to 1 / log(n)
    #'
    #' @return synthetic control weights


    syn_data <- data_out$synth_data
    ## create simplex constraints
    ## const %*% weights - val >= 0
    n <- dim(syn_data$Z0)[2]
    t <- dim(syn_data$Z0)[1]
    ## non-negative constraint
    const <- diag(rep(1, n))
    val <- rep(0,n)
    ## sum to one
    const <- rbind(const, rep(-1, n), rep(1, n))
    val <- c(val, -1.01, .99)

    ## inequality constraint function
    hin <- function(w) return(w)
    hin.jac <- function(w) return(matrix(1, length(w), length(w)))

    ## equality constraint function
    heq <- function(w) return(sum(w) - 1)
    heq.jac <- function(w) return(rep(1, length(w)))

    ## set regualrization parameter to 1 / log(n) if it isn't set
    if(is.null(lam)) {
        lam <- 1 / log(n)
    }
    
    ## objective function    
    y <- syn_data$Z1
    x <- syn_data$Z0
    obj <- function(w) {
        return(sum((y - x %*% w)^2) / (2 * n) +
               lam * sum(w * log(w)))
    }

    ## gradient
    grad <- function(w) {
        return(t(x) %*% (x %*% w - y) / n + lam * (log(w) + 1))
    }


    ## initialize at unregularized synthetic control weights
    syn_w <- fit_synth_formatted(data_out)$weights
    init_w <- rep(1/n, n)

    weights <- auglag(syn_w, obj, grad, hin, hin.jac, heq, heq.jac,
                  control.outer=list(trace=FALSE, kkt2.check=FALSE))$par
    ## solve the optimization problem
    #out <- list()
    #out$weights <- constrOptim(init_w, obj, grad, const, val, method="CG",
    #                   control=list(trace=2))$par

    return(list(weights=weights,
                controls=syn_data$Y0plot,
                is_treated=data_out$is_treated))
    
}


fit_entropy <- function(outcomes, metadata, trt_unit=1, lam=NULL) {
    #' Fit entropy regularized synthetic controls on outcomes
    #' Wrapper around fit_synth_formatted
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param lam Regularization parameter, defaults to 1 / log(n)
    #'
    #' @return Weights for synthetic controls, control outcomes as matrix,
    #'         and whether the unit is actually treated

    ## get the data into the right format
    data_out <- format_synth(outcomes, metadata, trt_unit)

    return(fit_entropy_formatted(data_out, lam))
}


get_entropy <- function(outcomes, metadata, trt_unit=1, lam=NULL) {
    #' Fit entropy regularized synthetic controls on outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param lam Regularization parameter, defaults to 1 / log(n)
    #'
    #' @return outcomes with additional synthetic control added and weights

    ## get the synthetic controls weights
    out <- fit_entropy(outcomes, metadata, trt_unit, lam)

    return(impute_controls(outcomes, out, trt_unit))
}
