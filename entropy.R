#############################################################
## Entropy regularized synthetic controls (experimental)
#############################################################
## helper log sum exp function
logsumexp <- function(x0) {
    m <- max(x0)
    return(log(sum(exp(x0 - m))) + m)
}


fit_entropy_formatted <- function(data_out, alpha=NULL) {
    #' Fit entropy regularized synthetic controls
    #' by solving the dual problem
    #' @param data_out formatted data from format_entropy
    #' @param alpha regularization parameter
    #'
    #' @return synthetic control weights

    syn_data <- data_out$synth_data

    ## data for objective and gradient
    y <- syn_data$Z1
    x <- t(syn_data$Z0)

    n <- dim(x)[1]
    t <- dim(x)[2]

    ## set alpha to n / 100 
    if(is.null(alpha)) {
        alpha <- n / 100
    }
    
    ## dual objective function
    obj <- function(lam) {
        obj1 <- logsumexp(-x %*% lam)
        obj2 <- t(y) %*% lam
        reg <- 1 / (4 * alpha) * sum(lam ^2) # regularization

        return(obj1 + obj2 + reg)
    }

    ## dual gradient
    grad <- function(lam) {
        eta <- x %*% lam
        ## compute weights for each value
        num <- colSums(as.numeric(exp(-eta)) * x)
        denom <- sum(exp(-eta))

        grad1 <- - num / denom
        grad2 <- y
        reg <- 1 / (2 * alpha) * lam #regularization

        return(grad1 + grad2 + reg)
    }

    ## initial value
    lam0 <- numeric(t)
    
    ## solve the optimization problem to get the dual variables
    out <- optim(lam0, obj, grad, method="L-BFGS-B")
    lam <- out$par

    ## get the primal weights from the dual variables
    eta <- x %*% lam
    weights <- exp(-eta) / sum(exp(-eta))

    return(list(weights=weights,
                dual=lam,
                controls=syn_data$Y0plot,
                is_treated=data_out$is_treated))
    
}


fit_entropy <- function(outcomes, metadata, trt_unit=1, alpha=NULL) {
    #' Fit entropy regularized synthetic controls on outcomes
    #' Wrapper around fit_synth_formatted
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha regularization parameter
    #'
    #' @return Weights for synthetic controls, control outcomes as matrix,
    #'         and whether the unit is actually treated

    ## get the data into the right format
    data_out <- format_synth(outcomes, metadata, trt_unit)

    return(fit_entropy_formatted2(data_out, alpha))
}


get_entropy <- function(outcomes, metadata, trt_unit=1, alpha=NULL) {
    #' Fit entropy regularized synthetic controls on outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha regularization parameter
    #'
    #' @return outcomes with additional synthetic control added and weights

    ## get the synthetic controls weights
    out <- fit_entropy2(outcomes, metadata, trt_unit, alpha)

    return(impute_controls(outcomes, out, trt_unit))
}


fit_l1_entropy_formatted <- function(data_out, eps=NULL) {
    #' Fit entropy regularized synthetic controls with l1 error bounds
    #' by solving the dual problem
    #' @param data_out formatted data from format_entropy
    #' @param eps width of error box around treated unit pre-treatment outcomes
    #'
    #' @return synthetic control weights

    syn_data <- data_out$synth_data

    ## data for objective and gradient
    y <- syn_data$Z1
    x <- t(syn_data$Z0)

    n <- dim(x)[1]
    t <- dim(x)[2]

    ## set eps to standard error of mean / t
    if(is.null(eps)) {
        eps <- sqrt(sum(colMeans(sweep(x, 2, colMeans(x))^2) / n)) 
    }    
    
    print(eps)
    
    ## dual objective function
    obj <- function(lam) {
        obj1 <- logsumexp(-x %*% lam)
        obj2 <- t(y) %*% lam
        reg <- sum(abs(lam) * eps)

        return(obj1 + obj2 + reg)
    }

    ## dual sub-gradient
    grad <- function(lam) {
        eta <- x %*% lam
        ## compute weights for each value
        num <- colSums(as.numeric(exp(-eta)) * x)
        denom <- sum(exp(-eta))

        grad1 <- - num / denom
        grad2 <- y
        reg <- sum(sign(lam) * eps) #regularization

        return(grad1 + grad2 + reg)
    }

    ## initial value
    lam0 <- numeric(t)
    
    ## solve the optimization problem to get the dual variables
    out <- optim(lam0, obj, grad, method="CG")
    lam <- out$par

    ## get the primal weights from the dual variables
    eta <- x %*% lam
    weights <- exp(-eta) / sum(exp(-eta))

    return(list(weights=weights,
                dual=lam,
                controls=syn_data$Y0plot,
                is_treated=data_out$is_treated))
    
}


fit_l1_entropy <- function(outcomes, metadata, trt_unit=1, eps=NULL) {
    #' Fit entropy regularized synthetic controls with l1 error bounds
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param eps width of error box around treated unit pre-treatment outcomes
    #'
    #' @return Weights for synthetic controls, control outcomes as matrix,
    #'         and whether the unit is actually treated

    ## get the data into the right format
    data_out <- format_synth(outcomes, metadata, trt_unit)

    return(fit_l1_entropy_formatted(data_out, eps))
}


get_l1_entropy <- function(outcomes, metadata, trt_unit=1, eps=NULL) {
    #' Fit entropy regularized synthetic controls on outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param eps width of error box around treated unit pre-treatment outcomes
    #'
    #' @return outcomes with additional synthetic control added and weights

    ## get the synthetic controls weights
    out <- fit_l1_entropy(outcomes, metadata, trt_unit, eps)
    ctrls <- impute_controls(outcomes, out, trt_unit)
    ctrls$dual <- out$dual
    return(ctrls)
}
