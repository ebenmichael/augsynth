#############################################################
## Entropy regularized synthetic controls
#############################################################

### helper functions

logsumexp <- function(x0) {
    #' Compute numerically stable logsumexp
    m <- max(x0)
    val <- log(sum(exp(x0 - m))) + m
    return(val)
}


logsumexp_grad <- function(eta, x) {
    #' Compute numerically stable logsumexp gradient with natural param eta
    #' and data x
    m <- max(-eta)
    num <- colSums(as.numeric(exp(-eta - m)) * x)
    denom <- sum(exp(-eta - m))
    return(-num / denom)

}


prox_l2 <- function(x, lam) {
    #' prox operator of lam * ||x||_2
    shrink <- max(0, 1 - lam / norm(x, type="2"))
    return(shrink * x)
}



prox_group <- function(x, lams, groups) {
    #' prox operator for group LASSO (generalization of prox_l2)
    #' @param x Vector to prox
    #' @param lams List of prox scalings for each group
    #' @param list of group indices

    ## initialize value
    proxx <- numeric(length(x))

    ## go through each group and shrink it
    for(i in 1:length(groups)) {
        g <- groups[[i]]
        gname <- names(groups)[i]
        proxx[g] <- prox_l2(x[g], lams[[gname]])
    }
    return(proxx)
}


fit_entropy_formatted <- function(data_out, eps) {
    #' Fit l2 entropy regularized synthetic controls
    #' by solving the dual problem
    #' @param data_out formatted data from format_entropy
    #' @param eps List of bounds on synthetic control differences for each
    #'            outcome, if only one outcome type then a scalar
    #'
    #' @return synthetic control weights

    syn_data <- data_out$synth_data

    ## data for objective and gradient
    y <- syn_data$Z1
    x <- t(syn_data$Z0)
    
    n <- dim(x)[1]
    t <- dim(x)[2]

    ## use prox gradient method to minimize f(x) + h(x)
    
    ## f(x)
    obj <- function(lam, ..) {
        obj1 <- logsumexp(x %*% lam)
        obj2 <- -t(y) %*% lam

        return(obj1 + obj2)
    }

    ## f'(x)
    grad <- function(lam, ...) {
        ## numerically stable gradient of log sum exp
        eta <- -x %*% lam
        grad1 <- -logsumexp_grad(eta, x)
        
        grad2 <- -y
        
        return(grad1 + grad2)
    }


    ## prox_h
    ## if there is no "groups" field in data_out, assume everything is in one group
    if(is.null(data_out$groups)) {
        groups <- list()
        groups$"1" <- 1:t
        epslist <- list()
        if(typeof(eps) == "list") {
            epslist$"1" <- eps[[1]]
        } else {
            epslist$"1" <- eps
        }
    } else {
        ## get the groups
        groups <- data_out$groups
        ## if eps is a scalar, then assume the same eps for each group
        if(length(eps) == 1) {
            epslist <- lapply(groups, function(x) eps)
        } else {
            epslist <- eps
        }
        ##groups <- groups[names(epslist)]
        epslist <- epslist[names(groups)]
    }
    prox <- function(lam, step, ...) {
        prox_group(lam, lapply(epslist, "*", step), groups)
    }


    ## initial value
    lam0 <- numeric(t)
    ## solve the dual problem with prx gradient
    
    out <- apg::apg(grad, prox, t, opts=list(MAX_ITERS=5000))
    lam <- out$x

    ## get the primal weights from the dual variables
    eta <- x %*% lam
    m <- max(eta)
    weights <- exp(eta - m) / sum(exp(eta - m))
    
    ## compute primal objective value
    primal_obj <- sqrt(sum((t(x) %*% weights - y)^2))
    
    ## primal objective value scaled by least squares difference for mean
    unif_primal_obj <- sqrt(sum((t(x) %*% rep(1/dim(x)[1], dim(x)[1]) - y)^2))
    scaled_primal_obj <- primal_obj / unif_primal_obj

    ## compute l2 error per group
    primal_group_obj <- lapply(groups,
                               function(g) sqrt(sum((t(x[,g]) %*% weights - y[g,])^2)))
    ## compute propensity scores
    pscores <- 1 / (1 + exp(-eta))

    ## get magnitude of vectors
    mags <- lapply(groups, function(g) sqrt(sum(y[g,]^2)))
    ## check for equality within 10^-3 * magnitude
    tol <- 10^-4
    equalfeasible <- mapply(function(ep, ob, mag)
        isTRUE(all.equal(ep, ob, tol * mag)),
        epslist, primal_group_obj, mags)
    lessfeasible <- mapply(function(ep, ob)
        ob < ep, 
        epslist, primal_group_obj)
    feasible <- mapply(function(a,b) a || b, equalfeasible, lessfeasible)
    feasible <- all(feasible)
    
    return(list(weights=weights,
                dual=lam,
                controls=syn_data$Y0plot,
                is_treated=data_out$is_treated,
                primal_obj=primal_obj,
                scaled_primal_obj=scaled_primal_obj,
                primal_group_obj=primal_group_obj,
                groups=groups,
                pscores=pscores,
                eta=eta,
                feasible=feasible))
    
}


fit_entropy <- function(outcomes, metadata, trt_unit=1, eps=NULL,
                           outcome_col=NULL) {
    #' Fit l2 entropy regularized synthetic controls on outcomes
    #' Wrapper around fit_l2_entropy_formatted
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param eps Bound on synthetic control differences
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #'
    #' @return Weights for synthetic controls, control outcomes as matrix,
    #'         and whether the unit is actually treated

    ## get the data into the right format
    data_out <- format_data(outcomes, metadata, trt_unit, outcome_col)

    return(fit_entropy_formatted(data_out, eps))
}


get_entropy <- function(outcomes, metadata, trt_unit=1, eps=NULL,
                           outcome_col=NULL) {
    #' Fit l2_entropy regularized synthetic controls on outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param eps Bound on synthetic control differences
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## get the synthetic controls weights
    data_out <- format_data(outcomes, metadata, trt_unit, outcome_col)
    out <- fit_entropy_formatted(data_out, eps)

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        outcomes <- outcomes %>% dplyr::arrange_(outcome_col)
    }
    
    ctrls <- impute_controls(data_out$outcomes, out, data_out$trt_unit)
    ctrls$dual <- out$dual
    ctrls$primal_obj <- out$primal_obj
    ctrls$pscores <- out$pscores
    ctrls$eta <- out$eta
    ctrls$groups <- out$groups
    ctrls$feasible <- out$feasible
    ctrls$primal_group_obj <- out$primal_group_obj
    ctrls$scaled_primal_obj <- out$scaled_primal_obj
    return(ctrls)
}
