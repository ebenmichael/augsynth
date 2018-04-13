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
    shrink <- max(0, 1 - lam / sqrt(sum(x^2)))
    return(shrink * x)
}

prox_l1 <- function(x, lam) {
    #' prox operator of sum(lam * abs(x))
    out <- (x - lam) * (x > lam) + (x + lam) * (x < -lam)
    return(out)
}



prox_group <- function(x, lams, groups) {
    #' prox operator for group LASSO (generalization of prox_l2)
    #' @param x Vector to prox
    #' @param lams List of prox scalings for each group
    #' @param list of group indices

    ## initialize value
    proxx <- numeric(length(x))

    ## ## go through each group and shrink it
    ## for(i in 1:length(groups)) {
    ##     g <- groups[[i]]
    ##     gname <- names(groups)[i]
    ##     proxx[g] <- prox_l2(x[g], lams[[gname]])
    ## }

    proxvals <- purrr::map2(groups, lams,
                            function(g, lam) {
                                px <- prox_l2(x[g], lam)
                                px
                            })
    proxx[unlist(groups)] <- unlist(proxvals)
    return(proxx)
}


fit_entropy_formatted <- function(data_out, eps, lasso=FALSE, max_iters=1000, tol=1e-8,
                                  init_theta=NULL) {
    #' Fit l2 entropy regularized synthetic controls
    #' by solving the dual problem
    #' @param data_out formatted data from format_entropy
    #' @param eps List of bounds on synthetic control differences for each
    #'            outcome, if only one outcome type then a scalar
    #' @param lasso Whether to do lasso (every covariate is separate)
    #' @param max_iters Maximum number of iterations
    #' @param tol Convergence tolerance
    #' @param init_theta Initial dual parameter to start at (for warm starts), default NULL
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
    if(is.null(data_out$groups)) {
        ## if doing lasso then set the groups to be every covariate
        if(lasso) {
            groups <- as.list(seq(1, dim(x)[2]))
            names(groups) <- paste(seq(1, dim(x)[2]))
            if(length(eps) == 1) {
                epslist <- lapply(groups, function(x) eps)
            } else {
                epslist <- as.list(eps)
                names(epslist) <- names(groups)
            }
        }
        else {
            ## if there is no "groups" field in data_out, assume everything is in one group
            groups <- list()
            groups$"1" <- 1:t
            epslist <- list()
            if(typeof(eps) == "list") {
                epslist$"1" <- eps[[1]]
            } else {
                epslist$"1" <- eps
            }
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

    ## if lasso (every group is separate)
    if(all(sapply(groups, function(g) length(g) == 1))) {
        epsvec <- unlist(epslist)
        prox <- function(lam, step, ...) {
            prox_l1(lam, epsvec * step)
        }
    }

    prox <- function(lam, step, ...) {
        prox_group(lam, lapply(epslist, "*", step), groups)
    }


    ## initial value
    if(is.null(init_theta)) {
        init_theta=numeric(t)
    }
    ## solve the dual problem with prx gradient
    out <- apg::apg(grad, prox, t,
                    opts=list(MAX_ITERS=max_iters, EPS=tol, X_INIT=init_theta))
    lam <- out$x

    ## get the primal weights from the dual variables
    eta <- x %*% lam
    m <- max(eta)
    weights <- exp(eta - m) / sum(exp(eta - m))
    
    ## compute primal objective value
    primal_obj <- sqrt(sum((t(x) %*% weights - y)^2))

    primal_obj <- sqrt(sum((syn_data$Z0 %*% weights - syn_data$Z1)^2))    
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
    tol <- 1e-1
    equalfeasible <- mapply(function(ep, ob, mag)
        if(ep > 0) {
            abs(ep - ob) / ep < tol
        } else {
            abs(ep - ob) < tol^3
        },
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
                           outcome_col=NULL, lasso=FALSE) {
    #' Fit l2 entropy regularized synthetic controls on outcomes
    #' Wrapper around fit_l2_entropy_formatted
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param eps Bound on synthetic control differences
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #' @param lasso Whether to do lasso (every covariate is separate)
    #'
    #' @return Weights for synthetic controls, control outcomes as matrix,
    #'         and whether the unit is actually treated

    ## get the data into the right format
    data_out <- format_data(outcomes, metadata, trt_unit, outcome_col)

    return(fit_entropy_formatted(data_out, eps, lasso))
}


get_entropy <- function(outcomes, metadata, trt_unit=1, eps=NULL,
                        outcome_col=NULL, lasso=FALSE,
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated"),
                        max_iters=1000, tol=1e-8) {
    #' Fit l2_entropy regularized synthetic controls on outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param eps Bound on synthetic control differences
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #' @param lasso Whether to do lasso (every covariate is separate)
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #' @param max_iters Maximum number of iterations
    #' @param tol Convergence tolerance
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## get the synthetic controls weights
    data_out <- format_data(outcomes, metadata, trt_unit, outcome_col, cols)
    out <- fit_entropy_formatted(data_out, eps, lasso, max_iters, tol=tol)

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        data_out$outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        data_out$outcomes <- data_out$outcomes %>% dplyr::arrange_(outcome_col)
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
    ctrls$controls <- out$controls
    return(ctrls)
}
