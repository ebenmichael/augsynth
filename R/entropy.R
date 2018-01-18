#############################################################
## Entropy regularized synthetic controls (experimental)
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
    ## if there is no "groups" field in data_out, assume everything is in one
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
        groups <- data_out$groups
        epslist <- eps
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
    primal_obj <- lapply(groups,
                         function(g) sqrt(sum((t(x[,g]) %*% weights - y[g,])^2)))
    ## compute propensity scores
    pscores <- 1 / (1 + exp(-eta))
    return(list(weights=weights,
                dual=lam,
                controls=syn_data$Y0plot,
                is_treated=data_out$is_treated,
                primal_obj=primal_obj,
                groups=groups,
                pscores=pscores,
                eta=eta))
    
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
    out <- fit_entropy(outcomes, metadata, trt_unit, eps, outcome_col)

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        outcomes <- outcomes %>% dplyr::arrange_(outcome_col)
    }
    ctrls <- impute_controls(outcomes, out, trt_unit)
    ctrls$dual <- out$dual
    ctrls$primal_obj <- out$primal_obj
    ctrls$pscores <- out$pscores
    ctrls$eta <- out$eta
    ctrls$groups <- out$groups
    return(ctrls)
}


### double robust estimation with synth weights

fit_ipw_formatted <- function(data_out, alpha_w) {
    #' Fit IPW weights with a logit propensity score model
    #' @param data_out formatted data from format_entropy
    #' @param alpha_w regularization parameter for weights
    #'
    #' @return inverse of predicted propensity scores
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated 
    
    syn_data <- data_out$synth_data
    ## get data into the right form
    x0 <- t(syn_data$Z0) # controls
    x1 <- t(syn_data$Z1) # treated
    x <- rbind(x0, x1)
    n <- dim(x0)[1]
    t <- dim(x0)[2]

    ## treatment indicator
    trt <- c(rep(0, n), 1)

    ## fit regression
    fit <- LiblineaR::LiblineaR(x, trt, 7, cost=alpha_w, bias=0)

    ## get predicted probabilities P(W=1|X)
    pscores <- predict(fit, x, proba=TRUE)$probabilities[1:n,2]
    weights <- pscores / (1-pscores)
    weights <- weights / sum(weights)
    primal_obj <- sum((weights %*% x0 - x1)^2)

    
    return(list(weights=weights,
                outparams=t(fit$W),
                controls=syn_data$Y0plot,
                treated=syn_data$Y1plot,
                primal_obj=primal_obj,
                is_treated=data_out$is_treated,
                dual=fit$W,
                pscores=pscores))
}


fit_ipw <- function(outcomes, metadata, trt_unit=1, alpha_w) {
    #' Fit IPW weights with a logit propensity score model
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha_w regularization parameter for weights
    #'
    #' @return inverse of predicted propensity scores
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated 

    ## format data
    data_out <- format_data(outcomes, metadata, trt_unit)

    ## get weights
    return(fit_ipw_formatted(data_out, alpha_w))
}

get_ipw <- function(outcomes, metadata, trt_unit=1, alpha_w=1) {
    #' Fit IPW weights with a logit propensity score model
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha_w regularization parameter for weights
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## get the synthetic controls weights
    out <- fit_ipw(outcomes, metadata, trt_unit, alpha_w)

    ctrls <- impute_controls(outcomes, out, trt_unit)
    ctrls$outparams <- out$outparams
    ctrls$primal_obj <- out$primal_obj
    ctrls$pscores <- out$pscores
    ctrls$dual <- out$dual
    return(ctrls)
}



fit_dr_formatted <- function(data_out, eps_w, alpha_o=NULL, syn=TRUE) {
    #' Fit a regularized outcome model and synthetic controls
    #' for a double robust estimator
    #' @param data_out formatted data from format_entropy
    #' @param eps_w List of bounds on synthetic control differences for each
    #'            outcome, if only one outcome type then a scalar
    #' @param alpha_o regularization parameter for outcome model, if NULL use CV
    #' @param syn whether to use synthetic control weights in DR estimate
    #'
    #' @return synthetic control weights,
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated 

    syn_data <- data_out$synth_data

    if(syn) {
        ## fit entropy regularized synthetic control weights
        ws <- fit_l2_entropy_formatted(data_out, eps_w)
    } else {
        ws <- fit_pscore_formatted(data_out, alpha_w)
    }

    ## fit regularized regression for outcomes for each post period
    x <- t(syn_data$Z0)
    n <- dim(x)[1]
    t <- dim(x)[2]
    ys <- t(syn_data$Y0[(t+1):dim(syn_data$Y0),])
    if(!is.null(alpha_o)) {
        outfit <- function(y) {
            fit <- glmnet::glmnet(x, y, alpha=0,
                                  lambda=alpha_o,
                                  intercept=FALSE)
            return(coef(fit)[-1,])
        }
    } else {
        outfit <- function(y) {
            alpha_o <- glmnet::cv.glmnet(x, y, alpha=0, intercept=FALSE)$lambda.min
            fit <- glmnet::glmnet(x, y, alpha=0,
                                  lambda=alpha_o,
                                  intercept=FALSE)
            return(coef(fit)[-1,])
        }
    }
    regweights <- apply(ys, 2,outfit)

    return(list(weights=ws$weights,
                dual=ws$dual,
                outparams=t(regweights),
                controls=syn_data$Y0plot,
                treated=syn_data$Y1plot,
                is_treated=data_out$is_treated,
                primal_obj=ws$primal_obj,
                pscores=ws$pscores))
    
    }



fit_dr <- function(outcomes, metadata, trt_unit=1, eps_w, alpha_o=NULL) {
    #' Fit a regularized outcome model and synthetic controls
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param eps_w List of bounds on synthetic control differences for each
    #'            outcome, if only one outcome type then a scalar
    #' @param alpha_o regularization parameter for outcome model, if NULL use CV
    #'
    #' @return synthetic control weights,
    #'         outcome regression parameters
    #'         control outcomes
    #'         boolean for treated
    
    ## get the data into the right format
    data_out <- format_synth(outcomes, metadata, trt_unit)

    return(fit_dr_formatted(data_out, eps_w, alpha_o))
    
}



get_dr <- function(outcomes, metadata, trt_unit=1, eps_w, alpha_o=NULL) {
    #' Fit a regularized outcome model and synthetic controls
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param eps_w List of bounds on synthetic control differences for each
    #'            outcome, if only one outcome type then a scalar
    #' @param alpha_o regularization parameter for outcome model, if NULL use CV
    #'
    #' @return synthetic control weights,
    #'         outcome regression parameters
    #'         control outcomes
    #'         boolean for treated
    #' @export
    
    ## get the data into the right format
    data_out <- format_synth(outcomes, metadata, trt_unit)

    ## fit outcome regression and weights
    fit <- fit_dr(outcomes, metadata, trt_unit, eps_w, alpha_o)
    
    ## compute the DR estimate to "impute" the controls
    ctrls <- impute_dr(outcomes, metadata, fit, trt_unit) 
    ctrls$dual <- fit$dual
    ctrls$primal_obj <- fit$primal_obj
    ctrls$pscores <- fit$pscores
    
    return(ctrls)
}



impute_dr <- function(outcomes, metadata, fit, trt_unit) {
    #' Impute the controls after fitting a dr estimator
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe with metadata, in particular a t_int column
    #' @param fit Output of fit_dr
    #' @param trt_unit Unit that is treated (target for regression)
    #'
    #' @return outcomes with additional synthetic control added,
    #'         synth weights
    #'         outcome regression weights
    


    ### weight the residuals

    ## separate out pre and post period controls
    t <- dim(fit$controls)[1]
    n <- dim(fit$contrsol)[2]
    
    t_int <- match(metadata$t_int, rownames(fit$controls))

    preC <- fit$controls[1:(t_int-1),]
    postC <- fit$control[(t_int):t,]

    ## and pre and post period treated
    preT <- fit$treated[1:(t_int - 1)]
    postT <- fit$treated[(t_int):t]
    
    ## find the residuals and weight them by synth weights
    resid <- postC - fit$outparams %*% preC
    wresid <- resid %*% fit$weights

    ## predict expected control outcome for treated
    mu0 <- fit$outparams %*% preT

    ## combine into DR estimate
    dr <- mu0 + wresid

    ## combine pre period with DR estimate into a "synthetic control"
    dr_ctrl <- c(preT, dr)
    
    ## replace true outcome with imputed value
    dr_outcomes <- outcomes %>%
        filter(unit == trt_unit,
        (potential_outcome == "Y(1)" & treated == TRUE) |
        (potential_outcome == "Y(0)" & treated == FALSE)) %>%
        mutate(outcome = dr_ctrl,
               synthetic = "DR",
               potential_outcome = "Y(0)")

    ## also include the synthetic control
    syn_ctrl <- fit$controls %*% fit$weights
    syn_outcomes <- outcomes %>%
        filter(unit == trt_unit,
        (potential_outcome == "Y(1)" & treated == TRUE) |
        (potential_outcome == "Y(0)" & treated == FALSE)) %>%
        mutate(outcome = syn_ctrl,
               synthetic = "Y",
               potential_outcome = "Y(0)")

    return(list(outcomes=rbind(outcomes, dr_outcomes, syn_outcomes),
                weights=fit$weights,
                dual=fit$dual,
                outparams=fit$outparams))
}


#### Choosing alpha

choosealpha <- function(outcomes, metadata, trt_unit=1, alphas) {
    #' Pick the hyperparameter which gives the best pre-period fit
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alphas regularization parameters to pick from
    #'
    #' @return Best alpha

    losses <- sapply(alphas,
                     function(alpha) fit_entropy(outcomes,
                                                 metadata,
                                                 trt_unit,
                                                 alpha)$primal_obj)
    return(alphas[which.min(losses)])
}
