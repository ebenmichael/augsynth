#############################################################
## Entropy regularized synthetic controls (experimental)
#############################################################
library(glmnet)
source("fit_synth.R")


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
        alpha <- 100 / n
    }
    
    ## dual objective function
    obj <- function(lam) {
        obj1 <- logsumexp(-x %*% lam)
        obj2 <- t(y) %*% lam
        reg <- alpha * sum(lam ^2) # regularization

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
        reg <- alpha * lam #regularization

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

    return(fit_entropy_formatted(data_out, alpha))
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
    out <- fit_entropy(outcomes, metadata, trt_unit, alpha)

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


### double robust estimation with synth weights

fit_dr_formatted <- function(data_out, alpha_w, alpha_o) {
    #' Fit a regularized outcome model and synthetic controls
    #' for a double robust estimator
    #' @param data_out formatted data from format_entropy
    #' @param alpha_w regularization parameter for weights
    #' @param alpha_o regularization parameter for outcome model
    #'
    #' @return synthetic control weights,
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated 

    syn_data <- data_out$synth_data

    ## fit entropy regularized synthetic control weights
    ent <- fit_entropy_formatted(data_out, alpha_w)

    ## fit regularized regression for outcomes for each post period
    x <- t(syn_data$Z0)
    n <- dim(x)[1]
    t <- dim(x)[2]
    ys <- t(syn_data$Y0[(t+1):dim(syn_data$Y0),])
    regweights <- apply(ys, 2,
                        function(y) {
                            fit <- glmnet(x, y, alpha=0,
                                          lambda=alpha_o, intercept=FALSE)
                            return(coef(fit)[-1,])
                        }
                        )

    return(list(weights=ent$weights,
                outparams=t(regweights),
                controls=syn_data$Y0plot,
                treated=syn_data$Y1plot,
                is_treated=data_out$is_treated))
    
    }



fit_dr <- function(outcomes, metadata, trt_unit=1, alpha_w, alpha_o) {
    #' Fit a regularized outcome model and synthetic controls
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha_w regularization parameter for weights
    #' @param alpha_o regularization parameter for outcome model
    #'
    #' @return synthetic control weights,
    #'         outcome regression parameters
    #'         control outcomes
    #'         boolean for treated
    
    ## get the data into the right format
    data_out <- format_synth(outcomes, metadata, trt_unit)

    return(fit_dr_formatted(data_out, alpha_w, alpha_o))
    
}


get_dr <- function(outcomes, metadata, trt_unit=1, alpha_w, alpha_o) {
    #' Fit a regularized outcome model and synthetic controls
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha_w regularization parameter for weights
    #' @param alpha_o regularization parameter for outcome model
    #'
    #' @return synthetic control weights,
    #'         outcome regression parameters
    #'         control outcomes
    #'         boolean for treated
    
    ## get the data into the right format
    data_out <- format_synth(outcomes, metadata, trt_unit)

    ## fit outcome regression and weights
    fit <- fit_dr(outcomes, metadata, trt_unit, alpha_w, alpha_o)

    ## compute the DR estimate to "impute" the controls
    return(impute_dr(outcomes, metadata, fit, trt_unit))
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

    ## predict expected c ontrol outcome for treated
    mu0 <- fit$outparams %*% preT

    ## combine pre period with DR estimate into a "synthetic control"
    dr_ctrl <- c(preT, mu0)
    
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
                outparams=fit$outparams))
}


##### Cross validation
#predfun <- function(Xtrain, Ytrain, Xtest, Ytest, ...) {
#    ## package train data into right format
#    data_out <- list(synth_data=list(Z0=Xtrain, Z1=Ytrain))
#    out <- fit_entropy_formatted(data_out, ...)
#
#    print(length(out$weights))
#    print(length(out$dual))
#    ## predict
#    print(dim(Xtest))
#    print(dim(matrix(out$weights)))
#    pred <- Xtest %*% matrix(out$weights, ncol=1)
#    return(sum((Ytest - pred)^2))
#}
