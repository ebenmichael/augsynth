################################################################################
## Implementations of other estimators (IPW, DR, Entropy Balancing)
################################################################################


fit_ipw_formatted <- function(data_out, alpha_w=NULL, net_param=1) {
    #' Fit IPW weights with a logit propensity score model
    #' @param data_out formatted data from format_ipw
    #' @param alpha_w regularization parameter for weights
    #' @param net_param elastic net parameter, 1: lasso, 0: ridge, default: 1
    #'
    #' @return inverse of predicted propensity scores
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated 
    
    ## get covariates and treatment indicator
    x <- data_out$X
    trt <- data_out$trt
    y <- data_out$y

    if(is.null(alpha_w)){
        cv <- glmnet::cv.glmnet(x, trt, family="binomial",
                                alpha=net_param,
                                intercept=FALSE)
        alpha_w <- cv$lambda.min
        
    }

    fit <- glmnet::glmnet(x, trt, family="binomial",
                          alpha=net_param,
                          lambda=alpha_w, intercept=FALSE)

    ## get predicted probabilities P(W=1|X)
    pscores <- predict(fit, x[trt==0,,drop=FALSE], type="response")
    weights <- pscores / (1-pscores)
    weights <- weights / sum(weights)
    x1 <- apply(x[trt == 1,,drop=FALSE], 2, mean)
    primal_obj <- sum((t(weights) %*% x[trt==0,,drop=FALSE] - x1)^2)

    unif_primal_obj <- sqrt(sum((t(x[trt==0,,drop=FALSE]) %*%
                                 rep(1/dim(x[trt==0,,drop=FALSE])[1],
                                     dim(x[trt==0,,drop=FALSE])[1]) - x1)^2))
    scaled_primal_obj <- primal_obj / unif_primal_obj
    
    return(list(weights=as.numeric(weights),
                controls=cbind(x[trt==0,,drop=FALSE],y[trt == 0,,drop=FALSE]),
                treated=cbind(x[trt == 1,,drop=FALSE],y[trt == 1,,drop=FALSE]),
                trt=trt,
                primal_obj=primal_obj,
                scaled_primal_obj=scaled_primal_obj,
                is_treated=data_out$is_treated,
                dual=fit$beta,
                pscores=pscores,
                reg_param=alpha_w,
                net_param=net_param))
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
    data_out <- format_ipw(outcomes, metadata)

    ## get weights
    return(fit_ipw_formatted(data_out, alpha_w))
}

get_ipw <- function(outcomes, metadata, alpha_w=NULL, net_param=1) {
    #' Fit IPW weights with a logit propensity score model
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha_w regularization parameter for weights
    #' @param net_param elastic net parameter, 1: lasso, 0: ridge, default: 1
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## get the synthetic controls weights
    data_out <- format_ipw(outcomes, metadata)
    out <- fit_ipw_formatted(data_out, alpha_w, net_param)

    ctrls <- impute_dr(data_out$outcomes, metadata, out)
    ctrls$outparams <- out$outparams
    ctrls$primal_obj <- out$primal_obj
    ctrls$scaled_primal_obj <- out$scaled_primal_obj
    ctrls$pscores <- out$pscores
    ctrls$dual <- out$dual
    ctrls$reg_param <- out$reg_param
    ctrls$net_param <- out$net_param
    return(ctrls)
}



fit_dr_formatted <- function(data_out, alpha_w=NULL, alpha_o=NULL, net_param=1) {
    #' Fit a regularized outcome model and synthetic controls
    #' for a double robust estimator
    #' @param data_out formatted data from format_ipw
    #' @param alpha_w regularization parameter for p score, if NULL use CV
    #' @param alpha_o regularization parameter for outcome model, if NULL use CV
    #' @param net_param elastic net parameter, 1: lasso, 0: ridge, default: 1
    #'
    #' @return synthetic control weights,
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated 


    ws <- fit_ipw_formatted(data_out, alpha_w, net_param)

    ## fit regularized regression for outcomes for each post period
    x <- data_out$X
    trt <- data_out$trt
    ys <- data_out$y
    
    if(!is.null(alpha_o)) {
        outfit <- function(y) {
            fit <- glmnet::glmnet(x, y, alpha=net_param,
                                  lambda=alpha_o,
                                  intercept=FALSE)
            return(coef(fit)[-1,])
        }
    } else {
        outfit <- function(y) {
            alpha_o <- glmnet::cv.glmnet(x, y, alpha=net_param,
                                         intercept=FALSE)$lambda.min
            fit <- glmnet::glmnet(x, y, alpha=0,
                                  lambda=alpha_o,
                                  intercept=FALSE)
            return(coef(fit)[-1,])
        }
    }
    regweights <- apply(ys, 2,outfit)

    return(list(weights=ws$weights,
                dual=ws$dual,
                outparams=regweights,
                controls=cbind(x[trt==0,,drop=FALSE],ys[trt == 0,,drop=FALSE]),
                treated=cbind(x[trt == 1,,drop=FALSE],ys[trt == 1,,drop=FALSE]),
                is_treated=data_out$is_treated,
                primal_obj=ws$primal_obj,
                scaled_primal_obj=ws$scaled_primal_obj,
                pscores=ws$pscores,
                p_reg_param=ws$reg_param,
                o_reg_param=alpha_o,
                net_param=net_param))
    
    }



fit_dr <- function(outcomes, metadata, trt_unit=1, alpha_w=NULL, alpha_o=NULL) {
    #' Fit a regularized outcome model and synthetic controls
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha_w regularization parameter for p score, if NULL use CV
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



get_dr <- function(outcomes, metadata, alpha_w=NULL, alpha_o=NULL,
                   net_param=1) {
    #' Fit a regularized outcome model and synthetic controls
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha_w regularization parameter for p score, if NULL use CV
    #' @param alpha_o regularization parameter for outcome model, if NULL use CV
    #' @param net_param elastic net parameter, 1: lasso, 0: ridge, default: 1
    #'
    #' @return synthetic control weights,
    #'         outcome regression parameters
    #'         control outcomes
    #'         boolean for treated
    #' @export
    
    ## get the data into the right format
    data_out <- format_ipw(outcomes, metadata)

    ## fit outcome regression and weights
    fit <- fit_dr_formatted(data_out, alpha_w, alpha_o, net_param)
    
    ## compute the DR estimate to "impute" the controls
    ctrls <- impute_dr(data_out$outcomes, metadata, fit) 
    ctrls$dual <- fit$dual
    ctrls$primal_obj <- fit$primal_obj
    ctrls$scaled_primal_obj <- fit$scaled_primal_obj
    ctrls$pscores <- fit$pscores
    ctrls$p_reg_param <- fit$p_reg_param
    ctrls$o_reg_param <- fit$o_reg_param
    ctrls$net_param <- fit$net_param
    
    return(ctrls)
}



impute_dr <- function(outcomes, metadata, fit) {
    #' Impute the controls after fitting a dr estimator
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe with metadata, in particular a t_int column
    #' @param fit Output of fit_dr
    #'
    #' @return outcomes with additional synthetic control added,
    #'         synth weights
    #'         outcome regression weights
    


    ### weight the residuals
    t <- dim(fit$controls)[2]
    ## separate out pre and post period controls
    
    t_int <- (metadata$t_int - min(outcomes$time) + 1)

    preC <- fit$controls[,1:(t_int-1)]
    postC <- fit$controls[,(t_int):t]

    ## and pre and post period treated
    preT <- fit$treated[,1:(t_int-1), drop=FALSE]
    postT <- fit$treated[,(t_int):t, drop=FALSE]

    ## find the residuals and weight them by synth weights
    if(is.null(fit$outparams)) {
        outparams <- matrix(0, nrow=dim(preC)[2], ncol=dim(postT)[2])
    } else {
        outparams <- fit$outparams
    }

    resid <- postC - preC %*% outparams

    wresid <- t(resid) %*% fit$weights

    ## predict expected control outcome for treated
    mu0 <- colMeans(preT %*% outparams)

    ## combine into DR estimate
    dr <- mu0 + wresid

    if(is.null(fit$outparams)) {
        ## if no outcome regression combine with weighted avg of pre-period
        dr_ctrl <- c(t(preC) %*% fit$weights, wresid)
    } else{
        ## combine pre period with DR estimate into a "synthetic control"
        dr_ctrl <- c(colMeans(preT), dr)
    }

    dr_ctrl <- c(t(preC) %*% fit$weights, wresid)
    ## replace true outcome with imputed value
    dr_outcomes <- outcomes %>%
        filter(unit == -1) %>%
        mutate(outcome = dr_ctrl,
               synthetic = "DR",
               potential_outcome = "Y(0)") %>% data.frame()

    ctrls <- outcomes %>% filter(!treated) %>% data.frame()
    avgs <- outcomes %>% filter(unit == -1) %>% data.frame()

    finalout <- bind_rows(ctrls, avgs, dr_outcomes)
    #finalout$outcome <- c(ctrls$outcome, avgs$outcome, dr_outcomes$outcome)
    return(list(outcomes=finalout,
                weights=fit$weights,
                dual=fit$dual,
                outparams=fit$outparams))
}


## entropy balancing
fit_ebal_formatted <- function(data_out) {
    #' Fit entropy balancing weights
    #' @param data_out formatted data from format_ipw
    #'
    #' @return entropy balancing weights
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated 
    
    ## get covariates and treatment indicator
    x <- data_out$X
    trt <- data_out$trt
    y <- data_out$y

    bal <- ebal::ebalance(trt, x, constraint.tolerance = 10^(-5))

    ## get predicted probabilities P(W=1|X)
    pscores <- bal$w / (1 + bal$w)
    weights <- bal$w / sum(bal$w)
    x1 <- apply(x[trt == 1,,drop=FALSE], 2, mean)
    primal_obj <- sum((t(weights) %*% x[trt==0,,drop=FALSE] - x1)^2)

    unif_primal_obj <- sqrt(sum((t(x[trt==0,,drop=FALSE]) %*%
                                 rep(1/dim(x[trt==0,,drop=FALSE])[1],
                                     dim(x[trt==0,,drop=FALSE])[1]) - x1)^2))
    scaled_primal_obj <- primal_obj / unif_primal_obj
    
    return(list(weights=as.numeric(weights),
                controls=cbind(x[trt==0,,drop=FALSE],y[trt == 0,,drop=FALSE]),
                treated=cbind(x[trt == 1,,drop=FALSE],y[trt == 1,,drop=FALSE]),
                trt=trt,
                primal_obj=primal_obj,
                scaled_primal_obj=scaled_primal_obj,
                is_treated=data_out$is_treated,
                dual=bal$coef,
                pscores=pscores))
}



get_ebal <- function(outcomes, metadata) {
    #' Fit entropy balancing weights
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #'
    #' @return outcomes with imputed control values, weights, dual variables
    #' @export

    ## get the synthetic controls weights
    data_out <- format_ipw(outcomes, metadata)
    out <- fit_ebal_formatted(data_out)

    ctrls <- impute_dr(data_out$outcomes, metadata, out)
    ctrls$outparams <- out$outparams
    ctrls$primal_obj <- out$primal_obj
    ctrls$scaled_primal_obj <- out$scaled_primal_obj
    ctrls$pscores <- out$pscores
    ctrls$dual <- out$dual
    return(ctrls)
}


fit_uniform_formatted <- function(data_out) {
    #' Use difference in means
    #' @param data_out Output of format_ipw
    #'
    #' @return uniform weights
    #' @export

    ## get covariates and treatment indicator
    x <- data_out$X
    trt <- data_out$trt
    y <- data_out$y


    ## get predicted probabilities P(W=1|X)
    weights <- rep(1/dim(x[trt==0,,drop=FALSE])[1], dim(x[trt==0,,drop=FALSE])[1])
    x1 <- apply(x[trt == 1,,drop=FALSE], 2, mean)
    primal_obj <- sum((t(weights) %*% x[trt==0,,drop=FALSE] - x1)^2)

    unif_primal_obj <- sqrt(sum((t(x[trt==0,,drop=FALSE]) %*%
                                 rep(1/dim(x[trt==0,,drop=FALSE])[1],
                                     dim(x[trt==0,,drop=FALSE])[1]) - x1)^2))
    scaled_primal_obj <- primal_obj / unif_primal_obj
    
    return(list(weights=as.numeric(weights),
                controls=cbind(x[trt==0,,drop=FALSE],y[trt == 0,,drop=FALSE]),
                treated=cbind(x[trt == 1,,drop=FALSE],y[trt == 1,,drop=FALSE]),
                trt=trt,
                primal_obj=primal_obj,
                scaled_primal_obj=scaled_primal_obj,
                is_treated=data_out$is_treated)
           )
}


get_uniform <- function(outcomes, metadata) {
    #' Use difference in means
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #'
    #' @return outcomes with imputed control values, weights, dual variables
    #' @export

    ## get the synthetic controls weights
    data_out <- format_ipw(outcomes, metadata)
    out <- fit_uniform_formatted(data_out)

    ctrls <- impute_dr(data_out$outcomes, metadata, out)
    ctrls$outparams <- out$outparams
    ctrls$primal_obj <- out$primal_obj
    ctrls$scaled_primal_obj <- out$scaled_primal_obj
    return(ctrls)
}
