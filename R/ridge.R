################################################################################
## Ridge-augmented SCM
################################################################################



fit_ridgeaug_formatted <- function(ipw_format, syn_format,
                                   lambda=NULL, scm=T) {
    #' Fit E[Y(0)|X] and for each post-period and balance these
    #'
    #' @param ipw_format Output of `format_ipw`
    #' @param syn_format Output of `syn_format`
    #' @param lambda Ridge hyper-parameter, if NULL use CV
    #' @param scm Include SCM or not
    #' 
    #' @return inverse of predicted propensity scores
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated

    X <- ipw_format$X
    y <- ipw_format$y
    trt <- ipw_format$trt

    ## center covariates
    X_cent <- apply(X, 2, function(x) x - mean(x[trt==0]))
    X_c <- X_cent[trt==0,,drop=FALSE]
    X_1 <- colMeans(X_cent[trt==1,,drop=FALSE])
    ## use CV to choose lambda if it's null
    if(is.null(lambda)) {
        if(ncol(y) > 1) {
            lambda <- glmnet::cv.glmnet(X_c, y[trt==0,,drop=FALSE], alpha=0, family="mgaussian")$lambda.min
        } else {
            lambda <- glmnet::cv.glmnet(X_c, y[trt==0], alpha=0, family="gaussian")$lambda.min
        }
    }

    ## if SCM fit scm
    if(scm) {
        syn <- fit_synth_formatted(syn_format)$weights
    } else {
        ## else sue uniform weights
        syn <- rep(1/sum(trt==0), sum(trt==0))
    }

    ## combine weights
    ridge_w <- t(X_1 - t(X_c) %*% syn) %*% solve(t(X_c) %*% X_c + lambda * diag(ncol(X_c))) %*% t(X_c)
    weights <- syn + t(ridge_w)

    data_out <- syn_format$synth_data
    

    primal_obj <- sqrt(sum((data_out$Z0 %*% weights - data_out$Z1)^2))
    ## primal objective value scaled by least squares difference for mean
    x <- t(data_out$Z0)
    y <- data_out$Z1
    unif_primal_obj <- sqrt(sum((t(x) %*% rep(1/dim(x)[1], dim(x)[1]) - y)^2))
    scaled_primal_obj <- primal_obj / unif_primal_obj        

    return(list(weights=weights,
                controls=data_out$Y0plot,
                primal_obj=primal_obj,
                scaled_primal_obj=scaled_primal_obj))
}



get_ridgeaug <- function(outcomes, metadata, trt_unit=1,
                         lambda=NULL, scm=T,
                         cols=list(unit="unit", time="time",
                                   outcome="outcome", treated="treated")) {
    #' Fit synthetic controls on outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param lambda Ridge hyper-parameter, if NULL use CV
    #' @param scm Include SCM or not
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## get the synthetic controls weights
    syn_data <- format_data(outcomes, metadata, trt_unit, cols=cols)
    ipw_data <- format_ipw(outcomes, metadata, cols=cols)
    out <- fit_ridgeaug_formatted(ipw_data, syn_data, lambda, scm)

    ctrls <- impute_controls(syn_data$outcomes, out, syn_data$trt_unit)
    ctrls$primal_obj <- out$primal_obj
    ctrls$scaled_primal_obj <- out$scaled_primal_obj
    
    return(ctrls)
}
