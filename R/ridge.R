################################################################################
## Ridge-augmented SCM
################################################################################

#' Ridge augmented weights (possibly with covariates)
#'
#' @param wide_data Output of `format_data`
#' @param synth_data Output of `format_synth`
#' @param Z Matrix of covariates, default is  NULL
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param ridge Include ridge or not
#' @param scm Include SCM or not
#' 
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate (zero in this case)}
#'          \item{"lambda"}{Value of the ridge hyperparameter}
#'          \item{"ridge_mhat"}{The ridge regression predictions (for estimating the bias)}
#'          \item{"synw"}{The synth weights(for estimating the bias)}
#' }
fit_ridgeaug_formatted <- function(wide_data, synth_data,
                                   Z=NULL, lambda=NULL, ridge=T, scm=T, lambda_ratio = 0.00001, lambda_steps = 40) {

    X <- wide_data$X
    y <- wide_data$y
    trt <- wide_data$trt


    ## center outcomes
    X_cent <- apply(X, 2, function(x) x - mean(x[trt==0]))
    X_c <- X_cent[trt==0,,drop=FALSE]
    X_1 <- matrix(colMeans(X_cent[trt==1,,drop=FALSE]), nrow=1)
    y_cent <- apply(y, 2, function(x) x - mean(x[trt==0]))
    y_c <- y_cent[trt==0,,drop=FALSE]

    ## if there are auxiliary covariates, use them
    if(!is.null(Z) & ridge) {
        ## center covariates
        Z_cent <- apply(Z, 2, function(x) x - mean(x[trt==0]))
        Z_c <- Z_cent[trt==0,,drop=FALSE]
        Z_1 <- matrix(colMeans(Z_cent[trt==1,,drop=FALSE]), nrow=1)
        

        ## use CV to choose lambda if it's null
        if(is.null(lambda)) {
            ## do CV with the intercept!
            if(ncol(y) > 1) {
                lambda <- glmnet::cv.glmnet(cbind(1, Z_c), cbind(X_c, y_cent[trt==0,,drop=FALSE]),
                                            alpha=0, family="mgaussian")$lambda.min
            } else {
                lambda <- glmnet::cv.glmnet(cbind(1, Z_c), cbind(X_c, y_cent[trt==0]),
                                            alpha=0, family="mgaussian")$lambda.min
            }        
        }

        ## regress out covariates
        Xc_hat <- Z_c %*% solve(t(Z_c) %*% Z_c + lambda * diag(ncol(Z_c))) %*% t(Z_c) %*% X_c
        X1_hat <- Z_1 %*% solve(t(Z_c) %*% Z_c + lambda * diag(ncol(Z_c))) %*% t(Z_c) %*% X_c

        res_t <- t(X_1  - X1_hat)
        res_c <- t(X_c - Xc_hat)

        new_synth_data <- synth_data
        new_synth_data$Z1 <- res_t
        new_synth_data$X1 <- res_t
        new_synth_data$Z0 <- res_c
        new_synth_data$X0 <- res_c
        ## if SCM fit scm
        if(scm) {
            syn <- fit_synth_formatted(new_synth_data)$weights
        } else {
            ## else use uniform weights
            syn <- rep(1/sum(trt==0), sum(trt==0))
        }
        ## get ridge weights
        ridge_w <- t(t(Z_1) - t(Z_c) %*% syn) %*% solve(t(Z_c) %*% Z_c + lambda * diag(ncol(Z_c))) %*% t(Z_c)
    
        weights <- syn + t(ridge_w)
    } else {
        ## use CV to choose lambda if it's null
        if(ridge) {
            if(is.null(lambda)) {
                if(ncol(y) > 1) {
                    # lambda <- glmnet::cv.glmnet(X_c, y[trt==0,,drop=FALSE], alpha=0, family="mgaussian")$lambda.min
                    lambda_max <- 1 + sqrt(sum((X_1 - t(as.matrix(apply(X_c, 2, mean))))^2))^2
                    scaler <- (lambda_ratio) ^ (1/lambda_steps)
                    lambdas <- lambda_max * (scaler ^ (seq(0:lambda_steps) - 1))
                    lambda_errors <- get_lambda_errors(lambdas, X_c, X_1, synth_data, trt, holdout_length=1)
                    lambda <- lambda_errors[which.min(lambda_errors)]
                    #plot(log(lambdas), get_lambda_errors(lambdas, X_c, X_1, synth_data, trt, holdout_length=10))
                } else {
                    # lambda <- glmnet::cv.glmnet(X_c, y[trt==0], alpha=0, family="gaussian")$lambda.min
                    lambda_max <- 1 + sqrt(sum((X_1 - t(as.matrix(apply(X_c, 2, mean))))^2))^2
                    scaler <- (lambda_ratio) ^ (1/lambda_steps)
                    lambdas <- lambda_max * (scaler ^ (seq(0:lambda_steps) - 1))
                    lambda_errors <- get_lambda_errors(lambdas, X_c, X_1, synth_data, trt, holdout_length=1)
                    lambda <- lambda_errors[which.min(lambda_errors)]
                }
            }
        }
        ## if SCM fit scm
        if(scm) {
            syn <- fit_synth_formatted(synth_data)$weights
        } else {
            ## else use uniform weights
            syn <- rep(1/sum(trt==0), sum(trt==0))
        }

        ## if ridge fit ridge
        if(ridge) {
            ridge_w <- t(t(X_1) - t(X_c) %*% syn) %*% solve(t(X_c) %*% X_c + lambda * diag(ncol(X_c))) %*% t(X_c)
        } else {
            ridge_w <- matrix(0, ncol=sum(trt==0), nrow=1)
        }
        ## combine weights
        weights <- syn + t(ridge_w)
    }

    l2_imbalance <- sqrt(sum((synth_data$X0 %*% weights - synth_data$X1)^2))
    
    ## primal objective value scaled by least squares difference for mean
    uni_w <- matrix(1/ncol(synth_data$X0), nrow=ncol(synth_data$X0), ncol=1)
    unif_l2_imbalance <- sqrt(sum((synth_data$X0 %*% uni_w - synth_data$X1)^2))
    scaled_l2_imabalance <- l2_imbalance / unif_l2_imbalance

    
    ## no outcome model
    mhat <- matrix(0, nrow=nrow(y), ncol=ncol(y))

    if(ridge) {
        if(is.null(Z)) {
            ridge_mhat <- X_cent %*% solve(t(X_c) %*% X_c +
                                           lambda * diag(ncol(X_c))) %*%
                t(X_c) %*% y_c
        } else {
            ridge_mhat <- Z_cent %*% solve(t(Z_c) %*% Z_c +
                                           lambda * diag(ncol(Z_c))) %*%
                t(Z_c) %*% y_c
        }
    } else {
        ridge_mhat <- NULL
    }
    return(list(weights=weights,
                l2_imbalance=l2_imbalance,
                scaled_l2_imbalance=scaled_l2_imabalance,
                mhat=mhat,
                lambda=lambda,
                ridge_mhat=ridge_mhat,
                synw=syn))
}
