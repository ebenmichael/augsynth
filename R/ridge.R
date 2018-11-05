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
#' }
fit_ridgeaug_formatted <- function(wide_data, synth_data,
                                   Z=NULL, lambda=NULL, ridge=T, scm=T) {

    X <- wide_data$X
    y <- wide_data$y
    trt <- wide_data$trt


    ## center outcomes
    X_cent <- apply(X, 2, function(x) x - mean(x[trt==0]))
    X_c <- X_cent[trt==0,,drop=FALSE]
    X_1 <- matrix(colMeans(X_cent[trt==1,,drop=FALSE]), nrow=1)
    y_cent <- apply(y, 2, function(x) x - mean(x[trt==0]))

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
        new_synth_data$Z0 <- res_c
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
                    lambda <- glmnet::cv.glmnet(X_c, y[trt==0,,drop=FALSE], alpha=0, family="mgaussian")$lambda.min
                } else {
                    lambda <- glmnet::cv.glmnet(X_c, y[trt==0], alpha=0, family="gaussian")$lambda.min
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

    l2_imbalance <- sqrt(sum((synth_data$Z0 %*% weights - synth_data$Z1)^2))
    
    ## primal objective value scaled by least squares difference for mean
    uni_w <- matrix(1/ncol(synth_data$Z0), nrow=ncol(synth_data$Z0), ncol=1)
    unif_l2_imbalance <- sqrt(sum((synth_data$Z0 %*% uni_w - synth_data$Z1)^2))
    scaled_l2_imabalance <- l2_imbalance / unif_l2_imbalance

    ## no outcome model
    mhat <- matrix(0, nrow=nrow(y), ncol=ncol(y))
    return(list(weights=weights,
                l2_imbalance=l2_imbalance,
                scaled_l2_imabalance=scaled_l2_imabalance,
                mhat=mhat,
                lambda=lambda))
}
