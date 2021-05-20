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
#' @param lambda_min_ratio Ratio of the smallest to largest lambda when tuning lambda values
#' @param n_lambda Number of lambdas to consider between the smallest and largest lambda value
#' @param lambda_max Initial (largest) lambda, if NULL sets it to be (1+norm(X_1-X_c))^2
#' @param holdout_length Length of conseuctive holdout period for when tuning lambdas 
#' @param min_1se If TRUE, chooses the maximum lambda within 1 standard error of the lambda that minimizes the CV error, if FALSE chooses the optimal lambda; default TRUE
#' @param V V matrix for synth, default NULL
#' @param residualize Whether to residualize auxiliary covariates or balance directly, default TRUE
#' @param ... optional arguments for outcome model
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate (zero in this case)}
#'          \item{"lambda"}{Value of the ridge hyperparameter}
#'          \item{"ridge_mhat"}{The ridge regression predictions (for estimating the bias)}
#'          \item{"synw"}{The synth weights(for estimating the bias)}
#'          \item{"lambdas"}{List of lambda values evaluated to tune ridge regression}
#'          \item{"lambda_errors"}{"The MSE associated with each lambda term in lambdas."}
#'          \item{"lambda_errors_se"}{"The SE of the MSE associated with each lambda term in lambdas."}
#' }
fit_ridgeaug_formatted <- function(wide_data, synth_data,
                                   Z=NULL, lambda=NULL, ridge=T, scm=T,
                                   lambda_min_ratio = 1e-8, n_lambda = 20,
                                   lambda_max = NULL,
                                   holdout_length = 1, min_1se = T,
                                   V = NULL,
                                   residualize = FALSE, ...) {
    extra_params = list(...)
    if (length(extra_params) > 0) {
        warning("Unused parameters in using ridge augmented weights: ", paste(names(extra_params), collapse = ", "))
    }

    X <- wide_data$X
    y <- wide_data$y
    trt <- wide_data$trt

    lambda_errors <- NULL
    lambda_errors_se <- NULL
    lambdas <- NULL

    ## center outcomes
    X_cent <- apply(X, 2, function(x) x - mean(x[trt==0]))
    X_c <- X_cent[trt==0,,drop=FALSE]
    X_1 <- matrix(colMeans(X_cent[trt==1,,drop=FALSE]), nrow=1)
    y_cent <- apply(y, 2, function(x) x - mean(x[trt==0]))
    y_c <- y_cent[trt==0,,drop=FALSE]

    t0 <- ncol(X_c)

    V <- make_V_matrix(t0, V)

    # apply V matrix transformation
    X_c <- X_c %*% V
    X_1 <- X_1 %*% V

    new_synth_data <- synth_data


    ## if there are auxiliary covariates, use them
    if(!is.null(Z)) {
        ## center covariates
        Z_cent <- apply(Z, 2, function(x) x - mean(x[trt==0]))
        Z_c <- Z_cent[trt==0,,drop=FALSE]
        Z_1 <- matrix(colMeans(Z_cent[trt==1,,drop=FALSE]), nrow=1)

        if(residualize) {
          ## regress out covariates
          Xc_hat <- Z_c %*% solve(t(Z_c) %*% Z_c) %*% t(Z_c) %*% X_c
          X1_hat <- Z_1 %*% solve(t(Z_c) %*% Z_c) %*% t(Z_c) %*% X_c

          # take residuals
          res_t <- X_1  - X1_hat
          res_c <- X_c - Xc_hat

          X_c <- res_c
          X_1 <- res_t

          X_cent[trt == 0,] <- res_c
          X_cent[trt == 1,] <- res_t


          new_synth_data$Z1 <- t(res_t)
          new_synth_data$X1 <- t(res_t)
          new_synth_data$Z0 <- t(res_c)
          new_synth_data$X0 <- t(res_c)
        } else {
            # standardize covariates to be on the same scale as the outcomes
            sdz <-  apply(Z_c, 2, sd)
            sdx <- sd(X_c)
            Z_c <- sdx * t(t(Z_c) / sdz)
            Z_1 <- sdx * Z_1 / sdz

          # concatenate
          X_c <- cbind(X_c, Z_c)
          X_1 <- cbind(X_1, Z_1)
          new_synth_data$Z1 <- t(X_1)
          new_synth_data$X1 <- t(X_1)
          new_synth_data$Z0 <- t(X_c)
          new_synth_data$X0 <- t(X_c)
          V <- diag(ncol(X_c))
        }
    } else {
        new_synth_data$Z1 <- t(X_1)
        new_synth_data$X1 <- t(X_1)
        new_synth_data$Z0 <- t(X_c)
        new_synth_data$X0 <- t(X_c)
    }
    out <- fit_ridgeaug_inner(X_c, X_1, trt, new_synth_data,
                               lambda, ridge, scm,
                               lambda_min_ratio, n_lambda,
                               lambda_max,
                               holdout_length, min_1se)

    weights <- out$weights
    synw <- out$synw
    lambda <- out$lambda
    lambdas <- out$lambdas
    lambda_errors <- out$lambda_errors
    lambda_errors_se <- out$lambda_errors_se

    # add back in covariate weights
    if(!is.null(Z)) {
        if(residualize) {
          no_cov_weights <- weights
          ridge_w <- t(t(Z_1) - t(Z_c) %*% weights) %*% 
                      solve(t(Z_c) %*% Z_c) %*% t(Z_c)
          weights <- weights + t(ridge_w)
        } else {
          no_cov_weights <- NULL
        }
    }

    l2_imbalance <- sqrt(sum((synth_data$X0 %*% weights - synth_data$X1)^2))

    ## primal objective value scaled by least squares difference for mean
    uni_w <- matrix(1/ncol(synth_data$X0), nrow=ncol(synth_data$X0), ncol=1)
    unif_l2_imbalance <- sqrt(sum((synth_data$X0 %*% uni_w - synth_data$X1)^2))
    scaled_l2_imabalance <- l2_imbalance / unif_l2_imbalance


    ## no outcome model
    mhat <- matrix(0, nrow=nrow(y), ncol=ncol(y))
    ridge_mhat <- mhat
    if(!is.null(Z)) {
      if(residualize) {
        ridge_mhat <- ridge_mhat + Z_cent %*% solve(t(Z_c) %*% Z_c) %*%
                        t(Z_c) %*% y_c

        ## regress out covariates for outcomes
        yc_hat <- ridge_mhat[trt == 0,, drop = F]
        # take residuals of outcomes
        y_c <- y_c - yc_hat
      } else {
        X_cent <- cbind(X_cent, Z_cent)
      }
    }

    if(ridge) {
        ridge_mhat <- ridge_mhat + X_cent %*% solve(t(X_c) %*% X_c +
                        lambda * diag(ncol(X_c))) %*%
                        t(X_c) %*% y_c
    }

    output <- list(weights = weights,
                l2_imbalance = l2_imbalance,
                scaled_l2_imbalance = scaled_l2_imabalance,
                mhat = mhat,
                lambda = lambda,
                ridge_mhat = ridge_mhat,
                synw = synw,
                lambdas = lambdas,
                lambda_errors = lambda_errors,
                lambda_errors_se = lambda_errors_se)

    if(!is.null(Z)) {
        output$no_cov_weights <- no_cov_weights

        z_l2_imbalance <- sqrt(sum((t(Z_c) %*% weights - t(Z_1))^2))
        z_unif_l2_imbalance <- sqrt(sum((t(Z_c) %*% uni_w - t(Z_1))^2))
        z_scaled_l2_imbalance <- z_l2_imbalance / z_unif_l2_imbalance

        output$covariate_l2_imbalance <- z_l2_imbalance
        output$scaled_covariate_l2_imbalance <- z_scaled_l2_imbalance

    }
    return(output)
}

#' Helper function to fit ridge ASCM
#' @param X_c Matrix of control lagged outcomes
#' @param X_1 Vector of treated leagged outcomes
#' @param trt Vector of treatment indicators
#' @param synth_data Output of `format_synth`
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param ridge Include ridge or not
#' @param scm Include SCM or not
#' @param lambda_min_ratio Ratio of the smallest to largest lambda when tuning lambda values
#' @param n_lambda Number of lambdas to consider between the smallest and largest lambda value
#' @param lambda_max Initial (largest) lambda, if NULL sets it to be (1+norm(X_1-X_c))^2
#' @param holdout_length Length of conseuctive holdout period for when tuning lambdas 
#' @param min_1se If TRUE, chooses the maximum lambda within 1 standard error of the lambda that minimizes the CV error, if FALSE chooses the optimal lambda; default TRUE
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"lambda"}{Value of the ridge hyperparameter}
#'          \item{"synw"}{The synth weights(for estimating the bias)}
#'          \item{"lambdas"}{List of lambda values evaluated to tune ridge regression}
#'          \item{"lambda_errors"}{"The MSE associated with each lambda term in lambdas."}
#'          \item{"lambda_errors_se"}{"The SE of the MSE associated with each lambda term in lambdas."}
#' }
fit_ridgeaug_inner <- function(X_c, X_1, trt, synth_data,
                               lambda, ridge, scm,
                               lambda_min_ratio, n_lambda,
                               lambda_max,
                               holdout_length, min_1se) {
    lambda_errors <- NULL
    lambda_errors_se <- NULL
    lambdas <- NULL

    ## if SCM fit scm
    if(scm) {
        syn <- fit_synth_formatted(synth_data)$weights
    } else {
        ## else use uniform weights
        syn <- rep(1 / sum(trt == 0), sum(trt == 0))
    }
    if(ridge) {
        if(is.null(lambda)) {
            cv_out <- cv_lambda(X_c, X_1, synth_data, trt, holdout_length, scm,
                      lambda_max, lambda_min_ratio, n_lambda, min_1se)

            lambda <- cv_out$lambda
            lambda_errors <- cv_out$lambda_errors
            lambda_errors_se <- cv_out$lambda_errors_se
            lambdas <- cv_out$lambdas
        }
        # get ridge weights
        ridge_w <- t(t(X_1) - t(X_c) %*% syn) %*%
                    solve(t(X_c) %*% X_c  + lambda * diag(ncol(X_c))) %*% t(X_c)
    } else {
        ridge_w <- matrix(0, ncol = sum(trt == 0), nrow=1)
    }
    ## combine weights
    weights <- syn + t(ridge_w)

    return(list(weights = weights,
                synw = syn,
                lambda = lambda,
                lambdas = lambdas,
                lambda_errors = lambda_errors,
                lambda_errors_se = lambda_errors_se))
}



#' Choose max lambda as largest eigenvalue of control X
#' @param X_c matrix of control lagged outcomes
#' @noRd
#' @return max lambda
get_lambda_max <- function(X_c) {
    svd(X_c)$d[1] ^ 2
}
#' Create list of lambdas
#' @param lambda_min_ratio Ratio of the smallest to largest lambda when tuning lambda values
#' @param n_lambda Number of lambdas to consider between the smallest and largest lambda value
#' @param lambda_max Initial (largest) lambda, if NULL sets it to be (1+norm(X_1-X_c))^2
#' @noRd
#' @return List of lambdas
create_lambda_list <- function(lambda_max, lambda_min_ratio, n_lambda) {
    scaler <- (lambda_min_ratio) ^ (1/n_lambda)
    lambdas <- lambda_max * (scaler ^ (seq(0:n_lambda) - 1))
    return(lambdas)
}

#' Choose either the lambda that minimizes CV MSE or largest lambda within 1 se of min
#' @param lambdas list of lambdas
#' @param lambda_errors The MSE associated with each lambda term in lambdas.
#' @param lambda_errors_se The SE of the MSE associated with each lambda
#' @param min_1se If TRUE, chooses the maximum lambda within 1 standard error of the lambda that minimizes the CV error, if FALSE chooses the optimal lambda; default TRUE
#' @noRd
#' @return optimal lambda
choose_lambda <- function(lambdas, lambda_errors, lambda_errors_se, min_1se) {
    # lambda with smallest error
    min_idx <- which.min(lambda_errors)
    min_error <- lambda_errors[min_idx]
    min_se <- lambda_errors_se[min_idx]
    lambda_min <- lambdas[min_idx]
    # max lambda with error within one se of min
    lambda_1se <- max(lambdas[lambda_errors <= min_error + min_se])
    return(if(min_1se) lambda_1se else lambda_min)
}

#' Choose best lambda with CV
#' @param X_c Matrix of control lagged outcomes
#' @param X_1 Vector of treated leagged outcomes
#' @param synth_data Output of `format_synth`
#' @param trt Vector of treatment indicators
#' @param holdout_length Length of conseuctive holdout period for when tuning lambdas 
#' @param scm Include SCM or not
#' @param lambda_max Initial (largest) lambda, if NULL sets it to be (1+norm(X_1-X_c))^2
#' @param lambda_min_ratio Ratio of the smallest to largest lambda when tuning lambda values
#' @param n_lambda Number of lambdas to consider between the smallest and largest lambda value
#' @param min_1se If TRUE, chooses the maximum lambda within 1 standard error of the lambda 
#' @noRd
#' @return \itemize{
#'          \item{"lambda"}{Value of the ridge hyperparameter}
#'          \item{"lambdas"}{List of lambda values evaluated to tune ridge regression}
#'          \item{"lambda_errors"}{"The MSE associated with each lambda term in lambdas."}
#'          \item{"lambda_errors_se"}{"The SE of the MSE associated with each lambda term}
#' }
cv_lambda <- function(X_c, X_1, synth_data, trt, holdout_length, scm,
                      lambda_max, lambda_min_ratio, n_lambda, min_1se) {
    if(is.null(lambda_max)) {
        lambda_max <- get_lambda_max(X_c) 
    }

    lambdas <- create_lambda_list(lambda_max, lambda_min_ratio, n_lambda)
    
    lambda_out <- get_lambda_errors(lambdas, X_c, X_1,
                                        synth_data, trt,
                                        holdout_length, scm)
    lambda_errors <- lambda_out$lambda_errors
    lambda_errors_se <- lambda_out$lambda_errors_se

    lambda <- choose_lambda(lambdas, lambda_errors, lambda_errors_se, min_1se)

    return(list(lambda = lambda, lambda_errors = lambda_errors,
                lambda_errors_se = lambda_errors_se, lambdas = lambdas))
}
