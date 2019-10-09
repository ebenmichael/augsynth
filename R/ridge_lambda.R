################################################################################
## Function to calculate error on different lambda values if using Ridge Augmented SCM
################################################################################

#' Get Lambda Errors
#'
#' @param lambdas Vector of lambda values to compute errors for
#' @param X_c Matrix of control group pre-treatment outcomes
#' @param X_t Matrix of treatment group pre-treatment outcomes
#' @param synth_data Output of `format_synth`
#' @param trt Boolean vector of treatment assignments
#' @param holdout_length Length of conseuctive holdout period for when tuning lambdas
#' @param scm Include SCM or not
#' 
#' @return List of lambda errors for each corresponding lambda in the lambdas parameter.
get_lambda_errors <- function(lambdas, X_c, X_t, synth_data, trt, holdout_length=1, scm=T) {
  lambda_errors = numeric(length(lambdas)) # vector that stores the sum MSE across all CV sets
  for (i in 1:(ncol(X_c) - holdout_length)) {
    X_0 <- X_c[,-(i:i+holdout_length-1)]
    X_1 <- matrix(X_t[-(i:i+holdout_length-1)])
    X_0v <- X_c[,i:i+holdout_length-1]
    X_1v <- matrix(X_t[i:i+holdout_length-1])
    new_synth_data <- synth_data
    new_synth_data$Z1 <- X_1
    new_synth_data$X1 <- X_1
    new_synth_data$Z0 <- t(X_0)
    new_synth_data$X0 <- t(X_0)
    if(scm) {
      syn <- fit_synth_formatted(new_synth_data)$weights
    } else {
      syn <- rep(1/sum(trt==0), sum(trt==0))
    }

    for (j in 1:length(lambdas)) {
      ridge_weights <- t(X_1 - t(X_0) %*% syn) %*% solve(t(X_0) %*% X_0 + lambdas[j] * diag(ncol(X_0))) %*% t(X_0)
      aug_weights <- syn + t(ridge_weights)
      error <- X_1v - t(X_0v) %*% aug_weights
      error <- sum(error ^ 2) # take sum of errors across the holdout time periods.
      lambda_errors[j] <- lambda_errors[j] + error
    }
  }
  return(lambda_errors)
}

