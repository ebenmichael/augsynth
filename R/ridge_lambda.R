################################################################################
## Function to calculate error on different lambda values if using Ridge Augmented SCM
################################################################################

#' Get Lambda Errors
#' @importFrom stats sd
#'
#' @param lambdas Vector of lambda values to compute errors for
#' @param X_c Matrix of control group pre-treatment outcomes
#' @param X_t Matrix of treatment group pre-treatment outcomes
#' @param synth_data Output of `format_synth`
#' @param trt Boolean vector of treatment assignments
#' @param holdout_length Length of conseuctive holdout period for when tuning lambdas
#' @param scm Include SCM or not
#' @noRd
#' @return List of lambda errors for each corresponding lambda in the lambdas parameter.
get_lambda_errors <- function(lambdas, X_c, X_t, synth_data, trt, holdout_length=1, scm=T) {
  # vector that stores the sum MSE across all CV sets
  errors <- matrix(0, nrow = ncol(X_c) - holdout_length, ncol = length(lambdas))
  lambda_errors = numeric(length(lambdas)) 
  lambda_errors_se = numeric(length(lambdas)) 

  for (i in 1:(ncol(X_c) - holdout_length)) {
    X_0 <- X_c[,-(i:(i + holdout_length - 1))]
    X_1 <- matrix(X_t[-(i:(i + holdout_length - 1))])
    X_0v <- X_c[,i:(i + holdout_length - 1)]
    X_1v <- matrix(X_t[i:(i + holdout_length - 1)], ncol = 1)
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
      # take sum of errors across the holdout time periods
      error <- sum(error ^ 2)
      errors[i, j] <-  error
      # lambda_errors[j] <- lambda_errors[j] + error
    }
  }
  lambda_errors <- apply(errors, 2, mean)
  lambda_errors_se <- apply(errors, 2, function(x) sd(x) / sqrt(length(x)))
  return(list(lambda_errors = lambda_errors, 
              lambda_errors_se = lambda_errors_se))
}