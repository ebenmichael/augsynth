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

