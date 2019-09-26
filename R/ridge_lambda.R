

get_lambda_errors <- function(lambdas, X_c, X_t, synth_t, trt, holdout_length=1, scm=T) {
  lambda_errors = c()
  for (lambda in lambdas) {
    errors = c()
    for (i in 1:(ncol(X_0) - holdout_length)) {
      X_0 <- X_c[,-(i:i+holdout_length-1)]
      X_1 <- X_t[-(i:i+holdout_length-1)]
      synth_t[0] <- synth_t[0][,-(i:i+holdout_length-1)] # TODO: validate synth_data format
      if(scm) {
        syn <- fit_synth_formatted(synth_t)$weights
      } else {
        ## else use uniform weights
        syn <- rep(1/sum(trt==0), sum(trt==0))
      }
    
      aug_weights = syn + t(t(X_1) - t(X_0) %*% syn) %*% solve(t(X_0) %*% X_0 + lambda * diag(ncol(X_0))) %*% t(X_0)
      X_0v <- X_c[,i:i+holdout_length-1]
      X_1v <- X_t[i:i+holdout_length-1]
      error <- X_1v - t(X_0v) %*% aug_weights
      error <- mean(error) # when holdout_length > 1, take the mean of errors across time periods
      errors <- c(errors, error)
    }
    mse <- mean(errors ^ 2)
    lambda_errors <- c(lambda_errors, mse) 
  }
  return lambda_errors
}
# TODO: validate trt?