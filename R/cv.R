drop_time_t <- function(wide_data, Z, t_drop) {
  new_wide_data <- list()
  new_wide_data$trt <- wide_data$trt
  
  if (is.list(wide_data$X)) {
    # TODO
  } else {
    new_wide_data$X <- wide_data$X[, -t_drop, drop = F]
    new_wide_data$y <- cbind(wide_data$X[, t_drop, drop = F], 
                             wide_data$y)
    
    X0 <- new_wide_data$X[new_wide_data$trt == 0,, drop = F]
    x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1,, drop = F]),
                 ncol=1)
    y0 <- new_wide_data$y[new_wide_data$trt == 0,, drop = F]
    y1 <- colMeans(new_wide_data$y[new_wide_data$trt == 1,, drop = F])
    
    new_synth_data <- list()
    new_synth_data$Z0 <- t(X0)
    new_synth_data$X0 <- t(X0)
    new_synth_data$Z1 <- x1
    new_synth_data$X1 <- x1
    
    return(list(wide_data = new_wide_data,
                synth_data = new_synth_data,
                Z = Z)) 
  }
}

drop_time_and_refit <- function(wide_data, Z, t_drop, progfunc, scm, fixedeff, ...) {
  new_data <- drop_time_t(wide_data, Z, t_drop)
  
  new_ascm <- do.call(fit_augsynth_internal,
                      c(list(wide = new_data$wide,
                             synth_data = new_data$synth_data,
                             Z = new_data$Z,
                             progfunc = progfunc,
                             scm = scm,
                             fixedeff = fixedeff, ...)))
  return(new_ascm)
}

cv <- function(lambdas, wide_data, Z, progfunc, scm, fixedeff, holdout_length, ...) {
    # TODO: instead of holdout_lenght, assume input parameter is a list of folds (each fold represented by the periods to drop)
    X <- wide_data$X
    
    lambda_error_vals <- vapply(lambdas, function(lambda){
      periods <- 1:(ncol(X) - holdout_length)
      print(periods)
      errors <- vapply(periods, function(t){
        t_drop <- t:(t+holdout_length - 1)
        print(t_drop)
        new_ascm <- drop_time_and_refit(wide_data, Z, t_drop, progfunc, scm, fixedeff, lambda = lambda, ...)
        err <- sum((predict(new_ascm, att = T)[(ncol(X)-holdout_length+1):ncol(X)])^2)
        err
      }, numeric(1))
      lambda_error <- mean(errors)
      lambda_error_se <- sd(errors) / sqrt(length(errors))
      c(lambda_error, lambda_error_se)
    }, numeric(2))
    return(list(lambda_errors = lambda_error_vals[1,], lambda_errors_se = lambda_error_vals[2,]))
}

cv2 <- function(lambdas, wide_data, Z, progfunc, scm, fixedeff, holdout_periods, ...) {
  X <- wide_data$X
  lambda_error_vals <- vapply(lambdas, function(lambda){
    errors <- apply(holdout_periods, 1, function(t_drop){
      print(t_drop)
      new_ascm <- drop_time_and_refit(wide_data, Z, t_drop, progfunc, scm, fixedeff, lambda = lambda, ...)
      err <- sum((predict(new_ascm, att = T)[(ncol(X)-length(t_drop)+1):ncol(X)])^2)
      err
    })
    lambda_error <- mean(errors)
    lambda_error_se <- sd(errors) / sqrt(length(errors))
    c(lambda_error, lambda_error_se)
  }, numeric(2))
  return(list(lambda_errors = lambda_error_vals[1,], lambda_errors_se = lambda_error_vals[2,]))
}