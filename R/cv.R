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
    # vector that stores the sum MSE across all CV sets
    errors <- matrix(0, nrow = ncol(X) - holdout_length, ncol = length(lambdas))
    lambda_errors = numeric(length(lambdas)) 
    lambda_errors_se = numeric(length(lambdas)) 
    
    # TODO: turn for loop into apply over lambda sequence, (vapply)
    for (i in 1:(ncol(X) - holdout_length)) {
      for (j in 1:length(lambdas)) {
        t_drop <- i:(i+holdout_length - 1)
        new_ascm <- drop_time_and_refit(wide_data, Z, t_drop, progfunc, scm, fixedeff, lambda = j, ...)
        err <- sum((predict(new_ascm, att = T)[(ncol(X)-holdout_length+1):ncol(X)])^2)
        errors[i, j] <-  err
      }
    }
    lambda_errors <- apply(errors, 2, mean)
    lambda_errors_se <- apply(errors, 2, function(x) sd(x) / sqrt(length(x)))
    return(list(lambda_errors = lambda_errors, 
                lambda_errors_se = lambda_errors_se))
  }