################################################################################
## Code for inference
################################################################################

#' Jackknife+ algorithm over time
#' @param ascm Fitted `augsynth` object
#' @param alpha Confidence level
#' @param conservative Whether to use the conservative jackknife+ procedure
#' @return List that contains:
#'         \itemize{
#'          \item{"att"}{Vector of ATT estimates}
#'          \item{"heldout_att"}{Vector of ATT estimates with the time period held out}
#'          \item{"se"}{Standard error, always NA but returned for compatibility}
#'          \item{"lb"}{Lower bound of 1 - alpha confidence interval}
#'          \item{"ub"}{Upper bound of 1 - alpha confidence interval}
#'          \item{"alpha"}{Level of confidence interval}
#'         }
time_jackknife_plus <- function(ascm, alpha = 0.05, conservative = F) {
    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    n_c <- dim(synth_data$Z0)[2]
    Z <- wide_data$Z

    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]

    jack_ests <- lapply(1:t0, 
        function(tdrop) {
            # drop unit i
            new_data <- drop_time_t(wide_data, Z, tdrop)
            # refit
            new_ascm <- do.call(fit_augsynth_internal,
                    c(list(wide = new_data$wide,
                            synth_data = new_data$synth_data,
                            Z = new_data$Z,
                            progfunc = ascm$progfunc,
                            scm = ascm$scm,
                            fixedeff = ascm$fixedeff),
                        ascm$extra_args))
            # get ATT estimates and held out error for time t
            # t0 is prediction for held out time
            est <- predict(new_ascm, att = F)[(t0 +1):t_final]
            est <- c(est, mean(est))
            err <- c(colMeans(wide_data$X[wide_data$trt == 1,
                                         tdrop,
                                         drop = F]) -
                    predict(new_ascm, att = F)[t0])
            list(err, rbind(est + abs(err), est - abs(err), est + err, est))
        })
    # get errors and jackknife distribution
    held_out_errs <- vapply(jack_ests, `[[`, numeric(1), 1)
    jack_dist <- vapply(jack_ests, `[[`,
                        matrix(0, nrow = 4, ncol = tpost + 1), 2)

    out <- list()
    att <- predict(ascm, att = T)
    out$att <- c(att, 
                 mean(att[(t0 + 1):t_final]))
    # held out ATT
    out$heldout_att <- c(held_out_errs, 
                          att[(t0 + 1):t_final], 
                          mean(att[(t0 + 1):t_final]))

    # out$se <- rep(NA, 10 + tpost)
    if(conservative) {
        qerr <- stats::quantile(abs(held_out_errs), 1 - alpha)
        out$lb <- c(rep(NA, t0), apply(jack_dist[4,,], 1, min) - qerr)
        out$ub <- c(rep(NA, t0), apply(jack_dist[4,,], 1, max) + qerr)
    } else {
        out$lb <- c(rep(NA, t0), apply(jack_dist[2,,], 1, stats::quantile, alpha / 2))
        out$ub <- c(rep(NA, t0), apply(jack_dist[1,,], 1, stats::quantile, 1 - alpha / 2))
    }
    # shift back to ATT scale
    y1 <- predict(ascm, att = F) + att
    y1 <-  c(y1, mean(y1[(t0 + 1):t_final]))
    shifted_lb <- y1 - out$ub
    shifted_ub <- y1 - out$lb
    out$lb <- shifted_lb
    out$ub <- shifted_ub
    out$alpha <- alpha


    return(out)
}

#' Drop time period from pre-treatment data
#' @param wide_data (X, y, trt)
#' @param Z Covariates matrix
#' @param t_drop Time to drop
#' @noRd
drop_time_t <- function(wide_data, Z, t_drop) {

        new_wide_data <- list()
        new_wide_data$trt <- wide_data$trt
        new_wide_data$X <- wide_data$X[, -t_drop, drop = F]
        new_wide_data$y <- cbind(wide_data$X[, t_drop, drop = F], 
                                 wide_data$y)

        X0 <- new_wide_data$X[new_wide_data$trt == 0,, drop = F]
        x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1,,
                                              drop = F]),
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

#' Conformal inference procedure to compute p-values and point-wise confidence intervals
#' @param ascm Fitted `augsynth` object
#' @param alpha Confidence level
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' @param grid_size Number of grid points to use when inverting the hypothesis test
#' @return List that contains:
#'         \itemize{
#'          \item{"att"}{Vector of ATT estimates}
#'          \item{"heldout_att"}{Vector of ATT estimates with the time period held out}
#'          \item{"se"}{Standard error, always NA but returned for compatibility}
#'          \item{"lb"}{Lower bound of 1 - alpha confidence interval}
#'          \item{"ub"}{Upper bound of 1 - alpha confidence interval}
#'          \item{"p_val"}{p-value for test of no post-treatment effect}
#'          \item{"alpha"}{Level of confidence interval}
#'         }
conformal_inf <- function(ascm, alpha = 0.05, type = "iid",
                          q = 1, ns = 1000, grid_size = 50) {
  wide_data <- ascm$data
  synth_data <- ascm$data$synth_data
  n <- nrow(wide_data$X)
  n_c <- dim(synth_data$Z0)[2]
  Z <- wide_data$Z

  t0 <- dim(synth_data$Z0)[1]
  tpost <- ncol(wide_data$y)
  t_final <- dim(synth_data$Y0plot)[1]

  # grid of nulls
  att <- predict(ascm, att = T)
  post_att <- att[(t0 +1):t_final]
  post_sd <- sqrt(mean(post_att ^ 2))
  # iterate over post-treatment periods to get pointwise CIs
  vapply(1:tpost,
         function(j) {
          # fit using t0 + j as a pre-treatment period and get reisduals
          new_wide_data <- wide_data
          new_wide_data$X <- cbind(wide_data$X, wide_data$y[, j, drop = TRUE])
          if(tpost > 1) {
            new_wide_data$y <- wide_data$y[, -j, drop = FALSE]
          } else {
            # set the post period has to be *something*
            new_wide_data$y <- matrix(1, nrow = n, ncol = 1)
          }


          # make a grid around the estimated ATT
          grid <- seq(att[t0 + j] - 2 * post_sd, att[t0 + j] + 2 * post_sd,
                      length.out = grid_size)
          compute_permute_ci(new_wide_data, ascm, grid, 1, alpha, "block", q, ns)
         },
         numeric(3)) -> cis

  # test a null post-treatment effect
  new_wide_data <- wide_data
  new_wide_data$X <- cbind(wide_data$X, wide_data$y)
  new_wide_data$y <- matrix(1, nrow = n, ncol = 1)
  null_p <- compute_permute_pval(new_wide_data, ascm, 0, ncol(wide_data$y), 
                                 type, q, ns)
  
  out <- list()
  att <- predict(ascm, att = T)
  out$att <- c(att, mean(att[(t0 + 1):t_final]))
  # out$se <- rep(NA, t_final)
  # out$sigma <- NA
  out$lb <- c(rep(NA, t0), cis[1, ], NA)
  out$ub <- c(rep(NA, t0), cis[2, ], NA)
  out$p_val <- c(rep(NA, t0), cis[3, ], null_p)
  out$alpha <- alpha
  return(out)
}

#' Compute conformal test statistics
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param h0 Null hypothesis to test
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' 
#' @return List that contains:
#'         \itemize{
#'          \item{"resids"}{Residuals after enforcing the null}
#'          \item{"test_stats"}{Permutation distribution of test statistics}
#'          \item{"stat_func"}{Test statistic function}
#'         }
#' @noRd
compute_permute_test_stats <- function(wide_data, ascm, h0,
                                       post_length, type,
                                       q, ns) {
  # format data
  new_wide_data <- wide_data
  t0 <- ncol(wide_data$X) - post_length
  tpost <- t0 + post_length
  # adjust outcomes for null
  new_wide_data$X[wide_data$trt == 1,(t0 + 1):tpost ] <- new_wide_data$X[wide_data$trt == 1,(t0 + 1):tpost] - h0
  X0 <- new_wide_data$X[new_wide_data$trt == 0,, drop = F]
  x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1,, drop = F]),
              ncol=1)

  new_synth_data <- list()
  new_synth_data$Z0 <- t(X0)
  new_synth_data$X0 <- t(X0)
  new_synth_data$Z1 <- x1
  new_synth_data$X1 <- x1

  # fit synth with adjusted data and get residuals
  new_ascm <- do.call(fit_augsynth_internal,
                    c(list(wide = new_wide_data,
                            synth_data = new_synth_data,
                            Z = wide_data$Z,
                            progfunc = ascm$progfunc,
                            scm = ascm$scm,
                            fixedeff = ascm$fixedeff),
                        ascm$extra_args))
  resids <- predict(new_ascm, att = T)[1:tpost]
  # permute residuals and compute test statistic
  if(post_length == 1) {
    stat_func <- function(x) (sum(abs(x) ^ q)  / sqrt(length(x))) ^ (1 / q)
  } else if(post_length > 0) {
    stat_func <- function(x) (abs(sum(x)) ^ q  / sqrt(length(x))) ^ (1 / q)
  }
  
  if(type == "iid") {
    test_stats <- sapply(1:ns, 
                        function(x) {
                          reorder <- sample(resids)
                          stat_func(reorder[(t0 + 1):tpost])
                        })
  } else {
    ## increment time by one step and wrap
    test_stats <- sapply(1:tpost,
                        function(j) {
                          reorder <- resids[(0:tpost -1 + j) %% tpost + 1]
                          stat_func(reorder[(t0 + 1):tpost])
                        })
  }
  
  return(list(resids = resids,
              test_stats = test_stats,
              stat_func = stat_func))
}


#' Compute conformal p-value
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param h0 Null hypothesis to test
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' 
#' @return Computed p-value
#' @noRd
compute_permute_pval <- function(wide_data, ascm, h0,
                                 post_length, type,
                                 q, ns) {
  t0 <- ncol(wide_data$X) - post_length
  tpost <- t0 + post_length
  out <- compute_permute_test_stats(wide_data, ascm, h0,
                                    post_length, type, q, ns)
  mean(out$stat_func(out$resids[(t0 + 1):tpost]) <= out$test_stats)
}

#' Compute conformal p-value
#' @param wide_data List containing pre- and post-treatment outcomes and outcome vector
#' @param ascm Fitted `augsynth` object
#' @param grid Set of null hypothesis to test for inversion
#' @param post_length Number of post-treatment periods
#' @param type Either "iid" for iid permutations or "block" for moving block permutations
#' @param q The norm for the test static `((sum(x ^ q))) ^ (1/q)`
#' @param ns Number of resamples for "iid" permutations
#' 
#' @return (lower bound of interval, upper bound of interval, p-value for null of 0 effect)
#' @noRd
compute_permute_ci <- function(wide_data, ascm, grid,
                               post_length, alpha, type,
                               q, ns) {
  # make sure 0 is in the grid
  grid <- c(grid, 0)
  ps <-sapply(grid, 
              function(x) {
                compute_permute_pval(wide_data, ascm, x, 
                                     post_length, type, q, ns)
              })
  c(min(grid[ps >= alpha]), max(grid[ps >= alpha]), ps[grid == 0])
}



#' Drop unit i from data
#' @param wide_data (X, y, trt)
#' @param Z Covariates matrix
#' @param i Unit to drop
#' @noRd
drop_unit_i <- function(wide_data, Z, i) {

        new_wide_data <- list()
        new_wide_data$trt <- wide_data$trt[-i]
        new_wide_data$X <- wide_data$X[-i,, drop = F]
        new_wide_data$y <- wide_data$y[-i,, drop = F]

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
        new_Z <- if(!is.null(Z)) Z[-i, , drop = F] else NULL

        return(list(wide_data = new_wide_data,
                    synth_data = new_synth_data,
                    Z = new_Z))
}

#' Drop unit i from data
#' @param wide_list (X, y, trt)
#' @param Z Covariates matrix
#' @param i Unit to drop
#' @noRd
drop_unit_i_multiout <- function(wide_list, Z, i) {

        new_wide_data <- list()
        new_wide_data$trt <- wide_list$trt[-i]
        new_wide_data$X <- lapply(wide_list$X, function(x) x[-i,, drop = F])
        new_wide_data$y <- lapply(wide_list$y, function(x) x[-i,, drop = F])
        new_Z <- if(!is.null(Z)) Z[-i, , drop = F] else NULL

        return(list(wide_list = new_wide_data,
                    Z = new_Z))
}


#' Estimate standard errors for single ASCM with the jackknife
#' Do this for ridge-augmented synth
#' @param ascm Fitted augsynth object
#' 
#' @return List that contains:
#'         \itemize{
#'          \item{"att"}{Vector of ATT estimates}
#'          \item{"se"}{Standard error estimate}
#'          \item{"lb"}{Lower bound of 1 - alpha confidence interval}
#'          \item{"ub"}{Upper bound of 1 - alpha confidence interval}
#'          \item{"alpha"}{Level of confidence interval}
#'         }
jackknife_se_single <- function(ascm) {

    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    n_c <- dim(synth_data$Z0)[2]
    Z <- wide_data$Z

    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)


    # only drop out control units with non-zero weights
    nnz_weights <- numeric(n)
    nnz_weights[wide_data$trt == 0] <- round(ascm$weights, 3) != 0
    # if more than one unit is treated, include them in the jackknife
    if(sum(wide_data$trt) > 1) {
      nnz_weights[wide_data$trt == 1] <- 1
    }

    trt_idxs <- (1:n)[as.logical(nnz_weights)]


    # jackknife estimates
    ests <- vapply(trt_idxs,
                   function(i) {
                       # drop unit i
                       new_data <- drop_unit_i(wide_data, Z, i)
                       # refit
                       new_ascm <- do.call(fit_augsynth_internal,
                                c(list(wide = new_data$wide,
                                       synth_data = new_data$synth_data,
                                       Z = new_data$Z,
                                       progfunc = ascm$progfunc,
                                       scm = ascm$scm,
                                       fixedeff = ascm$fixedeff),
                                  ascm$extra_args))
                       # get ATT estimates
                       est <- predict(new_ascm, att = T)[(t0 + 1):t_final]
                       c(est, mean(est))
                   },
                   numeric(tpost + 1))
    # convert to matrix
    ests <- matrix(ests, nrow = tpost + 1, ncol = length(trt_idxs))
    ## standard errors
    se2 <- apply(ests, 1,
                 function(x) (n - 1) / n * sum((x - mean(x, na.rm = T)) ^ 2))
    se <- sqrt(se2)

    out <- list()
    att <- predict(ascm, att = T)
    out$att <- c(att, mean(att[(t0 + 1):t_final]))

    out$se <- c(rep(NA, t0), se)
    # out$sigma <- NA
    return(out)
}


#' Compute standard errors using the jackknife
#' @param multisynth fitted multisynth object
#' @param relative Whether to compute effects according to relative time
#' @noRd
jackknife_se_multi <- function(multisynth, relative=NULL, alpha = 0.05) {
    ## get info from the multisynth object
    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    n_leads <- multisynth$n_leads
    n <- nrow(multisynth$data$X)
    att <- predict(multisynth, att=T)
    outddim <- nrow(att)

    J <- length(multisynth$grps)

    ## drop each unit and estimate overall treatment effect
    jack_est <- vapply(1:n,
                       function(i) {
                           msyn_i <- drop_unit_i_multi(multisynth, i)
                           pred <- predict(msyn_i[[1]], relative=relative, att=T)
                           if(length(msyn_i[[2]]) != 0) {
                               out <- matrix(NA, nrow=nrow(pred), ncol=(J+1))
                               out[,-(msyn_i[[2]]+1)] <- pred
                           } else {
                               out <- pred
                           }
                           out
                       },
                       matrix(0, nrow=outddim,ncol=(J+1)))

    se2 <- apply(jack_est, c(1,2),
                function(x) (n-1) / n * sum((x - mean(x,na.rm=T))^2, na.rm=T))
    lower_bound <- att - qnorm(1 - alpha / 2) * sqrt(se2)
    upper_bound <- att + qnorm(1 - alpha / 2) * sqrt(se2)
    return(list(att = att, se = sqrt(se2),
                lower_bound = lower_bound, upper_bound = upper_bound))

}

#' Helper function to drop unit i and refit
#' @param msyn multisynth_object
#' @param i Unit to drop
#' @noRd
drop_unit_i_multi <- function(msyn, i) {

    n <- nrow(msyn$data$X)
    time_cohort <- msyn$time_cohort
    which_t <- (1:n)[is.finite(msyn$data$trt)]

    not_miss_j <- which_t %in% setdiff(which_t, i)

    # drop unit i from data
    drop_i <- msyn$data
    drop_i$X <- msyn$data$X[-i, , drop = F]
    drop_i$y <- msyn$data$y[-i, , drop = F]
    drop_i$trt <- msyn$data$trt[-i]
    drop_i$mask <- msyn$data$mask[not_miss_j,, drop = F]

    if(!is.null(msyn$data$Z)) {
      drop_i$Z <- msyn$data$Z[-i, , drop = F]
    } else {
      drop_i$Z <- NULL
    }

    long_df <- msyn$long_df
    unit <- colnames(long_df)[1]
    # make alphabetical, because the ith unit is the index in alphabetical ordering
    long_df <- long_df[order(long_df[, unit, drop = TRUE]),]
    ith_unit <- unique(long_df[,unit, drop = TRUE])[i]
    long_df <- long_df[long_df[,unit, drop = TRUE] != ith_unit,]

    # re-fit everything
    args_list <- list(wide = drop_i, relative = msyn$relative,
                      n_leads = msyn$n_leads, n_lags = msyn$n_lags,
                      nu = msyn$nu, lambda = msyn$lambda,
                      V = msyn$V,
                      force = msyn$force, n_factors = msyn$n_factors,
                      scm = msyn$scm, time_w = msyn$time_w,
                      lambda_t = msyn$lambda_t,
                      fit_resids = msyn$fit_resids,
                      time_cohort = msyn$time_cohort, long_df = long_df,
                      how_match = msyn$how_match)
    msyn_i <- do.call(multisynth_formatted, c(args_list, msyn$extra_pars))

    # check for dropped treated units/time periods
    if(time_cohort) {
        dropped <- which(!msyn$grps %in% msyn_i$grps)
    } else {
        dropped <- which(!not_miss_j)
    }
    return(list(msyn_i,
                dropped))
}


#' Estimate standard errors for multi outcome ascm with jackknife
#' @param ascm Fitted augsynth object
#' @noRd
jackknife_se_multiout <- function(ascm) {

    wide_data <- ascm$data
    wide_list <- ascm$data_list
    n <- nrow(wide_data$X)
    Z <- wide_data$Z


    # only drop out control units with non-zero weights
    nnz_weights <- numeric(n)
    nnz_weights[wide_data$trt == 0] <- round(ascm$weights, 3) != 0

    trt_idxs <- (1:n)[as.logical(nnz_weights)]

    # jackknife estimates
    ests <- lapply(trt_idxs,
                   function(i) {
                       # drop unit i
                       new_data <- drop_unit_i_multiout(wide_list, Z, i)
                       # refit
                       new_ascm <- do.call(fit_augsynth_multiout_internal,
                                c(list(wide = new_data$wide,
                                       combine_method = ascm$combine_method,
                                       Z = new_data$Z,
                                       progfunc = ascm$progfunc,
                                       scm = ascm$scm,
                                       fixedeff = ascm$fixedeff),
                                  ascm$extra_args))
                        new_ascm$outcomes <- ascm$outcomes
                        new_ascm$data_list <- ascm$data_list
                        new_ascm$data$time <- ascm$data$time
                       # get ATT estimates
                       est <- predict(new_ascm, att = T)
                       est <- est[as.numeric(rownames(est)) >= ascm$t_int,, drop = F]
                       rbind(est, colMeans(est, na.rm = T))
                   })
    ests <- simplify2array(ests)
    ## standard errors
    se2 <- apply(ests, c(1, 2),
                 function(x) (n - 1) / n * sum((x - mean(x, na.rm = T)) ^ 2))
    se <- sqrt(se2)
    out <- list()
    att <- predict(ascm, att = T)
    att_post <- colMeans(att[as.numeric(rownames(att)) >= ascm$t_int,, drop = F],
                         na.rm = T)
    out$att <- rbind(att, att_post)
    t0 <- sum(as.numeric(rownames(att)) < ascm$t_int)
    out$se <- rbind(matrix(NA, t0, ncol(se)), se)
    out$sigma <- NA
    return(out)
}



#' Compute the weighted bootstrap distribution
#' @param multisynth fitted multisynth object
#' @param rweight Function to draw random weights as a function of n (e.g rweight(n))
#' @param relative Whether to compute effects according to relative time
#' @noRd
weighted_bootstrap_multi <- function(multisynth,
                                    rweight = rwild_b,
                                    n_boot = 1000,
                                    alpha = 0.05,
                                    relative=NULL) {
  ## get info from the multisynth object
  if(is.null(relative)) {
      relative <- multisynth$relative
  }

  n <- nrow(multisynth$data$X)
  att <- predict(multisynth, att=T)
  outddim <- nrow(att)
  n1 <- sum(is.finite(multisynth$data$trt))
  J <- length(multisynth$grps)


  # draw random weights to get bootstrap distribution
  bs_est <- vapply(1:n_boot,
                      function(i) {
                        Z <- rweight(n)# / sqrt(n1)

                        predict(multisynth, att=T, bs_weight = Z) - sum(Z) / n1 * att
                      },
                      matrix(0, nrow=outddim,ncol=(J+1)))

  se2 <- apply(bs_est, c(1,2),
              function(x) mean((x - mean(x))^2, na.rm=T))
  bias <- apply(bs_est, c(1,2),
              function(x) mean(x, na.rm=T))
  upper_bound <- att - apply(bs_est, c(1,2),
              function(x) quantile(x, alpha / 2, na.rm = T))
  
  lower_bound <- att - apply(bs_est, c(1,2),
              function(x) quantile(x, 1 - alpha / 2, na.rm = T))

  return(list(att = att,
              bias = bias,
              se = sqrt(se2),
              upper_bound = upper_bound,
              lower_bound = lower_bound))

}

#' Bayesian bootstrap
#' @param n Number of units
#' @export
rdirichlet_b <- function(n) {
  Z <- as.numeric(rgamma(n, 1, 1))
  return(Z / sum(Z) * n)
}

#' Non-parametric bootstrap
#' @param n Number of units
#' @export
rmultinom_b <- function(n) as.numeric(rmultinom(1, n, rep(1 / n, n)))

#' Wild bootstrap (Mammen 1993)
#' @param n Number of units
#' @export
rwild_b <- function(n) {
  sample(c(-(sqrt(5) - 1) / 2, (sqrt(5) + 1) / 2 ), n,
         replace = TRUE,
         prob = c((sqrt(5) + 1)/ (2 * sqrt(5)), (sqrt(5) - 1) / (2 * sqrt(5))))
}