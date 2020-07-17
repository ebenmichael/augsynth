################################################################################
## Code for inference
################################################################################

#' Estimate standard errors for single ASCM with the jackknife
#' Do this for ridge-augmented synth
#' @param ascm Fitted augsynth object
placebo_se_single <- function(ascm, homoskedastic = FALSE, 
                              include_treated = FALSE, unweighted = FALSE, ...) {

    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    n_c <- dim(synth_data$Z0)[2]
    Z <- wide_data$Z

    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)

    # LOO placebo estimates
    errs <- vapply((1:n)[wide_data$trt != 1],
                   function(i) {
                       # drop unit i
                       new_data <- placebo_unit_i(wide_data, Z, i)
                       # refit
                       new_ascm <- do.call(fit_augsynth_internal,
                                c(list(wide = new_data$wide,
                                       synth_data = new_data$synth_data,
                                       Z = new_data$Z,
                                       progfunc = ascm$progfunc,
                                       scm = ascm$scm,
                                       fixedeff = ascm$fixedeff),
                                  ascm$extra_args))
                       # get ATT estimates, truth is zero
                       est <- predict(new_ascm, att = T)[(t0 + 1):t_final]
                       c(est, mean(est))
                    #    est <- predict(new_ascm, att = F)[(t0 + 1):t_final]
                    #    est <- c(est, mean(est))
                    #    y0 <- c(wide_data$y[i, ], mean(wide_data$y[i, ]))
                    #    est - y0
                   },
                   numeric(tpost + 1))

    # convert to matrix
    errs <- matrix(errs, nrow = tpost + 1, ncol = (n - 1))
    ## standard errors
    if(homoskedastic) {
        se2 <- apply(errs, 1,
                function(x) {
                    trt_factor <- if(include_treated) 1 / sum(wide_data$trt == 1) else 0
                    sum(x ^ 2) * (trt_factor + sum(ascm$weights ^ 2))
                })
    } else if(!unweighted) {
        se2 <- apply(errs, 1, 
                    function(x) {
                        sum(ascm$weights ^ 2 * x ^ 2) / sum(ascm$weights) ^ 2
                    })
    } else {
        se2 <- apply(errs, 1, 
                    function(x) {
                        mean(x ^ 2) / length(x) 
                    })
    }
    # sig2 <- apply(errs ^ 2, 2, mean) ## estimate of variance

    se <- sqrt(se2)


    out <- list()
    att <- predict(ascm, att = T)
    out$att <- c(att, mean(att[(t0 + 1):t_final]))

    out$se <- c(rep(NA, t0), se)
    out$sigma <- NA
    return(out)
}


#' Use unit i as the treated unit
#' @param wide_data (X, y, trt)
#' @param Z Covariates matrix
#' @param i Unit to consider treated
placebo_unit_i <- function(wide_data, Z, i) {

        new_wide_data <- list()
        trt <- numeric(length(wide_data$trt))
        trt[i] <- 1
        new_wide_data$trt <- trt[wide_data$trt != 1]
        new_wide_data$X <- wide_data$X[wide_data$trt != 1, , drop = F]
        new_wide_data$y <- wide_data$y[wide_data$trt != 1, , drop = F]

        X0 <- new_wide_data$X[new_wide_data$trt == 0,, drop = F]
        x1 <- matrix(colMeans(new_wide_data$X[new_wide_data$trt == 1,, drop = F]),
                     ncol = 1)

        new_synth_data <- list()
        new_synth_data$Z0 <- t(X0)
        new_synth_data$X0 <- t(X0)
        new_synth_data$Z1 <- x1
        new_synth_data$X1 <- x1
        new_Z <- if(!is.null(Z)) Z[wide_data$trt != 1, , drop = F] else NULL

        return(list(wide_data = new_wide_data,
                    synth_data = new_synth_data,
                    Z = new_Z))
}


#' Estimate standard errors for single ASCM with residual bootstrap
#' Do this for ridge-augmented synth
#' @param ascm Fitted augsynth object
residual_bs_se_single <- function(ascm, b=1000, ...) {

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

    # trt_idxs <- (1:n)[as.logical(nnz_weights)]
    # n_jack <- length(trt_idxs)

    # LOO fitted values and residuals
    ests <- vapply((1:n)[wide_data$trt != 1],
                   function(i) {
                       # drop unit i
                       new_data <- placebo_unit_i(wide_data, Z, i)
                       # refit
                       new_ascm <- do.call(fit_augsynth_internal,
                                c(list(wide = new_data$wide,
                                       synth_data = new_data$synth_data,
                                       Z = new_data$Z,
                                       progfunc = ascm$progfunc,
                                       scm = ascm$scm,
                                       fixedeff = ascm$fixedeff),
                                  ascm$extra_args))
                       # get ATT estimates, truth is zero
                       est <- predict(new_ascm, att = F)[(t0 + 1):t_final]
                       est
                   },
                   numeric(tpost))
    # get residuals
    errs <- t(wide_data$y[wide_data$trt == 0, , drop = FALSE]) - ests


    bs_ests <- vapply(1:b, 
                      function(x) {
                        # resample residuals
                        smpl <- sample(1:ncol(errs), ncol(errs), 
                                        replace = TRUE)
                        new_y <- ests + errs[, smpl]
                        # get estimate
                        new_y %*% ascm$weights
                      },
                      numeric(tpost))
    # add in means
    bs_ests <- as.matrix(bs_ests, ncol = tpost)
    bs_ests <- rbind(bs_ests, colMeans(bs_ests))
    ## standard errors
    se2 <- apply(bs_ests, 1, var)
    se <- sqrt(se2)


    out <- list()
    att <- predict(ascm, att = T)
    out$att <- c(att, mean(att[(t0 + 1):t_final]))

    out$se <- c(rep(NA, t0), se)
    out$sigma <- NA
    return(out)
}

#' Jackknife+ algorithm over time
#' @param alpha Confidence level
time_jackknife_plus <- function(ascm, alpha = 0.05, conservative = F, ...) {
    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    n_c <- dim(synth_data$Z0)[2]
    Z <- wide_data$Z

    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)

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
    # use leave one out ATT estimates
    # out$att <- c(held_out_errs, 
    #              att[(t0 + 1):t_final], 
    #              mean(att[(t0 + 1):t_final]))
    out$att <- c(att, 
                 mean(att[(t0 + 1):t_final]))
    # held out ATT
    out$heldout_att <- c(held_out_errs, 
                          att[(t0 + 1):t_final], 
                          mean(att[(t0 + 1):t_final]))
    # se <- apply(jack_dist[3,,], 1, sd)
    # out$se <- c(rep(NA, t0), se)
    out$se <- rep(NA, 10 + tpost)
    out$sigma <- NA
    if(conservative) {
        qerr <- quantile(abs(held_out_errs), 1 - alpha)
        out$lb <- c(rep(NA, t0), apply(jack_dist[4,,], 1, min) - qerr)
        out$ub <- c(rep(NA, t0), apply(jack_dist[4,,], 1, max) + qerr)
    } else {
        out$lb <- c(rep(NA, t0), apply(jack_dist[2,,], 1, quantile, alpha / 2))
        out$ub <- c(rep(NA, t0), apply(jack_dist[1,,], 1, quantile, 1 - alpha / 2))
    }
    # shift back to ATT scale
    y1 <- predict(ascm, att = F) + att
    y1 <-  c(y1, mean(y1[(t0 + 1):t_final]))
    shifted_lb <- y1 - out$ub
    shifted_ub <- y1 - out$lb
    out$lb <- shifted_lb
    out$ub <- shifted_ub


    return(out)
}

#' Chernozhukov conformal t method
#' @param K Number of hold out periods
#' @param alpha Confidence level
chernozhukov_t <- function(ascm, K = 3, alpha = 0.05, ...) {
    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    n_c <- dim(synth_data$Z0)[2]
    Z <- wide_data$Z

    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    
    grps <- split(1:t0, sort(1:t0 %% K))

    holdouts <- lapply(grps, 
        function(tdrops) {
            # drop unit i
            new_data <- drop_time_t(wide_data, Z, tdrops)
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
            newt0 <- t0 - length(tdrops)
            est <- predict(new_ascm, att = T)[(newt0 + length(tdrops) + 1):t_final]
            est <- c(est, mean(est))
            err <- c(colMeans(wide_data$X[wide_data$trt == 1,
                                         tdrops,
                                         drop = F]) -
                    predict(new_ascm, att = F)[(newt0 + 1):(newt0 + length(tdrops))])
            list(err, est - mean(err))
        })

    # get errors and jackknife distribution
    centered_atts <- vapply(holdouts, `[[`,
                            numeric(tpost + 1), 2)
    mean_holdout <- rowMeans(centered_atts)
    out <- list()
    att <- predict(ascm, att = T)
    # use centered ATT estimates
    out$att <- c(att[1:t0], mean_holdout)

    se <- sqrt(rowSums((centered_atts - mean_holdout) ^ 2) / (K - 1))
    se <- se / sqrt(K) *
        sqrt(1 + K * pmin(floor(t0 / K), c(rep(1, tpost), tpost)) / tpost)

    out$se <- c(rep(NA, t0), se)
    out$sigma <- NA
    out$lb <- c(rep(NA, t0), mean_holdout - qt(1 - alpha / 2, K - 1) * se)
    out$ub <- c(rep(NA, t0), mean_holdout + qt(1 - alpha / 2, K - 1) * se)

    return(out)
}


#' Drop unit i from data
#' @param wide_data (X, y, trt)
#' @param Z Covariates matrix
#' @param t_drop Time to drop
drop_time_t <- function(wide_data, Z, t_drop) {

        new_wide_data <- list()
        new_wide_data$trt <- wide_data$trt
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

conformal_inf <- function(ascm, alpha = 0.05, type = "block", 
                          ns = 1000, grid_size = 100, ...) {
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
  grid <- c(seq( - 3 * max(abs(att)), 3 * max(abs(att)), length.out = grid_size))

  # iterate over post-treatment periods to get pointwise CIs
  vapply(1:tpost,
         function(j) {
          # fit using t0 + j as a pre-treatment period and get reisduals
          new_wide_data <- wide_data
          new_wide_data$X <- cbind(wide_data$X, wide_data$y[, j, drop = TRUE])
          if(tpost > 1) {
            new_wide_data$y <- wide_data$y[, -j, drop = FALSE]
          } else {
            # TODO: this is a stupid hack because the post period has to be *something*
            new_wide_data$y <- matrix(1, nrow = n, ncol = 1)
          }

          

          compute_permute_ci(new_wide_data, ascm, grid, alpha, type, ns)
         },
         numeric(3)) -> cis
  # get residuals for the average using the sliding average technique
  # if(ncol(wide_data$y) < ncol(wide_data$X)) {
  #   avg_data <- get_sliding_average(wide_data)
  #   avg_synth_data <- format_synth(avg_data$X, avg_data$trt, avg_data$y)
  #   avg_ci <- compute_permute_ci(avg_data, ascm, grid, alpha, type, ns)
  # } else {
    avg_ci <- c(NA, NA, NA)
  # }

  out <- list()
  att <- predict(ascm, att = T)
  out$att <- c(att, mean(att[(t0 + 1):t_final]))
  out$se <- rep(NA, t_final)
  out$sigma <- NA
  out$lb <- c(rep(NA, t0), cis[1, ], avg_ci[1])
  out$ub <- c(rep(NA, t0), cis[2, ], avg_ci[2])
  out$p_val <- c(rep(NA, t0), cis[3, ], avg_ci[3])
  
  return(out)
}


compute_permute_test_stats <- function(wide_data, ascm, h0, type, ns = 1000) {


  # format data
  new_wide_data <- wide_data
  t0 <- ncol(wide_data$X) - 1
  # adjust outcomes for null
  new_wide_data$X[wide_data$trt == 1,t0 + 1] <- new_wide_data$X[wide_data$trt == 1,t0 + 1] - h0
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

  resids <- predict(new_ascm, att = T)[1:(t0 + 1)]
  tpost <- t0 + 1
  # permute residuals and compute test statistic
  if(type == "iid") {
    test_stats <- sapply(1:ns, 
                        function(x) mean(abs(sample(resids)[(t0 + 1):tpost])))
  } else {
    ## increment time by one step and wrap
    test_stats <- sapply(1:tpost,
                        function(j) {
                          reorder <- resids[(1:tpost + j) %% tpost + 1]
                          mean(abs(reorder[(t0 + 1):tpost]))
                        })
  }
  
  return(list(resids = resids, 
              test_stats = test_stats))
}

compute_permute_pval <- function(wide_data, ascm, h0, type, ns = 1000) {
  t0 <- ncol(wide_data$X) - 1
  tpost <- t0 + 1
  out <- compute_permute_test_stats(wide_data, ascm, h0, type, ns)
  mean(mean(abs(out$resids[(t0 + 1):tpost])) <= out$test_stats)
}


compute_permute_ci <- function(wide_data, ascm, grid, alpha, type, ns = 1000) {
  # make sure 0 is in the grid
  grid <- c(grid, 0)
  ps <-sapply(grid, 
              function(x) {
                compute_permute_pval(wide_data, ascm, x, type, ns)
              })
  c(min(grid[ps >= alpha]), max(grid[ps >= alpha]), ps[grid == 0])
}


get_sliding_average <- function(wide_data) {
  tpost <- ncol(wide_data$y)
  t0 <- ncol(wide_data$X)
  # average the data in intervals of tpost
  new_wide_data <- wide_data
  new_wide_data$X <- t(apply(wide_data$X, 1,
                           function(x) {
                             stats::filter(x, 
                                           rep(1 / tpost, tpost), 
                                           sides = 1)[seq(tpost, t0, tpost)]
                           }))
  new_wide_data$y <- as.matrix(rowMeans(wide_data$y), ncol = 1)
  return(new_wide_data)
}

#' Estimate standard errors for single ASCM with residual bootstrap
#' Do this for ridge-augmented synth
#' @param ascm Fitted augsynth object
#' @param refit Whether to refit the weights
bs_se_single <- function(ascm, b=1000, alpha = 0.025, refit = TRUE, ...) {

    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    n_c <- dim(synth_data$Z0)[2]
    Z <- wide_data$Z

    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]


    bs_ests <- vapply(1:b,
                      function(x) {
                        # resample control units
                        n0 <- sum(wide_data$trt == 0)
                        smpl <- sample(1:n0, n0, replace = T)
                        new_data <- resample_controls(smpl, wide_data, Z)
                        # refit
                        if(refit) {
                            new_ascm <- do.call(fit_augsynth_internal,
                                    c(list(wide = new_data$wide,
                                        synth_data = new_data$synth_data,
                                        Z = new_data$Z,
                                        progfunc = ascm$progfunc,
                                        scm = ascm$scm,
                                        fixedeff = ascm$fixedeff),
                                    ascm$extra_args))
                        } else {
                            new_ascm <- ascm
                            new_ascm$weights <- ascm$weights[smpl]
                            new_ascm$weights <- new_ascm$weights / sum(new_ascm$weights)
                            new_ascm$data <- new_data$wide
                            new_ascm$data$synth_data <- new_data$synth_data
                            new_ascm$Z <- new_data$Z
                        }
                       # get ATT estimates, truth is zero
                       est <- predict(new_ascm, att = T)[(t0 + 1):t_final]
                       c(est, mean(est))
                      },
                      numeric(tpost + 1))

    out <- list()
    att <- predict(ascm, att = T)

    se <- apply(bs_ests, 1, sd)
    out$se <- c(rep(NA, t0), se)
    out$sigma <- NA
    out$lb <- c(rep(NA, t0), apply(bs_ests, 1, quantile, alpha ))
    out$ub <- c(rep(NA, t0), apply(bs_ests, 1, quantile, 1 - alpha))


    out$att <- c(att, mean(att[(t0 + 1):t_final]))

    out$se <- c(rep(NA, t0), se)
    out$sigma <- NA
    return(out)
}

#' Drop unit i from data
#' @param wide_data (X, y, trt)
#' @param Z Covariates matrix
resample_controls <- function(smpl, wide_data, Z) {

        new_wide_data <- list()
        new_wide_data$trt <- wide_data$trt
        new_wide_data$X <- wide_data$X
        new_wide_data$y <- wide_data$y

        # sample controls
        X0 <- wide_data$X[wide_data$trt == 0,, drop = F]
        y0 <- wide_data$y[wide_data$trt == 0,, drop = F]
        new_wide_data$X[wide_data$trt == 0,] <- X0[smpl,]

        new_wide_data$y[wide_data$trt == 0, ] <-y0[smpl,] 

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
        new_Z <- Z
        if(!is.null(Z)) {
            Z0 <- Z[trt == 0, , drop = F]
            new_Z[new_wide_data$trt == 0, , drop = F] <- Z0[smpl,, drop = F]
        }

        return(list(wide_data = new_wide_data,
                    synth_data = new_synth_data,
                    Z = new_Z))
}


# #' Estimate standard errors for single ASCM with Cattaneo, Feng, Titiunik
# #' Do this for ridge-augmented synth
# #' @param ascm Fitted augsynth object
# cft_se_single <- function(ascm, b=1000, ...) {

#     wide_data <- ascm$data
#     synth_data <- ascm$data$synth_data
#     n <- nrow(wide_data$X)
#     n_c <- dim(synth_data$Z0)[2]
#     Z <- wide_data$Z

#     t0 <- dim(synth_data$Z0)[1]
#     tpost <- ncol(wide_data$y)
#     t_final <- dim(synth_data$Y0plot)[1]

#     # compute pre-treatment residuals and sufficient stats
#     X0 <- ascm$data$X[ascm$data$trt == 0, ]
#     res <- predict(ascm, att=T)[1:t0]
#     gram_mat <- X0 %*% t(X0) / ncol(X0)
#     corr_vec <- t(X0) * res
#     corr_mat <- corr_vec - rowMeans(corr_vec)
    

#     # montecarlo to get quantiles
#     mc_bounds <- vapply(1:b, 
#                       function(x) {
#                         # sample noise
#                         noise_b <- rnorm(length(t0))
#                         # compute sup and inf over noise
#                         supinf <- compute_cft_bounds(ascm, gram_mat, corr_vec, noise)
#                         rbind(supinf$inf, supinf$sup)
#                       },
#                       matrix(NA, ncol = tpost, nrow = 2))
#     return(mc_bounds)
#     # add in means
#     bs_ests <- rbind(bs_ests, colMeans(bs_ests))
#     ## standard errors
#     se2 <- apply(bs_ests, 1, var)
#     se <- sqrt(se2)


#     out <- list()
#     att <- predict(ascm, att = T)
#     out$att <- c(att, mean(att[(t0 + 1):t_final]))

#     out$se <- c(rep(NA, t0), se)
#     out$sigma <- NA
#     return(out)
# }

# #' Compute sup and inf of errors for prediction quantiles
# compute_cft_bounds <- function(ascm, gram_mat, corr_vec, noise) {

#     t0 <- ncol(ascm$data$X)
#     n0 <- sum(ascm$data$trt == 0)
#     # montecarlo vector in constraint set
#     mc_vec <- colMeans(corr_vec * noise)

#     # make constraint set
#     if(ascm$progfunc == "None") {
#         # if just SCM, threshold weights
#         # FROM CFT_2020_JASA replication materia;s
#         eta <- sqrt(mean(res ^ 2) * log(n0) / t0) / min(apply(X0, 1, sd))
#         new_w <- ascm$weights
#         new_w[ascm$weights < eta] <- 0

#         # quadratic constraint for the errors
#         Qmat <- gram_mat
#         avec <- -2 * mc_vec - 2 * c(t(new_w) %*% gram_mat)
#         dscal <- 2 * sum(mc_vec * new_w) + new_w %*% gram_mat %*% w
#         quad_con <- optiSolve::quadcon(Q = Qmat,
#                                        a = avec,
#                                        d = dscal,
#                                        val= 1 / t0 ^ (7 / 6))
#         # linear constraint for the simplex constraint
#         lin_con <- optiSolve::lincon(A = matrix(1, nrow = n0), val = sum(new_w))
#     }

#     # optimize for lower and upper bounds

# }

#' Use leave out one estimates (placebo gaps) to estimate unit-level variance
#' Do this for ridge-augmented synth
#' @param wide_data Data formatted from format_data
#' @param synth_data Data formatted from foramt_synth
#' @param Z Matrix of auxiliary covariates
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param ridge Include ridge or not
#' @param scm Include SCM or not
#' @noRd
#' @return att estimates, test statistics, p-values
loo_se_ridgeaug <- function(wide_data, synth_data, Z=NULL,
                            lambda=NULL,
                            ridge=T, scm=T) {

    n_c <- dim(synth_data$Z0)[2]

    t0 <- dim(synth_data$Z0)[1]
    t_final <- dim(synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)

    ## att on actual sample
    aug_t <- fit_ridgeaug_formatted(wide_data, synth_data, Z, lambda, ridge, scm)
    att <- as.numeric(synth_data$Y1plot -
            synth_data$Y0plot %*% aug_t$weights)
    lam <- aug_t$lambda

    new_wide_data <- wide_data
    new_wide_data$X <- wide_data$X[wide_data$trt==0,,drop=F]
    new_wide_data$y <- wide_data$y[wide_data$trt==0,,drop=F]


    new_Z <- Z
    if(!is.null(new_Z)) {
        new_Z <- Z[wide_data$trt==0,,drop=F]
    }
    
    
    ## iterate over control units
    for(i in 1:n_c) {

        ## reset synth data to make a control a treated
        new_synth_data <- synth_data
        new_synth_data$Z0 <- synth_data$Z0[, -i]
        new_synth_data$X0 <- synth_data$X0[, -i]        
        new_synth_data$Y0plot <- synth_data$Y0plot[, -i]
        new_synth_data$Z1 <- synth_data$Z0[, i, drop=FALSE]
        new_synth_data$X1 <- synth_data$X0[, i, drop=FALSE]        
        new_synth_data$Y1plot <- synth_data$Y0plot[, i, drop=FALSE]

        ## reset ipw data to change treatment assignment
        new_wide_data$trt <- numeric(nrow(new_wide_data$X))
        new_wide_data$trt[i] <- 1
        
        aug <- fit_ridgeaug_formatted(new_wide_data, new_synth_data, new_Z, lam, ridge, scm)

        ## estimate satt
        errs[i,] <- new_synth_data$Y1plot[(t0+1):t_final,] -
            new_synth_data$Y0plot[(t0+1):t_final,] %*% aug$weights
    }

    ## standard errors

    sig2 <- apply(errs^2, 2, mean) ## estimate of variance

    se2 <- (#1 / sum(wide_data$trt==1) + ## contribution from treated unit
           sum(aug_t$weights^2)) * ## contribution from weights
        sig2

    # se2 <- t(errs) %*% aug_t$weights^2

    se <- sqrt(se2)
    ## se <- (1 / sqrt(sum(wide_data$trt==1)) + 
    ##        sqrt(sum(aug_t$weights^2))) * 
    ##     sig

    out <- list()
    out$att <- att

    out$se <- c(rep(NA, t0), se)
    out$sigma <- sqrt(sig2)
    return(out)
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
        # new_synth_data$Y1plot <- synth_data$Y1plot
        # # get the control unit index
        # yplot <- matrix(0, nrow = nrow(wide_data$X), 
        #                 ncol = ncol(wide_data$X) + ncol(wide_data$y))
        # yplot[wide_data$trt == 0,] <- t(synth_data$Y0plot)
        # yplot[wide_data$trt == 1,] <- synth_data$Y1plot
        # new_synth_data$Y0plot <- t(yplot[-i,][new_wide_data$trt == 0,, drop = F])
        new_Z <- if(!is.null(Z)) Z[-i, , drop = F] else NULL

        return(list(wide_data = new_wide_data,
                    synth_data = new_synth_data,
                    Z = new_Z)) 
}

#' Estimate standard errors with the jackknife
#' Do this for ridge-augmented synth
#' @param wide_data Data formatted from format_data
#' @param synth_data Data formatted from foramt_synth
#' @param Z Matrix of auxiliary covariates
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param ridge Include ridge or not
#' @param scm Include SCM or not
#' @param fixedeff Take out fixed effects or not
#' @noRd
#' @return att estimates, test statistics, p-values
jackknife_se_ridgeaug <- function(wide_data, synth_data, Z=NULL,
                            lambda=NULL,
                            ridge = T, scm = T, fixedeff = F) {
    
    n <- nrow(wide_data$X)
    n_c <- dim(synth_data$Z0)[2]

    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)

    ## att on actual sample
    aug_t <- fit_ridgeaug_formatted(wide_data, synth_data, Z, 
                                    lambda, ridge, scm)
    att <- as.numeric(synth_data$Y1plot -
            synth_data$Y0plot %*% aug_t$weights)
    lam <- aug_t$lambda

    # only drop out control units with non-zero weights
    nnz_weights <- numeric(n)
    nnz_weights[wide_data$trt == 0] <- round(aug_t$weights, 3) != 0

    trt_idxs <- (1:n)[as.logical(nnz_weights)]
    n_jack <- length(trt_idxs)
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
#' @noRd
jackknife_se_single <- function(ascm) {
    
    wide_data <- ascm$data
    # synth_data <- format_synth(wide_data$X, wide_data$trt, wide_data$y)
    synth_data <- ascm$data$synth_data
    n <- nrow(wide_data$X)
    n_c <- dim(synth_data$Z0)[2]
    Z <- wide_data$Z

    t0 <- dim(synth_data$Z0)[1]
    tpost <- ncol(wide_data$y)
    t_final <- dim(synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)

    ## att on actual sample
    # aug_t <- fit_ridgeaug_formatted(wide_data, synth_data, Z, 
    #                                 lambda, ridge, scm)
    # att <- as.numeric(synth_data$Y1plot -
    #         synth_data$Y0plot %*% aug_t$weights)
    # lam <- aug_t$lambda

    # only drop out control units with non-zero weights
    nnz_weights <- numeric(n)
    nnz_weights[wide_data$trt == 0] <- round(ascm$weights, 3) != 0

    trt_idxs <- (1:n)[as.logical(nnz_weights)]
    n_jack <- length(trt_idxs)

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
    out$sigma <- NA
    return(out)
}


#' Compute standard errors using the jackknife
#' @param multisynth fitted multisynth object
#' @param relative Whether to compute effects according to relative time
#' @noRd
jackknife_se_multi <- function(multisynth, relative=NULL) {
    ## get info from the multisynth object
    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    n_leads <- multisynth$n_leads
    n <- nrow(multisynth$data$X)
    outddim <- nrow(predict(multisynth, att=T))

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
    ## return(jack_est)
    se2 <- apply(jack_est, c(1,2),
                function(x) (n-1) / n * sum((x - mean(x,na.rm=T))^2, na.rm=T))

    return(sqrt(se2))

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
    drop_i <- list()
    drop_i$X <- msyn$data$X[-i, , drop = F]
    drop_i$y <- msyn$data$y[-i, , drop = F]
    drop_i$trt <- msyn$data$trt[-i]
    drop_i$mask <- msyn$data$mask[not_miss_j,, drop = F]
    
    long_df <- msyn$long_df
    unit <- colnames(long_df)[1]
    # make alphabetical, because the ith unit is the index in alphabetical ordering
    long_df <- long_df[order(long_df[, unit, drop = TRUE]),]
    ith_unit <- unique(long_df[,unit, drop = TRUE])[i]
    long_df <- long_df[long_df[,unit] != ith_unit,]

    # re-fit everything
    args_list <- list(wide = drop_i, relative = msyn$relative,
                      n_leads = msyn$n_leads, n_lags = msyn$n_lags,
                      nu = msyn$nu, lambda = msyn$lambda,
                      force = msyn$force, n_factors = msyn$n_factors,
                      scm = msyn$scm, time_w = msyn$time_w,
                      lambda_t = msyn$lambda_t,
                      fit_resids = msyn$fit_resids,
                      time_cohort = msyn$time_cohort, long_df = long_df)
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