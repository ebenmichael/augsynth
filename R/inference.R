################################################################################
## Code for inference
################################################################################

#' Use leave out one estimates (placebo gaps) to estimate unit-level variance
#' Do this for ridge-augmented synth
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param ridge Include ridge or not
#' @param scm Include SCM or not
#' 
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
#' param Z Covariates matrix
#' @param i Unit to drop
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

#' Estimate standard errors with the jackknife
#' Do this for ridge-augmented synth
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param ridge Include ridge or not
#' @param scm Include SCM or not
#' @param fixedeff Take out fixed effects or not
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

    # jackknife estimates
    ests <- vapply(trt_idxs, 
                   function(i) {
                       new_data <- drop_unit_i(wide_data, Z, i)
                       if(fixedeff) {
                           demeaned <- demean_data(new_data$wide, 
                                                   new_data$synth_data)
                            new_data$wide_data <- demeaned$wide
                            new_data$synth_data <- demeaned$synth_data
                            mhat <- demeaned$mhat[, (t0 + 1):t_final, drop = F]
                       } else {
                           mhat <- matrix(0, n - 1, tpost)
                       }
                       
                       new_y <- new_data$wide_data$y
                       new_trt <- new_data$wide_data$trt
                       aug <- do.call(fit_ridgeaug_formatted,
                                c(new_data,
                                list(lambda = lam, ridge = ridge, scm = scm)))
                       c(colMeans(mhat[new_trt == 1,, drop = FALSE]) +
                         t(new_y[new_trt == 0,, drop = F]) %*% aug$weights)
                   },
                   numeric(tpost))
    # convert to matrix
    ests <- matrix(ests, nrow = tpost, ncol = length(trt_idxs))
    ## standard errors
    se2 <- apply(ests, 1,
                 function(x) (n - 1) / n * sum((x - mean(x, na.rm = T)) ^ 2))
    se <- sqrt(se2)

    out <- list()
    out$att <- att

    out$se <- c(rep(NA, t0), se)
    out$sigma <- NA
    return(out)
}