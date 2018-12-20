################################################################################
## Code for inference
################################################################################

#' Use leave out one estimates (placebo gaps) to estimate unit-level variance
#' Do this for ridge-augmented synth
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param ridge Include ridge or not
#' @param scm Include SCM or not
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
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

    se2 <- (1 / sum(wide_data$trt==1) + ## contribution from treated unit
           sum(aug_t$weights^2)) * ## contribution from weights
        sig2

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
