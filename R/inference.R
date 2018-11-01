################################################################################
## Code for inference
################################################################################

#' Use leave out one estimates (placebo gaps) to estimate unit-level variance
#' Do this for ridge-augmented synth
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param covs Dataframe of covariates if not null then uses
#'             covariates in outcome model; default NULL
#' @param trt_unit Treated unit
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param scm Include SCM or not
#' @param ridge Include ridge or not
#' @param use_weights Whether to use weights in se estimate, default: FALSE
#' @param hc Type of HC variance estimate (-1 for homoskedastic)
#' @param trt_effect Whether to compute standard errors for treatment effect
#'                   or counterfactual mean
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' 
#' @return att estimates, test statistics, p-values
#' @export
loo_se_ridgeaug <- function(outcomes, metadata, covs=NULL,
                            trt_unit=1, lambda=NULL,
                            scm=T, ridge=T, use_weights=T, hc=-1,
                            trt_effect=T, 
                            cols=list(unit="unit", time="time",
                                      outcome="outcome", treated="treated")) {


    ## format data once
    data_out <- format_data(outcomes, metadata, trt_unit, cols=cols)
    ipw_dat <- format_ipw(outcomes, metadata, cols=cols)

    if(!is.null(covs)) {
        Z <- as.matrix(covs)
    } else {
        Z <- NULL
    }
    
    n_c <- dim(data_out$synth_data$Z0)[2]

    t0 <- dim(data_out$synth_data$Z0)[1]
    t_final <- dim(data_out$synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)

    ## att on actual sample
    if(is.null(Z)) {
        aug_t <- fit_ridgeaug_formatted(ipw_dat, data_out, lambda, scm, ridge)
    } else {
        aug_t <- fit_ridgeaug_cov_formatted(ipw_dat, data_out, Z, lambda, scm)
    }
    att <- as.numeric(data_out$synth_data$Y1plot -
            data_out$synth_data$Y0plot %*% aug_t$weights)
    lam <- aug_t$lambda
    
    ## iterate over control units
    for(i in 1:n_c) {

        ## reset synth data to make a control a treated
        new_data_out <- data_out
        new_data_out$synth_data$Z0 <- data_out$synth_data$Z0[, -i]
        new_data_out$synth_data$Y0plot <- data_out$synth_data$Y0plot[, -i]

        new_data_out$synth_data$Z1 <- data_out$synth_data$Z0[, i, drop=FALSE]
        new_data_out$synth_data$Y1plot <- data_out$synth_data$Y0plot[, i, drop=FALSE]

        ## reset ipw data to change treatment assignment
        new_ipw_dat <- ipw_dat
        new_ipw_dat$X <- ipw_dat$X[ipw_dat$trt==0,,drop=F]
        new_ipw_dat$y <- ipw_dat$y[ipw_dat$trt==0,,drop=F]
        new_ipw_dat$trt <- numeric(nrow(new_ipw_dat$X))
        new_ipw_dat$trt[i] <- 1

        if(is.null(Z)) {
            ## get ridge_aug weights
            aug <- fit_ridgeaug_formatted(new_ipw_dat, new_data_out, lam, scm, ridge)
        } else {
            new_Z <- Z[ipw_dat$trt==0,,drop=F]
            aug <- fit_ridgeaug_cov_formatted(new_ipw_dat, new_data_out, new_Z, lam, scm)
        }


        

        ## estimate satt
        errs[i,] <- new_data_out$synth_data$Y1plot[(t0+1):t_final,] -
            new_data_out$synth_data$Y0plot[(t0+1):t_final,] %*% aug$weights
    }

    ## standard errors
    if(use_weights) {
        if(hc == 0) {
            se <- sqrt(t(errs^2) %*% aug_t$weights^2)
        } else if(hc ==-1) {
            if(trt_effect) {
                se <- (1 / sqrt(sum(ipw_dat$trt==1)) + ## contribution from treated unit
                       sqrt(apply(errs^2, 2, mean))) * ## contribution from weights
                    sqrt(sum(aug_t$weights^2)) ## estimate of variance
            } else {
                se <- sqrt(apply(errs^2, 2, mean)) * sqrt(sum(aug_t$weights^2))
            }
            
        } else if(hc == "scm") {
            ## use synth weights on DI
            sc <- fit_synth_formatted(data_out)
            se <- sqrt(t(errs^2) %*% sc$weights)
        }
    } else {
        se <- sqrt(apply(errs^2, 2, mean))
    }

    ## combine into dataframe
    out <- outcomes %>% distinct(time)
    out$yhat <- as.numeric(data_out$synth_data$Y0plot %*% aug_t$weights)
    out$att <- att

    out$se <- c(rep(NA, t0), se)
    
    return(out)
}
