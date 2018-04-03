################################################################################
## Code for high dimensional options for Synth
## 1. LASSO as a covariate screen
## 2. Fitting E[Y(0)|X] and inputting into synth/maxent
## 3. DR approach: Fit E[Y(0)|X] and use synth/maxent to balance the residuals
################################################################################

#### Fitting and balancing the prognostic score

fit_prog_reg <- function(X, y, trt, opts=list(alpha=1)) {
    #' Use a separate regularized regression for each post period
    #' to fit E[Y(0)|X]
    #'
    #' @param X Matrix of covariates/lagged outcomes
    #' @param y Matrix of post-period outcomes
    #' @param trt Vector of treatment indicator
    #' @param opts List of options for glmnet; choice of alpha, default: 1
    #'
    #'
    #' @return \itemize{
    #'           \item{y0hat }{Predicted outcome under control}
    #'           \item{params }{Regression parameters}}

    ## helper function to fit regression with CV
    outfit <- function(x, y) {
            lam <- glmnet::cv.glmnet(x, y, alpha=opts$alpha)$lambda.min
            fit <- glmnet::glmnet(x, y, alpha=opts$alpha,
                                  lambda=lam)
            return(as.matrix(coef(fit)))
    }
    
    ## fit regressions
    regweights <- apply(y, 2,
                        function(yt) outfit(X[trt==0,],
                                            yt[trt==0]))
    
    ## Get predicted values
    y0hat <- cbind(rep(1, dim(X)[1]),
                   X) %*% regweights

    return(list(y0hat=y0hat,
                params=regweights))
    
}



fit_prog_rf <- function(X, y, trt, opts=NULL) {
    #' Use a separate random forest regression for each post period
    #' to fit E[Y(0)|X]
    #'
    #' @param X Matrix of covariates/lagged outcomes
    #' @param y Matrix of post-period outcomes
    #' @param trt Vector of treatment indicator
    #' @param opts List of options for randomForest
    #'
    #'
    #' @return \itemize{
    #'           \item{y0hat }{Predicted outcome under control}
    #'           \item{params }{Regression parameters}}

    ## helper function to fit RF
    outfit <- function(x, y) {
            fit <- randomForest::randomForest(x, y)
            return(fit)
    }
    
    ## fit regressions
    fits <- apply(y, 2,
                  function(yt) outfit(X[trt==0,],
                                      yt[trt==0]))


    ## fit synth with predicted values
    y0hat <- lapply(fits, function(fit) predict(fit,
                                                X)) %>%
        bind_rows() %>%
        as.matrix()


    ## keep feature importances
    imports <- lapply(fits, function(fit) importance(fit)) %>%
        bind_rows() %>%
        as.matrix()

    return(list(y0hat=y0hat,
                params=imports))
    
}



fit_prog_gsynth <- function(X, y, trt, opts=NULL) {
    #' Use gsynth to fit factor model for E[Y(0)|X]
    #'
    #' @param X Matrix of covariates/lagged outcomes
    #' @param y Matrix of post-period outcomes
    #' @param trt Vector of treatment indicator
    #' @param opts List of options for randomForest
    #'
    #'
    #' @return \itemize{
    #'           \item{y0hat }{Predicted outcome under control}
    #'           \item{params }{Regression parameters}}

    ## matrix with start of treatment
    t0 <- dim(X)[2]
    t_final <- t0 + dim(y)[2]
    n <- dim(X)[1]
    
    trtmat <- matrix(0, ncol=n, nrow=t_final)
    trtmat[t0:t_final, trt == 1] <- 1

    ## observed matrix
    I <- matrix(1, t_final, n)

    ## combine pre and post periods
    comb <- t(cbind(X, y))
    
    ## use internal gsynth function
    capture.output(gsyn <- gsynth:::synth.core(comb, NULL, trtmat, I, force=3, r.end=5, tol=0.001))

    ## get predicted outcomes
    y0hat <- matrix(0, nrow=n, ncol=(t_final-t0+1))
    y0hat[trt==0,]  <- t(gsyn$est.co$residuals[t0:t_final,] + gsyn$Y.co[t0:t_final, ])
    y0hat[trt==1,] <- gsyn$Y.ct[t0:t_final,]

    return(list(y0hat=y0hat,
                params=gsyn$est.co))
    
}


fit_progsyn_formatted <- function(ipw_format, syn_format,
                                  fit_progscore, fit_weights,
                                  opts.prog=NULL, opts.weights=NULL) {
    #' Fit E[Y(0)|X] and for each post-period and balance these
    #'
    #' @param ipw_format Output of `format_ipw`
    #' @param syn_format Output of `syn_format`
    #' @param fit_progscore Function to fit prognostic score
    #' @param fit_weights Function to fit synth weights
    #' @param opts.prog Optional options for fitting prognostic score
    #' @param opts.weights Optional options for fitting synth weights
    #' 
    #' @return inverse of predicted propensity scores
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated

    ## fit prognostic scores
    if(is.null(opts.prog)) {
        fitout <- fit_progscore(ipw_format$X, ipw_format$y, ipw_format$trt)
    } else {
        fitout <- fit_progscore(ipw_format$X, ipw_format$y, ipw_format$trt, opts.prog)
    }

    y0hat <- fitout$y0hat
    
    ## replace outcomes with fitted prognostic scores
    syn_format$synth_data$Z0 <- t(as.matrix(y0hat[ipw_format$trt == 0,]))
    syn_format$synth_data$Z1 <- as.matrix(y0hat[ipw_format$trt == 1,])

    ## fit synth/maxent weights
    syn <- fit_weights(syn_format)

    syn$params <- fitout$params
    return(syn)
}

get_progsyn <- function(outcomes, metadata, trt_unit=1,
                        progfunc=c("EN", "RF", "GSYN"),
                        weightfunc=c("SC","ENT"),
                        opts.prog = NULL,
                        opts.weights = NULL,
                        outcome_col=NULL,
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated")) {
    #' Fit synthetic controls on estimated outcomes under control
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param progfunc What function to use to impute control outcomes
    #'                 EN=Elastic Net, RF=Random Forest, GSYN=gSynth
    #' @param weightfunc What function to use to fit weights
    #'                   SC=Vanilla Synthetic Controls, ENT=Maximum Entropy
    #' @param opts.prog Optional options for fitting prognostic score
    #' @param opts.weights Optional options for fitting synth weights    
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## prognostic score and weight functions to use
    if(progfunc == "EN") {
        progf <- fit_prog_reg
    } else if(progfunc == "RF") {
        progf <- fit_prog_rf
    } else if(progfunc == "GSYN"){
        progf <- fit_prog_gsynth
    } else {
        stop("progfunc must be one of 'EN', 'RF', 'GSYN'")
    }

    if(weightfunc == "SC") {
        weightf <- fit_synth_formatted
    } else if(weightfunc == "ENT") {
        weightf <- fit_entropy_formatted
    }
    
    ## format data
    ipw_format <- format_ipw(outcomes, metadata, outcome_col, cols)
    syn_format <- format_data(outcomes, metadata, trt_unit, outcome_col, cols)

    ## fit weights
    out <- fit_progsyn_formatted(ipw_format, syn_format,
                                 progf, weightf,
                                 opts.prog, opt.weights)
                                 

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        data_out$outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        data_out$outcomes <- data_out$outcomes %>% dplyr::arrange_(outcome_col)
    }


    ctrls <- impute_controls(syn_format$outcomes, out, trt_unit)

    ctrls$params <- out$params
    ctrls$dual <- out$dual
    ctrls$primal_obj <- out$primal_obj
    ctrls$pscores <- out$pscores
    ctrls$eta <- out$eta
    ctrls$groups <- out$groups
    ctrls$feasible <- out$feasible
    ctrls$primal_group_obj <- out$primal_group_obj
    ctrls$scaled_primal_obj <- out$scaled_primal_obj
    ctrls$controls <- out$controls
    return(ctrls)
}
    
