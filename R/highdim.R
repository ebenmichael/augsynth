################################################################################
## Code for high dimensional options for Synth
## 1. LASSO as a covariate screen
## 2. Fitting E[Y(0)|X] and inputting into synth/maxent
## 3. DR approach: Fit E[Y(0)|X] and use synth/maxent to balance the residuals
################################################################################

#### Fitting and balancing the prognostic score

fit_prog_reg <- function(X, y, trt, alpha=1, lambda=NULL,
                         poly_order=1, type="sep") {
    #' Use a separate regularized regression for each post period
    #' to fit E[Y(0)|X]
    #'
    #' @param X Matrix of covariates/lagged outcomes
    #' @param y Matrix of post-period outcomes
    #' @param trt Vector of treatment indicator
    #' @param alpha Mixing between L1 and L2, default: 1 (LASSO)
    #' @param lambda Regularization hyperparameter, if null then CV
    #' @param poly_order Order of polynomial to fit, default 1
    #' @param type How to fit outcome model(s)
    #'             \itemize{
    #'              \item{sep }{Separate outcome models}
    #'              \item{avg }{Average responses into 1 outcome}
    #'              \item{multi }{Use multi response regression in glmnet}}
    #'
    #' @return \itemize{
    #'           \item{y0hat }{Predicted outcome under control}
    #'           \item{params }{Regression parameters}}

    X <- matrix(poly(matrix(X),poly_order), nrow=dim(X)[1])

    ## helper function to fit regression with CV
    outfit <- function(x, y) {
        if(is.null(lambda)) {
            lam <- glmnet::cv.glmnet(x, y, alpha=alpha)$lambda.min
        } else {
            lam <- lambda
        }
        fit <- glmnet::glmnet(x, y, alpha=alpha,
                              lambda=lam)
        
        return(as.matrix(coef(fit)))
    }

    if(type=="avg") {
        ## if fitting the average post period value, stack post periods together
        stacky <- c(y)
        stackx <- do.call(rbind,
                          lapply(1:dim(y)[2],
                                 function(x) X))
        stacktrt <- rep(trt, dim(y)[2])
        regweights <- outfit(stackx[stacktrt==0,],
                             stacky[stacktrt==0])
    } else if(type=="sep"){
        ## fit separate regressions for each post period
        regweights <- apply(as.matrix(y), 2,
                            function(yt) outfit(X[trt==0,],
                                                yt[trt==0]))
    } else {
        ## fit multi response regression
        lam <- glmnet::cv.glmnet(X, y, family="mgaussian",
                                 alpha=alpha)$lambda.min
        fit <- glmnet::glmnet(X, y, family="mgaussian",
                              alpha=alpha,
                              lambda=lam)
        regweights <- as.matrix(do.call(cbind, coef(fit)))
    }


    ## Get predicted values
    y0hat <- cbind(rep(1, dim(X)[1]),
                   X) %*% regweights

    return(list(y0hat=y0hat,
                params=regweights))
    
}




fit_prog_rf <- function(X, y, trt, avg=FALSE) {
    #' Use a separate random forest regression for each post period
    #' to fit E[Y(0)|X]
    #'
    #' @param X Matrix of covariates/lagged outcomes
    #' @param y Matrix of post-period outcomes
    #' @param trt Vector of treatment indicator
    #' @param opts List of options for randomForest
    #'             \itemize{\item{avg }{Fit the average post-period rather than time periods separately}}
    #'
    #' @return \itemize{
    #'           \item{y0hat }{Predicted outcome under control}
    #'           \item{params }{Regression parameters}}

    ## helper function to fit RF
    outfit <- function(x, y) {
            fit <- randomForest::randomForest(x, y)
            return(fit)
    }


    if(avg | dim(y)[2] == 1) {
        ## if fitting the average post period value, stack post periods together
        stacky <- c(y)
        stackx <- do.call(rbind,
                          lapply(1:dim(y)[2],
                                 function(x) X))
        stacktrt <- rep(trt, dim(y)[2])
        fit <- outfit(stackx[stacktrt==0,],
                      stacky[stacktrt==0])

        ## predict outcome
        y0hat <- matrix(predict(fit, X), ncol=1)

        
        ## keep feature importances
        imports <- randomForest::importance(fit)

        
    } else {
        ## fit separate regressions for each post period
        fits <- apply(as.matrix(y), 2,
                      function(yt) outfit(X[trt==0,],
                                          yt[trt==0]))
        
        ## predict outcome
        y0hat <- lapply(fits, function(fit) as.matrix(predict(fit,X))) %>%
            bind_rows() %>%
            as.matrix()

        
        ## keep feature importances
        imports <- lapply(fits, function(fit) randomForest::importance(fit)) %>%
            bind_rows() %>%
            as.matrix()

    }


    return(list(y0hat=y0hat,
                params=imports))
    
}



fit_prog_gsynth <- function(X, y, trt, r=0, r.end=5, force=3, CV=1) {
    #' Use gsynth to fit factor model for E[Y(0)|X]
    #'
    #' @param X Matrix of covariates/lagged outcomes
    #' @param y Matrix of post-period outcomes
    #' @param trt Vector of treatment indicator
    #' @param r Number of factors to use (or start with if CV==1)
    #' @param r.end Max number of factors to consider if CV==1
    #' @param force=c(0,1,2,3) Fixed effects (0=none, 1=unit, 2=time, 3=two-way)
    #' @param CV Whether to do CV (0=no CV, 1=yes CV)
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
    capture.output(gsyn <- gsynth:::synth.core(comb, NULL, trtmat, I,
                                               r=r, r.end=r.end,
                                               force=force, CV=CV,
                                               tol=0.001))
    ## get predicted outcomes
    y0hat <- matrix(0, nrow=n, ncol=(t_final-t0))
    y0hat[trt==0,]  <- t(gsyn$Y.co[(t0+1):t_final,,drop=FALSE] -
                         gsyn$est.co$residuals[(t0+1):t_final,,drop=FALSE])

    y0hat[trt==1,] <- gsyn$Y.ct[(t0+1):t_final,]

    ## add treated prediction for whole pre-period
    gsyn$est.co$Y.ct <- gsyn$Y.ct

    ## control and treated residuals
    gsyn$est.co$ctrl_resids <- gsyn$est.co$residuals
    gsyn$est.co$trt_resids <- colMeans(cbind(X[trt==1,,drop=FALSE],
                                            y[trt==1,,drop=FALSE])) -
        rowMeans(gsyn$est.co$Y.ct)
    
    return(list(y0hat=y0hat,
                params=gsyn$est.co))
    
}



fit_prog_complete <- function(X, y, trt, rank.max=5, lambda=0, type="svd") {
    #' Use nuclear norm matrix completion to fit outcome model
    #'
    #' @param X Matrix of covariates/lagged outcomes
    #' @param y Matrix of post-period outcomes
    #' @param trt Vector of treatment indicator
    #' @param rank.max Max rank of the solution
    #' @param lambdaNuclear norm regularization parameter
    #' @param type "svd" is soft-thresholded SVD and "als" is alternating ridge
    #'
    #' @return \itemize{
    #'           \item{y0hat }{Predicted outcome under control}
    #'           \item{params }{Regression parameters}}

    t0 <- dim(X)[2]
    t_final <- t0 + dim(y)[2]
    n <- dim(X)[1]    

    ## construct matrix
    mismat <- matrix(y, nrow=n)
    mismat[trt==1 ,] <- NA
    mismat <- cbind(X,mismat)
    ## fit matrix completion
    fit_comp <- softImpute::softImpute(softImpute::biScale(mismat),
                                       rank.max, lambda, type)

    ## impute matrix
    imp_mat <- softImpute::complete(matrix(NA, ncol=t_final, nrow=n,),
                        fit_comp)
    
    
    trtmat <- matrix(0, ncol=n, nrow=t_final)
    trtmat[t0:t_final, trt == 1] <- 1

    ## get predicted outcomes
    y0hat <- imp_mat[,(t0+1):t_final,drop=FALSE]
    params <- fit_comp

    params$trt_resids <- colMeans(cbind(X[trt==1,,drop=FALSE],
                                        y[trt==1,,drop=FALSE])) -
        rowMeans(imp_mat[trt==1,,drop=FALSE])

    params$ctrl_resids <- t(cbind(X[trt==0,,drop=FALSE],
                                y[trt==0,,drop=FALSE]) - imp_mat[trt==0,,drop=FALSE])
    params$Y.ct <- t(imp_mat[trt==1,,drop=FALSE])
    return(list(y0hat=y0hat,
                params=params))
    
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

    X <- ipw_format$X
    y <- ipw_format$y
    trt <- ipw_format$trt
    
    ## fit prognostic scores
    if(is.null(opts.prog)) {
        fitout <- fit_progscore(X, y, trt)
    } else {
        fitout <- do.call(fit_progscore,
                          c(list(X=X, y=y, trt=trt), opts.prog))
    }

    y0hat <- fitout$y0hat

    ## replace outcomes with fitted prognostic scores
    syn_format$synth_data$Z0 <- t(as.matrix(y0hat[ipw_format$trt == 0,,drop=FALSE]))
    syn_format$synth_data$Z1 <- as.matrix(colMeans(as.matrix(y0hat[ipw_format$trt == 1,,drop=FALSE])))

    ## fit synth/maxent weights
    if(is.null(opts.weights)) {
        syn <- fit_weights(syn_format)
    }
    syn <- do.call(fit_weights, c(list(data_out=syn_format), opts.weights))

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
    } else if(progfunc == "COMP"){
        progf <- fit_prog_complete
    } else {
        stop("progfunc must be one of 'EN', 'RF', 'GSYN', 'COMP'")
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
                                 opts.prog, opts.weights)
                                 

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        data_out$outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        data_out$outcomes <- data_out$outcomes %>% dplyr::arrange_(outcome_col)
    }


    ctrls <- impute_controls(syn_format$outcomes, out, syn_format$trt_unit)

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
    
####### Apply a covariate screen then fit synth

lasso_screen <- function(ipw_format, syn_format, alpha=1, type="sep") {
    #' Screen covariates for the outcome process
    #'
    #' @param ipw_format Output of `format_ipw`
    #' @param syn_format Output of `syn_format`
    #' @param alpha Elastic Net parameter
    #' @param type How to fit outcome model(s)
    #'             \itemize{
    #'              \item{sep }{Separate outcome models}
    #'              \item{avg }{Average responses into 1 outcome}
    #'              \item{multi }{Use multi response regression in glmnet}}
    #'
    #' @return \itemize{
    #'           \item{selX }{Selected covariates}
    #'           \item{params }{Regression parameters}}

    X <- ipw_format$X
    y <- ipw_format$y
    trt <- ipw_format$trt
    
    ## helper function to fit regression with CV
    outfit <- function(x, y) {
            lam <- glmnet::cv.glmnet(x, y, alpha=alpha, intercept=FALSE)$lambda.min
            fit <- glmnet::glmnet(x, y, alpha=alpha,
                                  lambda=lam, intercept=FALSE)
            return(as.matrix(coef(fit))[-1,])
    }

    if(type=="avg") {
        ## if fitting the average post period value, stack post periods together
        stacky <- c(y)
        stackx <- do.call(rbind,
                          lapply(1:dim(y)[2],
                                 function(x) X))
        stacktrt <- rep(trt, dim(y)[2])
        regweights <- outfit(stackx[stacktrt==0,],
                             stacky[stacktrt==0])
    } else if(type=="sep"){
        ## fit separate regressions for each post period
        regweights <- apply(as.matrix(y), 2,
                            function(yt) outfit(X[trt==0,],
                                                yt[trt==0]))
    } else {
        ## fit multi response regression
        lam <- glmnet::cv.glmnet(X, y, family="mgaussian",
                                 alpha=alpha)$lambda.min
        fit <- glmnet::glmnet(X, y, family="mgaussian",
                              alpha=alpha,
                              lambda=lam)
        regweights <- as.matrix(do.call(cbind, coef(fit))[-1,])
    }
    

    ## get covariates with non-zero regression weight
    selected <- apply(as.matrix(regweights), 1, function(beta) 1 * (sum(abs(beta)) > 0))
    ## only return those covariates
    selX <- X[, selected == 1]
    
    return(list(selX=selX,
                params=list(regparams=regweights,
                            selected=selected)))
}



double_screen <- function(ipw_format, syn_format, alpha=1, type="sep", mine=0, by=1) {
    #' Screen covariates for the outcome process with LASSO and selection with
    #' SC with infinity norm
    #'
    #' @param ipw_format Output of `format_ipw`
    #' @param syn_format Output of `syn_format`
    #' @param alpha Elastic net parameter
    #' @param type How to fit outcome model(s)
    #'             \itemize{
    #'              \item{sep }{Separate outcome models}
    #'              \item{avg }{Average responses into 1 outcome}
    #'              \item{multi }{Use multi response regression in glmnet}}
    #' @param mine Smallest imbalance to consider, default 0
    #' @param by Step size for binary search of minimal L infinity error
    #'
    #' @return \itemize{
    #'           \item{selX }{Selected covariates}
    #'           \item{params }{Regression parameters}}

    ## screen covariates for outcome process
    lasout <- lasso_screen(ipw_format, syn_format, alpha, type)

    ## screen using L infinity

    
    ## create the feasibility function by changing the LASSO hyper parameter
    feasfunc <- function(ep) {
        suppressMessages(feas <- fit_entropy_formatted(syn_format, ep, lasso=TRUE)$feasible)
        return(feas)
    }

    
    ## find the best epsilon
    minep <- bin_search(mine, 10 * max(ipw_format$X), by, feasfunc)

    ## if it failed, then stop everything
    if(minep < 0) {
        stop("Failed to find a synthetic control with good enough balance")
    }

    ## fit with minep
    suppressMessages(linfent <- fit_entropy_formatted(syn_format, minep, lasso=TRUE))
    
    ## get covariates with non-zero dual value
    selected_p <- 1 * (abs(linfent$dual) > 0)

    ## take union of outcome and selection covariates
    selected <- 1 * ((selected_p + lasout$params$selected) > 0)

    ## only return those covariates
    selX <- ipw_format$X[, selected == 1]
    
    return(list(selX=selX,
                params=list(regparams=lasout$params$regparams,
                            selparams=linfent$dual,
                            minep=minep,
                            selected=selected)))
}


fit_screensyn_formatted <- function(ipw_format, syn_format,
                                    screen_x, fit_weights,
                                    opts.screen=NULL, opts.weights=NULL) {
    #' Select covariates for E[Y(0)|X], then balance those
    #'
    #' @param ipw_format Output of `format_ipw`
    #' @param syn_format Output of `syn_format`
    #' @param fit_progscore Function to fit prognostic score
    #' @param fit_weights Function to fit synth weights
    #' @param opts.screen Optional options for covariate screening
    #' @param opts.weights Optional options for fitting synth weights
    #' 
    #' @return inverse of predicted propensity scores
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated

    ## fit prognostic scores
    if(is.null(opts.screen)) {
        fitout <- screen_x(ipw_format, syn_format)
    } else {
        fitout <- do.call(screen_x, c(list(ipw_format=ipw_format,
                                         syn_format=syn_format),
                                      opts.screen))
    }

    selX <- fitout$selX
    
    ## replace outcomes with fitted prognostic scores
    syn_format$synth_data$Z0 <- t(as.matrix(selX[ipw_format$trt == 0,]))
    syn_format$synth_data$Z1 <- as.matrix(colMeans(as.matrix(selX[ipw_format$trt == 1,,drop=FALSE])))

    ## fit synth/maxent weights
    if(is.null(opts.weights)) {
        syn <- fit_weights(syn_format)
    } else {
        syn <- do.call(fit_weights,
                       c(list(data_out=syn_format),
                         opts.weights))
    }
        
    

    syn$params <- fitout$params
    return(syn)
}


get_screensyn <- function(outcomes, metadata, trt_unit=1,
                        screenfunc=c("LAS", "2"),
                        weightfunc=c("SC","ENT"),
                        opts.screen = NULL,
                        opts.weights = NULL,
                        outcome_col=NULL,
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated")) {
    #' Fit synthetic controls on estimated outcomes under control
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param screenfunc What function to use to impute control outcomes
    #'                 LAS=LASSO on outcome model,
    #'                 2=LASSO on outcome and L infinity dual on selection
    #' @param weightfunc What function to use to fit weights
    #'                   SC=Vanilla Synthetic Controls, ENT=Maximum Entropy
    #' @param opts.screen Optional options for fitting prognostic score
    #' @param opts.weights Optional options for fitting synth weights    
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## prognostic score and weight functions to use
    if(screenfunc == "LAS") {
        screenf <- lasso_screen
    } else if(screenfunc == "2") {
        screenf <- double_screen
    } else {
        stop("screen must be one of 'LAS', '2'")
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
    out <- fit_screensyn_formatted(ipw_format, syn_format,
                                 screenf, weightf,
                                 opts.screen, opts.weights)
                                 

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        data_out$outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        data_out$outcomes <- data_out$outcomes %>% dplyr::arrange_(outcome_col)
    }


    ctrls <- impute_controls(syn_format$outcomes, out, syn_format$trt_unit)

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


##### Doubly Robust estimation combining and outcome model and selection model

fit_augsyn_formatted <- function(ipw_format, syn_format,
                                fit_progscore, fit_weights,
                                opts.prog=NULL, opts.weights=NULL) {
    #' Fit E[Y(0)|X] and for each post-period and balance pre-period
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

    X <- ipw_format$X
    y <- ipw_format$y
    trt <- ipw_format$trt
    
    ## fit prognostic scores
    if(is.null(opts.prog)) {
        fitout <- fit_progscore(X, y, trt)
    } else {
        fitout <- do.call(fit_progscore,
                          c(list(X=X, y=y, trt=trt),
                            opts.prog))
    }

    y0hat <- fitout$y0hat
    
    ## fit synth/maxent weights
    if(is.null(opts.weights)) {        
        syn <- fit_weights(syn_format)
    } else {
        syn <- do.call(fit_weights,
                       c(list(data_out=syn_format),
                         opts.weights))
    }

    syn$params <- fitout$params

    ## return predicted values for treatment and control
    syn$y0hat_c <- y0hat[ipw_format$trt == 0,]
    syn$y0hat_t <- colMeans(y0hat[ipw_format$trt == 1,,drop=FALSE])

    ## residuals for controls

    syn$resid <- ipw_format$y[ipw_format$trt == 0,] - y0hat[ipw_format$trt == 0,]

    ## difference between observed treated and predicted control
    syn$tauhat <- ipw_format$y[ipw_format$trt == 1,] - y0hat[ipw_format$trt == 1,]

    ## and treated pre outcomes
    syn$treatout <- colMeans(X[trt ==1,,drop=FALSE])

    ## and control pre-outcomes
    syn$pre_ctrls <- ipw_format$X[ipw_format$trt == 0,]
    
    return(syn)
}



impute_synaug <- function(outcomes, metadata, fit, trt_unit) {
    #' Impute the controls after fitting a dr estimator
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe with metadata, in particular a t_int column
    #' @param fit Output of fit_dr
    #'
    #' @return outcomes with additional synthetic control added,
    #'         synth weights
    #'         outcome regression weights

    ## weight control residuals
    wresid <- t(fit$resid) %*% fit$weights

    ## combine weighted residuals and predicted value into DR estimate
    dr <- fit$y0hat_t + wresid

    
    ## combine weighted pre-period controls with
    ## augmented estimate into a "synthetic control"
    dr_ctrl <- c(t(fit$pre_ctrls) %*% fit$weights, dr)

    ## replace true outcome with imputed value
    dr_outcomes <- outcomes %>%
        filter(unit == trt_unit) %>%
        mutate(outcome = dr_ctrl,
               synthetic = "Y",
               potential_outcome = "Y(0)") %>% data.frame()

    ctrls <- outcomes %>% filter(!treated) %>% data.frame()
    avgs <- outcomes %>% filter(unit == trt_unit) %>% data.frame()

    finalout <- bind_rows(ctrls, avgs, dr_outcomes)
    #finalout$outcome <- c(ctrls$outcome, avgs$outcome, dr_outcomes$outcome)
    return(list(outcomes=finalout,
                weights=fit$weights,
                dual=fit$dual,
                outparams=fit$params))
}


get_augsyn <- function(outcomes, metadata, trt_unit=1,
                        progfunc=c("EN", "RF", "GSYN", "COMP"),
                        weightfunc=c("SC","ENT"),
                        opts.prog = NULL,
                        opts.weights = NULL,
                        outcome_col=NULL,
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated")) {
    #' Fit outcome model and balance residuals
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param progfunc What function to use to impute control outcomes
    #'                 EN=Elastic Net, RF=Random Forest, GSYN=gSynth,
    #'                 Comp=softImpute
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
    } else if(progfunc == "COMP"){
        progf <- fit_prog_complete
    } else {
        stop("progfunc must be one of 'EN', 'RF', 'GSYN', 'COMP'")
    }

    if(weightfunc == "SC") {
        weightf <- fit_synth_formatted
    } else if(weightfunc == "ENT") {
        weightf <- fit_entropy_formatted
    } else if(weightfunc == "NONE") {
        ## still fit synth even if none
        ## TODO: This is a dumb wasteful hack
        weightf <- fit_synth_formatted
    } else {
        stop("weightfunc must be one of `SC`, `ENT`, `NONE`")
    }
    
    ## format data
    ipw_format <- format_ipw(outcomes, metadata, outcome_col, cols)
    syn_format <- format_data(outcomes, metadata, trt_unit, outcome_col, cols)

    ## fit outcomes and weights
    out <- fit_augsyn_formatted(ipw_format, syn_format,
                                 progf, weightf,
                                 opts.prog, opts.weights)
                                 

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        data_out$outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        data_out$outcomes <- data_out$outcomes %>% dplyr::arrange_(outcome_col)
    }

    ## if weightfunc is none, set weights to zero
    if(weightfunc == "NONE") {
        out$weights <- rep(0, length(out$weights))
    }

    ctrls <- impute_synaug(syn_format$outcomes, metadata, out, syn_format$trt_unit)

    ## outcome model estimate
    ctrls$outest <- out$tauhat
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
    ctrls$tauhat <- out$tauhat
    ctrls$y0hat_t <- out$y0hat_t
    ctrls$resid <- out$resid
    
    return(ctrls)
}


### Combine synth and gsynth by balancing pre-period residuals


fit_residaug_formatted <- function(ipw_format, syn_format,
                                  fit_progscore, fit_weights,
                                  opts.prog=NULL, opts.weights=NULL) {
    #' Fit E[Y(0)|X] and for each post-period and balance pre-period
    #'
    #' @param ipw_format Output of `format_ipw`
    #' @param syn_format Output of `syn_format`
    #' @param fit_progscore Function to fit prognostic score
    #' @param fit_weights Function to fit synth weights
    #' @param opts.gsyn Optional options for gsynth
    #' @param opts.weights Optional options for fitting synth weights
    #' 
    #' @return inverse of predicted propensity scores
    #'         outcome regression parameters
    #'         control outcomes
    #'         treated outcomes
    #'         boolean for treated

    X <- ipw_format$X
    y <- ipw_format$y
    trt <- ipw_format$trt

    ## fit prognostic scores
    if(is.null(opts.prog)) {
        fitout <- fit_progscore(X, y, trt)
    } else {
        fitout <- do.call(fit_progscore,
                          c(list(X=X, y=y, trt=trt),
                            opts.prog))
    }

    
    ## ## fit prognostic scores
    ## if(is.null(opts.gsyn)) {
    ##     gsyn <- fit_prog_gsynth(X, y, trt)
    ## } else {
    ##     gsyn <- do.call(fit_prog_gsynth,
    ##                       c(list(X=X, y=y, trt=trt),
    ##                         opts.gsyn))
    ## }
    
    y0hat <- fitout$y0hat

    ## get residuals
    ctrl_resids <- fitout$params$ctrl_resids
    trt_resids <- fitout$params$trt_resids

    ## trt_resids <- colMeans(cbind(X[trt==1,,drop=FALSE],
    ##                              y[trt==1,,drop=FALSE])) -
    ##     rowMeans(gsyn$params$Y.ct)
    
    ## replace outcomes with gsynth pre-period residuals
    t0 <- dim(X)[2]

    syn_format$synth_data$Z0 <- ctrl_resids[1:t0, ]
    syn_format$synth_data$Z1 <- as.matrix(trt_resids[1:t0])
    
    ## fit synth/maxent weights
    if(is.null(opts.weights)) {
        syn <- fit_weights(syn_format)
    } else {
        syn <- do.call(fit_weights,
                       c(list(data_out=syn_format),
                         opts.weights))
    }

    syn$params <- fitout$params    

    ## return predicted values for treatment and control
    syn$y0hat_c <- y0hat[ipw_format$trt == 0,]
    syn$y0hat_t <- y0hat[ipw_format$trt == 1,]

    ## residuals for controls
    
    ## difference between observed treated and predicted control
    syn$tauhat <- trt_resids

    ## and treated pre outcomes
    syn$treatout <- colMeans(cbind(X[trt==1,,drop=FALSE], y[trt==1,,drop=FALSE]))

    
    return(syn)
}



impute_residaug <- function(outcomes, metadata, fit, trt_unit) {
    #' Impute the controls after fitting gynsth and reweighting residuals
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe with metadata, in particular a t_int column
    #' @param fit Output of fit_gsynaug_formatted
    #'
    #' @return outcomes with additional synthetic control added,
    #'         synth weights
    #'         outcome regression weights

    ## reweight residuals
    wresid <- fit$params$ctrl_resids %*% fit$weights

    ## combine weighted residuals and predicted value into augmented estimate

    aug_ctrl <- rowMeans(fit$params$Y.ct) + wresid

    ## keep track of difference
    tauhat <- fit$treatout - aug_ctrl

    ## replace true outcome with imputed value
    aug_outcomes <- outcomes %>%
        filter(unit == trt_unit) %>%
        mutate(outcome = aug_ctrl,
               synthetic = "Y",
               potential_outcome = "Y(0)") %>% data.frame()

    ctrls <- outcomes %>% filter(!treated) %>% data.frame()
    avgs <- outcomes %>% filter(unit == trt_unit) %>% data.frame()

    finalout <- bind_rows(ctrls, avgs, aug_outcomes)
    #finalout$outcome <- c(ctrls$outcome, avgs$outcome, dr_outcomes$outcome)
    return(list(outcomes=finalout,
                weights=fit$weights,
                dual=fit$dual,
                outparams=fit$params,
                tauhat_aug=tauhat))
}



get_residaug <- function(outcomes, metadata, trt_unit=1,
                        progfunc=c("GSYN", "COMP"),
                        weightfunc=c("SC","ENT","NONE"),
                        opts.prog = NULL,
                        opts.weights = NULL,
                        outcome_col=NULL,
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated")) {
    #' Fit outcome model and balance residuals
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param progfunc What function to use to impute control outcomes
    #'                 GSYN=gSynth, COMP=Matrix completion
    #' @param weightfunc What function to use to fit weights
    #'                   SC=Vanilla Synthetic Controls, ENT=Maximum Entropy
    #'                   NONE=No reweighting, just gsynth
    #' @param opts.prog Optional options for gsynth
    #' @param opts.weights Optional options for fitting synth weights    
    #' @param outcome_col Column name which identifies outcomes, if NULL then
    #'                    assume only one outcome
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return outcomes with additional synthetic control added and weights
    #' @export

    ## prognostic score and weight functions to use
    if(progfunc == "GSYN"){
        progf <- fit_prog_gsynth
    } else if(progfunc == "COMP"){
        progf <- fit_prog_complete
    } else {
        stop("progfunc must be one of 'GSYN', 'COMP'")
    }

    
    ## weight function to use
    if(weightfunc == "SC") {
        weightf <- fit_synth_formatted
    } else if(weightfunc == "ENT") {
        weightf <- fit_entropy_formatted
    } else if(weightfunc == "NONE") {
        ## still fit synth even if none
        ## TODO: This is a dumb wasteful hack
        weightf <- fit_synth_formatted
    } else {
        stop("weightfunc must be one of 'SC', 'ENT', 'NONE'")
    }
    
    ## format data
    ipw_format <- format_ipw(outcomes, metadata, outcome_col, cols)
    syn_format <- format_data(outcomes, metadata, trt_unit, outcome_col, cols)

    ## fit outcomes and weights
    out <- fit_residaug_formatted(ipw_format, syn_format,
                                 progf, weightf,
                                 opts.prog, opts.weights)
                                 

    ## match outcome types to synthetic controls
    if(!is.null(outcome_col)) {
        data_out$outcomes[[outcome_col]] <- factor(outcomes[[outcome_col]],
                                          levels = names(out$groups))
        data_out$outcomes <- data_out$outcomes %>% dplyr::arrange_(outcome_col)
    }

    ## if no weighting, set weights to 0
    if(weightfunc == "NONE") {
        out$weights <- rep(0, length(out$weights))
    }
    
    ctrls <- impute_residaug(syn_format$outcomes, metadata, out, syn_format$trt_unit)

    ## outcome model estimate
    ctrls$outest <- out$tauhat
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
