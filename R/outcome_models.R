################################################################################
## Code to fit various outcome models
################################################################################

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
fit_prog_reg <- function(X, y, trt, alpha=1, lambda=NULL,
                         poly_order=1, type="sep", ...) {
    if(!requireNamespace("glmnet", quietly = TRUE)) {
        stop("In order to fit an elastic net outcome model, you must install the glmnet package.")
    }
    
    X <- matrix(poly(matrix(X),degree=poly_order), nrow=dim(X)[1])

    ## helper function to fit regression with CV
    outfit <- function(x, y) {
        if(is.null(lambda)) {
            lam <- glmnet::cv.glmnet(x, y, alpha=alpha, ...)$lambda.min
        } else {
            lam <- lambda
        }
        fit <- glmnet::glmnet(x, y, alpha=alpha,
                              lambda=lam, ...)
        
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



#' Use a separate random forest regression for each post period
#' to fit E[Y(0)|X]
#'
#' @param X Matrix of covariates/lagged outcomes
#' @param y Matrix of post-period outcomes
#' @param trt Vector of treatment indicator
#' @param avg Predict the average post-treatment outcome
#' @param opts List of options for randomForest
#'             \itemize{\item{avg }{Fit the average post-period rather than time periods separately}}
#'
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_prog_rf <- function(X, y, trt, avg=FALSE, ...) {

    if(!requireNamespace("randomForest", quietly = TRUE)) {
        stop("In order to fit a random forest outcome model, you must install the randomForest package.")
    }

    
    ## helper function to fit RF
    outfit <- function(x, y) {
            fit <- randomForest::randomForest(x, y, ...)
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
fit_prog_gsynth <- function(X, y, trt, r=0, r.end=5, force=3, CV=1, ...) {

    if(!requireNamespace("gsynth", quietly = TRUE)) {
        stop("In order to fit generalized synthetic controls, you must install the gsynth package.")
    }

    
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


#' Use Athey (2017) matrix completion panel data code
#'
#' @param X Matrix of covariates/lagged outcomes
#' @param y Matrix of post-period outcomes
#' @param trt Vector of treatment indicator
#' @param unit_fixed Whether to estimate unit fixed effects
#' @param time_fixed Whether to estimate time fixed effects
#'
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_prog_mcpanel <- function(X, y, trt, unit_fixed=1, time_fixed=1, ...) {


    if(!requireNamespace("MCPanel", quietly = TRUE)) {
        stop("In order to fit matrix completion, you must install the MCPanel package.")
    }

    
    ## create matrix and missingness matrix

    t0 <- dim(X)[2]
    t_final <- t0 + dim(y)[2]
    n <- dim(X)[1]    
    
    fullmat <- cbind(X, y)
    maskmat <- matrix(1, nrow=nrow(fullmat), ncol=ncol(fullmat))
    maskmat[trt==1, (t0+1):t_final] <- 0

    ## estimate matrix
    mcp <- MCPanel::mcnnm_cv(fullmat, maskmat,
                             to_estimate_u=unit_fixed, to_estimate_v=time_fixed)
    
    ## impute matrix
    imp_mat <- mcp$L +
        sweep(matrix(0, nrow=nrow(fullmat), ncol=ncol(fullmat)), 1, mcp$u, "+") + # unit fixed
        sweep(matrix(0, nrow=nrow(fullmat), ncol=ncol(fullmat)), 2, mcp$v, "+") # time fixed
    
    
    trtmat <- matrix(0, ncol=n, nrow=t_final)
    trtmat[t0:t_final, trt == 1] <- 1

    ## get predicted outcomes
    y0hat <- imp_mat[,(t0+1):t_final,drop=FALSE]
    params <- mcp

    params$trt_resids <- colMeans(cbind(X[trt==1,,drop=FALSE],
                                        y[trt==1,,drop=FALSE])) -
        rowMeans(imp_mat[trt==1,,drop=FALSE])

    params$ctrl_resids <- t(cbind(X[trt==0,,drop=FALSE],
                                y[trt==0,,drop=FALSE]) - imp_mat[trt==0,,drop=FALSE])
    params$Y.ct <- t(imp_mat[trt==1,,drop=FALSE])
    return(list(y0hat=y0hat,
                params=params))
    
}


#' Fit a Comparitive interupted time series
#' to fit E[Y(0)|X]
#'
#' @param X Matrix of covariates/lagged outcomes
#' @param y Matrix of post-period outcomes
#' @param trt Vector of treatment indicator
#' @param poly_order Order of time trend polynomial to fit, default 1
#' @param weights Weights to use in WLS, default is no weights
#'
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_prog_cits <- function(X, y, trt, poly_order=1, weights=NULL, ...) {

    ## combine back into a panel structure
    ids <- 1:nrow(X)
    t0 <- dim(X)[2]
    t_final <- t0 + dim(y)[2]
    n <- nrow(X)


    if(is.null(weights)) {
        weights <- rep(1, n)
    }
    
    pnl1 <- data.frame(X)
    colnames(pnl1) <- 1:(t0)

    pnl1 <- pnl1 %>% mutate(trt=trt, post=0, id=ids, weight=weights) %>%
        gather(time, val, -trt, -post, -id, -weight) %>%
        mutate(time=as.numeric(time))

    pnl2 <- data.frame(y)
    colnames(pnl2) <- (t0+1):t_final
    pnl2 <- pnl2 %>% mutate(trt=trt, post=1, id=ids, weight=weights) %>%
        gather(time, val, -trt, -post, -id, -weight) %>%
        mutate(time=as.numeric(time))
    
    
    pnl <- bind_rows(pnl1, pnl2)
    
    ## fit regression
    if(poly_order == "fixed") {
        fit <- pnl %>%
            filter(!((post==1) & (trt==1))) %>% ## filter out post-period treated outcomes
            lm(val ~  as.factor(id) + as.factor(time),
              .,
              weights = .$weight 
              )
    } else if(poly_order > 0) {
        fit <- pnl %>%
            filter(!((post==1) & (trt==1))) %>% ## filter out post-period treated outcomes
        lm(val ~ poly(time, poly_order) + post + trt + poly(time * trt, poly_order),
              ., 
              weights = .$weight
              )
    } else {

        fit <- pnl %>%
            filter(!((post==1) & (trt==1))) %>% ## filter out post-period treated outcomes
            lm(val ~  post + trt,
              .,
              weights = .$weight 
              )
    }

    
    ## get predicted post-period outcomes
    
    y0hat <- matrix(0, nrow=n, ncol=(t_final-t0))
    y0hat[trt==0,]  <- matrix(predict(fit,
                                      pnl %>% filter(post==1 & trt==0)),
                              ncol=ncol(y))

    y0hat[trt==1,] <- matrix(predict(fit,
                                     pnl %>% filter(post==1 & trt==1)),
                             ncol=ncol(y))


    params <- list()

    
    ## add treated prediction for whole pre-period
    params$Y.ct <- matrix(predict(fit,
                                  pnl %>% filter(trt==1),
                                  ncol=(ncol(X) + ncol(y))))

    ## and control prediction
    ctrl_pred <- matrix(predict(fit,
                                pnl %>% filter(trt==0)),
                                ncol=(ncol(X) + ncol(y)))

    ## control and treated residuals
    params$ctrl_resids <- t(cbind(X[trt==0,,drop=FALSE],
                                y[trt==0,,drop=FALSE])) - 
        t(ctrl_pred)
    params$trt_resids <- colMeans(cbind(X[trt==1,,drop=FALSE],
                                            y[trt==1,,drop=FALSE])) -
        rowMeans(params$Y.ct)
    
    return(list(y0hat=y0hat,
                params=params))
    
}




#' Fit a bayesian structural time series
#' to fit E[Y(0)|X]
#'
#' @param X Matrix of covariates/lagged outcomes
#' @param y Matrix of post-period outcomes
#' @param trt Vector of treatment indicator
#'
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Model parameters}}
fit_prog_causalimpact <- function(X, y, trt, ...) {


    if(!requireNamespace("CausalImpact", quietly = TRUE)) {
        stop("In order to fit bayesian structural time series, you must install the CausalImpact package.")
    }

    ## structure data accordingly
    ids <- 1:nrow(X)
    t0 <- dim(X)[2]
    t_final <- t0 + dim(y)[2]
    n <- nrow(X)

    comb <- cbind(X, y)

    imp_dat <- t(rbind(colMeans(comb[trt==1,,drop=F]), comb[trt==0,,drop=F]))

    
    ## get predicted post-period outcomes
    ## TODO: is this the way to use CausalImpact??
    ci_func <- function(i) {
        ## fit causal impact using controls
        CausalImpact::CausalImpact(t(rbind(comb[i,], comb[-i,][trt[-i]==0,])),
                                   pre.period=c(1, t0), post.period=c(t0+1, t_final)
                                   )$series$point.pred
        
    }

    y0hat <- t(sapply(1:n, ci_func))

    params <- list()

    
    ## add treated prediction for whole pre-period
    params$Y.ct <- t(y0hat[trt==1,,drop=F])

    ## and control prediction
    ctrl_pred <- y0hat[trt==0,,drop=F]

    ## control and treated residuals
    params$ctrl_resids <- t(cbind(X[trt==0,,drop=FALSE],
                                y[trt==0,,drop=FALSE])) - 
        t(ctrl_pred)
    
    params$trt_resids <- colMeans(cbind(X[trt==1,,drop=FALSE],
                                            y[trt==1,,drop=FALSE])) -
        rowMeans(params$Y.ct)
    return(list(y0hat=y0hat[,(t0+1):t_final, drop=F],
                params=params))
    
}




#' Fit a seq2seq model with a feedforward net
#' to fit E[Y(0)|X]
#'
#' @param X Matrix of covariates/lagged outcomes
#' @param y Matrix of post-period outcomes
#' @param trt Vector of treatment indicator
#' @param layers List of (n_hidden_units, activation function) pairs to define layers
#' @param epochs Number of epochs for training
#' @param patience Number of epochs to wait before early stopping
#' @param val_split Proportion of control units to use for validation
#' @param verbose Whether to print training progress
#'
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Model parameters}}
fit_prog_seq2seq <- function(X, y, trt,
                             layers=list(c(50, "relu"), c(5, "relu")),
                             epochs=500,
                             patience=5,
                             val_split=0.2,
                             verbose=F, ...) {

    if(!requireNamespace("keras", quietly = TRUE)) {
        stop("In order to fit a neural network, you must install the keras package.")
    }
    
    ## structure data accordingly
    ids <- 1:nrow(X)
    t0 <- dim(X)[2]
    t_final <- t0 + dim(y)[2]
    n <- nrow(X)


    Xctrl <- X[trt==0,,drop=F]
    yctrl <- y[trt==0,,drop=F]

    ## create first layer
    model <- keras::keras_model_sequential() %>%
        keras::layer_dense(units = layers[[1]][1], activation = layers[[1]][2],
                    input_shape = ncol(Xctrl))

    ## add layers
    for(layer in layers[-1]) {
        model %>% keras::layer_dense(units = layer[1], activation = layer[2])
    }

    ## output lyaer
    model %>% keras::layer_dense(units=ncol(yctrl))

    ## compile
    model %>% keras::compile(optimizer="rmsprop", loss="mse", metrics=c("mae")) 

    ## fit model
    learn <- model %>%
        keras::fit(x=Xctrl, y=yctrl,
            epochs=epochs,
            batch_size=nrow(Xctrl),
            validation_split=val_split,
            callbacks=list(keras::callback_early_stopping(patience=patience)),
            verbose=verbose)

    ## predict for everything
    y0hat <- model %>% predict(X)
    params=list(model=model, learn=learn)
    
    return(list(y0hat=y0hat,
                params=params))
}



