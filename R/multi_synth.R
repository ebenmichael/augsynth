################################################################################
## Synth with multiple time periods
################################################################################

#' Internal function to fit synth with staggered adoption
#' @param X Matrix of pre-final intervention outcomes
#' @param trt Vector of treatment levels/times
#' @param relative Whether to re-index time according to treatment date, default T
#' @param gap Number of time periods after treatment to impute control values.
#'            For units treated at time T_j, all units treated after T_j + gap
#'            will be used as control values. If larger than the number of periods,
#'            only never never treated units (pure controls) will be used as comparison units
#' @param link Link function/dispersion function. Default is logit, for internal use only
#' @param regularizer Form of pooling. "multi_ridge" pools parameters towards the group mean.
#'                    "ridge" pools parameters towards the group mean and
#'                    "nuc" imposes a low-rank restriction on the parameters.
#' @param lambda Regularization hyper-parameter. If NULL then solutions will be computed for a
#'               range of lambdas equally spaced on the log scale.
#' @param nlambda Number of values of lambda to consider if lambda is NULL
#' @param lambda.min.ratio Ratio between largest and smallest lambda values.
#' @param alpha Hyper-parameter that controls trade-off between overall and individual balance.
#'              Larger values of alpha place more emphasis on individual balance.
#'              Regularization is
#'                (1-alpha) * lambda ||global|| + alpha * lambda ||individual||
#' @param opts Additional options for optimization
#' 
multisynth_ <- function(X, trt, mask, relative, gap,
                       link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                       regularizer=c("nuc", "ridge"),
                       lambda=NULL, nlambda=20, lambda.min.ratio=1e-3,
                       alpha=1,
                       opts=list()) {

    ## map string args to actual params
    params <- map_to_param(X, link, regularizer, alpha)
    weightfunc <- params[[1]]
    weightptr <- params[[2]]
    proxfunc <- params[[3]]
    balancefunc <- params[[4]]
    prox_opts <- params[[5]]


    if(!relative) {
        out <- multisynth_absolute_(X, trt, mask, gap, weightfunc, weightptr,
                                    proxfunc, balancefunc, lambda,
                                    nlambda, lambda.min.ratio,
                                    opts, prox_opts)
    } else {
       out <- multisynth_relative_(X, trt, mask, gap, weightfunc, weightptr,
                                    proxfunc, balancefunc, lambda,
                                    nlambda, lambda.min.ratio,
                                    opts, prox_opts)
    }
    return(out)
}

#' Internal function to fit synth with staggered adoption according to absolute time
multisynth_absolute_ <- function(X, trt, mask, gap, weightfunc, weightfunc_ptr,
                                 proxfunc, balancefunc, lambda=NULL,
                                 nlambda=20, lambda.min.ratio=1e-3,
                                 opts=list(),
                                 prox_opts=list()) {

    n <- dim(X)[1]

    d <- dim(X)[2]
    
    ## average over treatment times/groups
    grps <- sort(unique(trt[is.finite(trt)]))
    J <- length(grps)
    x_t <- vapply(1:J,
                  function(j) colMeans(X[(trt ==grps[j]), , drop=FALSE]) * mask[j,],
                  numeric(d))
    x_t <- as.matrix(x_t)



    ## global balance
    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))

    ## compute global average according to calendar time
    avg <- apply(x_t, 1, function(x) sum(x * n1)) / sum(n1)
    x_t <- cbind(avg, x_t)

    
    loss_opts = list(X=X,
                     Xt=x_t,
                     mask=mask,
                     trt=trt,
                     unique_trt=grps,
                     gap=gap,
                     weight_func=weightfunc_ptr,
                     n1=n1
                     )
        

    ## initialize at 0 if no initialization is passed
    init <- matrix(0, nrow=d, ncol=ncol(x_t))

    ## combine opts with defaults
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))
    
    
    ## if hyperparam is NULL, start from reference weights and decrease
    if(is.null(lambda)) {
        lam0 <- balancefunc(balancing_grad_multisynth(init, loss_opts))
        lam1 <- lam0 * lambda.min.ratio
        ## decrease on log scale
        lambda <- exp(seq(log(lam0), log(lam1), length.out=nlambda))
    }


    ## collect results
    out <- list()
    out$theta <- matrix(,nrow=d, ncol=length(lambda))
    out$imbalance <- matrix(,nrow=d, ncol=length(lambda))    
    out$weights <- matrix(0, nrow=n, ncol=length(lambda))
    out$weightfunc <- weightfunc


    ## with multiple hyperparameters do warm starts        
    prox_opts = c(prox_opts,
                  list(lam=1))
    
    apgout <- balancer:::apg_warmstart(make_grad_multisynth_absolute(),
                            proxfunc, loss_opts, prox_opts,
                            lambda,
                            opts$x, opts$max_it, opts$eps,
                            opts$alpha, opts$beta, opts$accel, opts$verbose)

    ## weights and theta
    out$theta <- apgout
    out$theta <- lapply(apgout,
                        function(th) cbind(th[,1], th[,1] + th[,-1]))
    weights <- lapply(
        out$theta,
        function(theta) {
            weights <- matrix(0, nrow=dim(X), ncol=(ncol(x_t)-1))
            for(j in 1:(ncol(x_t)-1)) {
                print(c(grps[j], grps[j] + gap, sum(trt > grps[j] + gap)))
                ## restrict to units treated later than T_j + gap
                Xmat <- X[trt > grps[j] + gap,,drop=F]
                weights[trt > grps[j] + gap,j] <- weightfunc(Xmat,
                                                             theta[,(j+1),drop=F])
            }
            weights
        })

    out$weights <- weights
    
    out$imbalance <- lapply(out$theta,
                           function(th) grad_multisynth_absolute(as.matrix(th), loss_opts))

    out$lambda <- lambda

    return(out)
}


#' Internal function to fit synth with staggered adoption according to time relative to treatment
multisynth_relative_ <- function(X, trt, mask, gap, weightfunc, weightfunc_ptr,
                                 proxfunc, balancefunc, lambda=NULL,
                                 nlambda=20, lambda.min.ratio=1e-3,
                                 opts=list(),
                                 prox_opts=list()) {

    n <- dim(X)[1]

    d <- dim(X)[2]
    
    ## average over treatment times/groups
    grps <- sort(unique(trt[is.finite(trt)]))
    J <- length(grps)

    x_t <- vapply(1:J,
                  function(j) {
                      vec <- colMeans(X[(trt ==grps[j]), ,
                                        drop=FALSE]) * mask[j,]
                      c(tail(vec, (d-grps[j])), head(vec, grps[j]))
                  },
                  numeric(d))
    x_t <- as.matrix(x_t)


    ## pure controls
    Xc <- X[!is.finite(trt),,drop=FALSE]

    ## global balance
    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))

    relative <- T

    ## compute glboal average according to time relative to treatment
    avg <- apply(x_t, 1, function(x) sum(x * n1)) / sum(n1)
    x_t <- cbind(avg, x_t)


    
    loss_opts = list(X=X,
                     Xt=x_t,
                     mask=mask,
                     trt=trt,
                     unique_trt=grps,
                     gap=gap,
                     weight_func=weightfunc_ptr,
                     n1=n1
                     )
    

    ## initialize at 0 if no initialization is passed
    init <- matrix(0, nrow=d, ncol=ncol(x_t))

    ## combine opts with defaults
    opts <- c(opts,
              list(max_it=5000,
                   eps=1e-8,
                   alpha=1.01, beta=.9,
                   accel=T,
                   x=init,
                   verbose=F))
    
    
    ## if hyperparam is NULL, start from reference weights and decrease
    if(is.null(lambda)) {
        lam0 <- balancefunc(grad_multisynth_relative(init, loss_opts))
        lam1 <- lam0 * lambda.min.ratio
        ## decrease on log scale
        lambda <- exp(seq(log(lam0), log(lam1), length.out=nlambda))
    }


    ## collect results
    out <- list()
    out$theta <- matrix(,nrow=d, ncol=length(lambda))
    out$imbalance <- matrix(,nrow=d, ncol=length(lambda))    
    out$weights <- matrix(0, nrow=n, ncol=length(lambda))
    out$weightfunc <- weightfunc


    ## with multiple hyperparameters do warm starts        
    prox_opts = c(prox_opts,
                  list(lam=1))
    
    apgout <- balancer:::apg_warmstart(make_grad_multisynth_relative(),
                            proxfunc, loss_opts, prox_opts,
                            lambda,
                            opts$x, opts$max_it, opts$eps,
                            opts$alpha, opts$beta, opts$accel, opts$verbose)

    ## weights and theta
    out$theta <- apgout
    out$theta <- lapply(apgout,
                        function(th) cbind(th[,1], th[,1] + th[,-1]))
    weights <- lapply(
        out$theta,
        function(theta) {
            weights <- matrix(0, nrow=dim(X), ncol=(ncol(x_t)-1))
            for(j in 1:(ncol(x_t)-1)) {
                ## shift the lagged outcomes
                Xmat <- t(apply(X[trt > grps[j] + gap,,drop=F],
                              1,
                              function(vec) c(tail(vec * mask[j,], (d-grps[j])),
                                              head(vec * mask[j,], grps[j]))))
                weights[trt > grps[j] + gap,j] <- weightfunc(Xmat,
                                                             theta[,(j+1),drop=F])
            }
            weights
        })

    out$weights <- weights
    
    out$imbalance <- lapply(out$theta,
                           function(th) grad_multisynth_relative(as.matrix(th), loss_opts))

    out$lambda <- lambda

    return(out)
}








map_to_param <- function(X, link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                         regularizer=c("nuc", "ridge"), alpha=1) {
    #' Map string choices to the proper parameters for the balancer sub-functions
    #' @param X n x d matrix of covariates
    #' @param type Find balancing weights for ATT, subgroup ATTs,
    #'             subgroup ATTs with multilevel p-score, multilevel observational studies,
    #'             ATT with missing outcomes, and heterogeneous effects
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param alpha Elastic net parameter \eqn{\frac{1-\alpha}{2}\|\beta\|_2^2 + \alpha\|\beta\|_1}, defaults to 1
    #'
    #' @return Parameters for balancer
    
    
    if(link == "logit") {
        ## if(normalized) {
        weightfunc <- balancer:::softmax_weights
        weightptr <- balancer:::make_softmax_weights()
        ## } else {
        ##     weightfunc <- exp_balancer:::weights_ipw
        ##     weightptr <- balancer:::make_exp_weights_ipw()
        ## }
    } else if(link == "linear") {
        weightfunc <- balancer:::lin_weights_ipw
        weightptr <- balancer:::make_lin_weights_ipw()
    } else if(link == "pos-linear") {
        weightfunc <- balancer:::pos_lin_weights_ipw
        weightptr <- balancer:::make_pos_lin_weights_ipw()        
    } else if(link == "enet") {
        stop("Elastic Net not impleneted")
    } else if(link == "pos-enet") {
        stop("Elastic Net not impleneted")
    } else {
        stop("link must be one of ('logit', 'linear', 'pos-linear')")
    }


    prox_opts <- list(alpha=alpha)
    
    if(regularizer == "ridge") {
        proxfunc <- balancer:::make_prox_multilevel_ridge()
        balancefunc <- function(x) (1 - alpha) * (1 + balancer:::l2(x[,1]))^2 + alpha * (1 + balancer:::l2(x[,1]))^2
    } else if(regularizer == "nuc") {
        proxfunc <- balancer:::make_prox_multilevel_ridge_nuc()
        balancefunc <- function(x) (1 - alpha) * (1 + balancer:::l2(x[,1]))^2 + alpha * balancer:::op_norm(x[,-1])
    }else {
        stop("regularizer must be one of ('ridge', 'nuc')")
    }


    return(list(weightfunc, weightptr, proxfunc, balancefunc, prox_opts))
}

