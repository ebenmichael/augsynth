################################################################################
## Synth with multiple time periods
################################################################################

multisynth <- function(X, trt, mask, y=NULL,
                       link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                       regularizer=c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc",
                                     "l1_all", "l1_nuc"),
                       lambda=NULL, nlambda=20, lambda.min.ratio=1e-3,
                       interact=F, normalized=TRUE, alpha=1, Q=NULL,
                       ipw_weights=NULL, mu0=NULL,
                       opts=list()) {

    ## map string args to actual params
    params <- map_to_param(X, link, regularizer, ipw_weights, normalized, Q, alpha)
    weightfunc <- params[[1]]
    weightptr <- params[[2]]
    proxfunc <- params[[3]]
    balancefunc <- params[[4]]
    prox_opts <- params[[5]]

    
    prep <- preprocess(X, trt, ipw_weights, "", link, normalized)
    X <- prep$X
    ipw_weights <- prep$ipw_weights

    if(normalized) {
        mask <- cbind(1, mask)
    }
          


    out <- multisynth_(X, trt, mask, weightfunc, weightptr,
                       proxfunc, balancefunc, lambda,
                       nlambda, lambda.min.ratio,
                       ipw_weights, mu0, opts, prox_opts)
    return(out)
}


multisynth_ <- function(X, trt, mask, weightfunc, weightfunc_ptr,
                        proxfunc, balancefunc, lambda=NULL,
                        nlambda=20, lambda.min.ratio=1e-3,
                        ipw_weights=NULL, mu0=NULL, opts=list(),
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


    ## pure controls
    Xc <- X[!is.finite(trt),,drop=FALSE]

    ## global balance
    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))
    avg <- apply(x_t, 1, function(x) sum(x * n1)) / sum(n1)
    x_t <- cbind(avg, x_t)

    ipw_weights <- ipw_weights[!is.finite(trt),,drop=F]
    
    loss_opts = list(Xc=Xc,
                     Xt=x_t,
                     mask=mask,
                     weight_func=weightfunc_ptr,
                     ipw_weights=ipw_weights,
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
    
    apgout <- balancer:::apg_warmstart(make_balancing_grad_multisynth(),
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
                weights[!is.finite(trt),j] <- weightfunc(X[!is.finite(trt),,drop=F],
                                                (theta[,(j+1),drop=F]) * mask[j,],
                                                ipw_weights)
            }
            weights
        })

    out$weights <- weights
    
    out$imbalance <- lapply(out$theta,
                           function(th) balancing_grad_multisynth(as.matrix(th), loss_opts))

    out$lambda <- lambda

    return(out)
}








map_to_param <- function(X, link=c("logit", "linear", "pos-linear", "pos-enet", "posenet"),
                         regularizer=c(NULL, "l1", "grpl1", "l2", "ridge", "linf", "nuc",
                                       "l1_all", "l1_nuc"),
                         ipw_weights=NULL,
                         normalized=F, Q=NULL, alpha=1) {
    #' Map string choices to the proper parameters for the balancer sub-functions
    #' @param X n x d matrix of covariates
    #' @param type Find balancing weights for ATT, subgroup ATTs,
    #'             subgroup ATTs with multilevel p-score, multilevel observational studies,
    #'             ATT with missing outcomes, and heterogeneous effects
    #' @param link Link function for weights
    #' @param regularizer Dual of balance criterion
    #' @param normalized Whether to normalize the weights
    #' @param Q Matrix for generalized ridge, if null then use inverse covariance matrix
    #' @param alpha Elastic net parameter \eqn{\frac{1-\alpha}{2}\|\beta\|_2^2 + \alpha\|\beta\|_1}, defaults to 1
    #'
    #' @return Parameters for balancer
    
    
    if(link == "logit") {
        if(normalized) {
            weightfunc <- balancer:::softmax_weights_ipw
        weightptr <- balancer:::make_softmax_weights_ipw()
        } else {
            weightfunc <- exp_balancer:::weights_ipw
            weightptr <- balancer:::make_exp_weights_ipw()
        }
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
    
    if(is.null(regularizer)) {
        proxfunc <- balancer:::make_no_prox()
    } else if(regularizer == "l1") {
        proxfunc <- if(normalized) balancer:::make_prox_l1_normalized() else balancer:::make_prox_l1()
        balancefunc <- balancer:::linf
    } else if(regularizer == "grpl1") {
        proxfunc <- if(normalized) balancer:::make_prox_l1_grp_normalized() else balancer:::make_prox_l1_grp()
        balancefunc <- balancer:::linf_grp
    } else if(regularizer == "l1grpl1") {
        proxfunc <- if(normalized) balancer:::make_prox_l1_grp_l1_normalized() else balancer:::make_prox_l1_grp_l1()
        ## double the covariate matrix to include two sets of parameters
        X <- cbind(X,X)
        balancefunc <- balancer:::linf_grp_linf
    } else if(regularizer == "l2") {
        proxfunc <- if(normalized) balancer:::make_prox_l2_normalized() else balancer:::make_prox_l2()
        balancefunc <- balancer:::l2
    } else if(regularizer == "ridge") {
        proxfunc <- if(normalized) balancer:::make_prox_l2_sq_normalized() else balancer:::make_prox_l2_sq()
        ## balancefunc <- l2sq
        balancefunc <- function(x) (1 + balancer:::l2(x))^2
    } else if (regularizer=="mahal") {
        proxfunc <- if(normalized) balancer:::make_prox_l2_sq_Q_normalized() else balancer:::make_prox_l2_Q_sq()
        if(is.null(Q)) {
            ## use inverse covariance matrix of Q is null
            ## just do svd once
            Xsvd <- svd(X)
            prox_opts$evec <- Xsvd$v
            prox_opts$eval <- 1 / (Xsvd$d^2 / nrow(X))
            Q <- prox_opts$evec %*% diag(prox_opts$eval) %*% t(prox_opts$evec)
        } else {
            ## eigendecomposition once
            Qsvd <- eigen(Q)
            prox_opts$evec <- Qsvd$vectors
            prox_opts$eval <- Qsvd$values
        }

        if(normalized) {
            Q <- rbind(0, cbind(0, Q))
        }
        balancefunc <- function(x) l2Q(x, Q)

    } else if(regularizer == "enet") {
        proxfunc <- if(normalized) balancer:::make_prox_enet_normalized() else balancer:::make_prox_enet()
        balancefunc <- function(x) (1 - alpha) * (1 + balancer:::l2(x))^2 + alpha * balancer:::linf(x)
    }else if(regularizer == "linf") {
        stop("L infinity regularization not implemented")
    } else if(regularizer == "nuc") {
        proxfunc <- if(normalized) balancer:::make_prox_nuc_normalized()else balancer:::make_prox_nuc()
        balancefunc <- balancer:::op_norm
    } else if(regularizer == "nucl1") {
        proxfunc <- if(normalized) balancer:::make_prox_nuc_l1_normalized() else balancer:::make_prox_nuc_l1()
        ## double the covariate matrix to include two sets of parameters
        X <- cbind(X,X)
        balancefunc <- balancer:::linf_op
    } else if(regularizer == "l1_all") {
        proxfunc <- balancer:::make_prox_l1_all()
    } else if(regularizer == "l1_nuc") {
        proxfunc <- balancer:::make_prox_l1_nuc()
    } else if(regularizer == "multi_ridge") {
        proxfunc <- if(normalized) balancer:::make_prox_multilevel_ridge_normalized() else balancer:::make_prox_multilevel_ridge()
        balancefunc <- function(x) (1 - alpha) * (1 + balancer:::l2(x[,1]))^2 + alpha * (1 + balancer:::l2(x[,1]))^2
    } else if(regularizer == "multi_ridge_nuc") {
        proxfunc <- if(normalized) balancer:::make_prox_multilevel_ridge_nuc_normalized() else balancer:::make_prox_multilevel_ridge_nuc()
        balancefunc <- function(x) (1 - alpha) * (1 + balancer:::l2(x[,1]))^2 + alpha * balancer:::op_norm(x[,-1])
    }else {
        stop("regularizer must be one of (NULL, 'l1', 'grpl1', 'l1grpl1', 'ridge', 'linf', 'nuc', 'nucl1')")
    }


    return(list(weightfunc, weightptr, proxfunc, balancefunc, prox_opts))
}

#' Preprocess covariates
#' 
#' @param X n x d matrix of covariates
#' @param trt Vector of treatment status indicators
#' @param ipw_weights Optional ipw weights to measure dispersion against
#' @param type Find balancing weights for ATT, subgroup ATTs,
#'             subgroup ATTs with multilevel p-score, multilevel observational studies,
#'             ATT with missing outcomes, and heterogeneous effects#'
#' @param link Link function for weights
#' @param normalized Whether to normalize the weights
#'
#' @return Processed covariate matrix
preprocess <- function(X, trt, ipw_weights, type, link, normalized) {

    ## ## center covariates by control means
    ## if(type == "att") {
    ##     X <- apply(X, 2, function(x) x - mean(x[trt==0]))
    ## }

    if(is.null(ipw_weights)) {
        ipw_weights = matrix(1, length(trt), 1)
    } else {
        ipw_weights = matrix(ipw_weights, length(trt), 1)    
    }
    ## add intercept
    if(normalized) {
        ## X <- cbind(sum(trt)/sum(1-trt), X)
        X <- cbind(1, X)
    }
    return(list(X=X, ipw_weights=ipw_weights))
    
}
