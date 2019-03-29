
#' Fit staggered synth
#' DEPRECATED
#' @param form outcome ~ treatment | auxillary covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param data Panel data as dataframe
#' @param relative Whether to compute balance by relative time
#' @param gap How long past treatment effects should be estimated for
#' @param lambda Regularization hyperparameter
#' @param alpha Fraction of balance for individual balance
#' @param force Include "none", "unit", "time", "two-way" fixed effects. Default: "two-way"
#' @param n_factors Number of factors for interactive fixed effects
#' @param opts_weights Optional options for fitting synth weights
#'
#' @return augsynth object that contains:
#'         \itemize{
#'          \item{"weights"}{weights}
#'          \item{"data"}{Panel data as matrices}
#'         }
multisynth_old <- function(form, unit, time, data,
                       relative=T, gap=NULL,
                       lambda=NULL, alpha=0.5,
                       force="two-way",
                       n_factors=NULL,
                       opts_weights=NULL) {
    
    call_name <- match.call()
    
    form <- Formula::Formula(form)
    unit <- enquo(unit)
    time <- enquo(time)
    
    ## format data
    outcome <- terms(formula(form, rhs=1))[[2]]
    trt <- terms(formula(form, rhs=1))[[3]]
    wide <- format_data_stag(outcome, trt, unit, time, data)

    
    
    ## if gap is NULL set it to be the size of X
    if(is.null(gap)) {
        gap <- ncol(wide$X) + 1
    } else if(gap > ncol(wide$X)) {
        gap <- ncol(wide$X) + 1
    }


    force <- case_when(force == "none" ~ 0,
                       force == "unit" ~ 1,
                       force == "time" ~ 2,
                       TRUE ~ 3)
    ## fit interactive fixed effects model
    if(is.null(n_factors)) {
        out <- fit_gsynth_multi(cbind(wide$X, wide$y), wide$trt, force=force)
        y0hat <- out$y0hat
        params <- out$params
        
    } else if(force == 0 & n_factors == 0) {
        ## if no fixed effects or factors, just do nothing
        y0hat <- matrix(0, nrow=nrow(wide$X), ncol=(ncol(wide$X) + ncol(wide$y)))
        params <- NULL
    } else {
        ## if number of factors is provided don't do CV
        out <- fit_gsynth_multi(cbind(wide$X, wide$y), wide$trt,
                                r=n_factors, r.end=n_factors,
                                CV=0, force=force)
        y0hat <- out$y0hat
        params <- out$params        

    }


    
    ## get residuals from outcome model
    residuals <- cbind(wide$X, wide$y) - y0hat

    
    
    ## fit multisynth
    opts_weights <- c(opts_weights,
                      list(link="logit",
                           regularizer="ridge",
                           nlambda=20, lambda.min.ratio=1e-2,
                           opts=NULL))

    ## If no lambda or multiple lambdas, search over possible lambdas and choose the one with best balance
    if(is.null(lambda) || length(lambda) > 1) {
        suppressWarnings(
            msynth <- multisynth_(X=residuals[,1:ncol(wide$X)],
                                  trt=wide$trt,
                                  mask=wide$mask, gap=gap,
                                  relative=relative, lambda=lambda, alpha=alpha,
                                  link=opts_weights$link, regularizer=opts_weights$regularizer,
                                  nlambda=opts_weights$nlambda,
                                  lambda.min.ratio=opts_weights$lambda.min.ratio,
                                  opts=opts_weights$opts)
        )
        ## Balance for aggregate estimate
        global_l2 <- lapply(msynth$imbalance,
                                   function(imbal) sqrt(sum(imbal[,1]^2)))

        ## balance for individual estimates
        ind_op <- 
            lapply(msynth$imbalance,
                   function(imbal) {
                       if(all(is.finite(imbal))) {
                           svd(imbal[,-1])$d[1]
                       } else {
                           Inf
                       }})

        ## l2 imbalance for individual estimates
        ind_l2 <- lapply(msynth$imbalance,
                         function(imbal)  {
                             sqrt(sum(imbal[,-1]^2))
                         })
        
        ## get the setting of lambda with the best weighted balance
        best <- which.min((1-alpha) * as.numeric(global_l2) +
                          alpha * as.numeric(ind_op))

        msynth$global_l2 <- global_l2[[best]]
        msynth$ind_l2 <- ind_l2[[best]]
        msynth$ind_op <- ind_op[[best]]
        msynth$weights <- msynth$weights[[best]]
        msynth$theta <- msynth$theta[[best]]
        msynth$imbalance <- msynth$imbalance[[best]]
        msynth$lambda <- msynth$lambda[best]

        
    } else {
        msynth <- multisynth_(X=residuals[,1:ncol(wide$X)],
                              trt=wide$trt,
                              mask=wide$mask, gap=gap,
                              relative=relative, lambda=lambda, alpha=alpha,
                              link=opts_weights$link, regularizer=opts_weights$regularizer,
                              nlambda=opts_weights$nlambda,
                              lambda.min.ratio=opts_weights$lambda.min.ratio,
                              opts=opts_weights$opts)

        ## Balance for aggregate estimate
        global_l2 <- lapply(msynth$imbalance,
                                   function(imbal) sqrt(sum(imbal[,1]^2)))

        ## balance for individual estimates
        ind_op <- lapply(msynth$imbalance,
                         function(imbal) svd(imbal[,-1])$d[1])

        ## l2 imbalance for individual estimates
        ind_l2 <- lapply(msynth$imbalance,
                         function(imbal)  {
                             sqrt(sum(imbal[,-1]^2))
                         })
        
        msynth$global_l2 <- global_l2[[1]]
        msynth$ind_l2 <- ind_l2[[1]]        
        msynth$ind_op <- ind_op[[1]]
        msynth$weights <- msynth$weights[[1]]
        msynth$theta <- msynth$theta[[1]]
        msynth$imbalance <- msynth$imbalance[[1]]


    }
    msynth$data <- wide
    msynth$data$time <- data %>% distinct(!!time) %>% pull(!!time)
    msynth$call <- call_name
    msynth$relative <- relative
    msynth$gap <- gap
    msynth$alpha <- alpha

    ## average together treatment groups
    grps <- unique(wide$trt)
    J <- length(grps)-1

    msynth$y0hat <- y0hat
    msynth$residuals <- residuals
    
    ## Get imbalance for uniform weights on raw data
    ## TODO: Get rid of this stupid hack of just fitting the weights again with zero steps
    unif <- multisynth_(X=wide$X, ## X=residuals[,1:ncol(wide$X)],
                        trt=wide$trt,
                        mask=wide$mask, gap=gap,
                        relative=relative, lambda=lambda, alpha=alpha,
                        link=opts_weights$link, regularizer=opts_weights$regularizer,
                        nlambda=opts_weights$nlambda,
                        lambda.min.ratio=opts_weights$lambda.min.ratio,
                        opts=list(max_it=0))
    ## Balance for aggregate estimate
    msynth$scaled_global_l2 <- msynth$global_l2  / sqrt(sum(unif$imbalance[[1]][,1]^2))

    ## balance for individual estimates
    msynth$scaled_ind_op <- msynth$ind_op / svd(unif$imbalance[[1]][,-1])$d[1]
    msynth$scaled_ind_l2 <- msynth$ind_l2  / sqrt(sum(unif$imbalance[[1]][,-1]^2))
    
    ## outcome model parameters
    msynth$params <- params
    
    ##format output
    class(msynth) <- "multisynth"
    return(msynth)
}
