################################################################################
## Fitting, plotting, summarizing staggered synth
################################################################################

#' Fit staggered synth
#' @param form outcome ~ treatment
#' @param unit Name of unit column
#' @param time Name of time column
#' @param data Panel data as dataframe
#' @param n_leads How long past treatment effects should be estimated for, default is number of post treatment periods for last treated unit
#' @param n_lags Number of pre-treatment periods to balance, default is to balance all periods
#' @param nu Fraction of balance for individual balance
#' @param lambda Regularization hyperparameter, default = 0
#' @param fixedeff Whether to include a unit fixed effect, default F 
#' @param n_factors Number of factors for interactive fixed effects, setting to NULL fits with CV, default is 0
#' @param scm Whether to fit scm weights
#' @param eps_abs Absolute error tolerance for osqp
#' @param eps_rel Relative error tolerance for osqp
#' @param verbose Whether to print logs for osqp
#' @param ... Extra arguments
#'
#' @return multisynth object that contains:
#'         \itemize{
#'          \item{"weights"}{weights matrix where each column is a set of weights for a treated unit}
#'          \item{"data"}{Panel data as matrices}
#'          \item{"imbalance"}{Matrix of treatment minus synthetic control for pre-treatment time periods, each column corresponds to a treated unit}
#'          \item{"global_l2"}{L2 imbalance for the pooled synthetic control}
#'          \item{"scaled_global_l2"}{L2 imbalance for the pooled synthetic control, scaled by the imbalance for unitform weights}
#'          \item{"ind_l2"}{Average L2 imbalance for the individual synthetic controls}
#'          \item{"scaled_ind_l2"}{Average L2 imbalance for the individual synthetic controls, scaled by the imbalance for unitform weights}
#'         \item{"n_leads", "n_lags"}{Number of post treatment outcomes (leads) and pre-treatment outcomes (lags) to include in the analysis}
#'          \item{"nu"}{Fraction of balance for individual balance}
#'          \item{"lambda"}{Regularization hyperparameter}
#'          \item{"scm"}{Whether to fit scm weights}
#'          \item{"grps"}{Time periods for treated units}
#'          \item{"y0hat"}{Pilot estimates of control outcomes}
#'          \item{"residuals"}{Difference between the observed outcomes and the pilot estimates}
#'          \item{"n_factors"}{Number of factors for interactive fixed effects}
#'         }
#' @export
multisynth <- function(form, unit, time, data,
                       n_leads=NULL, n_lags=NULL,
                       nu=NULL, lambda=0,
                       fixedeff = FALSE,
                       n_factors=0,
                       scm=T, 
                       eps_abs = 1e-4,
                       eps_rel = 1e-4,
                       verbose = FALSE, ...) {

    call_name <- match.call()

    form <- Formula::Formula(form)
    unit <- enquo(unit)
    time <- enquo(time)

    ## format data
    outcome <- terms(formula(form, rhs=1))[[2]]
    trt <- terms(formula(form, rhs=1))[[3]]
    wide <- format_data_stag(outcome, trt, unit, time, data)

    force <- if(fixedeff) 3 else 2

    # if n_leads is NULL set it to be the largest possible number of leads
    # for the last treated unit
    if(is.null(n_leads)) {
        n_leads <- ncol(wide$y)
    } else if(n_leads > max(apply(1-wide$mask, 1, sum)) + ncol(wide$y)) {
        n_leads <- max(apply(1-wide$mask, 1, sum)) + ncol(wide$y)
    }

    ## if n_lags is NULL set it to the largets number of pre-treatment periods
    if(is.null(n_lags)) {
        n_lags <- ncol(wide$X)
    } else if(n_lags > ncol(wide$X)) {
        n_lags <- ncol(wide$X)
    }

    msynth <- multisynth_formatted(wide = wide, relative = T,
                                n_leads = n_leads, n_lags = n_lags,
                                nu = nu, lambda = lambda,
                                force = force, n_factors = n_factors,
                                scm = scm, time_w = F,
                                lambda_t = 0,
                                fit_resids = TRUE, eps_abs = eps_abs,
                                eps_rel = eps_rel, verbose = verbose, ...)

    if(scm) {
        ## Get imbalance for uniform weights on raw data
        ## TODO: Get rid of this stupid hack of just fitting the weights again with big lambda
        unif <- multisynth_qp(X=wide$X, ## X=residuals[,1:ncol(wide$X)],
                            trt=wide$trt,
                            mask=wide$mask,
                            n_leads=n_leads,
                            n_lags=n_lags,
                            relative=T,
                            nu=0, lambda=1e10,
                            eps_rel = eps_rel, 
                            eps_abs = eps_abs,
                            verbose = verbose)
        ## scaled global balance
        ## msynth$scaled_global_l2 <- msynth$global_l2  / sqrt(sum(unif$imbalance[,1]^2))
        msynth$scaled_global_l2 <- msynth$global_l2  / unif$global_l2

        ## balance for individual estimates
        ## msynth$scaled_ind_l2 <- msynth$ind_l2  / sqrt(sum(unif$imbalance[,-1]^2))
        msynth$scaled_ind_l2 <- msynth$ind_l2  / unif$ind_l2
    }

    msynth$call <- call_name

    return(msynth)

}


#' Internal funciton to fit staggered synth with formatted data
#' @param wide List containing data elements
#' @param relative Whether to compute balance by relative time
#' @param n_leads How long past treatment effects should be estimated for
#' @param n_lags Number of pre-treatment periods to balance, default is to balance all periods
#' @param nu Fraction of balance for individual balance
#' @param lambda Regularization hyperparameter, default = 0
#' @param force c(0,1,2,3) what type of fixed effects to include
#' @param n_factors Number of factors for interactive fixed effects, default does CV
#' @param scm Whether to fit scm weights
#' @param time_w Whether to fit time weights
#' @param fit_resids Whether to fit SCM on the residuals or not
#' @param eps_abs Absolute error tolerance for osqp
#' @param eps_rel Relative error tolerance for osqp
#' @param verbose Whether to print logs for osqp
#'
#' @return multisynth object
multisynth_formatted <- function(wide, relative=T, n_leads, n_lags,
                       nu, lambda,
                       force,
                       n_factors,
                       scm, time_w,
                       lambda_t,
                       fit_resids,
                       eps_abs, eps_rel,
                       verbose, ...) {

    ## average together treatment groups
    ## grps <- unique(wide$trt) %>% sort()
    grps <- wide$trt[is.finite(wide$trt)]
    J <- length(grps)

    ## fit outcome models
    if(time_w) {
        # Autoregressive model
        out <- fit_time_reg(cbind(wide$X, wide$y), wide$trt,
                            n_leads, lambda_t, ...)
        y0hat <- out$y0hat
        residuals <- out$residuals
        params <- out$time_weights
    } else if(is.null(n_factors)) {
        out <- fit_gsynth_multi(cbind(wide$X, wide$y), wide$trt, force=force)
        y0hat <- out$y0hat
        params <- out$params
        n_factors <- ncol(params$factor)
        ## get residuals from outcome model
        residuals <- cbind(wide$X, wide$y) - y0hat
        
        
    } else if(force == 0 & n_factors == 0) {
        # if no fixed effects or factors, just take out 
        # control averages at each time point
        # time fixed effects from pure controls
        pure_ctrl <- cbind(wide$X, wide$y)[!is.finite(wide$trt), , drop = F]
        y0hat <- matrix(colMeans(pure_ctrl),
                          nrow = nrow(wide$X), ncol = ncol(pure_ctrl), 
                          byrow = T)
        residuals <- cbind(wide$X, wide$y) - y0hat
        params <- NULL
    } else if(force != 0) {
        ## take out pre-treatment averages
        fullmask <- cbind(wide$mask, matrix(0, nrow=nrow(wide$mask),
                                            ncol=ncol(wide$y)))
        out <- fit_feff(cbind(wide$X, wide$y), wide$trt, fullmask, force)
        y0hat <- out$y0hat
        residuals <- out$residuals
        params <- NULL
    } else {
        ## if number of factors is provided don't do CV
        out <- fit_gsynth_multi(cbind(wide$X, wide$y), wide$trt,
                                r=n_factors, r.end=n_factors,
                                CV=0, force=force)
        y0hat <- out$y0hat
        params <- out$params        

        ## get residuals from outcome model
        residuals <- cbind(wide$X, wide$y) - y0hat

    }

    ## balance the residuals
    if(fit_resids) {
        if(time_w) {
            # fit scm on residuals after taking out unit fixed effects
            fullmask <- cbind(wide$mask, matrix(0, nrow=nrow(wide$mask),
                                            ncol=ncol(wide$y)))
            out <- fit_feff(cbind(wide$X, wide$y), wide$trt, fullmask, force)
            bal_mat <- lapply(out$residuals, function(x) x[,1:ncol(wide$X)])
        } else if(typeof(residuals) == "list") {
            bal_mat <- lapply(residuals, function(x) x[,1:ncol(wide$X)])
        } else {
            bal_mat <- residuals[,1:ncol(wide$X)]
        }
    } else {
        # if not balancing residuals, then take out control averages
        # for each time
        ctrl_avg <- matrix(colMeans(wide$X[!is.finite(wide$trt), , drop = F]),
                          nrow = nrow(wide$X), ncol = ncol(wide$X), byrow = T)
        bal_mat <- wide$X - ctrl_avg
        bal_mat <- wide$X
    }
    

    if(scm) {
        ## if no nu value is provided, use default based on
        ## global and individual imbalance for no-pooling estimator
        if(is.null(nu)) {
            ## fit with nu = 0
            nu_fit <- multisynth_qp(X=bal_mat,
                                    trt=wide$trt,
                                    mask=wide$mask,
                                    n_leads=n_leads,
                                    n_lags=n_lags,
                                    relative=relative,
                                    nu=0, lambda=lambda,
                                    eps_rel = eps_rel,
                                    eps_abs = eps_abs,
                                    verbose = verbose)
            ## select nu by triangle inequality ratio
            glbl <- sqrt(sum(nu_fit$imbalance[,1]^2))
            ind <- sum(apply(nu_fit$imbalance[,-1], 2, function(x) sqrt(sum(x^2))))
            nu <- glbl / ind

        }

        msynth <- multisynth_qp(X=bal_mat,
                                trt=wide$trt,
                                mask=wide$mask,
                                n_leads=n_leads,
                                n_lags=n_lags,
                                relative=relative,
                                nu=nu, lambda=lambda,
                                eps_rel = eps_rel,
                                eps_abs = eps_abs,
                                verbose = verbose)
    } else {
        msynth <- list(weights = matrix(0, nrow = nrow(wide$X), ncol = J),
                       imbalance=NA,
                       global_l2=NA,
                       ind_l2=NA)
    }

    ## put in data and hyperparams
    msynth$data <- wide
    msynth$relative <- relative
    msynth$n_leads <- n_leads
    msynth$n_lags <- n_lags
    msynth$nu <- nu
    msynth$lambda <- lambda
    msynth$scm <- scm


    msynth$grps <- grps
    msynth$y0hat <- y0hat
    msynth$residuals <- residuals

    msynth$n_factors <- n_factors
    msynth$force <- force


    ## outcome model parameters
    msynth$params <- params

    # more arguments
    msynth$scm <- scm
    msynth$time_w <- time_w
    msynth$lambda_t <- lambda_t
    msynth$fit_resids <- fit_resids
    msynth$extra_pars <- c(list(eps_abs = eps_abs, 
                                eps_rel = eps_rel, 
                                verbose = verbose), 
                           list(...))

    ##format output
    class(msynth) <- "multisynth"
    return(msynth)
}






#' Get prediction of average outcome under control or ATT
#' @param multisynth Fit multisynth object
#' @param relative Whether to aggregate estimates according to calendar or relative time
#' @param att Whether to estimate the ATT or the missing counterfactual
#'
#' @return Matrix of predicted post-treatment control outcomes for each treated unit
#' @export
predict.multisynth <- function(multisynth, relative=NULL, att=F) {


    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    n_leads <- multisynth$n_leads
    d <- ncol(multisynth$data$X)
    n <- nrow(multisynth$data$X)
    fulldat <- cbind(multisynth$data$X, multisynth$data$y)
    ttot <- ncol(fulldat)
    grps <- multisynth$grps
    J <- length(grps)
    which_t <- (1:n)[is.finite(multisynth$data$trt)]    

    n1 <- sapply(1:J, function(j) 1)
    fullmask <- cbind(multisynth$data$mask, matrix(0, nrow=J, ncol=(ttot-d)))
    

    ## estimate the post-treatment values to get att estimates
    mu1hat <- vapply(1:J,
                     function(j) colMeans(fulldat[which_t[j],
                                                , drop=FALSE]),
                     numeric(ttot))



    ## get average outcome model estimates and reweight residuals
    if(typeof(multisynth$y0hat) == "list") {
        mu0hat <- vapply(1:J,
                        function(j) {
                            y0hat <- colMeans(multisynth$y0hat[[j]][which_t[j],
                                                                  , drop=FALSE])
                            if(!all(multisynth$weights == 0)) {
                                y0hat + t(multisynth$residuals[[j]]) %*%
                                    multisynth$weights[,j] / 
                                    sum(multisynth$weights[,j])
                            } else {
                                y0hat
                            }
                        }
                       , numeric(ttot)
                        )
    } else {
        mu0hat <- vapply(1:J,
                        function(j) {
                            y0hat <- colMeans(multisynth$y0hat[which_t[j],
                                                                  , drop=FALSE])
                            if(!all(multisynth$weights == 0)) {
                                y0hat + t(multisynth$residuals) %*%
                                    multisynth$weights[,j] / 
                                    sum(multisynth$weights[,j])
                            } else {
                                y0hat
                            }
                        }
                       , numeric(ttot)
                        )
    }

    tauhat <- mu1hat - mu0hat

    ## re-index time if relative to treatment
    if(relative) {
        total_len <- min(d + n_leads, ttot + d - min(grps)) ## total length of predictions
        mu0hat <- vapply(1:J,
                         function(j) {
                             vec <- c(rep(NA, d-grps[j]),
                                      mu0hat[1:grps[j],j],
                                      mu0hat[(grps[j]+1):(min(grps[j] + n_leads, ttot)), j])
                             ## last row is post-treatment average
                             c(vec,
                               rep(NA, total_len - length(vec)),
                               mean(mu0hat[(grps[j]+1):(min(grps[j] + n_leads, ttot)), j]))
                               
                         },
                         numeric(total_len +1
                                 ))
        
        tauhat <- vapply(1:J,
                         function(j) {
                             vec <- c(rep(NA, d-grps[j]),
                                      tauhat[1:grps[j],j],
                                      tauhat[(grps[j]+1):(min(grps[j] + n_leads, ttot)), j])
                             ## last row is post-treatment average
                             c(vec,
                               rep(NA, total_len - length(vec)),
                               mean(tauhat[(grps[j]+1):(min(grps[j] + n_leads, ttot)), j]))
                         },
                         numeric(total_len +1
                                 ))
        ## get the overall average estimate
        avg <- apply(mu0hat, 1, function(z) sum(n1 * z, na.rm=T) / sum(n1 * !is.na(z)))
        mu0hat <- cbind(avg, mu0hat)

        avg <- apply(tauhat, 1, function(z) sum(n1 * z, na.rm=T) / sum(n1 * !is.na(z)))
        tauhat <- cbind(avg, tauhat)
        
    } else {

        ## remove all estimates for t > T_j + n_leads
        vapply(1:J,
               function(j) c(mu0hat[1:min(grps[j]+n_leads, ttot),j],
                             rep(NA, max(0, ttot-(grps[j] + n_leads)))),
               numeric(ttot)) -> mu0hat

        vapply(1:J,
               function(j) c(tauhat[1:min(grps[j]+n_leads, ttot),j],
                             rep(NA, max(0, ttot-(grps[j] + n_leads)))),
               numeric(ttot)) -> tauhat

        
        ## only average currently treated units
        avg1 <- rowSums(t(fullmask) *  mu0hat * n1) /
                rowSums(t(fullmask) *  n1)
        avg2 <- rowSums(t(1-fullmask) *  mu0hat * n1) /
            rowSums(t(1-fullmask) *  n1)
        avg <- replace_na(avg1, 0) * apply(fullmask, 2, min) +
            replace_na(avg2,0) * apply(1-fullmask, 2, max)
        cbind(avg, mu0hat) -> mu0hat

        ## only average currently treated units
        avg1 <- rowSums(t(fullmask) *  tauhat * n1) /
            rowSums(t(fullmask) *  n1)
        avg2 <- rowSums(t(1-fullmask) *  tauhat * n1) /
            rowSums(t(1-fullmask) *  n1)
        avg <- replace_na(avg1, 0) * apply(fullmask, 2, min) +
            replace_na(avg2,0) * apply(1 - fullmask, 2, max)
        cbind(avg, tauhat) -> tauhat
    }
    

    if(att) {
        return(tauhat)
    } else {
        return(mu0hat)
    }
}


#' Print function for multisynth
#' @export
print.multisynth <- function(multisynth) {
    ## straight from lm
    cat("\nCall:\n", paste(deparse(multisynth$call), 
        sep="\n", collapse="\n"), "\n\n", sep="")

    # print att estimates
    att_post <- predict(multisynth, att=T)[,1]
    att_post <- att_post[length(att_post)]

    cat(paste("Average ATT Estimate: ",
              format(round(mean(att_post),3), nsmall = 3), "\n\n", sep=""))
}



#' Plot function for multisynth
#' @param levels Treatment levels to plot for, default plots for everything
#' @param se Whether to plot standard errors
#' @param jackknife Whether to compute jackknife standard errors, default T
#' @export
plot.multisynth <- function(multisynth, levels=NULL, 
                            se=T, jackknife=T) {
    plot(summary(multisynth, jackknife=jackknife), levels, se)
}

compute_se <- function(multisynth, relative=NULL) {


    ## get info from the multisynth object
    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    n_leads <- multisynth$n_leads
    d <- ncol(multisynth$data$X)
    fulldat <- cbind(multisynth$data$X, multisynth$data$y)
    ttot <- ncol(fulldat)
    J <- length(multisynth$grps)
    n1 <- multisynth$data$trt[is.finite(multisynth$data$trt)] %>%
        table() %>% as.numeric()
    grps <- multisynth$grps
    fullmask <- cbind(multisynth$data$mask, matrix(0, nrow=J, ncol=(ttot-d)))
    
    

    ## use weighted control residuals to estimate variance for treated units
    if(typeof(multisynth$residuals) == "list") {
        trt_var <- vapply(1:J,
                          function(j) {
                              colSums(multisynth$residuals[[j]]^2 * multisynth$weights[,j]) / n1[j]
                          },
                          numeric(ttot))

        ## standard error estimate of imputed counterfactual mean
        ## from control residuals and weights
        ctrl_var <- vapply(1:J,
                           function(j) colSums(multisynth$residuals[[j]]^2 * multisynth$weights[,j]^2),
                           numeric(ttot))

    } else {
        trt_var <- vapply(1:J,
                          function(j) {
                              colSums(multisynth$residuals^2 * multisynth$weights[,j]) / n1[j]
                          },
                          numeric(ttot))

        ## standard error estimate of imputed counterfactual mean
        ## from control residuals and weights
        ctrl_var <- vapply(1:J,
                           function(j) colSums(multisynth$residuals^2 * multisynth$weights[,j]^2),
                           numeric(ttot))
        
    }
    
    

    ## standard error
    se <- sqrt(trt_var + ctrl_var)

    ## re-index time if relative to treatment
    if(relative) {
        total_len <- min(d + n_leads, ttot + d - min(grps)) ## total length of predictions
        
        se <- vapply(1:J,
                     function(j) {
                         vec <- c(rep(NA, d-grps[j]),
                                  se[1:grps[j],j],
                                  se[(grps[j]+1):(min(grps[j] + n_leads, ttot)), j])
                         c(vec, rep(NA, total_len - length(vec)))
                         },
                         numeric(total_len))
        ## get the overall standard error estimate
        avg_se <- apply(se, 1, function(z) sqrt(sum(n1^2 * z^2, na.rm=T)) / sum(n1 * !is.na(z)))
        se <- cbind(avg_se, se)
        
    } else {

        ## remove all estimates for t > T_j + n_leads
        vapply(1:J,
               function(j) c(se[1:min(grps[j]+n_leads, ttot),j],
                             rep(NA, max(0, ttot-(grps[j] + n_leads)))),
               numeric(ttot)) -> tauhat

        
        ## only average currently treated units
        avg1 <- sqrt(rowSums(t(fullmask) *  se^2 * n1^2)) /
                rowSums(t(fullmask) *  n1)
        avg2 <- sqrt(rowSums(t(1-fullmask) *  se^2 * n1^2)) /
            rowSums(t(1-fullmask) *  n1)
        avg_se <- replace_na(avg1, 0) * apply(fullmask, 2, min) +
            replace_na(avg2,0) * apply(1-fullmask, 2, max)
        se <- cbind(avg_se, se)

    }
    
    
    return(se)
}




#' Summary function for multisynth
#' @param relative Whether to estimate effects for time relative to treatment
#' @param jackknife Whether to use jackknice standard errors, default T
#' 
#' @return summary.multisynth object that contains:
#'         \itemize{
#'          \item{"att"}{Dataframe with ATT estimates, standard errors for each treated unit}
#'          \item{"global_l2"}{L2 imbalance for the pooled synthetic control}
#'          \item{"scaled_global_l2"}{L2 imbalance for the pooled synthetic control, scaled by the imbalance for unitform weights}
#'          \item{"ind_l2"}{Average L2 imbalance for the individual synthetic controls}
#'          \item{"scaled_ind_l2"}{Average L2 imbalance for the individual synthetic controls, scaled by the imbalance for unitform weights}
#'         \item{"n_leads", "n_lags"}{Number of post treatment outcomes (leads) and pre-treatment outcomes (lags) to include in the analysis}
#'         }
#' @export
summary.multisynth <- function(multisynth, jackknife=T) {
    relative <- T

    
    grps <- unique(multisynth$data$trt) %>% sort()
    J <- length(grps)
    n_leads <- multisynth$n_leads
    d <- ncol(multisynth$data$X)
    n <- nrow(multisynth$data$X)
    ttot <- d + ncol(multisynth$data$y)
    which_t <- (1:n)[is.finite(multisynth$data$trt)]
    times <- multisynth$data$time
    
    summ <- list()
    ## post treatment estimate for each group and overall
    att <- predict(multisynth, relative, att=T)
    
    if(jackknife) {
        se <- jackknife_se_multi(multisynth, relative)
    } else {
        se <- compute_se(multisynth, relative)
    }
    

    if(relative) {
        att <- data.frame(cbind(c(-(d-1):min(n_leads, ttot-min(grps)), NA),
                                att))
        names(att) <- c("Time", "Average", 
                        as.character(multisynth$data$units[which_t]))

        att %>% gather(Level, Estimate, -Time) %>%
            rename("Time"=Time) %>%
            mutate(Time=Time-1) -> att

        se <- data.frame(cbind(c(-(d-1):min(n_leads, ttot-min(grps)), NA),
                               se))
        names(se) <- c("Time", "Average", 
                       as.character(multisynth$data$units[which_t]))
        se %>% gather(Level, Std.Error, -Time) %>%
            rename("Time"=Time) %>%
            mutate(Time=Time-1)-> se

    } else {
        att <- data.frame(cbind(times, att))
        names(att) <- c("Time", "Average", times[grps[1:J]])        
        att %>% gather(Level, Estimate, -Time) -> att

        se <- data.frame(cbind(times, se))
        names(se) <- c("Time", "Average", times[grps[1:J]])        
        se %>% gather(Level, Std.Error, -Time) -> se
    }

    summ$att <- inner_join(att, se, by = c("Time", "Level"))


    summ$relative <- relative
    summ$grps <- grps
    summ$call <- multisynth$call
    summ$global_l2 <- multisynth$global_l2
    summ$scaled_global_l2 <- multisynth$scaled_global_l2

    summ$ind_l2 <- multisynth$ind_l2
    summ$scaled_ind_l2 <- multisynth$scaled_ind_l2

    summ$n_leads <- multisynth$n_leads
    summ$n_lags <- multisynth$n_lags

    class(summ) <- "summary.multisynth"
    return(summ)
}

#' Print function for summary function for multisynth
#' @param level Which treatment level to show summary for, default is average
#' @export
print.summary.multisynth <- function(summ, level = NULL) {
    ## straight from lm
    cat("\nCall:\n", paste(deparse(summ$call), sep="\n", collapse="\n"), "\n\n", sep="")

    if(is.null(level)) level <- "Average"

    first_lvl <- summ$att %>% filter(Level != "Average") %>% pull(Level) %>% min()
    
    ## get ATT estimates for treatment level, post treatment
    if(summ$relative) {
        summ$att %>%
            filter(Time >= 0, Level==level) %>%
            rename("Time Since Treatment"=Time) -> att_est
    } else if(level == "Average") {
        summ$att %>% filter(Time > first_lvl, Level=="Average") -> att_est
    } else {
        summ$att %>% filter(Time > level, Level==level) -> att_est
    }


    cat(paste("Global L2 Imbalance (Scaled): ",
              format(round(summ$global_l2,3), nsmall=3), "  (",
              format(round(summ$scaled_global_l2,3), nsmall=3), ")\n\n",
              "Individual L2 Imbalance (Scaled): ",
              format(round(summ$ind_l2,3), nsmall=3), "  (",
              format(round(summ$scaled_ind_l2,3), nsmall=3), ")\t",
              "\n\n",
              sep=""))

    cat(paste("Average ATT Estimate (Std. Error): ",
              summ$att %>%
              filter(Level == level, is.na(Time)) %>%
              pull(Estimate) %>%
              round(3) %>% format(nsmall=3),
              "  (",
              summ$att %>%
              filter(Level == level, is.na(Time)) %>%
              pull(Std.Error) %>%
              round(3) %>% format(nsmall=3),
              ")\n\n", sep=""))

    print(att_est, row.names=F)

}

#' Plot function for summary function for multisynth
#' @param levels Treatment levels to plot for, default plots for everything
#' @param se Whether to plot standard errors
#' @export
plot.summary.multisynth <- function(summ, levels=NULL, se=T) {

    ## get the last time period for each level
    summ$att %>%
        filter(!is.na(Estimate),
               Time >= -summ$n_lags,
               Time <= summ$n_leads) %>%
        group_by(Level) %>%
        summarise(last_time=max(Time)) -> last_times

    if(is.null(levels)) levels <- unique(summ$att$Level)

    summ$att %>% inner_join(last_times) %>%
        filter(Level %in% levels) %>%
        mutate(label=ifelse(Time == last_time, Level, NA),
               is_avg = ifelse(("Average" %in% levels) * (Level == "Average"),
                               "A", "B")) %>%
        ggplot2::ggplot(ggplot2::aes(x=Time, y=Estimate,
                                     group=Level,
                                     color=is_avg,
                                     alpha=is_avg)) +
            ggplot2::geom_line(size=1) +
            ggplot2::geom_point(size=1) + 
            ggrepel::geom_label_repel(ggplot2::aes(label=label),
                                      nudge_x=1, na.rm=T) + 
            ggplot2::geom_hline(yintercept=0, lty=2) -> p

    if(summ$relative) {
        p <- p + ggplot2::geom_vline(xintercept=0, lty=2) +
            ggplot2::xlab("Time Relative to Treatment")
    } else {
        p <- p + ggplot2::geom_vline(aes(xintercept=as.numeric(Level)),
                                     lty=2, alpha=0.5,
                                     summ$att %>% filter(Level != "Average"))
    }

    ## add ses
    if(se) {
        if("Average" %in% levels) {
            p <- p + ggplot2::geom_ribbon(
                ggplot2::aes(ymin=Estimate-2*Std.Error,
                             ymax=Estimate+2*Std.Error),
                alpha=0.2, color=NA,
                data = summ$att %>% 
                        filter(Level == "Average",
                               Time >= 0))
            
        } else {
            p <- p + ggplot2::geom_ribbon(
                ggplot2::aes(ymin=Estimate-2*Std.Error,
                             ymax=Estimate+2*Std.Error),
                             data = . %>% filter(Time >= 0),
                alpha=0.2, color=NA)
        }
    }

    p <- p + ggplot2::scale_alpha_manual(values=c(1, 0.5)) +
        ggplot2::scale_color_manual(values=c("#333333", "#818181")) +
        ggplot2::guides(alpha=F, color=F) + 
        ggplot2::theme_bw()
    return(p)

}
