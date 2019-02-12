################################################################################
## Fitting, plotting, summarizing staggered synth
################################################################################





#' Fit staggered synth
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
#' @export
multisynth <- function(form, unit, time, data,
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
                           regularizer="nuc",
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

        
        ## get the setting of lambda with the best weighted balance
        best <- which.min((1-alpha) * as.numeric(global_l2) +
                          alpha * as.numeric(ind_op))

        msynth$global_l2 <- global_l2[[best]]
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

        
        msynth$global_l2 <- global_l2[[1]]
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
    
    ## outcome model parameters
    msynth$params <- params
    
    ##format output
    class(msynth) <- "multisynth"
    return(msynth)
}

#' Get prediction of average outcome under control
#' @param multisynth Fit multisynth object
#' @param relative Whether to aggregate estimates according to calendar or relative time
#' @param att Whether to estimate the ATT or the missing counterfactual
#'
#' @return Vector of predicted post-treatment control averages for each treatment group
predict.multisynth <- function(multisynth, relative=NULL, att=F) {


    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    gap <- multisynth$gap
    d <- ncol(multisynth$data$X)
    fulldat <- cbind(multisynth$data$X, multisynth$data$y)
    ttot <- ncol(fulldat)
    grps <- unique(multisynth$data$trt) %>% sort()
    J <- length(grps) - 1
    n1 <- multisynth$data$trt[is.finite(multisynth$data$trt)] %>%
        table() %>% as.numeric()
    fullmask <- cbind(multisynth$data$mask, matrix(0, nrow=J, ncol=(ttot-d)))
    

    ## estimate the post-treatment values to gett att estimates
    mu1hat <- vapply(1:J,
                     function(j) colMeans(fulldat[multisynth$data$trt ==grps[j],
                                                , drop=FALSE]),
                     numeric(ttot))
    ## average outcome model estimates for treatment groups
    y0hat_avg <- vapply(1:J,
                     function(j) colMeans(multisynth$y0hat[multisynth$data$trt ==grps[j],
                                                , drop=FALSE]),
                     numeric(ttot))

    
    ## residuals for treated units
    trt_r <- vapply(1:J,
                    function(j) colMeans(multisynth$residuals[multisynth$data$trt ==grps[j],
                                                , drop=FALSE]),
                    numeric(ttot))

    ## tauhat <- trt_r - t(multisynth$residuals) %*% multisynth$weights
    ## mu0hat <- tauhat - mu1hat
    mu0hat <- y0hat_avg + t(multisynth$residuals) %*% multisynth$weights
    tauhat <- mu1hat - mu0hat
    
    ## re-index time if relative to treatment
    if(relative) {
        total_len <- min(d + gap, ttot + d - grps[1]) ## total length of predictions
        mu0hat <- vapply(1:J,
                         function(j) {
                             vec <- c(rep(NA, d-grps[j]),
                                      mu0hat[1:grps[j],j],
                                      mu0hat[(grps[j]+1):(min(grps[j] + gap, ttot)), j])
                             c(vec, rep(NA, total_len - length(vec)))
                         },
                         numeric(total_len))
        
        tauhat <- vapply(1:J,
                         function(j) {
                             vec <- c(rep(NA, d-grps[j]),
                                      tauhat[1:grps[j],j],
                                      tauhat[(grps[j]+1):(min(grps[j] + gap, ttot)), j])
                             c(vec, rep(NA, total_len - length(vec)))
                         },
                         numeric(total_len))
        ## get the overall average estimate
        avg <- apply(mu0hat, 1, function(z) sum(n1 * z, na.rm=T) / sum(n1 * !is.na(z)))
        mu0hat <- cbind(avg, mu0hat)

        avg <- apply(tauhat, 1, function(z) sum(n1 * z, na.rm=T) / sum(n1 * !is.na(z)))
        tauhat <- cbind(avg, tauhat)
        
    } else {

        ## remove all estimates for t > T_j + gap
        vapply(1:J,
               function(j) c(mu0hat[1:min(grps[j]+gap, ttot),j],
                             rep(NA, max(0, ttot-(grps[j] + gap)))),
               numeric(ttot)) -> mu0hat

        vapply(1:J,
               function(j) c(tauhat[1:min(grps[j]+gap, ttot),j],
                             rep(NA, max(0, ttot-(grps[j] + gap)))),
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
            replace_na(avg2,0) * apply(1-fullmask, 2, max)
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
    cat("\nCall:\n", paste(deparse(multisynth$call), sep="\n", collapse="\n"), "\n\n", sep="")

    ## ## print att estimates
    ## att_post <- colMeans(augsynth$data$y[augsynth$data$trt == 1,,drop=F]) -
    ##     predict(augsynth)

    ## cat(paste("Average ATT Estimate: ",
    ##           format(round(mean(att_post),3), nsmall = 3), "\n\n", sep=""))
}



#' Plot function for multisynth
#' @param relative Whether to estimate effects for time relative to treatment
#' @param levels Treatment levels to plot for, default plots for everything
#' @export
plot.multisynth <- function(multisynth, relative=NULL, levels=NULL) {
    plot(summary(multisynth, relative), levels)
}

compute_se <- function(multisynth, relative=NULL) {


    ## get info from the multisynth object
    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    gap <- multisynth$gap
    d <- ncol(multisynth$data$X)
    fulldat <- cbind(multisynth$data$X, multisynth$data$y)
    ttot <- ncol(fulldat)
    grps <- unique(multisynth$data$trt) %>% sort()
    J <- length(grps) - 1
    n1 <- multisynth$data$trt[is.finite(multisynth$data$trt)] %>%
        table() %>% as.numeric()
    fullmask <- cbind(multisynth$data$mask, matrix(0, nrow=J, ncol=(ttot-d)))
    
    
    ## standard error estimate of treated means from treated residuals
    ## heteroskedastic errors between treated/untreated
    ## trt_var <- vapply(1:J,
    ##                   function(j) {
    ##                       apply(multisynth$residuals[is.finite(multisynth$data$trt),
    ##                                                , drop=FALSE], 2, var) / n1[j]
    ##                   },
    ##                   numeric(ttot))

    ## use weighted control residuals to estimate variance for treated units
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
    
    ## standard error
    se <- sqrt(trt_var + ctrl_var)


    ## re-index time if relative to treatment
    if(relative) {
        total_len <- min(d + gap, ttot + d - grps[1]) ## total length of predictions
        
        se <- vapply(1:J,
                     function(j) {
                         vec <- c(rep(NA, d-grps[j]),
                                  se[1:grps[j],j],
                                  se[(grps[j]+1):(min(grps[j] + gap, ttot)), j])
                         c(vec, rep(NA, total_len - length(vec)))
                         },
                         numeric(total_len))
        ## get the overall standard error estimate
        avg_se <- apply(se, 1, function(z) sqrt(sum(n1^2 * z^2, na.rm=T)) / sum(n1 * !is.na(z)))
        se <- cbind(avg_se, se)
        
    } else {

        ## remove all estimates for t > T_j + gap
        vapply(1:J,
               function(j) c(se[1:min(grps[j]+gap, ttot),j],
                             rep(NA, max(0, ttot-(grps[j] + gap)))),
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
#' @export
summary.multisynth <- function(multisynth, relative=NULL, level=NULL) {

    if(is.null(relative)) {
        relative <- multisynth$relative
    }

    if(is.null(level)) level <- "Average"

    
    grps <- unique(multisynth$data$trt) %>% sort()
    J <- length(grps) - 1    
    gap <- multisynth$gap
    d <- ncol(multisynth$data$X)
    ttot <- d + ncol(multisynth$data$y)

    times <- multisynth$data$time
    
    summ <- list()
    ## post treatment estimate for each group and overall
    att <- predict(multisynth, relative, att=T)
    se <- compute_se(multisynth, relative)
    

    if(relative) {
        att <- data.frame(cbind(-(d-1):min(gap, ttot-grps[1]), att))
        names(att) <- c("Time", "Average", times[grps[1:J]])
        att %>% gather(Level, Estimate, -Time) %>%
            rename("Time"=Time) -> att

        se <- data.frame(cbind(-(d-1):min(gap, ttot-grps[1]), se))
        names(se) <- c("Time", "Average", times[grps[1:J]])
        se %>% gather(Level, Std.Error, -Time) %>%
            rename("Time"=Time) -> se
        
        
    } else {
        att <- data.frame(cbind(times, att))
        names(att) <- c("Time", "Average", times[grps[1:J]])        
        att %>% gather(Level, Estimate, -Time) -> att

        se <- data.frame(cbind(times, se))
        names(se) <- c("Time", "Average", times[grps[1:J]])        
        se %>% gather(Level, Std.Error, -Time) -> se
    }

    summ$att <- inner_join(att, se)


    summ$relative <- relative
    summ$grps <- grps
    summ$call <- multisynth$call
    summ$global_l2 <- multisynth$global_l2
    summ$scaled_global_l2 <- multisynth$scaled_global_l2

    summ$ind_op <- multisynth$ind_op
    summ$scaled_ind_op <- multisynth$scaled_ind_op

    summ$level <- level
    
    class(summ) <- "summary.multisynth"
    return(summ)
}

#' Print function for summary function for multisynth
#' @param level Which treatment level to show summary for, default is average
#' @export
print.summary.multisynth <- function(summ) {
    ## straight from lm
    cat("\nCall:\n", paste(deparse(summ$call), sep="\n", collapse="\n"), "\n\n", sep="")

    level <- summ$level

    first_lvl <- summ$att %>% filter(Level != "Average") %>% pull(Level) %>% min()
    
    ## get ATT estimates for treatment level, post treatment
    if(summ$relative) {
        summ$att %>%
            filter(Time > 0, Level==level) %>%
            rename("Time Since Treatment"=Time) -> att_est
    } else if(level == "Average") {
        summ$att %>% filter(Time > first_lvl, Level=="Average") -> att_est
    } else {
        summ$att %>% filter(Time > level, Level==level) -> att_est
    }


    cat(paste("Global L2 Imbalance (Scaled): ",
              format(round(summ$global_l2,3), nsmall=3), "  (",
              format(round(summ$scaled_global_l2,3), nsmall=3), ")\n\n",
              "Individual Operator Imbalance (Scaled): ",
              format(round(summ$ind_op,3), nsmall=3), "  (",
              format(round(summ$scaled_ind_op,3), nsmall=3), ")\t",
              "\n\n",
              sep=""))


    print(att_est, row.names=F)
    
}

#' Plot function for summary function for multisynth
#' @param levels Treatment levels to plot for, default plots for everything
#' @export
plot.summary.multisynth <- function(summ, levels=NULL) {

    ## get the last time period for each level
    summ$att %>% na.omit() %>%
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
            ggplot2::geom_line() +
            ggrepel::geom_label_repel(ggplot2::aes(label=label),
                                      nudge_x=1, na.rm=T) + 
            ggplot2::geom_hline(yintercept=0, lty=2) -> p

    if(summ$relative) {
        p <- p + ggplot2::geom_vline(xintercept=0, lty=2)
    } else {
        p <- p + ggplot2::geom_vline(aes(xintercept=as.numeric(Level)),
                                     lty=2, alpha=0.5,
                                     summ$att %>% filter(Level != "Average"))
    }

    ## add ses
    if("Average" %in% levels) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=Estimate-2*Std.Error,
                                                   ymax=Estimate+2*Std.Error),
                                      alpha=0.2, color=NA,
                                      data=summ$att %>% filter(Level == "Average"))
    } else {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=Estimate-2*Std.Error,
                                                   ymax=Estimate+2*Std.Error),
                                      alpha=0.2, color=NA)
    }
    
    p <- p + ggplot2::scale_alpha_manual(values=c(1, 0.5)) +
        ggplot2::scale_color_manual(values=c("#333333", "#818181")) +
        ggplot2::guides(alpha=F, color=F) + 
        ggplot2::theme_bw()
    return(p)
    
}
