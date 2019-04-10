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
#' @param alpha Fraction of balance for individual balance
#' @param lambda Regularization hyperparameter, default = 0
#' @param force Include "none", "unit", "time", "two-way" fixed effects. Default: "two-way"
#' @param n_factors Number of factors for interactive fixed effects, default does CV
#'
#' @return augsynth object that contains:
#'         \itemize{
#'          \item{"weights"}{weights}
#'          \item{"data"}{Panel data as matrices}
#'         }
multisynth <- function(form, unit, time, data,
                       relative=T, gap=NULL,
                       alpha=NULL, lambda=0,
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
        ## get residuals from outcome model
        residuals <- cbind(wide$X, wide$y) - y0hat
        
        
    } else if(force == 0 & n_factors == 0) {
        ## if no fixed effects or factors, just do nothing
        y0hat <- matrix(0, nrow=nrow(wide$X), ncol=(ncol(wide$X) + ncol(wide$y)))
        params <- NULL
        ## get residuals from outcome model
        residuals <- cbind(wide$X, wide$y) - y0hat
        
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
    if(typeof(residuals) == "list") {
        bal_mat <- lapply(residuals, function(x) x[,1:ncol(wide$X)])
    } else {
        bal_mat <- residuals[,1:ncol(wide$X)]
    }


    ## if no alpha value is provided, use default based on
    ## global and individual imbalance for no-pooling estimator
    if(is.null(alpha)) {
        ## fit with alpha = 0
        alpha_fit <- multisynth_qp(X=bal_mat,
                                   trt=wide$trt,
                                   mask=wide$mask, gap=gap,
                                   relative=relative,
                                   alpha=0, lambda=lambda)

        alpha <- alpha_fit$ind_l2^2 / (alpha_fit$global_l2^2 + alpha_fit$ind_l2^2)
        
    }
    
    msynth <- multisynth_qp(X=bal_mat,
                            trt=wide$trt,
                            mask=wide$mask, gap=gap,
                            relative=relative,
                            alpha=alpha, lambda=lambda)

    ## put in data and hyperparams
    msynth$data <- wide
    msynth$data$time <- data %>% distinct(!!time) %>% pull(!!time)
    msynth$call <- call_name
    msynth$relative <- relative
    msynth$gap <- gap
    msynth$alpha <- alpha

    ## average together treatment groups
    grps <- unique(wide$trt) %>% sort()
    J <- length(grps)-1
    msynth$grps <- grps
    msynth$y0hat <- y0hat
    msynth$residuals <- residuals

    msynth$n_factors <- n_factors
    msynth$force <- force
    
    ## Get imbalance for uniform weights on raw data
    ## TODO: Get rid of this stupid hack of just fitting the weights again with big lambda
    unif <- multisynth_qp(X=wide$X, ## X=residuals[,1:ncol(wide$X)],
                          trt=wide$trt,
                          mask=wide$mask, gap=gap,
                          relative=relative,
                          alpha=0, lambda=1e10)
    ## scaled global balance
    ## msynth$scaled_global_l2 <- msynth$global_l2  / sqrt(sum(unif$imbalance[,1]^2))
    msynth$scaled_global_l2 <- msynth$global_l2  / unif$global_l2

    ## balance for individual estimates
    ## msynth$scaled_ind_l2 <- msynth$ind_l2  / sqrt(sum(unif$imbalance[,-1]^2))
    msynth$scaled_ind_l2 <- msynth$ind_l2  / unif$ind_l2

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
#' @export
predict.multisynth <- function(multisynth, relative=NULL, att=F) {


    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    gap <- multisynth$gap
    d <- ncol(multisynth$data$X)
    fulldat <- cbind(multisynth$data$X, multisynth$data$y)
    ttot <- ncol(fulldat)
    grps <- multisynth$grps
    J <- length(grps) - 1

    n1 <- sapply(1:J, function(j) sum(multisynth$data$trt == grps[j]))

    fullmask <- cbind(multisynth$data$mask, matrix(0, nrow=J, ncol=(ttot-d)))
    

    ## estimate the post-treatment values to get att estimates
    mu1hat <- vapply(1:J,
                     function(j) colMeans(fulldat[multisynth$data$trt ==grps[j],
                                                , drop=FALSE]),
                     numeric(ttot))


    ## get average outcome model estimates and reweight residuals
    if(typeof(multisynth$y0hat) == "list") {
        mu0hat <- vapply(1:J,
                        function(j) {
                            y0hat <- colMeans(multisynth$y0hat[[j]][multisynth$data$trt ==grps[j],
                                                                  , drop=FALSE])
                            y0hat + t(multisynth$residuals[[j]]) %*%
                                multisynth$weights[,j] / sum(multisynth$weights[,j])
                        }
                       , numeric(ttot)
                        )
    } else {
        mu0hat <- vapply(1:J,
                        function(j) {
                            y0hat <- colMeans(multisynth$y0hat[multisynth$data$trt ==grps[j],
                                                                  , drop=FALSE])
                            y0hat + t(multisynth$residuals) %*%
                                multisynth$weights[,j] / sum(multisynth$weights[,j])
                        }
                       , numeric(ttot)
                        )
    }
    
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
#' @param se Whether to plot standard errors
#' @param jackknice Whether to compute jackknife standard errors, default T
#' @export
plot.multisynth <- function(multisynth, relative=NULL, levels=NULL, se=T, jackknife=T) {
    plot(summary(multisynth, relative, jackknife=jackknife), levels, se)
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


#' Compute standard errors using the jackknife
#' @param multisynth fitted multisynth object
#' @param relative Whether to compute effects according to relative time
jackknife <- function(multisynth, relative=NULL) {
    ## get info from the multisynth object
    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    gap <- multisynth$gap
    n <- nrow(multisynth$data$X)
    outddim <- nrow(predict(multisynth, att=T))
    J <- length(multisynth$grps) - 1
    ## drop each unit and estimate overall treatment effect   
    jack_est <- vapply(1:n,
                       function(i) {
                           msyn_i <- drop_unit_i_(multisynth, i)
                           predict(msyn_i, relative=relative, att=T)[,1]
                       },
                       numeric(outddim))



    
    se <- apply(jack_est, 1, sd, na.rm=T)
    ## se[1:(length(se)-gap-1)] <- NA
    ## pad with NA
    ## TODO: this is a dumb hack to work with existing code easily, fix this
    se <- cbind(se, matrix(NA, nrow=length(se), ncol=J))
    return(se)

}

#' Helper function to drop unit i from the multisynth object
drop_unit_i_ <- function(msyn, i) {

    msyn_i <- msyn
    msyn_i$data$X <- msyn$data$X[-i,]
    msyn_i$data$y <- msyn$data$y[-i,]
    msyn_i$data$trt <- msyn$data$trt[-i]
    if(typeof(msyn$residuals) == "list") {
        msyn_i$residuals <- lapply(msyn$residuals, function(x) x[-i,])
        msyn_i$y0hat <- lapply(msyn$y0hat, function(x) x[-i,])

    } else {
        msyn_i$residuals <- msyn$residuals[-i,]
        msyn_i$y0hat <- msyn$y0hat[-i,]
    }

    msyn_i$weights <- msyn$weights[-i,]

    ## ## refit weights
    ## if(typeof(msyn_i$residuals) == "list") {
    ##     bal_mat <- lapply(msyn_i$residuals, function(x) x[,1:ncol(msyn_i$data$X)])
    ## } else {
    ##     bal_mat <- msyn_i$residuals[,1:ncol(msyn_i$data$X)]
    ## }
    ## msyn_i$grps <- unique(msyn_i$data$trt) %>% sort()
    ## not_miss_j <- msyn$grps[is.finite(msyn$grps)] %in% msyn_i$grps[is.finite(msyn_i$grps)]

    ## msyn_i$data$mask <- msyn$data$mask[not_miss_j,]
    ## msyn_i$weights <- multisynth_qp(X=bal_mat,
    ##                                 trt=msyn_i$data$trt,
    ##                                 mask=msyn_i$data$mask,
    ##                                 gap=msyn_i$gap,
    ##                                 relative=msyn_i$relative,
    ##                                 alpha=msyn_i$alpha, lambda=msyn_i$lambda)$weights
    
    
    return(msyn_i)
}

#' Summary function for multisynth
#' @param relative Whether to estimate effects for time relative to treatment
#' @param level What treatment level to report
#' @param jackknife Whether to use jackknice standard errors, default T
#' @export
summary.multisynth <- function(multisynth, relative=NULL, level=NULL, jackknife=T) {
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
    
    if(jackknife) {
        se <- jackknife(multisynth, relative)
    } else {
        se <- compute_se(multisynth, relative)
    }
    

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

    summ$ind_l2 <- multisynth$ind_l2
    summ$scaled_ind_l2 <- multisynth$scaled_ind_l2
    
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
              "Individual L2 Imbalance (Scaled): ",
              format(round(summ$ind_l2,3), nsmall=3), "  (",
              format(round(summ$scaled_ind_l2,3), nsmall=3), ")\t",
              "\n\n",
              sep=""))


    print(att_est, row.names=F)
    
}

#' Plot function for summary function for multisynth
#' @param levels Treatment levels to plot for, default plots for everything
#' @param se Whether to plot standard errors
#' @export
plot.summary.multisynth <- function(summ, levels=NULL, se=T) {

    ## get the last time period for each level
    summ$att %>%
        filter(!is.na(Estimate)) %>%
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
            ggplot2::geom_point(size=1) + 
            ggrepel::geom_label_repel(ggplot2::aes(label=label),
                                      nudge_x=1, na.rm=T) + 
            ggplot2::geom_hline(yintercept=0, lty=2) -> p

    if(summ$relative) {
        p <- p + ggplot2::geom_vline(xintercept=0, lty=2) +
            xlab("Time Relative to Treatment")
    } else {
        p <- p + ggplot2::geom_vline(aes(xintercept=as.numeric(Level)),
                                     lty=2, alpha=0.5,
                                     summ$att %>% filter(Level != "Average"))
    }

    ## add ses
    if(se) {
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
    }
    
    p <- p + ggplot2::scale_alpha_manual(values=c(1, 0.5)) +
        ggplot2::scale_color_manual(values=c("#333333", "#818181")) +
        ggplot2::guides(alpha=F, color=F) + 
        ggplot2::theme_bw()
    return(p)
    
}
