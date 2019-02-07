################################################################################
## Fitting, plotting, summarizing staggered synth
################################################################################





#' Fit staggered synth
#' @param form outcome ~ treatment | auxillary covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param data Panel data as dataframe
#' @param opts_weights Optional options for fitting synth weights
#'
#' @return augsynth object that contains:
#'         \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#'          \item{"data"}{Panel data as matrices}
#'         }
#' @export
multisynth <- function(form, unit, time, data,
                       relative=T, gap=NULL,
                       lambda=NULL, alpha=0.5,
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
    
    ## fit multisynth
    opts_weights <- c(opts_weights,
                      list(link="logit",
                           regularizer="nuc",
                           nlambda=20, lambda.min.ratio=1e-3,
                           opts=NULL))
    
    msynth <- multisynth_(X=wide$X, trt=wide$trt,
                          mask=wide$mask, gap=gap,
                          relative=relative, lambda=lambda, alpha=alpha,
                          link=opts_weights$link, regularizer=opts_weights$regularizer,
                          nlambda=opts_weights$nlambda,
                          lambda.min.ratio=opts_weights$lambda.min.ratio,
                          opts=opts_weights$opts)
    msynth$data <- wide
    msynth$data$time <- data %>% distinct(!!time) %>% pull(!!time)
    msynth$call <- call_name
    msynth$relative <- relative
    msynth$gap <- gap
    msynth$alpha <- alpha

    ## Get imbalance for uniform weights
    ## TODO: Get rid of this stupid hack of just fitting the weights again with zero steps
    unif <- multisynth_(X=wide$X, trt=wide$trt,
                          mask=wide$mask, gap=gap,
                          relative=relative, lambda=lambda, alpha=alpha,
                          link=opts_weights$link, regularizer=opts_weights$regularizer,
                          nlambda=opts_weights$nlambda,
                          lambda.min.ratio=opts_weights$lambda.min.ratio,
                          opts=list(max_it=0))
    
    ## Balance for aggregate estimate
    msynth$global_l2 <- sqrt(sum(msynth$imbalance[[1]][,1]^2))
    msynth$scaled_global_l2 <- msynth$global_l2 / sqrt(sum(unif$imbalance[[1]][,1]^2))

    ## balance for individual estimates
    msynth$ind_op <- svd(msynth$imbalance[[1]][,-1])$d[1]
    msynth$scaled_ind_op <- msynth$ind_op / svd(unif$imbalance[[1]][,-1])$d[1]
    
    
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
    gap <- tmp$gap
    d <- ncol(multisynth$data$X)
    fulldat <- cbind(multisynth$data$X, multisynth$data$y)
    ttot <- ncol(fulldat)

    grps <- unique(multisynth$data$trt) %>% sort()
    J <- length(grps) - 1
    n1 <- multisynth$data$trt[is.finite(multisynth$data$trt)] %>%
        table() %>% as.numeric()
    fullmask <- cbind(multisynth$data$mask, matrix(0, nrow=J, ncol=(ttot-d)))
    

    ## get weighted average estimates
    mu0hat <- lapply(multisynth$weights,
                     function(w) t(fulldat) %*% w)

    ## estimate the post-treatment values to gett att estimates
    mu1hat <- vapply(1:J,
                     function(j) colMeans(fulldat[multisynth$data$trt ==grps[j],
                                                , drop=FALSE]),
                     numeric(ttot))

    tauhat <- lapply(mu0hat, function(x) mu1hat - x)
    
    ## re-index time if relative to treatment
    if(relative) {
        total_len <- min(d + gap, ttot + d - grps[1]) ## total length of predictions
        mu0hat <- lapply(mu0hat,
                         function(x) {
                             vapply(1:J,
                                    function(j) {
                                        vec <- c(rep(NA, d-grps[j]),
                                          x[1:grps[j],j],
                                          x[(grps[j]+1):(min(grps[j] + gap, ttot)), j])
                                        c(vec, rep(NA, total_len - length(vec)))
                                    },
                                    numeric(total_len))
                         })

        tauhat <- lapply(tauhat,
                         function(x) {
                             vapply(1:J,
                                    function(j) {
                                        vec <- c(rep(NA, d-grps[j]),
                                                 x[1:grps[j],j],
                                                 x[(grps[j]+1):(min(grps[j] + gap, ttot)), j])
                                        c(vec, rep(NA, total_len - length(vec)))
                                    },
                                    numeric(total_len))
                         })
    ## get the overall average estimate
    lapply(mu0hat,
           function(x) {
               avg <- apply(x, 1, function(z) sum(n1 * z, na.rm=T) / sum(n1 * !is.na(z)))
               cbind(avg, x)
           }) -> mu0hat

    lapply(tauhat,
           function(x) {
               avg <- apply(x, 1, function(z) sum(n1 * z, na.rm=T) / sum(n1 * !is.na(z)))
               cbind(avg, x)
           }) -> tauhat
        
    } else {

        ## remove all estimates for t > T_j + gap
        lapply(mu0hat, function(x) {
            vapply(1:J,
                   function(j) c(x[1:min(grps[j]+gap, ttot),j],
                                 rep(NA, max(0, ttot-(grps[j] + gap)))),
                   numeric(ttot))
        }) -> mu0hat

        lapply(tauhat, function(x) {
            vapply(1:J,
                   function(j) c(x[1:min(grps[j]+gap, ttot),j],
                                 rep(NA, max(0, ttot-(grps[j] + gap)))),
                   numeric(ttot))
        }) -> tauhat

        
        ## only average currently treated units
        lapply(mu0hat, function(x) {
            avg1 <- rowSums(t(fullmask) *  x * n1) /
                rowSums(t(fullmask) *  n1)
            avg2 <- rowSums(t(1-fullmask) *  x * n1) /
                rowSums(t(1-fullmask) *  n1)
            avg <- replace_na(avg1, 0) * apply(fullmask, 2, min) +
                replace_na(avg2,0) * apply(1-fullmask, 2, max)
            cbind(avg, x)
        }) -> mu0hat

        ## only average currently treated units
        lapply(tauhat, function(x) {
            avg1 <- rowSums(t(fullmask) *  x * n1) /
                rowSums(t(fullmask) *  n1)
            avg2 <- rowSums(t(1-fullmask) *  x * n1) /
                rowSums(t(1-fullmask) *  n1)
            avg <- replace_na(avg1, 0) * apply(fullmask, 2, min) +
                replace_na(avg2,0) * apply(1-fullmask, 2, max)
            cbind(avg, x)
        }) -> tauhat        
    }
    

    if(att) {
        return(tauhat[[1]])
    } else {
        return(mu0hat[[1]])
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
#' @export
plot.multisynth <- function(multisynth) {
    plot(summary(multisynth))
}


#' Summary function for multisynth
#' @param relative Whether to estimate effects for time relative to treatment
#' @export
summary.multisynth <- function(multisynth, relative=NULL) {

    if(is.null(relative)) {
        relative <- multisynth$relative
    }

    grps <- unique(multisynth$data$trt) %>% sort()
    J <- length(grps) - 1    
    gap <- tmp$gap
    d <- ncol(multisynth$data$X)
    ttot <- d + ncol(multisynth$data$y)

    times <- multisynth$data$time
    
    summ <- list()
    ## post treatment estimate for each group and overall
    att <- predict(multisynth, relative, att=T)


    if(relative) {
        att <- data.frame(cbind(-d:min(gap-1, ttot-grps[1]-1), att))
        names(att) <- c("Time", "Average", times[grps[1:J]])
        att %>% gather(Level, Estimate, -Time) %>%
            rename("Time"=Time) -> att
    } else {
        att <- data.frame(cbind(times, att))
        names(att) <- c("Time", "Average", times[grps[1:J]])        
        att %>% gather(Level, Estimate, -Time) -> att
    }

    summ$att <- att


    summ$relative <- relative
    summ$grps <- grps
    summ$call <- multisynth$call
    summ$global_l2 <- multisynth$global_l2
    summ$scaled_global_l2 <- multisynth$scaled_global_l2

    summ$ind_op <- multisynth$ind_op
    summ$scaled_ind_op <- multisynth$scaled_ind_op
    ## get estimated bias

    
    class(summ) <- "summary.multisynth"
    return(summ)
}

#' Print function for summary function for multisynth
#' @param level Which treatment level to show summary for, default is average
#' @export
print.summary.multisynth <- function(summ, level=NULL) {
    ## straight from lm
    cat("\nCall:\n", paste(deparse(summ$call), sep="\n", collapse="\n"), "\n\n", sep="")

    if(is.null(level)) level <- "Average"
    
    ## get ATT estimates for treatment level, post treatment
    att_est <- summ$att %>%
        filter(Level == level)
    if(summ$relative) {
        att_est %>% filter(Time >= 0) -> att_est
    } else if(level == "Average") {
        att_est %>% filter(Time >= min(Level)) -> att_est
    } else {
        att_est %>% filter(Time >= level) -> att_est
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
#' @export
plot.summary.multisynth <- function(summ) {

    summ$att %>%
        ggplot2::ggplot(ggplot2::aes(x=Time, y=Estimate)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin=Estimate-2*Std.Error,
                        ymax=Estimate+2*Std.Error),
                    alpha=0.2) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=summ$t_int, lty=2) +
        ggplot2::geom_hline(yintercept=0, lty=2) + 
        ggplot2::theme_bw()
    
}
