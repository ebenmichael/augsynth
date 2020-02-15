################################################################################
## Main functions for single-period treatment augmented synthetic controls Method
################################################################################


#' Fit Augmented SCM
#' @importFrom stats terms
#' @importFrom stats formula
#' @importFrom stats update 
#' @importFrom stats delete.response 
#' @importFrom stats model.matrix 
#' @importFrom stats model.frame 
#' 
#' @param form outcome ~ treatment | auxillary covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param t_int Time of intervention
#' @param data Panel data as dataframe
#' @param progfunc What function to use to impute control outcomes
#'                 Ridge=Ridge regression (allows for standard errors),
#'                 None=No outcome model,
#'                 EN=Elastic Net, RF=Random Forest, GSYN=gSynth,
#'                 MCP=MCPanel, CITS=CITS
#'                 CausalImpact=Bayesian structural time series with CausalImpact
#'                 seq2seq=Sequence to sequence learning with feedforward nets
#' @param scm Whether the SCM weighting function is used
#' @param fixedeff Whether to include a unit fixed effect, default F 
#' @param cov_agg Covariate aggregation functions, if NULL then use mean with NAs omitted
#' @param ... optional arguments for outcome model
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
single_augsynth <- function(form, unit, time, t_int, data,
                     progfunc=c("Ridge", "None", "EN", "RF", "GSYN", "MCP",
                                "CITS", "CausalImpact", "seq2seq"),
                     scm=T,
                     fixedeff = FALSE,
                     cov_agg=NULL, ...) {
    call_name <- match.call()

    form <- Formula::Formula(form)
    unit <- enquo(unit)
    time <- enquo(time)

    ## format data
    outcome <- terms(formula(form, rhs=1))[[2]]
    trt <- terms(formula(form, rhs=1))[[3]]

    wide <- format_data(outcome, trt, unit, time, t_int, data)

    synth_data <- do.call(format_synth, wide)

    

    ## add covariates
    if(length(form)[2] == 2) {

        ## if no aggregation functions, use the mean (omitting NAs)
        cov_agg <- c(function(x) mean(x, na.rm=T))

        cov_form <- update(formula(delete.response(terms(form, rhs=2, data=data))),
                           ~. - 1) ## ensure that there is no intercept

        ## pull out relevant covariates and aggregate
        pre_data <- data %>% 
            filter(!! (time) < t_int)
        model.matrix(cov_form,
                     model.frame(cov_form, pre_data,
                                 na.action=NULL) ) %>%
            data.frame() %>%
            mutate(unit=pull(pre_data, !!unit)) %>%
            group_by(unit) %>%
            summarise_all(cov_agg) %>%
            select(-unit) %>%
            as.matrix() -> Z
    } else {
        Z <- NULL
    }

    # fit augmented SCM
    augsynth <- fit_augsynth_internal(wide, synth_data, Z, progfunc, 
                                      scm, fixedeff, ...)

    # add some extra data
    augsynth$data$time <- data %>% distinct(!!time) %>% pull(!!time)
    augsynth$call <- call_name
    augsynth$t_int <- t_int 

    return(augsynth)
}

#' Internal function to fit augmented SCM
#' @param wide Data formatted from format_data
#' @param synth_data Data formatted from foramt_synth
#' @param Z Matrix of auxiliary covariates
#' @param progfunc outcome model to use
#' @param scm Whether to fit SCM
#' @param fixedeff Whether to de-mean synth
#' @param ... Extra args for outcome model
fit_augsynth_internal <- function(wide, synth_data, Z, progfunc,
                                  scm, fixedeff, ...) {

    n <- nrow(wide$X)
    t0 <- ncol(wide$X)
    ttot <- t0 + ncol(wide$y)
    if(fixedeff) {
        demeaned <- demean_data(wide, synth_data)
        fit_wide <- demeaned$wide
        fit_synth_data <- demeaned$synth_data
        mhat <- demeaned$mhat
    } else {
        fit_wide <- wide
        fit_synth_data <- synth_data
        mhat <- matrix(0, n, ttot)
    }
    ## fit augsynth
    if(progfunc == "Ridge") {
        if(scm) {
            ## Ridge ASCM
            augsynth <- do.call(fit_ridgeaug_formatted,
                                list(wide_data = fit_wide, 
                                   synth_data = fit_synth_data, 
                                   Z = Z, ...))
        } else {
            ## Just ridge regression
            augsynth <- do.call(fit_ridgeaug_formatted, list(wide_data = fit_wide, 
                                   synth_data = fit_synth_data,
                                   Z = Z, ridge = T, scm = F, ...))
        }
    } else if(progfunc == "None") {
        ## Just SCM
        augsynth <- do.call(fit_ridgeaug_formatted,
                        c(list(wide_data = fit_wide, 
                               synth_data = fit_synth_data,
                               Z = Z, ridge = F, scm = T)))
    } else {
        ## Other outcome models
        augsynth <- fit_augsyn(fit_wide, fit_synth_data, 
                               progfunc, scm, ...)
    }

    augsynth$mhat <- mhat + cbind(matrix(0, nrow = n, ncol = t0), 
                                  augsynth$mhat)
    augsynth$data <- wide
    augsynth$data$Z <- Z
    augsynth$progfunc <- progfunc
    augsynth$scm <- scm
    augsynth$fixedeff <- fixedeff
    augsynth$extra_args <- list(...)
    ##format output
    class(augsynth) <- "augsynth"
    return(augsynth)
}

#' Get prediction of ATT or average outcome under control
#' @param object augsynth object
#' @param ... Optional arguments
#'
#' @return Vector of predicted post-treatment control averages
#' @export
predict.augsynth <- function(object, ...) {
    if ("att" %in% names(list(...))) {
        att <- list(...)$att
    } else {
        att <- F
    }
    augsynth <- object
    
    X <- augsynth$data$X
    y <- augsynth$data$y
    comb <- cbind(X, y)
    trt <- augsynth$data$trt
    mhat <- augsynth$mhat
    
    m1 <- colMeans(mhat[trt==1,,drop=F])

    resid <- (comb[trt==0,,drop=F] - mhat[trt==0,drop=F])

    y0 <- m1 + t(resid) %*% augsynth$weights
    if(att) {
        return(colMeans(comb[trt == 1,, drop = F]) - c(y0))
    } else {
        return(y0)
    }
}


#' Print function for augsynth
#' @param x augsynth object
#' @param ... Optional arguments
#' @export
print.augsynth <- function(x, ...) {
    augsynth <- x
    
    ## straight from lm
    cat("\nCall:\n", paste(deparse(augsynth$call), sep="\n", collapse="\n"), "\n\n", sep="")

    ## print att estimates
    tint <- ncol(augsynth$data$X)
    ttotal <- tint + ncol(augsynth$data$y)
    att_post <- predict(augsynth, att = T)[(tint + 1):ttotal]

    cat(paste("Average ATT Estimate: ",
              format(round(mean(att_post),3), nsmall = 3), "\n\n", sep=""))
}


#' Plot function for augsynth
#' @importFrom graphics plot
#' 
#' @param x Augsynth object to be plotted
#' @param ... Optional arguments
#' @export
plot.augsynth <- function(x, ...) {
    if ("se" %in% names(list(...))) {
        se <- list(...)$se
    } else {
        se <- T
    }
    augsynth <- x
    plot(summary(augsynth, se), se = se)
}


#' Summary function for augsynth
#' @param object augsynth object
#' @param ... Optional arguments
#' @export
summary.augsynth <- function(object, ...) {
    augsynth <- object
    if ("se" %in% names(list(...))) {
        se <- list(...)$se
    } else {
        se <- T
    }
    
    
    summ <- list()

    t0 <- ncol(augsynth$data$X)
    t_final <- t0 + ncol(augsynth$data$y)

    if(se) {
        att_se <- jackknife_se_single(augsynth)
        att <- data.frame(Time = augsynth$data$time,
                          Estimate = att_se$att[1:t_final],
                          Std.Error = att_se$se[1:t_final])
        att_avg <- att_se$att[t_final + 1]
        att_avg_se <- att_se$se[t_final + 1]

    } else {
        t0 <- ncol(augsynth$data$X)
        t_final <- t0 + ncol(augsynth$data$y)
        att_est <- predict(augsynth, att = T)
        att <- data.frame(Time = augsynth$data$time,
                          Estimate = att_est)
        att$Std.Error <- NA
        att_avg <- mean(att_est[(t0 + 1):t_final])
        att_avg_se <- NA
    }

    summ$att <- att
    summ$average_att <- data.frame(Estimate = att_avg, Std.Error = att_avg_se)
    summ$t_int <- augsynth$t_int
    summ$call <- augsynth$call
    summ$l2_imbalance <- augsynth$l2_imbalance
    summ$scaled_l2_imbalance <- augsynth$scaled_l2_imbalance
    ## get estimated bias

    if(augsynth$progfunc == "Ridge") {
        mhat <- augsynth$ridge_mhat
        w <- augsynth$synw
    } else {
        mhat <- augsynth$mhat
        w <- augsynth$weights
    }
    trt <- augsynth$data$trt
    m1 <- colMeans(mhat[trt==1,,drop=F])

    summ$bias_est <- m1 - t(mhat[trt==0,,drop=F]) %*% w
    if(augsynth$progfunc == "None" | (!augsynth$scm)) {
        summ$bias_est <- NA
    }
    
    class(summ) <- "summary.augsynth"
    return(summ)
}

#' Print function for summary function for augsynth
#' @param x summary object
#' @param ... Optional arguments
#' @export
print.summary.augsynth <- function(x, ...) {
    summ <- x
    
    ## straight from lm
    cat("\nCall:\n", paste(deparse(summ$call), sep="\n", collapse="\n"), "\n\n", sep="")

    t_final <- nrow(summ$att)

    ## distinction between pre and post treatment
    att_est <- summ$att$Estimate
    t_total <- length(att_est)
    t_int <- summ$att %>% filter(Time <= summ$t_int) %>% nrow()
    
    att_pre <- att_est[1:(t_int-1)]
    att_post <- att_est[t_int:t_total]

    ## pool the standard error estimates to summarise it
    se_est <- summ$att$Std.Error

    se_pool <- sqrt(mean(se_est[t_int:t_total]^2))

    att_post <- summ$average_att$Estimate
    se_pool <- summ$average_att$Std.Error
    
    cat(paste("Average ATT Estimate (Std. Error): ",
              format(round(att_post,3), nsmall=3), "  (",
              format(round(se_pool,3)), ")\n",
              "L2 Imbalance (Scaled): ",
              format(round(summ$l2_imbalance,3), nsmall=3), "  (",
              format(round(summ$scaled_l2_imbalance,3), nsmall=3), ")\t",
              "Avg Estimated Bias: ",
              format(round(mean(summ$bias_est), 3),nsmall=3), "\n\n",
              sep=""))


    print(summ$att[t_int:t_final,], row.names=F)
    
}

#' Plot function for summary function for augsynth
#' @param x Summary object
#' @param ... Optional arguments
#' @export
plot.summary.augsynth <- function(x, ...) {
    summ <- x
    if ("se" %in% names(list(...))) {
        se <- list(...)$se
    } else {
        se <- T
    }
    
    p <- summ$att %>%
        ggplot2::ggplot(ggplot2::aes(x=Time, y=Estimate))
    if(se) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=Estimate-2*Std.Error,
                        ymax=Estimate+2*Std.Error),
                    alpha=0.2)
    }
    p + ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=summ$t_int, lty=2) +
        ggplot2::geom_hline(yintercept=0, lty=2) + 
        ggplot2::theme_bw()
    
}



#' augsynth: A package implementing the Augmented Synthetic Controls Method
#' @docType package
#' @name augsynth-package
#' @importFrom magrittr "%>%"
#' @importFrom purrr reduce
#' @import dplyr
#' @import LowRankQP
#' @import tidyr
NULL
