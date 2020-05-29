#' Fit Augmented SCM with multiple outcomes
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
augsynth_multiout <- function(form, unit, time, t_int, data,
                              progfunc=c("Ridge", "None", "EN", "RF",
                                          "GSYN", "MCP",
                                          "CITS", "CausalImpact",
                                          "seq2seq"),
                              scm=T,
                              fixedeff = FALSE,
                              cov_agg=NULL,
                              ...) {
    call_name <- match.call()

    form <- Formula::Formula(form)
    unit <- enquo(unit)
    time <- enquo(time)
    
    ## format data
    outcome <- terms(formula(form, rhs=1))[[2]]
    trt <- terms(formula(form, rhs=1))[[3]]

    outcomes_str <- all.vars(outcome)
    outcomes <- sapply(outcomes_str, quo)
    # get outcomes as a list
    wide_list <- format_data_multi(outcomes, trt, unit, time, t_int, data)
    

    

    ## add covariates
    if(length(form)[2] == 2) {
        Z <- extract_covariates(form, unit, time, t_int, data, cov_agg)
    } else {
        Z <- NULL
    }

    # combine method is just concatenation
    # TODO: Include other options
    combine_method <- "concat"

    # fit augmented SCM
    augsynth <- fit_augsynth_multiout_internal(wide_list, combine_method, Z,
                                               progfunc, scm,
                                               fixedeff, ...)

    # add some extra data
    augsynth$data_list <- wide_list
    augsynth$data$time <- data %>% distinct(!!time) %>% pull(!!time)
    augsynth$call <- call_name
    augsynth$t_int <- t_int 
    augsynth$outcomes <- outcomes_str
    augsynth$combine_method <- combine_method
    return(augsynth)
}


#' Internal function to fit augmented SCM with multiple outcomes
#' @param wide_list List of matrices for each outcome formatted from format_data
#' @param combine_method How to combine outcomes
#' @param Z Matrix of auxiliary covariates
#' @param progfunc outcome model to use
#' @param scm Whether to fit SCM
#' @param fixedeff Whether to de-mean synth
#' @param ... Extra args for outcome model
#' @noRd
fit_augsynth_multiout_internal <- function(wide_list, combine_method, Z,
                                           progfunc, scm, fixedeff, ...) {


    # combine into a matrix for fitting and balancing
    out <- combine_outcomes(wide_list, combine_method, fixedeff, ...)
    wide_bal <- out$wide_bal
    mhat <- out$mhat
    V <- out$V
    synth_data <- do.call(format_synth, wide_bal)

    # set Y1 and Y0plot to be raw concatenated outcomes
    X <- do.call(cbind, wide_list$X)
    y <- do.call(cbind, wide_list$y)
    trt <- wide_list$trt
    synth_data$Y0plot <- t(cbind(X, y)[trt == 0,, drop = F])
    synth_data$Y1plot <- colMeans(cbind(X, y)[trt == 1,, drop = F])


    augsynth <- fit_augsynth_internal(wide_bal, synth_data, Z, progfunc, 
                                      scm, fixedeff, V = V, ...)

    # potentially add back in fixed effects
    augsynth$mhat <- mhat + augsynth$mhat

    augsynth$data = list(X = X, trt = trt, y = y)
    ##format output
    class(augsynth) <- c("augsynth_multiout", "augsynth")
    return(augsynth)
}

#' Helper function to combine multiple outcomes into a single balance matrix
#' @param wide_list List of lists of pre/post treatment data for each outcome
#' @param combine_method How to combine outcomes
#' @param fixedeff Whether to take out unit fixed effects or not
#' @param k Number of principal directions to keep, default all
#' @param ... Extra argumemnts for combination
#' @noRd
#' @return \itemize{
#'          \item{"X"}{Matrix of combined pre-treatment outcomes}
#'          \item{"trt"}{Vector of treatment assignments}
#'          \item{"y"}{Matrix of combined post-treatment outcomes}
#'         }
combine_outcomes <- function(wide_list, combine_method, fixedeff,
                             k= NULL, ...) {
    print(class(wide_list$X))
    print(is.list(wide_list$X))
    print(class(wide_list$X$lngdpcapita))
    n_outs <- length(wide_list$X)
    total_pre <- Map(ncol, wide_list$X) %>% Reduce(`+`, .)
    total_post <- Map(ncol, wide_list$y) %>% Reduce(`+`, .)
    total_dim <- total_pre + total_post
    n_units <- Map(nrow, wide_list$X) %>% Reduce(max, .)
    # take out unit fixed effects
    demean_j <- function(j) {
        means <- rowMeans(wide_list$X[[j]])

        new_wide_data <- list()
        new_X <- wide_list$X[[j]] - means
        new_y <- wide_list$y[[j]] - means

        new_wide_data$X <- new_X
        new_wide_data$y <- new_y
        new_wide_data$mhat_pre <- replicate(ncol(wide_list$X[[j]]),
                                            means)
        new_wide_data$mhat_post <- replicate(ncol(wide_list$y[[j]]),
                                            means)
        return(new_wide_data)
    }
    if(fixedeff) {
        new_wide_list <- lapply(1:n_outs, demean_j)
        wide_list$X <- lapply(new_wide_list, function(x) x$X)
        wide_list$y <- lapply(new_wide_list, function(x) x$y)
        mhat_pre <- do.call(cbind, lapply(new_wide_list, function(x) x$mhat_pre))
        mhat_post <- do.call(cbind, lapply(new_wide_list, function(x) x$mhat_post))
        mhat <- cbind(mhat_pre, mhat_post)
    } else {
        mhat <- matrix(0, nrow = n_units, ncol = total_dim)
    }

    # combine outcomes
    if(combine_method == "concat") {
        # center X and scale by overall variance for outcome
        # X <- lapply(wide_list$X, function(x) t(t(x) - colMeans(x)) / sd(x))
        wide_bal <- list(X = do.call(cbind, wide_list$X),
                         y = do.call(cbind, wide_list$y),
                         trt = wide_list$trt)

        # V matrix scales by inverse variance for outcome and number of periods
        V <- do.call(c, 
            lapply(wide_list$X, 
                function(x) rep(1 / (sqrt(ncol(x)) * 
                        sd(x[wide_list$trt == 0, , drop = F], na.rm=T)), 
                        ncol(x))))
    } else if(combine_method == "svd") {
        wide_bal <- list(X = do.call(cbind, wide_list$X),
                         y = do.call(cbind, wide_list$y),
                         trt = wide_list$trt)

        # first get the standard deviations of the outcomes to put on the same scale
        sds <- do.call(c, 
            lapply(wide_list$X, 
                function(x) rep((sqrt(ncol(x)) * sd(x, na.rm=T)), ncol(x))))

        # do an SVD on centered and scaled outcomes
        X0 <- wide_bal$X[wide_bal$trt == 0, , drop = FALSE]
        X0 <- t((t(X0) - colMeans(X0)) / sds)
        k <- if(is.null(k)) ncol(X0) else k
        V <- diag(1 / sds) %*% svd(X0)$v[, 1:k, drop = FALSE]

    }else {
        stop(paste("combine_method should be one of ('concat'),", 
            combine_method, " is not a valid combining option"))
    }

    return(list(wide_bal = wide_bal, mhat = mhat, V = V))

}

#' Get prediction of ATT or average outcome under control
#' @param object augsynth_multiout object
#' @param ... Optional arguments, including \itemize{\item{"att"}{Whether to return the ATT or average outcome under control}}
#'
#' @return Vector of predicted post-treatment control averages
#' @export
predict.augsynth_multiout <- function(object, ...) {
    if ("att" %in% names(list(...))) {
        att <- list(...)$att
    } else {
        att <- F
    }

    # call augsynth predict
    pred <- NextMethod()

    if(att) {
        pred_names <- names(pred)
    } else {
        pred_names <- rownames(pred)
    }
    pred <- c(pred)
    names(pred) <- pred_names

    # separate out by outcome
    n_outs <- length(object$outcomes)
    pred_reshape <- matrix(NA, ncol = n_outs, 
                               nrow = length(object$data$time))
    rownames(pred_reshape) <- object$data$time
    colnames(pred_reshape) <- object$outcomes
    # get outcome names for predictions
    pre_outs <- do.call(c, 
                        sapply(1:n_outs, 
                               function(j) {
                                   rep(object$outcomes[j],
                                       ncol(object$data_list$X[[j]]))
                               }, simplify = FALSE))

    post_outs <- do.call(c,
                         sapply(1:n_outs, 
                                function(j) {
                                    rep(object$outcomes[j],
                                        ncol(object$data_list$y[[j]]))
                               }, simplify = FALSE))
    pred_reshape[cbind(names(pred), c(pre_outs, post_outs))] <- pred

    return(pred_reshape)
}


#' Print function for augsynth
#' @param x augsynth_multiout object
#' @param ... Optional arguments
#' @export
print.augsynth_multiout <- function(x, ...) {
    ## straight from lm
    cat("\nCall:\n", paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")

    ## print att estimates
    att <- predict(x, att = T)
    att_post <- data.frame(
        colMeans(att[as.numeric(rownames(att)) >= x$t_int,, drop = F]))
    names(att_post) <- c("")
    cat("Average ATT Estimate:\n")
    print(att_post)
    cat("\n\n")
}

#' Summary function for augsynth
#' @param object augsynth_multiout object
#' @param ... Optional arguments, including \itemize{\item{"se"}{Whether to plot standard error}}
#' @export
summary.augsynth_multiout <- function(object, ...) {
    if ("se" %in% names(list(...))) {
        se <- list(...)$se
    } else {
        se <- T
    }
    summ <- list()

    if(se) {
        att_se <- jackknife_se_multiout(object)
        t_final <- nrow(att_se$att)

        att_df <- data.frame(att_se$att[1:(t_final - 1),, drop=F])
        names(att_df) <- object$outcomes
        att_df$Time <- object$data$time
        att_df <- att_df %>% gather(Outcome, Estimate, -Time)

        se_df <- data.frame(att_se$se[1:(t_final - 1),, drop=F])
        names(se_df) <- object$outcomes
        se_df$Time <- object$data$time
        se_df <- se_df %>% gather(Outcome, Std.Error, -Time)

        att <- inner_join(att_df, se_df, by = c("Time", "Outcome"))


        att_avg <- data.frame(att_se$att[t_final,, drop = F])
        names(att_avg) <- object$outcomes
        att_avg <- gather(att_avg, Outcome, Estimate)
        att_avg_se <- data.frame(att_se$se[t_final,, drop = F])
        names(att_avg_se) <- object$outcomes
        att_avg_se <- gather(att_avg_se, Outcome, Std.Error)
        average_att <- inner_join(att_avg, att_avg_se, by="Outcome")

    } else {
        att_est <- predict(object, att = T)
        att_df <- data.frame(att_est)
        names(att_df) <- object$outcomes
        att_df$Time <- object$data$time
        att <- att_df %>% gather(Outcome, Estimate, -Time)
        att$Std.Error <- NA
        att_avg <- colMeans(
            att[as.numeric(rownames(att)) >= object$t_int,, drop = F])
        names(att_avg) <- object$outcomes
        average_att <- gather(att_avg, Outcome, Estimate)
        average_att$Std.Error <- NA
    }

    summ$att <- att
    summ$average_att <- average_att
    summ$t_int <- object$t_int
    summ$call <- object$call
    summ$l2_imbalance <- object$l2_imbalance
    summ$scaled_l2_imbalance <- object$scaled_l2_imbalance
    ## get estimated bias

    if(object$progfunc == "Ridge") {
        mhat <- object$ridge_mhat
        w <- object$synw
    } else {
        mhat <- object$mhat
        w <- object$weights
    }
    trt <- object$data$trt
    m1 <- colMeans(mhat[trt==1,,drop=F])

    summ$bias_est <- m1 - t(mhat[trt==0,,drop=F]) %*% w
    if(object$progfunc == "None" | (!object$scm)) {
        summ$bias_est <- NA
    }
    
    class(summ) <- "summary.augsynth_multiout"
    return(summ)
}


#' Print function for summary function for augsynth
#' @param x summary.augsynth_multiout object
#' @param ... Optional arguments
#' @export
print.summary.augsynth_multiout <- function(x, ...) {
    ## straight from lm
    cat("\nCall:\n", paste(deparse(x$call), sep="\n", collapse="\n"), "\n\n", sep="")
    
    att_est <- x$att$Estimate
    ## get pre-treatment fit by outcome
    imbal <- x$att %>% 
        filter(Time < x$t_int) %>%
        group_by(Outcome) %>%
        summarise(Pre.RMSE = sqrt(mean(Estimate ^ 2)))

    cat(paste("Overall L2 Imbalance (Scaled):",
              format(round(x$l2_imbalance,3), nsmall=3), "  (",
              format(round(x$scaled_l2_imbalance,3), nsmall=3), ")\n\n",
            #   "Avg Estimated Bias: ",
            #   format(round(mean(summ$bias_est), 3),nsmall=3), "\n\n",
              sep=""))
    cat("Average ATT Estimate:\n")
    print(inner_join(x$average_att, imbal, by = "Outcome"))
    cat("\n\n")
}


#' Plot function for summary function for augsynth
#' @param x summary.augsynth_multiout object
#' @param ... Optional arguments, including \itemize{\item{"se"}{Whether to plot standard error}}
#' 
#' @export
plot.summary.augsynth_multiout <- function(x, ...) {
    if ("se" %in% names(list(...))) {
        se <- list(...)$se
    } else {
        se <- T
    }
    p <- x$att %>%
        ggplot2::ggplot(ggplot2::aes(x=Time, y=Estimate))
    if(se) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=Estimate-2*Std.Error,
                        ymax=Estimate+2*Std.Error),
                    alpha=0.2)
    }
    p + ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=x$t_int, lty=2) +
        ggplot2::geom_hline(yintercept=0, lty=2) +
        ggplot2::facet_wrap(~ Outcome, scales = "free_y") +
        ggplot2::theme_bw()
    
}