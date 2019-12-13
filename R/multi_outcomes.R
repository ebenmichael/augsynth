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
#' @param ... optional arguments for outcome model
#' @param cov_agg Covariate aggregation functions, if NULL then use mean with NAs omitted
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
                                    combine_method = c("concat"),
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
#' @param wide Data formatted from format_data
#' @param synth_data Data formatted from foramt_synth
#' @param Z Matrix of auxiliary covariates
#' @param progfunc outcome model to use
#' @param scm Whether to fit SCM
#' @param fixed_eff Whether to de-mean synth
#' @param ... Extra args for outcome model
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
#' @param fixed_eff Whether to take out unit fixed effects or not
#' @param k Number of principal directions to keep, default all
#' 
#' @return \itemize{
#'          \item{"X"}{Matrix of combined pre-treatment outcomes}
#'          \item{"trt"}{Vector of treatment assignments}
#'          \item{"y"}{Matrix of combined post-treatment outcomes}
#'         }
combine_outcomes <- function(wide_list, combine_method, fixedeff, 
                             k= NULL, ...) {

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
#' @param augsynth augsynth_multiout object
#' @param att Whether to return the ATT or average outcome under control
#'
#' @return Vector of predicted post-treatment control averages
#' @export
predict.augsynth_multiout <- function(augsynth, att = F) {
    

    # X <- augsynth$data$X
    # y <- augsynth$data$y
    # Y0 <- t(augsynth$data$synth_data$Y0plot)
    # Y1 <- augsynth$data$synth_data$Y1plot

    # trt <- augsynth$data$trt
    # mhat <- augsynth$mhat
    
    # m1 <- colMeans(mhat[trt==1,,drop=F])

    # resid <- (Y0 - mhat[trt==0,drop=F])
    # y0 <- m1 + t(resid) %*% augsynth$weights
    # if(att) {
    #     pred <- Y1 - c(y0)
    # } else {
    #     pred <- y0
    # }
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
    n_outs <- length(augsynth$outcomes)
    pred_reshape <- matrix(NA, ncol = n_outs, 
                               nrow = length(augsynth$data$time))
    rownames(pred_reshape) <- augsynth$data$time
    colnames(pred_reshape) <- augsynth$outcomes
    # get outcome names for predictions
    pre_outs <- do.call(c, 
                        sapply(1:n_outs, 
                               function(j) {
                                   rep(augsynth$outcomes[j],
                                       ncol(augsynth$data_list$X[[j]]))
                               }, simplify = FALSE))

    post_outs <- do.call(c,
                         sapply(1:n_outs, 
                                function(j) {
                                    rep(augsynth$outcomes[j],
                                        ncol(augsynth$data_list$y[[j]]))
                               }, simplify = FALSE))
    pred_reshape[cbind(names(pred), c(pre_outs, post_outs))] <- pred

    return(pred_reshape)
}


#' Print function for augsynth
#' @export
print.augsynth_multiout <- function(augsynth) {
    ## straight from lm
    cat("\nCall:\n", paste(deparse(augsynth$call), sep="\n", collapse="\n"), "\n\n", sep="")

    ## print att estimates
    att <- predict(augsynth, att = T)
    att_post <- data.frame(
        colMeans(att[as.numeric(rownames(att)) >= augsynth$t_int,, drop = F]))
    names(att_post) <- c("")
    cat("Average ATT Estimate:\n")
    print(att_post)
    cat("\n\n")
}

#' Summary function for augsynth
#' @param se Whether to compute standard errors, default TRUE
#' @export
summary.augsynth_multiout <- function(augsynth, se = T) {

    summ <- list()

    if(se) {
        att_se <- jackknife_se_multiout(augsynth)
        t_final <- nrow(att_se$att)

        att_df <- data.frame(att_se$att[1:(t_final - 1),, drop=F])
        names(att_df) <- augsynth$outcomes
        att_df$Time <- augsynth$data$time
        att_df <- att_df %>% gather(Outcome, Estimate, -Time)

        se_df <- data.frame(att_se$se[1:(t_final - 1),, drop=F])
        names(se_df) <- augsynth$outcomes
        se_df$Time <- augsynth$data$time
        se_df <- se_df %>% gather(Outcome, Std.Error, -Time)

        att <- inner_join(att_df, se_df, by = c("Time", "Outcome"))


        att_avg <- data.frame(att_se$att[t_final,, drop = F])
        names(att_avg) <- augsynth$outcomes
        att_avg <- gather(att_avg, Outcome, Estimate)
        att_avg_se <- data.frame(att_se$se[t_final,, drop = F])
        names(att_avg_se) <- augsynth$outcomes
        att_avg_se <- gather(att_avg_se, Outcome, Std.Error)
        average_att <- inner_join(att_avg, att_avg_se, by="Outcome")

    } else {
        att_est <- predict(augsynth, att = T)
        att_df <- data.frame(att_est)
        names(att_df) <- augsynth$outcomes
        att_df$Time <- augsynth$data$time
        att <- att_df %>% gather(Outcome, Estimate, -Time)
        att$Std.Error <- NA
        att_avg <- colMeans(
            att[as.numeric(rownames(att)) >= augsynth$t_int,, drop = F])
        names(att_avg) <- augsynth$outcomes
        average_att <- gather(att_avg, Outcome, Estimate)
        average_att$Std.Error <- NA
    }

    summ$att <- att
    summ$average_att <- average_att
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
    
    class(summ) <- "summary.augsynth_multiout"
    return(summ)
}


#' Print function for summary function for augsynth
#' @export
print.summary.augsynth_multiout <- function(summ) {
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

    ## get pre-treatment fit by outcome
    imbal <- summ$att %>% 
        filter(Time < summ$t_int) %>%
        group_by(Outcome) %>%
        summarise(Pre.RMSE = sqrt(mean(Estimate ^ 2)))

    att_post <- summ$average_att$Estimate
    se_pool <- summ$average_att$Std.Error
    cat(paste("Overall L2 Imbalance (Scaled):",
              format(round(summ$l2_imbalance,3), nsmall=3), "  (",
              format(round(summ$scaled_l2_imbalance,3), nsmall=3), ")\n\n",
            #   "Avg Estimated Bias: ",
            #   format(round(mean(summ$bias_est), 3),nsmall=3), "\n\n",
              sep=""))
    cat("Average ATT Estimate:\n")
    print(inner_join(summ$average_att, imbal, by = "Outcome"))
    cat("\n\n")
}


#' Plot function for summary function for augsynth
#' @param se Whether to plot standard error
#' @export
plot.summary.augsynth_multiout <- function(summ, se = T) {

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
        ggplot2::facet_wrap(~ Outcome, scales = "free_y") +
        ggplot2::theme_bw()
    
}