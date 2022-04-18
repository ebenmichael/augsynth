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
#' @param combine_method How to combine outcomes: `concat` concatenates outcomes and `avg` averages them, default: 'avg'
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
                              progfunc=c("Ridge", "None"),
                              scm=T,
                              fixedeff = FALSE,
                              cov_agg=NULL,
                              combine_method = "avg",
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
        cov_form <- paste(deparse(terms(formula(form, rhs = 2))[[3]]), collapse = "")
        new_form <- as.formula(paste("~", cov_form))
        Z <- extract_covariates(new_form, unit, time, t_int, data, cov_agg)
    } else {
        Z <- NULL
    }

    # only allow ridge augmentation
    if(! tolower(progfunc) %in% c("none", "ridge")) {
      stop(paste(progfunc, "is not a valid augmentation function with multiple outcomes. Only `none` or `ridge` are allowable options for `prog_func`"))
    }

    # fit augmented SCM
    augsynth <- fit_augsynth_multiout_internal(wide_list, combine_method, Z,
                                               progfunc, scm,
                                               fixedeff, outcomes_str, ...)

    # add some extra data
    augsynth$data$time <- data %>% distinct(!!time) %>% pull(!!time)
    augsynth$call <- call_name
    augsynth$t_int <- t_int 
    augsynth$combine_method <- combine_method

    treated_units <- data %>% filter(!!trt == 1) %>% distinct(!!unit) %>% pull(!!unit)
    control_units <- data %>% filter(!(!!unit %in% treated_units)) %>% 
                        distinct(!!unit) %>% pull(!!unit)
    augsynth$weights <- matrix(augsynth$weights)
    rownames(augsynth$weights) <- control_units

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
                                           progfunc, scm, fixedeff, 
                                           outcomes_str, ...) {


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
    augsynth$mhat <- mhat# + augsynth$mhat

    augsynth$data = list(X = X, trt = trt, y = y, Z = Z)
    augsynth$data_list <- wide_list
    augsynth$outcomes <- outcomes_str
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
        mhat_pre <- lapply(new_wide_list, function(x) x$mhat_pre)
        mhat_post <- lapply(new_wide_list, function(x) x$mhat_post)
    } else {
        mhat_pre <- lapply(
          1:n_outs,
          function(j) matrix(0, nrow = n_units, ncol = ncol(wide_list$X[[j]])))
        mhat_post <- lapply(
          1:n_outs,
          function(j) matrix(0, nrow = n_units, ncol = ncol(wide_list$y[[j]])))
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
    # } else if(combine_method == "svd") {
    #     wide_bal <- list(X = do.call(cbind, wide_list$X),
    #                      y = do.call(cbind, wide_list$y),
    #                      trt = wide_list$trt)

    #     # first get the standard deviations of the outcomes to put on the same scale
    #     sds <- do.call(c, 
    #         lapply(wide_list$X, 
    #             function(x) rep((sqrt(ncol(x)) * sd(x, na.rm=T)), ncol(x))))

    #     # do an SVD on centered and scaled outcomes
    #     X0 <- wide_bal$X[wide_bal$trt == 0, , drop = FALSE]
    #     X0 <- t((t(X0) - colMeans(X0)) / sds)
    #     k <- if(is.null(k)) ncol(X0) else k
    #     V <- diag(1 / sds) %*% svd(X0)$v[, 1:k, drop = FALSE]
    } else if(combine_method == "avg") {
        # average pre-treatment outcomes, dividing by standard deviation
      wide_bal <- list(X = Reduce(`+`, lapply(wide_list$X, function(x) x / sd(x[wide_list$trt == 0,]))),
                       #X = Reduce(`+`, wide_list$X),
                       y = Reduce(`+`, wide_list$y),
                       trt = wide_list$trt)

      V <- diag(ncol(wide_bal$X))
      # mhat_pre <- Reduce(`+`, mhat_pre)
      # mhat_post <- Reduce(`+`, mhat_post)
      # mhat <- cbind(mhat_pre, mhat_post)



    } else {
        stop(paste("combine_method should be one of ('avg', 'concat'),", 
            combine_method, " is not a valid combining option"))
    }

  mhat_pre <- do.call(cbind, mhat_pre)
  mhat_post <- do.call(cbind, mhat_post)
  mhat <- cbind(mhat_pre, mhat_post)

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

    # separate out by outcome
    n_outs <- length(object$data_list$X)
    max_t <- max(sapply(1:n_outs, 
      function(k) ncol(object$data_list$X[[k]]) + ncol(object$data_list$y[[k]])))
    pred_reshape <- matrix(NA, ncol = n_outs, 
                               nrow = max_t)
    colnames <- lapply(1:n_outs, 
      function(k) colnames(cbind(object$data_list$X[[k]], 
                                 object$data_list$y[[k]])))
    rownames(pred_reshape) <- colnames[[which.max(sapply(colnames, length))]]
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
    # print(pred)
    # print(cbind(names(pred), c(pre_outs, post_outs)))
    
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
summary.augsynth_multiout <- function(object, inf = T, inf_type = "jackknife", ...) {
    

    summ <- list()

    if(inf) {
        if(inf_type == "jackknife") {
            att_se <- jackknife_se_multiout(object)
        } else if(inf_type == "jackknife+") {
          att_se <- time_jackknife_plus_multiout(object, ...)
        } else if(inf_type == "conformal") {
          att_se <- conformal_inf_multiout(object, ...)
        } else {
            stop(paste(inf_type, "is not a valid choice of 'inf_type'"))
        }

        t_final <- nrow(att_se$att)

        att_df <- data.frame(att_se$att[1:(t_final - 1),, drop=F])
        names(att_df) <- object$outcomes
        att_df$Time <- object$data$time
        att_df <- att_df %>% gather(Outcome, Estimate, -Time)

        if(inf_type == "jackknife") {
          se_df <- data.frame(att_se$se[1:(t_final - 1),, drop=F])
          names(se_df) <- object$outcomes
          se_df$Time <- object$data$time
          se_df <- se_df %>% gather(Outcome, Std.Error, -Time)

          att <- inner_join(att_df, se_df, by = c("Time", "Outcome"))
        } else if(inf_type %in% c("conformal", "jackknife+")) {
          
          lb_df <- data.frame(att_se$lb[1:(t_final - 1),, drop=F])
          names(lb_df) <- object$outcomes
          lb_df$Time <- object$data$time
          lb_df <- lb_df %>% gather(Outcome, lower_bound, -Time)

          ub_df <- data.frame(att_se$ub[1:(t_final - 1),, drop=F])
          names(ub_df) <- object$outcomes
          ub_df$Time <- object$data$time
          ub_df <- ub_df %>% gather(Outcome, upper_bound, -Time)

          att <- inner_join(att_df, lb_df, by = c("Time", "Outcome")) %>%
              inner_join(ub_df, by = c("Time", "Outcome")) 
          if(inf_type == "conformal") {

            pval_df <- data.frame(att_se$p_val[1:(t_final - 1),, drop=F])
            names(pval_df) <- object$outcomes
            pval_df$Time <- object$data$time
            pval_df <- pval_df %>% gather(Outcome, p_val, -Time)
            att <- inner_join(att, pval_df, by = c("Time", "Outcome")) 
          }
        }

        att_avg <- data.frame(att_se$att[t_final,, drop = F])
        names(att_avg) <- object$outcomes
        att_avg <- gather(att_avg, Outcome, Estimate)

        if(inf_type == "jackknife") {
          att_avg_se <- data.frame(att_se$se[t_final,, drop = F])
          names(att_avg_se) <- object$outcomes
          att_avg_se <- gather(att_avg_se, Outcome, Std.Error)
          average_att <- inner_join(att_avg, att_avg_se, by="Outcome")
        } else if(inf_type %in% c("conformal", "jackknife+")){
          att_avg_lb <- data.frame(att_se$lb[t_final,, drop = F])
          names(att_avg_lb) <- object$outcomes
          att_avg_lb <- gather(att_avg_lb, Outcome, lower_bound)

          att_avg_ub <- data.frame(att_se$ub[t_final,, drop = F])
          names(att_avg_ub) <- object$outcomes
          att_avg_ub <- gather(att_avg_ub, Outcome, upper_bound)
          

          average_att <- inner_join(att_avg, att_avg_lb, by="Outcome") %>%
              inner_join(att_avg_ub, by = "Outcome")
          
          if(inf_type == "conformal") {
            att_avg_pval <- data.frame(att_se$p_val[t_final,, drop = F])
            names(att_avg_pval) <- object$outcomes
            att_avg_pval <- gather(att_avg_pval, Outcome, p_val)

             average_att <- inner_join(average_att, att_avg_pval, by = "Outcome")
          }
        } else {
          average_att <- gather(att_avg, Outcome, Estimate)
        }
        

    } else {
        att_est <- predict(object, att = T)
        att_df <- data.frame(att_est)
        names(att_df) <- object$outcomes
        att_df$Time <- object$data$time
        att <- att_df %>% gather(Outcome, Estimate, -Time)
        att$Std.Error <- NA
        t_int <- min(sapply(object$data_list$X, ncol))
        att_avg <- data.frame(colMeans(
            att_est[as.numeric(rownames(att)) >= t_int,, drop = F]))
        names(att_avg) <- object$outcomes
        average_att <- gather(att_avg, Outcome, Estimate)
        average_att$Std.Error <- NA
    }

      # get average of all outcomes
    sds <- data.frame(Outcome = object$outcomes,
                      sdo = sapply(object$data_list$X, sd))
    att %>%
      inner_join(sds, by = "Outcome") %>%
      mutate(Estimate = Estimate / sdo) %>%
      group_by(Time) %>%
      summarise(Estimate = mean(Estimate)) %>%
      mutate(Outcome = "Average") %>%
      bind_rows(att, .) -> att

    summ$att <- att
    summ$average_att <- average_att
    summ$t_int <- object$t_int
    summ$call <- object$call
    summ$l2_imbalance <- object$l2_imbalance
    summ$scaled_l2_imbalance <- object$scaled_l2_imbalance
    summ$inf_type <- inf_type
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
#' @importFrom graphics plot
#' @param x summary.augsynth_multiout object
#' @param inf Boolean, whether to plot uncertainty intervals, default TRUE
#' @param plt_avg Boolean, whether to plot the average of the outcomes, default FALSE
#' @param ... Optional arguments for summary function
#' 
#' @export
plot.augsynth_multiout  <- function(x, inf = T, plt_avg = F, ...) {
  plot(summary(x, ...), inf =  inf, plt_avg = plt_avg)
}

#' Plot function for summary function for augsynth
#' @param x summary.augsynth_multiout object
#' @param inf Boolean, whether to plot uncertainty intervals, default TRUE
#' @param plt_avg Boolean, whether to plot the average of the outcomes, default FALSE
#' 
#' @export
plot.summary.augsynth_multiout <- function(x, inf = T, plt_avg = F, ...) {
    if(plt_avg) {
      p <- x$att %>%
        ggplot2::ggplot(ggplot2::aes(x=Time, y=Estimate))
    } else {
      p <- x$att %>%
        filter(Outcome != "Average") %>% 
        ggplot2::ggplot(ggplot2::aes(x=Time, y=Estimate))
    }
    if(inf) {
      if(x$inf_type == "jackknife") {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=Estimate-2*Std.Error,
                        ymax=Estimate+2*Std.Error),
                    alpha=0.2, data = . %>% filter(Outcome != "Average"))
      } else if(x$inf_type %in% c("conformal", "jackknife+")) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=lower_bound,
                        ymax=upper_bound),
                    alpha=0.2, data =  . %>% filter(Outcome != "Average"))
      }

    }
    p + ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=x$t_int, lty=2) +
        ggplot2::geom_hline(yintercept=0, lty=2) +
        ggplot2::facet_wrap(~ Outcome, scales = "free_y") +
        ggplot2::theme_bw()

}