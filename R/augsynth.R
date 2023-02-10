################################################################################
## Main functions for single-period treatment augmented synthetic controls Method
################################################################################


#' Fit Augmented SCM
#' 
#' @param form outcome ~ treatment | auxillary covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param t_int Time of intervention
#' @param data Panel data as dataframe
#' @param progfunc What function to use to impute control outcomes
#'                 ridge=Ridge regression (allows for standard errors),
#'                 none=No outcome model,
#'                 en=Elastic Net, RF=Random Forest, GSYN=gSynth,
#'                 mcp=MCPanel, 
#'                 cits=Comparitive Interuppted Time Series
#'                 causalimpact=Bayesian structural time series with CausalImpact
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
                     progfunc = "ridge",
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
    
    treated_units <- data %>% filter(!!trt == 1) %>% distinct(!!unit) %>% pull(!!unit)
    control_units <- data %>% filter(!(!!unit %in% treated_units)) %>% 
                        distinct(!!unit) %>% arrange(!!unit) %>% pull(!!unit)
        ## add covariates
    if(length(form)[2] == 2) {
        Z <- extract_covariates(form, unit, time, t_int, data, cov_agg)
    } else {
        Z <- NULL
    }
    
    # fit augmented SCM
    augsynth <- fit_augsynth_internal(wide, synth_data, Z, progfunc, 
                                      scm, fixedeff, ...)
    
    # add some extra data
    augsynth$data$time <- data %>% distinct(!!time) %>%
                                   arrange(!!time) %>% pull(!!time)
    augsynth$call <- call_name
    augsynth$t_int <- t_int 
    
    augsynth$weights <- matrix(augsynth$weights)
    rownames(augsynth$weights) <- control_units

    return(augsynth)
}


#' Internal function to fit augmented SCM
#' @param wide Data formatted from format_data
#' @param synth_data Data formatted from foramt_synth
#' @param Z Matrix of auxiliary covariates
#' @param progfunc outcome model to use
#' @param scm Whether to fit SCM
#' @param fixedeff Whether to de-mean synth
#' @param V V matrix for Synth, default NULL
#' @param ... Extra args for outcome model
#' 
#' @noRd
#' 
fit_augsynth_internal <- function(wide, synth_data, Z, progfunc,
                                  scm, fixedeff, V = NULL, ...) {

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
    if (is.null(progfunc)) {
        progfunc = "none"
    }
    progfunc = tolower(progfunc)
    ## fit augsynth
    if(progfunc == "ridge") {
        # Ridge ASCM
        augsynth <- do.call(fit_ridgeaug_formatted,
                            list(wide_data = fit_wide,
                                 synth_data = fit_synth_data,
                                 Z = Z, V = V, scm = scm, ...))
    } else if(progfunc == "none") {
        ## Just SCM
        augsynth <- do.call(fit_ridgeaug_formatted,
                        c(list(wide_data = fit_wide, 
                               synth_data = fit_synth_data,
                               Z = Z, ridge = F, scm = T, V = V, ...)))
    } else {
        ## Other outcome models
        progfuncs = c("ridge", "none", "en", "rf", "gsyn", "mcp",
                      "cits", "causalimpact", "seq2seq")
        if (progfunc %in% progfuncs) {
            augsynth <- fit_augsyn(fit_wide, fit_synth_data, 
                                   progfunc, scm, ...)
        } else {
            stop("progfunc must be one of 'EN', 'RF', 'GSYN', 'MCP', 'CITS', 'CausalImpact', 'seq2seq', 'None'")
        }
        
    }

    augsynth$mhat <- mhat + cbind(matrix(0, nrow = n, ncol = t0), 
                                  augsynth$mhat)
    augsynth$data <- wide
    augsynth$data$Z <- Z
    augsynth$data$synth_data <- synth_data
    augsynth$progfunc <- progfunc
    augsynth$scm <- scm
    augsynth$fixedeff <- fixedeff
    augsynth$extra_args <- list(...)
    if(progfunc == "ridge") {
        augsynth$extra_args$lambda <- augsynth$lambda
    } else if(progfunc == "gsyn") {
        augsynth$extra_args$r <- ncol(augsynth$params$factor)
        augsynth$extra_args$CV <- 0
    }
    ##format output
    class(augsynth) <- "augsynth"
    return(augsynth)
}

#' Get prediction of ATT or average outcome under control
#' @param object augsynth object
#' @param att If TRUE, return the ATT, if FALSE, return imputed counterfactual
#' @param ... Optional arguments
#'
#' @return Vector of predicted post-treatment control averages
#' @export
predict.augsynth <- function(object, att = F, ...) {
    # if ("att" %in% names(list(...))) {
    #     att <- list(...)$att
    # } else {
    #     att <- F
    # }
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
        rnames <- rownames(y0)
        y0_vec <- c(y0)
        names(y0_vec) <- rnames
        return(y0_vec)
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
#' @param inf Boolean, whether to get confidence intervals around the point estimates
#' @param cv If True, plot cross validation MSE against hyper-parameter, otherwise plot effects
#' @param ... Optional arguments
#' @export
plot.augsynth <- function(x, inf = T, cv = F, ...) {
    # if ("se" %in% names(list(...))) {
    #     se <- list(...)$se
    # } else {
    #     se <- T
    # }

    augsynth <- x
    
    if (cv == T) {
        errors = data.frame(lambdas = augsynth$lambdas,
                            errors = augsynth$lambda_errors,
                            errors_se = augsynth$lambda_errors_se)
        p <- ggplot2::ggplot(errors, ggplot2::aes(x = lambdas, y = errors)) +
              ggplot2::geom_point(size = 2) + 
              ggplot2::geom_errorbar(
                ggplot2::aes(ymin = errors,
                             ymax = errors + errors_se),
                width=0.2, size = 0.5) 
        p <- p + ggplot2::labs(title = bquote("Cross Validation MSE over " ~ lambda),
                              x = expression(lambda), y = "Cross Validation MSE", 
                              parse = TRUE)
        p <- p + ggplot2::scale_x_log10()
        
        # find minimum and min + 1se lambda to plot
        min_lambda <- choose_lambda(augsynth$lambdas,
                                   augsynth$lambda_errors,
                                   augsynth$lambda_errors_se,
                                   F)
        min_1se_lambda <- choose_lambda(augsynth$lambdas,
                                       augsynth$lambda_errors,
                                       augsynth$lambda_errors_se,
                                       T)
        min_lambda_index <- which(augsynth$lambdas == min_lambda)
        min_1se_lambda_index <- which(augsynth$lambdas == min_1se_lambda)

        p <- p + ggplot2::geom_point(
            ggplot2::aes(x = min_lambda, 
                         y = augsynth$lambda_errors[min_lambda_index]),
            color = "gold")
        p + ggplot2::geom_point(
              ggplot2::aes(x = min_1se_lambda,
                           y = augsynth$lambda_errors[min_1se_lambda_index]),
              color = "gold") +
            ggplot2::theme_bw()
    } else {
        plot(summary(augsynth, ...), inf = inf)
    }
}


#' Summary function for augsynth
#' @param object augsynth object
#' @param inf Boolean, whether to get confidence intervals around the point estimates
#' @param inf_type Type of inference algorithm. Options are
#'         \itemize{
#'          \item{"conformal"}{Conformal inference (default)}
#'          \item{"jackknife+"}{Jackknife+ algorithm over time periods}
#'          \item{"jackknife"}{Jackknife over units}
#'         }
#' @param ... Optional arguments for inference, for more details for each `inf_type` see
#'         \itemize{
#'          \item{"conformal"}{`conformal_inf`}
#'          \item{"jackknife+"}{`time_jackknife_plus`}
#'          \item{"jackknife"}{`jackknife_se_single`}
#'         }
#' @export
summary.augsynth <- function(object, inf = T, inf_type = "conformal", ...) {
    augsynth <- object
    # if ("inf" %in% names(list(...))) {
    #     inf <- list(...)$inf
    # } else {
    #     inf <- T
    # }
    # if ("inf_type" %in% names(list(...))) {
    #     inf_type <- list(...)$inf_type
    # } else {
    #     inf_type <- "conformal"
    # }
    
    
    summ <- list()

    t0 <- ncol(augsynth$data$X)
    t_final <- t0 + ncol(augsynth$data$y)

    if(inf) {
        if(inf_type == "jackknife") {
            att_se <- jackknife_se_single(augsynth)
        } else if(inf_type == "jackknife+") {
            att_se <- time_jackknife_plus(augsynth, ...)
        } else if(inf_type == "conformal") {
          att_se <- conformal_inf(augsynth, ...)
        } else {
            stop(paste(inf_type, "is not a valid choice of 'inf_type'"))
        }

        att <- data.frame(Time = augsynth$data$time,
                          Estimate = att_se$att[1:t_final])
        if(inf_type == "jackknife") {
          att$Std.Error <- att_se$se[1:t_final]
          att_avg_se <- att_se$se[t_final + 1]
        } else {
          att_avg_se <- NA
        }
        att_avg <- att_se$att[t_final + 1]
        if(inf_type %in% c("jackknife+", "nonpar_bs", "t_dist", "conformal")) {
            att$lower_bound <- att_se$lb[1:t_final]
            att$upper_bound <- att_se$ub[1:t_final]
        }
        if(inf_type == "conformal") {
          att$p_val <- att_se$p_val[1:t_final]
        }

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
    if(inf) {
      if(inf_type %in% c("jackknife+", "conformal")) {
          summ$average_att$lower_bound <- att_se$lb[t_final + 1]
          summ$average_att$upper_bound <- att_se$ub[t_final + 1]
          summ$alpha <-  att_se$alpha
      }
      if(inf_type == "conformal") {
        summ$average_att$p_val <- att_se$p_val[t_final + 1]
      }
    }
    summ$t_int <- augsynth$t_int
    summ$call <- augsynth$call
    summ$l2_imbalance <- augsynth$l2_imbalance
    summ$scaled_l2_imbalance <- augsynth$scaled_l2_imbalance
    if(!is.null(augsynth$covariate_l2_imbalance)) {
      summ$covariate_l2_imbalance <- augsynth$covariate_l2_imbalance
      summ$scaled_covariate_l2_imbalance <- augsynth$scaled_covariate_l2_imbalance
    }
    ## get estimated bias

    if(tolower(augsynth$progfunc) == "ridge") {
        mhat <- augsynth$ridge_mhat
        w <- augsynth$synw
    } else {
        mhat <- augsynth$mhat
        w <- augsynth$weights
    }
    trt <- augsynth$data$trt
    m1 <- colMeans(mhat[trt==1,,drop=F])

    if(tolower(augsynth$progfunc) == "none" | (!augsynth$scm)) {
        summ$bias_est <- NA
    } else {
      summ$bias_est <- m1 - t(mhat[trt==0,,drop=F]) %*% w
    }
    
    
    summ$inf_type <- if(inf) inf_type else "None"
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


    out_msg <- ""


    # print out average post treatment estimate
    att_post <- summ$average_att$Estimate
    se_est <- summ$att$Std.Error
    if(summ$inf_type == "jackknife") {
      se_avg <- summ$average_att$Std.Error

      out_msg <- paste("Average ATT Estimate (Jackknife Std. Error): ",
                        format(round(att_post,3), nsmall=3), 
                        "  (",
                        format(round(se_avg,3)), ")\n")
      inf_type <- "Jackknife over units"
    } else if(summ$inf_type == "conformal") {
      p_val <- summ$average_att$p_val
      out_msg <- paste("Average ATT Estimate (p Value for Joint Null): ",
                        format(round(att_post,3), nsmall=3), 
                        "  (",
                        format(round(p_val,3)), ")\n")
      inf_type <- "Conformal inference"
    } else if(summ$inf_type == "jackknife+") {
      out_msg <- paste("Average ATT Estimate: ",
                        format(round(att_post,3), nsmall=3), "\n")
      inf_type <- "Jackknife+ over time periods"
    } else {
      out_msg <- paste("Average ATT Estimate: ",
                        format(round(att_post,3), nsmall=3), "\n")
      inf_type <- "None"
    }


    out_msg <- paste(out_msg, 
              "L2 Imbalance: ",
              format(round(summ$l2_imbalance,3), nsmall=3), "\n",
              "Percent improvement from uniform weights: ",
              format(round(1 - summ$scaled_l2_imbalance,3)*100), "%\n\n",
              sep="")
  if(!is.null(summ$covariate_l2_imbalance)) {

    out_msg <- paste(out_msg,
                     "Covariate L2 Imbalance: ",
                     format(round(summ$covariate_l2_imbalance,3), 
                                  nsmall=3),
                    "\n",
                     "Percent improvement from uniform weights: ",
                     format(round(1 - summ$scaled_covariate_l2_imbalance,3)*100), 
                     "%\n\n",
                     sep="")

  }
  out_msg <- paste(out_msg, 
              "Avg Estimated Bias: ",
              format(round(mean(summ$bias_est), 3),nsmall=3), "\n\n",
              "Inference type: ",
              inf_type,
              "\n\n",
              sep="")
  cat(out_msg)

    if(summ$inf_type == "jackknife") {
      out_att <- summ$att[t_int:t_final,] %>% 
              select(Time, Estimate, Std.Error)
    } else if(summ$inf_type == "conformal") {
      out_att <- summ$att[t_int:t_final,] %>% 
              select(Time, Estimate, lower_bound, upper_bound, p_val)
      names(out_att) <- c("Time", "Estimate", 
                          paste0((1 - summ$alpha) * 100, "% CI Lower Bound"),
                          paste0((1 - summ$alpha) * 100, "% CI Upper Bound"),
                          paste0("p Value"))
    } else if(summ$inf_type == "jackknife+") {
      out_att <- summ$att[t_int:t_final,] %>% 
              select(Time, Estimate, lower_bound, upper_bound)
      names(out_att) <- c("Time", "Estimate", 
                          paste0((1 - summ$alpha) * 100, "% CI Lower Bound"),
                          paste0((1 - summ$alpha) * 100, "% CI Upper Bound"))
    } else {
      out_att <- summ$att[t_int:t_final,] %>% 
              select(Time, Estimate)
    }
    out_att %>%
      mutate_at(vars(-Time), ~ round(., 3)) %>%
      print(row.names = F)

    
}

#' Plot function for summary function for augsynth
#' @param x Summary object
#' @param inf Boolean, whether to plot confidence intervals
#' @param ... Optional arguments
#' @export
plot.summary.augsynth <- function(x, inf = T, ...) {
    summ <- x
    # if ("inf" %in% names(list(...))) {
    #     inf <- list(...)$inf
    # } else {
    #     inf <- T
    # }
    
    p <- summ$att %>%
        ggplot2::ggplot(ggplot2::aes(x=Time, y=Estimate))
    if(inf) {
        if(all(is.na(summ$att$lower_bound))) {
            p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=Estimate-2*Std.Error,
                        ymax=Estimate+2*Std.Error),
                    alpha=0.2)
        } else {
            p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=lower_bound,
                        ymax=upper_bound),
                    alpha=0.2)
        }

    }
    p + ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept=summ$t_int, lty=2) +
        ggplot2::geom_hline(yintercept=0, lty=2) + 
        ggplot2::theme_bw()

}



#' augsynth
#' 
#' @description A package implementing the Augmented Synthetic Controls Method
#' @docType package
#' @name augsynth-package
#' @importFrom magrittr "%>%"
#' @importFrom purrr reduce
#' @import dplyr
#' @import LowRankQP
#' @import tidyr
#' @importFrom stats terms
#' @importFrom stats formula
#' @importFrom stats update 
#' @importFrom stats delete.response 
#' @importFrom stats model.matrix 
#' @importFrom stats model.frame 
#' @importFrom stats na.omit
NULL
