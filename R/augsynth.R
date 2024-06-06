
#### Main functions for single-period treatment augmented synthetic controls Method ####



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
#' @param plot Whether or not to return a plot of the augsynth model
#' @param inf_type Type of inference algorithm. Options are
#'         \itemize{
#'          \item{"conformal"}{Conformal inference (default)}
#'          \item{"jackknife+"}{Jackknife+ algorithm over time periods}
#'          \item{"jackknife"}{Jackknife over units}
#'          \item{"permutation"}{Placebo permutation, raw ATT}
#'          \item{"permutation_rstat"}{Placebo permutation, RMSPE adjusted ATT}
#'         }
#' @param ... Optional arguments for inference, for more details for each `inf_type` see
#'         \itemize{
#'          \item{"conformal"}{`conformal_inf`}
#'          \item{"jackknife+"}{`time_jackknife_plus`}
#'          \item{"jackknife"}{`jackknife_se_single`}
#'          \item{"permutation"}{`permutation_inf`}
#'         }
#' @return augsynth object that contains:
#'         \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#'          \item{"data"}{Panel data as matrices}
#'         }
#' @param ... optional arguments for outcome model
#' @export
single_augsynth <- function(form, unit, time, t_int, data,
                            progfunc = "ridge",
                            scm=T,
                            fixedeff = FALSE,
                            cov_agg=NULL,
                            inf_type= 'none', #c( 'none', "conformal", "jackknife+", "jackknife", "permutation", "permutation_rstat"),
                            ...) {

    call_name <- match.call()
    #inf_type = match.arg(inf_type)

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

    # TODO update similar attribute for multi <- if want to use for plotting later
    augsynth$trt_unit <- data %>% filter(!!as.name(trt) == 1) %>%
        pull(quo_name(unit)) %>% unique()
    augsynth$time_var <- quo_name(time)
    augsynth$unit_var <- quo_name(unit)
    augsynth$raw_data <- data
    augsynth$form <- form
    augsynth$cov_agg <- cov_agg

    if (tolower(inf_type) != 'none') {
        augsynth <- add_inference(augsynth, inf_type = inf_type)
    }

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
        progfuncs = c("ridge", "None", "en", "rf", "gsyn", "mcp",
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




#' Function to add inference to augsynth object
#' @param object augsynth object
#' @param inf_type Type of inference algorithm. Options are
#'         \itemize{
#'          \item{"conformal"}{Conformal inference (default)}
#'          \item{"jackknife+"}{Jackknife+ algorithm over time periods}
#'          \item{"jackknife"}{Jackknife over units}
#'          \item{"permutation"}{Placebo permutation, raw ATT}
#'          \item{"permutation_rstat"}{Placebo permutation, RMSPE adjusted ATT}
#'         }
#' @param linear_effect Boolean, whether to invert the conformal inference hypothesis test to get confidence intervals for a linear-in-time treatment effect: intercept + slope * time
#' @param ... Optional arguments for inference, for more details for each `inf_type` see
#'         \itemize{
#'          \item{"conformal"}{`conformal_inf`}
#'          \item{"jackknife+"}{`time_jackknife_plus`}
#'          \item{"jackknife"}{`jackknife_se_single`}
#'          \item{"permutation"}{`permutation_inf`}
#'         }
#' @export

add_inference <- function(object, inf_type = "conformal", linear_effect = F, ...) {
    augsynth <- object

    t0 <- ncol(augsynth$data$X)
    t_final <- t0 + ncol(augsynth$data$y)

    if (tolower(inf_type) != "none") {

        if(inf_type == "jackknife") {
            att_se <- jackknife_se_single(augsynth)
        } else if(inf_type == "jackknife+") {
            att_se <- time_jackknife_plus(augsynth, ...)
        } else if(inf_type == "conformal") {
            att_se <- conformal_inf(augsynth, ...)

            if(linear_effect) {
                att_linear <- conformal_inf_linear(augsynth, ...)
            }

        } else if (inf_type %in% c('permutation', 'permutation_rstat')) {
            if (is.null(augsynth$results$permutations)) {
                augsynth <- add_placebo_distribution(augsynth)
            }
            att_se <- permutation_inf(augsynth, inf_type)
        } else {
            stop(paste(inf_type, "is not a valid choice of 'inf_type'"))
        }

        att <- data.frame(Time = augsynth$data$time,
                          Estimate = att_se$att[1:t_final])
        rownames(att) <- att$Time

        if(inf_type == "jackknife") {
            att$Std.Error <- att_se$se[1:t_final]
            att_avg_se <- att_se$se[t_final + 1]
        } else {
            att_avg_se <- NA
        }
        if( length( att_se$att ) > t_final ) {
            att_avg <- att_se$att[t_final + 1]
        } else {
            att_avg <- mean(att$Estimate[(t0 + 1):t_final])
        }
        if(inf_type %in% c("jackknife+", "nonpar_bs", "t_dist", "conformal", "permutation", "permutation_rstat")) {
            att$lower_bound <- att_se$lb[1:t_final]
            att$upper_bound <- att_se$ub[1:t_final]
        }
        if(inf_type %in% c("conformal", "permutation", "permutation_rstat")) {
            att$p_val <- att_se$p_val[1:t_final]
        }
    } else {
        # No inference, make table of estimates and NAs for SEs, etc.
        att_est <- predict(augsynth, att = TRUE)
        att <- data.frame(Time = augsynth$data$time,
                          Estimate = att_est)
        att$Std.Error <- NA
        att_avg <- mean(att_est[(t0 + 1):t_final])
        att_avg_se <- NA
    }

    augsynth$results$att <- att
    augsynth$results$average_att <- data.frame(Value = "Average Post-Treatment Effect",
                                               Estimate = att_avg, Std.Error = att_avg_se)


    if(inf_type %in% c("jackknife+", "conformal", "permutation", "permutation_rstat")) {
        augsynth$results$average_att$lower_bound <- att_se$lb[t_final + 1]
        augsynth$results$average_att$upper_bound <- att_se$ub[t_final + 1]
        augsynth$results$alpha <- att_se$alpha

        if (inf_type == 'conformal') {
            if(linear_effect) {
                augsynth$results$average_att <- data.frame(
                    Value = c("Average Post-Treatment Effect",
                              "Treatment Effect Intercept",
                              "Treatment Effect Slope"),
                    Estimate = c(att_avg, att_linear$est_int,
                                 att_linear$est_slope),
                    Std.Error = c(att_avg_se, NA, NA),
                    p_val = c(att_se$p_val[t_final + 1], NA, NA),
                    lower_bound = c(att_se$lb[t_final + 1],
                                    att_linear$ci_int[1],
                                    att_linear$ci_slope[1]),
                    upper_bound =  c(att_se$ub[t_final + 1],
                                     att_linear$ci_int[2],
                                     att_linear$ci_slope[2])
                )
            } else {
                augsynth$results$average_att <- data.frame(
                    Value = c("Average Post-Treatment Effect"),
                    Estimate = att_avg,
                    Std.Error = att_avg_se,
                    p_val = att_se$p_val[t_final + 1],
                    lower_bound = att_se$lb[t_final + 1],
                    upper_bound =  att_se$ub[t_final + 1]
                )

            }
        }
    }
    if(inf_type %in% c("conformal", "permutation", "permutation_rstat")) {
        augsynth$results$average_att$p_val <- att_se$p_val[t_final + 1]
    }

    augsynth$results$inf_type <- inf_type

    return(augsynth)

}


#### Summary methods ####


#' Summary function for augsynth
#'
#' Summary summarizes an augsynth result by (usually) adding an
#' inferential result, if that has not been calculated already, and
#' calculating a few other summary statistics such as estimated bias.
#' This method does this via `add_inference()`, if inference is
#' needed.
#'
#' @param object augsynth object
#'
#' @param inf_type Type of inference algorithm. If left NULL, inherits
#'   inf_type from `object` or otherwise defaults to "conformal."
#'   Options are
#'         \itemize{
#'          \item{"conformal"}{Conformal inference (default)}
#'          \item{"jackknife+"}{Jackknife+ algorithm over time periods}
#'          \item{"jackknife"}{Jackknife over units}
#'          \item{"permutation"}{Placebo permutation, raw ATT}
#'          \item{"permutation_rstat"}{Placebo permutation, RMSPE adjusted ATT}
#'         }
#' @param ... Optional arguments for inference, for more details for
#'   each `inf_type` see
#'         \itemize{
#'          \item{"conformal"}{`conformal_inf`}
#'          \item{"jackknife+"}{`time_jackknife_plus`}
#'          \item{"jackknife"}{`jackknife_se_single`}
#'          \item{"permutation", "permutation_rstat"}{`permutation_inf`}
#'         }
#' @export
summary.augsynth <- function(object, inf_type = NULL, ...) {
    augsynth <- object

    t0 <- ncol(augsynth$data$X)
    t_final <- t0 + ncol(augsynth$data$y)

    if (is.null(inf_type)) {
        if (!is.null(augsynth$results)) {
            inf_type <- augsynth$results$inf_type
        } else {
            inf_type <- 'conformal'
        }
    }
    if (is.null(augsynth$results)) {
        augsynth <- add_inference(augsynth, inf_type = inf_type)
    } else if (augsynth$results$inf_type != inf_type) {
        augsynth <- add_inference(augsynth, inf_type = inf_type)
    }

    # Copy over all of OG object except for data
    nms = names(augsynth)
    nms = nms[!nms %in% c("data", "raw_data", "results")]
    summ <- augsynth$results
    summ <- c( summ, augsynth[nms] )


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

    summ$n_unit = n_unit(augsynth)
    summ$n_time = n_time(augsynth)
    summ$n_tx = n_treated(augsynth)[1]
    summ$time_tx = t0
    summ$donor_table = donor_table( augsynth )

    summ$treated_table = treated_table( augsynth )

    class(summ) <- "summary.augsynth"
    return(summ)
}




#' Methods for accessing details of summary.augsynth object
#'
#' @param x summary.augsynth result object
#'
#' @rdname summary.augsynth_class
#'
#' @return is.summary.augsynth: TRUE if object is a augsynth object.
#'
#' @export
is.summary.augsynth <- function(x) {
    inherits(x, "summary.augsynth")
}



#' Print function for summary function for augsynth
#'
#' @param x summary object
#' @param ... Optional arguments
#' @export
print.summary.augsynth <- function(x, ...) {
    summ <- x

    ## straight from lm
    cat("\nCall:\n", paste(deparse(summ$call), sep="\n", collapse="\n"), "\n", sep="")

    t_final <- nrow(summ$att)
    cat( sprintf( "    Fit to %d units and %d+%d = %d time points; %g treated at %s %g.\n",
                  summ$n_unit, summ$time_tx, t_final - summ$time_tx, t_final,
                  summ$n_tx,
                  summ$time_var,
                  summ$att$Time[[summ$time_tx+1]]) )
    cat( "\n" )

    ## distinction between pre and post treatment
    att_est <- summ$att$Estimate
    t_total <- length(att_est)
    t_int <- summ$att %>% filter(Time <= summ$t_int) %>% nrow()

    att_pre <- att_est[1:(t_int-1)]
    att_post <- att_est[t_int:t_total]


    out_msg <- ""


    # print out average post treatment estimate
    att_post <- summ$average_att$Estimate[1]
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
                         format(att_post, digits = 3),
                         "  (",
                         format(p_val, digits = 2), ")\n")
        inf_type <- "Conformal inference"

      if("Treatment Effect Slope" %in% summ$average_att$Value) {
        lowers <- summ$average_att$lower_bound[2:3]
        uppers <- summ$average_att$upper_bound[2:3]
        out_msg_line2 <- paste0("Confidence intervals for linear-in-time treatment effects (Intercept + Slope * Time)\n",
        "\tIntercept: [", format(lowers[1], digits = 3), ",",
        format(uppers[1], digits = 3), "]\n",
        "\tSlope: [", format(lowers[2], digits = 3), ",",
        format(uppers[2], digits = 3), "]\n")
        out_msg <- paste0(out_msg, out_msg_line2)
      }

    } else if(summ$inf_type == "jackknife+") {
        out_msg <- paste("Average ATT Estimate: ",
                         format(round(att_post,3), nsmall=3), "\n")
        inf_type <- "Jackknife+ over time periods"
    } else if (summ$inf_type %in% c('permutation', "permutation_rstat")) {
        out_msg <- paste("Average ATT Estimate: ",
                         format(round(att_post,3), nsmall=3), "\n")
        inf_type <- ifelse(summ$inf_type == 'permutation',
                           "Permutation inference",
                           "Permutation inference (RMSPE-adjusted)")
        out_msg <-paste0( out_msg, "\n",
                         ( sprintf( "Donor RMSPE range from %.2f to %.2f\n",
                                    min( summ$donor_table$RMSPE ), max( summ$donor_table$RMSPE ) ) ) )

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

    rng = range( summ$donor_table$weight[ summ$donor_table$weight > 1/(1000*summ$n_unit) ] )
    cat( sprintf( "%d donor units used with weights of %.3f to %.3f\n",
                  sum( summ$donor_table$weight > 1/(1000*summ$n_unit) ), rng[[1]], rng[[2]] ) )
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
    } else if (summ$inf_type %in% c('permutation', "permutation_rstat")) {
        out_att <- summ$att[t_int:t_final, ] %>%
            select(Time, Estimate, lower_bound, upper_bound, p_val)
        names(out_att) <- c("Time", "Estimate",
                            paste0((1 - summ$alpha) * 100, "% CI Lower Bound"),
                            paste0((1 - summ$alpha) * 100, "% CI Upper Bound"),
                            paste0('p Value'))
    } else {
        out_att <- summ$att[t_int:t_final,] %>%
            select(Time, Estimate)
    }
    out_att %>%
        mutate_at(vars(-Time), ~ round(., 3)) %>%
        print(row.names = F)

    invisible( summ )
}



#' Plot function for summary function for augsynth
#'
#' @param x summary.augsynth object
#' @inheritParams plot_augsynth_results
#'
#' @export
plot.summary.augsynth <- function(x,
                                  plot_type = 'estimate',
                                  ...) {

    summ <- x

    plot_augsynth_results( summ, plot_type = plot_type, ... )

}





#### Package documentation ####


#' augsynth
#'
#' @description A package implementing the Augmented Synthetic Controls Method
#' @name augsynth-package
#' @importFrom magrittr "%>%"
#' @importFrom purrr reduce
#' @import dplyr
#' @import tidyr
#' @importFrom stats terms
#' @importFrom stats formula
#' @importFrom stats update
#' @importFrom stats delete.response
#' @importFrom stats model.matrix
#' @importFrom stats model.frame
#' @importFrom stats na.omit
"_PACKAGE"

