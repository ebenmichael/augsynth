

#' Plot function for augsynth or summary.augsynth objects
#'
#' @importFrom graphics plot
#'
#' @param augsynth Augsynth or summary.augsynth object to be plotted
#' @param plot_type The stylized plot type to be returned. Options include
#'        \itemize{
#'          \item{"estimate"}{The ATT and 95\% confidence interval}
#'          \item{"estimate only"}{The ATT without a confidence interval}
#'          \item{"outcomes"}{The level of the outcome variable for the treated and synthetic control units.}
#'          \item{"outcomes raw average"}{The level of the outcome variable for the treated and synthetic control units, along with the raw average of the donor units.}
#'          \item{"placebo"}{The ATTs resulting from placebo tests on the donor units.}  }
#' @param cv If True, plot cross validation MSE against hyper-parameter, otherwise plot effects
#' @param inf_type Type of inference algorithm. Inherits inf_type from `object` or otherwise defaults to "conformal". Options are
#'         \itemize{
#'          \item{"conformal"}{Conformal inference (default)}
#'          \item{"jackknife+"}{Jackknife+ algorithm over time periods}
#'          \item{"jackknife"}{Jackknife over units}
#'          \item{"permutation"}{Placebo permutation, raw ATT}
#'          \item{"permutation_rstat"}{Placebo permutation, RMSPE adjusted ATT}
#'          \item{"None"}{Return ATT Estimate only}
#'         }
#' @param ... Optional arguments for inference, for more details for each `inf_type` see
#'         \itemize{
#'          \item{"conformal"}{`conformal_inf`}
#'          \item{"jackknife+"}{`time_jackknife_plus`}
#'          \item{"jackknife"}{`jackknife_se_single`}
#'          \item{"permutation"}{`permutation_inf`}
#'         }
#' @param ... Optional arguments
plot_augsynth_results <- function( augsynth,
                          cv = FALSE,
                          plot_type = 'estimate',
                          inf_type = NULL, ...) {

    stopifnot(tolower(inf_type) %in% c('conformal', 'jackknife', 'jackknife+', 'permutation', 'permutation_rstat', 'none'))

    # Summarize object if needed.
    if ( is.augsynth(augsynth) ) {
        augsynth = summary(augsynth, inf_type=inf_type)
    }

    if (cv == TRUE) {
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
        return(p)
    } else if (plot_type == 'estimate only') {
        p <- augsynth_plot_from_results(augsynth, ci = FALSE)
    } else if (plot_type == 'estimate') {
        p <- augsynth_plot_from_results(augsynth, ci = TRUE)
    } else if (grepl('placebo', plot_type)) {
        p <- permutation_plot(augsynth, inf_type = augsynth$inf_type)
    } else if (plot_type == 'outcomes') {
        p <- augsynth_outcomes_plot(augsynth, measure = 'synth')
    } else if (plot_type == 'outcomes raw average') {
        p <- augsynth_outcomes_plot(augsynth, measure = c('synth', 'average'))
    }
    return(p)
}






#' Plot function for summary function for augsynth
#'
#' @param augsynth Summary object
#' @param ... Optional arguments
#'
#' @noRd
augsynth_plot_from_results <- function(augsynth,
                                       ci = TRUE,
                                       ...) {

    pdat = NA
    if ( is.augsynth(augsynth) ) {
        pdat <- augsynth$results$att
    } else {
        # Summary object
        pdat <- augsynth$att
    }

    p <- pdat %>%
        ggplot2::ggplot(ggplot2::aes(x = Time, y = Estimate))

    if (ci) {
        if(all(is.na(pdat$lower_bound))) {
            p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = Estimate - 2 * Std.Error,
                                                       ymax = Estimate + 2 * Std.Error),
                                          alpha = 0.2)
        } else {
            p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower_bound,
                                                       ymax = upper_bound),
                                          alpha = 0.2)
        }

    }
    p <- p + ggplot2::geom_line() +
        ggplot2::geom_vline(xintercept = augsynth$t_int, lty = 2) +
        ggplot2::geom_hline(yintercept = 0, lty = 2) +
        ggplot2::labs(x = augsynth$time_var) +
        ggplot2::theme_bw()

    return(p)

}








#' Plot the original level of the outcome variable for the treated
#' unit and its synthetic counterfactual
#'
#' @param augsynth Augsynth object or augsynth summary object to be plotted
#' @param  measure Whether to plot the synthetic counterfactual or the
#'   raw average of donor units.  Can list both if desired.
#'
#' @noRd
augsynth_outcomes_plot <- function(augsynth, measure = c("synth", "average")) {

    if (!is.summary.augsynth(augsynth)) {
        augsynth <- summary(augsynth)
    }

    series = augsynth$treated_table

    all_y = c(series$Yobs, series$Yhat, series$raw_average)
    max_y <- max(all_y)
    min_y <- min(all_y)

    pt = which( series$tx == 1 )[[1]]
    cut_time = (series$time[pt] + series$time[pt + 1]) / 2

    p <- ggplot2::ggplot( series ) +
        ggplot2::geom_line(aes(x = time, y = Yobs, linetype = as.character(augsynth$trt_unit)))

    if ('synth' %in% measure) {
        p <- p +
            ggplot2::geom_line(aes(x = time, y = Yhat, linetype = 'Synthetic counterfactual'))
    }

    if ('average' %in% measure) {
        p <- p +
            ggplot2::geom_line(aes(x = time, y = raw_average, linetype = 'Donor raw average'))
    }

    p <- p +
        ggplot2::labs(linetype = NULL,
                      x = augsynth$time_var,
                      y = 'Outcome') +
        ggplot2::ylim(min_y, max_y) +
        ggplot2::theme_bw() +
        ggplot2::geom_vline(xintercept = cut_time, linetype = 'dashed') +
        ggplot2::theme(legend.position = 'bottom',
                       legend.key = ggplot2::element_rect(fill = scales::alpha("white", 0.5)),
                       legend.background = ggplot2::element_rect(fill = scales::alpha("white", 0))
                       )

    return(p)
}


