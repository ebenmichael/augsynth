

# Methods for interacting with the augsynth object






#' Print function for augsynth
#' @param x augsynth object
#' @param ... Optional arguments
#' @export
print.augsynth <- function(x, ...) {
    augsynth <- x

    ## straight from lm
    cat("\nCall:\n", paste(deparse(augsynth$call), sep="\n", collapse="\n"), "\n", sep="")
    cat( sprintf( "    Fit to %d units and %d time points\n\n", n_unit(augsynth), n_time(augsynth) ) )
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
#'
#' @param augsynth Augsynth object to be plotted
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
#' @export
plot.augsynth <- function(augsynth,
                          cv = FALSE, # ND note — not sure what this does?
                          plot_type = 'estimate',
                          inf_type = NULL, ...) {

    if (is.null(inf_type) & !is.null(augsynth$results)) {
        inf_type = augsynth$results$inf_type
    } else if (is.null(inf_type) & is.null(augsynth$results)) {
        # if no inf_type given for a basic (backwards compatible) augsynth object, set inf_type to conformal
        inf_type = 'conformal'
    }

    # if inf_type is set to "none", then only return a raw treatment estimate or an outcomes plot (treated/synth trajectories)
    if ((inf_type %in% c('None', 'none')) & (!grepl('outcomes', plot_type))) {
        plot_type = 'estimate only'
    }

    # if the user specifies the "placebo" plot type without accompanying inference, default to placebo and show message
    if ((plot_type == 'placebo') & (!inf_type %in% c('permutation', 'permutation_rstat'))) {
        message('Placebo plots are only available for permutation-based inference. The plot shows results from "inf_type = "permutation""')
        return(permutation_plot(augsynth, inf_type = 'permutation'))
    }

    if (is.null(inf_type)) {
        if (!is.null(augsynth$results)) {
            inf_type <- augsynth$results$inf_type
        } else {
            inf_type <- 'conformal'
        }
    }

    if (is.null(augsynth$results)) {
        augsynth <- add_inference(augsynth, inf_type = inf_type)
    }

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
        return(p)
    } else if (plot_type == 'estimate only') {
        p <- augsynth_plot_from_results(augsynth, inf_type = 'none')
    } else if (plot_type == 'estimate') {
        p <- augsynth_plot_from_results(augsynth, inf_type = inf_type)
    } else if (grepl('placebo', plot_type)) {
        p <- permutation_plot(augsynth, inf_type = inf_type)
    } else if (plot_type == 'outcomes') {
        p <- augsynth_outcomes_plot(augsynth, measure = 'synth')
    } else if (plot_type == 'outcomes raw average') {
        p <- augsynth_outcomes_plot(augsynth, measure = c('synth', 'average'))
    }
    return(p)
}




#' Methods for accessing details of augsynth result object (of class augsynth)
#'
#'
#' @param x augsynth result object
#'
#' @rdname augsynth_class
#'
#' @return is.augsynth: TRUE if object is a augsynth object.
#'
#' @export
is.augsynth <- function(x) {
    inherits(x, "augsynth")
}



#'
#' @return dim: Dimension of data as pair of (# units, # time points).
#'
#' @rdname augsynth_class
#' @export
#'
dim.augsynth <- function(x, ... ) {
    n_unit = length( unique( x$raw_data[[ x$unit_var ]] ) )
    n_time = length( unique( x$raw_data[[ x$time_var ]] ) )
    return( c( n_unit, n_time ) )
}




#'
#' @return Single number (of unique units).
#'
#' @rdname synth_class
#' @export
#'
n_unit <- function(x, ... ) {
    UseMethod( "n_unit" )
}

#' @title Number of time points in fit data
#'
#' @rdname synth_class
#' @return Single number (of unique time points).
#' @export
n_time <- function(x, ... ) {
    UseMethod( "n_time" )
}



#' @title Number of treated units in fit data
#'
#' @rdname synth_class
#' @return Single number (of number of treated units).
#' @export
n_treated <- function(x, ... ) {
    UseMethod( "n_treated" )
}



#'
#' @return Single number (of unique units).
#'
#' @rdname augsynth_class
#'
#' @export
#'
n_unit.augsynth <- function(x, ... ) {
    dim.augsynth(x)[[1]]
}

#' @title Number of time points in augsynth
#'
#' @rdname augsynth_class
#'
#' @return Single number (of unique time points).
#' @export
n_time.augsynth <- function(x, ... ) {
    dim.augsynth(x)[[2]]
}


#'
#' @rdname augsynth_class
#'
#' @return Number of treated units (always 1 for augsynth)
#' @export
n_treated.augsynth <- function(x, ... ) {
    return( 1 )
}

