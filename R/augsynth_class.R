

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
#' @export
plot.augsynth <- function(augsynth,
                          cv = FALSE,
                          plot_type = 'estimate',
                          inf_type = NULL, ...) {

    message("Plotting augsynth objects may slow execution time. For faster results, plot from an augsynth summary object using plot.summary.augsynth()")

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

    plot_augsynth_results( augsynth=augsynth, cv=cv, plot_type=plot_type, inf_type=inf_type, ... )

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


#' RMSPE for treated unit
#'
#' @param augsynth Augsynth object
#' @return RMSPE (Root mean squared predictive error) for the treated unit in pre-treatment era
#'
#' @export
RMSPE <- function( augsynth ) {
    stopifnot( is.augsynth(augsynth) )

    pd = predict( augsynth, att = TRUE )
    sqrt( mean( pd[1:ncol(augsynth$data$X)]^2 ) )
}





#' Get placebo distribution
#'
#' @param augsynth Augsynth object or summery object with permutation inference of some sort.
#'
#' @return Data frame holding the placebo distribution, one row per placebo unit and time point.
#'
#' @export
placebo_distribution <- function( augsynth ) {
    inf_type = NA
    if ( is.summary.augsynth(augsynth) ) {
        inf_type = augsynth$inf_type
    } else if ( is.augsynth(augsynth) ) {
        inf_type = augsynth$results$inf_type
    } else {
        stop( "Object must be an Augsynth object or summary object" )
    }
    if ( !is.null( inf_type ) && inf_type %in% c( "permutation", "permutation_rstat" ) ) {
        if ( is.augsynth(augsynth) ) {
            return( augsynth$results$permutations$placebo_dist )
        } else {
            return( augsynth$permutations$placebo_dist )
        }
    } else {
        stop( "Placebo distribution only available for permutation inference" )
    }
}







#' Return a summary data frame for the treated unit
#'
#' @param augsynth Augsynth object of interest
#'
#' @return Dataframe of information about the treated unit, one row
#'   per time point.  This includes the measured outcome, predicted
#'   outcome from the synthetic unit, the average of all donor units
#'   (as reference, called `raw_average`), and the estimated impact
#'   (`ATT`), and the r-statistic (ATT divided by RMSPE).
#'
#' @seealso [donor_table()]
#' @export
treated_table <- function(augsynth) {

    if ( is.summary.augsynth( augsynth ) ) {
        return( augsynth$treated_table )
    }

    # Calculate the time series of the treated, the synthetic control,
    # and the overall donor pool average
    trt_index <- which(augsynth$data$trt == 1)
    df <- bind_cols(augsynth$data$X, augsynth$data$y)
    # synth_unit <- t(df[-trt_index, ]) %*% augsynth$weights
    synth_unit <- predict(augsynth)
    average_unit <- df[-trt_index, ] %>% colMeans()
    treated_unit <- t(df[trt_index, ])
    lvls = tibble(
        time = as.numeric( colnames(df) ),
        Yobs = as.numeric( treated_unit ),
        Yhat = as.numeric( synth_unit ),
        raw_average = as.numeric( average_unit )
    )

    #lvls <- df %>%
    #    group_by( !!sym(augsynth$time_var ), ever_Tx) %>%
    #    summarise( Yavg = mean( Yobs ), .groups="drop" ) %>%
    #    pivot_wider( names_from = ever_Tx, values_from = Yavg )
    #colnames(lvls)[2:3] <- c("raw_average", "Yobs")

    t0 <- ncol(augsynth$data$X)
    tpost <- ncol(augsynth$data$y)
    lvls$tx = rep( c(0,1), c( t0, tpost ) )
    #lvls$Yhat = predict( augsynth )
    lvls$ATT = lvls$Yobs - lvls$Yhat
    lvls$rstat = lvls$ATT / sqrt( mean( lvls$ATT[ lvls$tx == 0 ]^2 ) )

    lvls <- dplyr::relocate( lvls,
                             time, tx, Yobs, Yhat, raw_average, ATT, rstat )

    return( lvls )
}


