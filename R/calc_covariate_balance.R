

#' Make covariate balance table.
#'
#' Make a table comparing means of covariates in the treated group,
#' the raw control group, and the new weighted control group (the
#' synthetic control)
#'
#' @param ascm An augsynth result object from single_augsynth
#' @param pre_period List of names of the pre-period timepoints to
#'   calculate balance for.  NULL means none.
#'
#' @export
covariate_balance_table = function( ascm, pre_period = NULL ) {

    stopifnot( is.augsynth( ascm ) )

    trt = ascm$data$trt
    weight = rep( 0, length( trt ) )
    weight[ trt == 0 ] = ascm$synw
    stopifnot( abs( sum( weight ) - 1 ) < 0.000001 )

    Z = ascm$data$Z
    if ( !is.null( pre_period ) ) {
        xx <- ascm$data$X[ , colnames(ascm$data$X) %in% pre_period, drop = FALSE]
        if ( ncol( xx ) > 0 ) {
            Z = cbind( Z, xx )
        }
    }

    Co_means = t( Z ) %*% weight

    # Means of the outcome at lagged time points
    #Co_means_2 = ascm$data$synth_data$Z0 %*% ascm$synw

    # Unform weighting
    n_donor = length(ascm$synw)
    unit_weight = rep( 1 / n_donor, nrow(Z) )
    unit_weight[ trt == 1 ] = 0
    raw_means = t( Z ) %*% unit_weight

    Tx_means = Z[ trt == 1, ]

    means = tibble( variable = names( Tx_means ),
                        Tx = as.numeric( Tx_means ),
                        Co = as.numeric( Co_means ),
                        Raw = as.numeric( raw_means ) )


    means
}
