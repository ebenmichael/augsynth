
# Utility function to make fake data for testing and whatnot



#' A synthetic data simulator
#'
#' This generates data with time shocks that interact with the unit
#' latent factors.  It also gives a unit fixed effect and a time fixed
#' effect.
#'
#' @param n_time Number of time periods
#' @param n_post Number of the time periods that are after tx onset
#' @param n_U Number of latent factors
#' @param N Number of units
#' @param N_tx number of treated units, out of total N
#' @param sd_time How much the time varying shocks matter.
#' @param sd_unit_fe How much the units have different shifts
#' @param sd_time_fe How much the time intercept shifts vary
#' @param tx_shift Add this shift to the distribution of the
#'   treatment latent factors to create differences in tx and co
#'   groups.
#'
#' @noRd
make_synth_data = function(n_U, N, n_time, N_tx = 1, n_post = 3,
                           long_form = FALSE,
                           tx_impact = 1:n_post,
                           sd_time = 0.3,
                           sd_unit_fe = 1,
                           sd_time_fe = 0.2,
                           sd_e = 1,
                           tx_shift = 0.5 ) {

    stopifnot( N_tx < N )
    stopifnot( n_post < n_time )

    # Make correlation structure for latent factors
    Sigma = matrix(0.15, nrow = n_U, ncol = n_U)
    diag(Sigma) = 1
    #solve(Sigma)
    U = MASS::mvrnorm(N, mu = rep(1, n_U), Sigma = Sigma)
    U = abs( U )
    U = cbind( U, 1 )
    U[1,] = sd_unit_fe * U[1,]
    #U
    #summary(U)
    U[1:N_tx,1:n_U] = U[1:N_tx,1:n_U] + tx_shift

    # The time varying component, with the first row being an intercept
    shocks = matrix(rnorm(n_U * n_time, sd=sd_time), nrow = n_U)
    shocks = abs( shocks )
    shocks = rbind( 1, shocks )
    #shocks[1, ] = 2 * sort(shocks[1, ])
    #shocks[2, ] = sort(shocks[2, ])
    #shocks[2, ] = 0.5 * rev(sort(shocks[2, ]))
    #browser()

    # Make a time drift by having the fixed effect time shock increase over time.
    shocks[n_U+1,] = sort( shocks[n_U+1,] )
    shocks[n_U+1,] = sd_time_fe * shocks[n_U+1,] / sd_time
    # Alt way to create time drift
    # shocks = shocks + rep( 1:n_time, each=n_U ) / n_time

    #shocks = shocks + (1:nrow(shocks))/nrow(shocks)
    #qplot(1:n_time, shocks[1, ])
    #dim(U)
    #dim(shocks)

    # Outcome is latent factors times time shocks, with extra noise
    # added.
    Y = U %*% shocks
    Y = Y + rnorm( length(Y), sd = sd_e )
    #dim(Y)

    dat = as.data.frame(Y)
    colnames(dat) = 1:n_time
    dat$ID = 1:N
    dat$Tx = 0
    dat$Tx[1:N_tx] = 1


    # Add in some covariates predictive of latent U (skipping
    # intercept)
    X = U[,1:n_U, drop=FALSE]
    for (i in 1:n_U) {
        X[, i] = X[, i] + rnorm(nrow(X))
    }
    X = round( X, digits = 1 )

    #head(dat)
    colnames(X) = paste0("X", 1:ncol(X))
    X = as.data.frame(X)
    dat = bind_cols(dat, X)

    # Add in a treatment impact!
    tx_index = (n_time-n_post+1):n_time
    imp = matrix( tx_impact,
                  ncol=n_post,
                  nrow=N_tx,
                  byrow = TRUE )
    dat[ dat$Tx == 1, tx_index ] = dat[ dat$Tx == 1, tx_index ] + imp


    # Final packaging of the data
    if ( long_form ) {

        ldat = pivot_longer(
            dat,
            cols = all_of(1:n_time),
            names_to = "time",
            values_to = "Y"
        )

        # Make Y look nice.
        ldat$Y = round( ldat$Y, digits = 1 ) #round( ldat$Y / 5, digits = 1 ) * 5

        ldat$time = as.numeric(ldat$time)
        ldat = mutate( ldat,
                       ever_Tx = Tx,
                       Tx = ifelse( ever_Tx & time > n_time - n_post, 1, 0 ) )
        ldat
    } else {
        dat
    }
}


