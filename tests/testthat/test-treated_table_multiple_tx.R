context("Testing treated_table() with multiple treated units")

library( augsynth )
set.seed( 4040440 )

n_units <- 12
n_time <- 10

# Single treated unit (regression guard) and multiple treated units (the bug).
dat1 <- augsynth:::make_synth_data( n_time = n_time, n_U = 5, N = n_units,
                                    N_tx = 1, tx_impact = 2, tx_shift = 1,
                                    long_form = TRUE )
dat3 <- augsynth:::make_synth_data( n_time = n_time, n_U = 5, N = n_units,
                                    N_tx = 3, tx_impact = 2, tx_shift = 1,
                                    long_form = TRUE )

syn1 <- augsynth( Y ~ Tx, unit = ID, time = time, data = dat1, progfunc = "none" )
syn3 <- augsynth( Y ~ Tx, unit = ID, time = time, data = dat3, progfunc = "none" )


test_that( "treated_table() has one row per time period with multiple treated units", {

    expect_equal( sum( syn3$data$trt == 1 ), 3 )

    tt <- treated_table( syn3 )

    # One row per time period, regardless of the number of treated units.
    expect_equal( nrow( tt ), n_time )
    expect_equal( length( tt$Yobs ), n_time )

    # Yobs is the average of the treated units' observed outcomes.
    df <- dplyr::bind_cols( syn3$data$X, syn3$data$y )
    trt_index <- which( syn3$data$trt == 1 )
    expect_equal( tt$Yobs, as.numeric( colMeans( df[ trt_index, , drop = FALSE ] ) ) )
})


test_that( "treated_table() ATT matches predict(att = TRUE) with multiple treated units", {

    tt <- treated_table( syn3 )

    # The averaged treated series is the correct comparison against the pooled
    # synthetic control, so Yobs - Yhat must equal augsynth's own ATT series.
    expect_equal( tt$ATT, as.numeric( predict( syn3, att = TRUE ) ),
                  tolerance = 1e-8 )
})


test_that( "summary.augsynth() works with multiple treated units", {

    # Prior to the fix this errored inside treated_table() with
    # "Tibble columns must have compatible sizes".
    expect_error( summary( syn3 ), NA )
    expect_error( print( summary( syn3 ) ), NA )
})


test_that( "treated_table() is unchanged for a single treated unit", {

    tt <- treated_table( syn1 )

    expect_equal( nrow( tt ), n_time )

    # With one treated unit Yobs is exactly that unit's observed series.
    df <- dplyr::bind_cols( syn1$data$X, syn1$data$y )
    trt_index <- which( syn1$data$trt == 1 )
    expect_equal( tt$Yobs, as.numeric( df[ trt_index, ] ) )
    expect_equal( tt$ATT, as.numeric( predict( syn1, att = TRUE ) ),
                  tolerance = 1e-8 )
})
