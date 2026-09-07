context("Testing permutation inference guard with multiple treated units")

library( augsynth )
set.seed( 4040440 )

n_units <- 12
n_time <- 10

dat1 <- augsynth:::make_synth_data( n_time = n_time, n_U = 5, N = n_units,
                                    N_tx = 1, tx_impact = 2, tx_shift = 1,
                                    long_form = TRUE )
dat3 <- augsynth:::make_synth_data( n_time = n_time, n_U = 5, N = n_units,
                                    N_tx = 3, tx_impact = 2, tx_shift = 1,
                                    long_form = TRUE )

syn1 <- augsynth( Y ~ Tx, unit = ID, time = time, data = dat1, progfunc = "none" )
syn3 <- augsynth( Y ~ Tx, unit = ID, time = time, data = dat3, progfunc = "none" )

msg <- "only supported for a single treated unit"


test_that( "permutation inference errors clearly with multiple treated units", {

    expect_equal( sum( syn3$data$trt == 1 ), 3 )

    # Previously these tripped the opaque internal assertion
    # "Two versions of estimated impacts do not correspond. Serious error."
    expect_error( summary( syn3, inf_type = "permutation" ), msg )
    expect_error( summary( syn3, inf_type = "permutation_rstat" ), msg )
    expect_error( augsynth:::add_placebo_distribution( syn3 ), msg )
    expect_error( augsynth:::get_placebo_gaps( syn3 ), msg )

    # donor_table() requests permutation inference to compute RMSPEs.
    expect_error( donor_table( syn3 ), msg )
    expect_error( donor_table( syn3, include_RMSPE = TRUE ), msg )

    # The message must not be the old cryptic assertion.
    err <- tryCatch( summary( syn3, inf_type = "permutation" ), error = function( e ) conditionMessage( e ) )
    expect_false( grepl( "Please contact package maintainers", err ) )
})


test_that( "permutation inference still works with a single treated unit", {

    expect_error( summary( syn1, inf_type = "permutation" ), NA )
    expect_error( summary( syn1, inf_type = "permutation_rstat" ), NA )

    s <- summary( syn1, inf_type = "permutation" )
    expect_equal( nrow( s$att ), n_time )
    expect_true( "p_val" %in% names( s$att ) )
})


test_that( "the guard only blocks the permutation path, not other code", {

    # Non-permutation paths must not be caught by the guard. donor_table()
    # without RMSPEs never requests permutation inference, and predict() and
    # single-treated-unit conformal inference are all unrelated to it.
    expect_error( donor_table( syn3, include_RMSPE = FALSE ), NA )
    expect_error( predict( syn3, att = TRUE ), NA )
    expect_error( summary( syn1, inf_type = "conformal" ), NA )
})
