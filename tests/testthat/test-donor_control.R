

test_that("donor control methods work", {

    set.seed(5344)

    dat = augsynth:::make_synth_data( n_time = 10, n_U = 5, N = 12, long_form = TRUE, tx_impact = 2, tx_shift = 1 )

    syn = augsynth( Y ~ Tx | X1, unit = ID, time = time, data=dat, fixedeff = TRUE, scm = FALSE )
    syn

    a1 = update_augsynth( syn, drop = 1000 )

    syn$call = a1$call=NULL
    expect_equal( syn, a1 )

    synt = donor_table( syn )
    a1 = update_augsynth( syn, drop = 1 )
    a1t = donor_table( a1 )
    expect_true( nrow( a1t ) < nrow( synt ) )

    syn = single_augsynth( Y ~ Tx, unit = ID, time = time, t_int = 8,
                           data=dat, scm = TRUE, progfunc = "none",
                           fixedeff = FALSE )
    syn

    sum <- summary( syn, inf_type = "permutation_rstat" )
    RMSPE( syn )



    synt = donor_table( syn )
    synt
    a1 = update_augsynth( syn, drop = c("2", "4", "5") )
    t1 = treated_table( syn )
    t1$raw_average = NULL
    t2 = treated_table( a1 )
    t2$raw_average = NULL
    expect_equal( t1, t2 )

    names( syn )
    names( a1 )
    syn$call = a1$call=NULL
    expect_equal( syn$l2_imbalance, a1$l2_imbalance )


})
