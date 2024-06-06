

test_that("summary.augsynth works", {

    set.seed( 4040440 )

    dat = augsynth:::make_synth_data( n_time = 10, n_U = 5, N = 12, long_form = TRUE, tx_impact = 2, tx_shift = 1 )
    head( dat )

    #ggplot( dat, aes( time, Y, color = as.factor( ever_Tx ), group=ID ) ) +
    #    geom_line()

    library( augsynth )
    syn = augsynth( Y ~ Tx | X1 + X2 + X3 + X4 + X5, unit = ID, time = time, data=dat, inf_type = "jackknife" )

    expect_output( print( syn ), "Fit to 12 units and 10 time points.*Average ATT Estimate" )

    sum = summary( syn )
    print( sum )
    expect_output( print( sum ), "5 donor units used with weights.*Avg Estimated Bias: .*Jackknife over units" )

    s2 = summary( syn, inf_type = "jackknife+" )
    print( s2 )
    expect_output( print( s2 ),
                   "5 donor units used with weights.*Avg Estimated Bias: .*Jackknife\\+ over time periods" )


    s2 = summary( syn, inf_type = "none" )
    print( s2 )
    expect_output( print( s2 ),
                   "Inference type: None" )

})
