library( augsynth )
set.seed( 4040440 )

n_units <- 12
n_time <- 10

dat = augsynth:::make_synth_data( n_time = n_time, n_U = 5, N = n_units, long_form = TRUE, tx_impact = 2, tx_shift = 1 )
syn = augsynth( Y ~ Tx | X1 + X2 + X3 + X4 + X5, unit = ID, time = time, data = dat, progfunc="none")

test_that( "covariate table works", {
    cov = covariate_balance_table( syn )
    cov
    expect_true( all( paste0("X", 1:5 ) %in% cov$variable ) )

    expect_output( print( cov ), "variable.*Tx.*Co.*Raw")

    c2 = covariate_balance_table( syn, pre_period = c("1", "2") )
    c2
    aa <- dat$Y[ dat$ever_Tx & dat$time %in% c("1", "2") ]
    expect_equal( c2$Tx[ c2$variable == "2" ], aa[[2]] )
    expect_equal( c2$Tx[ c2$variable == "1" ], aa[[1]] )
})


test_that("summary.augsynth works", {

    expect_output( print( syn ), glue::glue("Fit to {n_units} units and {n_time} time points.*Average ATT Estimate") )

    sum = summary( syn )

    # get donor table, add up number of units without weights and then check that the number in the summary is the same as that
    n_donor <- sum$donor_table %>% filter(weight > 0) %>% nrow()
    expect_output( print( sum ), glue::glue("{n_donor} donor units used with weights."))

    s2 = summary( syn, inf_type = "jackknife+" )
    n_donor2 <- s2$donor_table %>% filter(weight > 0) %>% nrow()
    expect_output( print( s2 ),
                   glue::glue("{n_donor2} donor units used with weights.*Avg Estimated Bias: .*Jackknife\\+ over time periods"))


    s3 = summary( syn, inf_type = "none" )
    expect_output( print( s3 ),
                   "Inference type: None" )

})
