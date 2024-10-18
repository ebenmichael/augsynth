

set.seed(7393)
dat = augsynth:::make_synth_data( n_time = 10, n_U = 5, N = 12, long_form = TRUE, tx_impact = 2, tx_shift = 1 )
# ggplot( dat, aes( time, Y, group=ID ) ) + geom_line()

syn = augsynth( Y ~ Tx | X1, unit = ID, time = time, data = dat, fixedeff = TRUE, progfunc = "none" )

test_that("Donor control - high RMSPE multiple doesn't drop any units", {
    syn2 <- update_augsynth(syn, drop = 1000)
    expect_equal(nrow(donor_table(syn2)), 11)
})



test_that("Donor control drops the correct units based on RMSPE", {

    ssyn = summary( syn, inf_type="permutation" )
    if ( FALSE ) {
        plot( ssyn, "outcomes raw average" )
    }
    treated_table(syn)
    ssyn

    dtable <- donor_table(syn)
    dtable
    nrow( dtable )
    trt_RMSPE <- add_inference(syn, inf_type = 'permutation')$results$permutations$placebo_dist %>%
        filter(time < syn$t_int) %>%
        filter(ID == 1) %>%
        pull(RMSPE) %>% unique()

    drop_factor = 1.0
    drop_units <- dtable %>%
        filter(RMSPE > trt_RMSPE * drop_factor) %>%
        pull(ID)
    drop_units

    syn3 <- update_augsynth(syn, drop = drop_factor)
    d3 <- donor_table(syn3)
    expect_true( all( !( drop_units %in% d3$ID ) ) )

    # TODO: Why should the new donor table have 8 rows?   I dropped
    # that check in the following:

    #criteria <- (nrow(donor_table(syn3)) == 8) & (drop_units %in% $ID %>% all() == FALSE)
    #expect_true(criteria)


    # Check that we get the same model parameters from dropping units manually
    dat_tmp <- dat %>% filter(!ID %in% drop_units)
    syn_manual <- augsynth( Y ~ Tx | X1, unit = ID, time = time, data = dat_tmp,
                            fixedeff = TRUE, progfunc = "none" )
    s_syn_manual <- summary(syn_manual)
    s_syn3 <- summary(syn3)

    expect_equal( s_syn_manual$att$Estimate, s_syn3$att$Estimate )

    expect_equal( nrow(donor_table(syn_manual)), nrow(donor_table(syn3)))


    # Check that dropping useless donors doesn't change anything
    dd = setdiff( as.character( 1:12 ), dtable$ID )
    synX = update_augsynth( syn, drop = dd )
    expect_equal( syn$att$Estimate, synX$att$Estimate )
    expect_equal( nrow(donor_table(syn)), nrow(donor_table(synX)))

})



test_that("Donor control drops the correct units based on unit names", {
    drop_units <- c('2', '3', '4')

    syn4 <- update_augsynth(syn, drop = drop_units)

    expect_true( all( ! (drop_units %in% donor_table(syn4)$ID ) ) )


    # Check that we get the same model parameters from dropping units manually
    dat_tmp <- dat %>% filter(!ID %in% drop_units)
    syn_manual <- augsynth( Y ~ Tx | X1, unit = ID, time = time,
                            data = dat_tmp, fixedeff = TRUE,
                            progfunc = "none" )
    s_syn_manual <- summary(syn_manual)
    s_syn4 <- summary(syn4)

    expect_equal( s_syn_manual$att$Estimate, s_syn4$att$Estimate )
    expect_equal( nrow(donor_table(syn_manual)), nrow(donor_table(syn4)))


})


test_that("`update_augsynth` returns the same basic SCM in cases when updates are not applied", {

    syn1 <- augsynth(lngdpcapita ~ treated, fips, year_qtr, kansas,
                     progfunc = "None", scm = T)
    syn2 <- update_augsynth(syn1, Inf) # dropping based on RMSPE multiple
    syn3 <- update_augsynth(syn1, "") # dropping based on unit names

    # test that weights are the same
    expect_equal(syn2$weights, syn3$weights)
    expect_equal(syn1$weights, syn2$weights)

    # test that ATTs are the same
    summ1 <- summary(syn1)
    summ2 <- summary(syn2)
    summ3 <- summary(syn3)
    expect_equal(summ2$att['Estimate'], summ3$att['Estimate'])
    expect_equal(summ1$att['Estimate'], summ2$att['Estimate'])
})




