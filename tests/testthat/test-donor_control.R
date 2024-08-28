

set.seed(7393)
dat = augsynth:::make_synth_data( n_time = 10, n_U = 5, N = 12, long_form = TRUE, tx_impact = 2, tx_shift = 1 )

syn = augsynth( Y ~ Tx | X1, unit = ID, time = time, data = dat, fixedeff = TRUE, scm = FALSE )

test_that("Donor control - high RMSPE multiple doesn't drop any units", {
    syn2 <- update_augsynth(syn, drop = 1000)
    expect_equal(nrow(donor_table(syn2)), 11)
})

test_that("Donor control drops the correct units based on 1.5x RMSPE", {

    dtable <- donor_table(syn)
    trt_RMSPE <- add_inference(syn, inf_type = 'permutation')$results$permutations$placebo_dist %>%
        filter(time < syn$t_int) %>%
        filter(ID == 1) %>%
        pull(RMSPE) %>% unique()

    drop_units <- dtable %>%
        filter(RMSPE > trt_RMSPE * 1.5) %>% # should drop units with ID 3, 8, and 11 based on 1.5x treated RMSPE
        pull(ID)

    syn3 <- update_augsynth(syn, drop = 1.5)
    criteria <- (nrow(donor_table(syn3)) == 8) & (drop_units %in% donor_table(syn3)$ID %>% all() == FALSE)
    expect_true(criteria)


    # Check that we get the same model parameters from dropping units manually
    dat_tmp <- dat %>% filter(!ID %in% drop_units)
    syn_manual <- augsynth( Y ~ Tx | X1, unit = ID, time = time, data = dat_tmp, fixedeff = TRUE, scm = FALSE )
    s_syn_manual <- summary(syn_manual)
    s_syn3 <- summary(syn3)

    criteria_1 <- (s_syn_manual$att$Estimate == s_syn3$att$Estimate) %>% all()
    criteria_2 <- (nrow(donor_table(syn_manual)) == nrow(donor_table(syn3)))

    expect_true(criteria_1 & criteria_2)

})

test_that("Donor control drops the correct units based on unit names", {
    drop_units <- c('2', '3', '4')

    syn4 <- update_augsynth(syn, drop = drop_units)
    criteria_1 <- nrow(donor_table(syn4)) == 8
    criteria_2 <- drop_units %in% donor_table(syn4)$ID %>% all() == FALSE

    expect_true(criteria_1 & criteria_2)


    # Check that we get the same model parameters from dropping units manually
    dat_tmp <- dat %>% filter(!ID %in% drop_units)
    syn_manual <- augsynth( Y ~ Tx | X1, unit = ID, time = time, data = dat_tmp, fixedeff = TRUE, scm = FALSE )
    s_syn_manual <- summary(syn_manual)
    s_syn4 <- summary(syn4)

    criteria_1 <- (s_syn_manual$att$Estimate == s_syn4$att$Estimate) %>% all()
    criteria_2 <- (nrow(donor_table(syn_manual)) == nrow(donor_table(syn4)))

    expect_true(criteria_1 & criteria_2)


})

