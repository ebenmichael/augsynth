context("Testing if hyperparameter selection methods work")

test_that("No errors in lexical time selection", {

    run_one <- function() {
        set.seed(12345)
        ## simulate data
        mo <- sim_factor_model(50, 50, 10, 10, 1)
        ## fit synthetic controls for one outcome
        out <- lexical_time(mo$outcomes %>% filter(outcome_id == 1),
                              mo$metadata)
    }
    expect_error(run_one(), NA)
    })


test_that("No errors in group time selection", {

    run_one <- function() {
        set.seed(12345)
        ## simulate data
        mo <- sim_factor_model(50, 50, 10, 10, 1)
        ## fit synthetic controls for one outcome
        out <- recent_group(mo$outcomes %>% filter(outcome_id == 1),
                              mo$metadata, 5)
    }
    expect_error(run_one(), NA)
})

test_that("No errors in LASSO time selection", {

    run_one <- function() {
        set.seed(12345)
        ## simulate data
        mo <- sim_factor_model(50, 50, 10, 10, 1)
        ## fit synthetic controls for one outcome
        out <- sep_lasso(mo$outcomes %>% filter(outcome_id == 1),
                         mo$metadata)
    }
    expect_error(run_one(), NA)
})

