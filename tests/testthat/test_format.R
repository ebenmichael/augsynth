context("Test data formatting")

test_that("format_synth creates matrices of the right dimensions", {
    
    ## simulate data
    mo <- sim_factor_model(20, 25, 15, 10, 1)
    data_out <- format_synth(mo$outcomes %>% filter(outcome_id == 1),
                             mo$metadata)
    test_dim <- function(obj, d) {
        expect_equivalent(dim(obj), d)
        }
    test_dim(data_out$synth_data$Z0,
                      c(14, 19))
    test_dim(data_out$synth_data$Z1,
                      c(14, 1))
    test_dim(data_out$synth_data$Y0plot,
                      c(25, 19))
    test_dim(data_out$synth_data$Y1plot,
             c(25, 1))
    }
    )

test_that("format_synth and format_synth_multi have the same output with one outcome", {
    ## simulate data
    mo <- sim_factor_model(20, 25, 15, 10, 1)
    data_out1 <- format_data(mo$outcomes %>% filter(outcome_id == 1),
                              mo$metadata)
    data_out2 <- format_data(mo$outcomes %>% filter(outcome_id == 1),
                                    mo$metadata,
                                    outcome_col = "outcome_id")

    expect_equivalent(data_out1$synth_data$Z0, data_out2$synth_data$Z0)
    expect_equivalent(data_out1$synth_data$Z1, data_out2$synth_data$Z1)
    expect_equivalent(data_out1$synth_data$Y0plot, data_out2$synth_data$Y0plot)
    expect_equivalent(data_out1$synth_data$Y1plot, data_out2$synth_data$Y1plot)
})
