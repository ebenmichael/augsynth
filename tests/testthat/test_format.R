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
