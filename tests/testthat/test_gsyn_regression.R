context("GSYN prognostic function on panels with 10+ time periods")

## time ids 1..30 sort differently as strings than as numbers; the basque
## panel's 4-digit years do not, which is what hid #107
test_that("progfunc = 'GSYN' runs on a panel with time ids 1..30", {
    skip_if_not_installed("gsynth")
    set.seed(42)
    n_t <- 30
    units <- paste0("u", 1:6)
    df <- expand.grid(region = units, period = 1:n_t, stringsAsFactors = FALSE)
    df$y <- 100 + 5 * sin(df$period / 3) +
        as.numeric(factor(df$region)) * 10 + rnorm(nrow(df))
    df$trt <- ifelse(df$region == "u1" & df$period > 25, 1, 0)

    asyn <- augsynth(y ~ trt, region, period, df, t_int = 26,
                     progfunc = "GSYN", scm = TRUE, fixedeff = TRUE,
                     CV = 0, r = 0)
    expect_true(is.finite(asyn$scaled_l2_imbalance))
    expect_true(all(is.finite(predict(asyn))))
    expect_equal(summary(asyn, inf_type = "none")$average_att$Estimate,
                 0.4796335, tolerance = 1e-4)
})
