context("Testing that we can load packaged data")

test_that("kansas data loads", {
    expect_error(data(kansas), NA)
})