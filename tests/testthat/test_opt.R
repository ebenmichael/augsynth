context("Test objective functions and gradients for optimization")

test_that("logsumexp is identity for length one vectors", {
    test_equal <- function(x) expect_equal(logsumexp(x), x)
    test_equal(1)
    test_equal(c(.5))
    test_equal(-.5)
})

test_that("logsumexp doesn't over/underflow", {
    expect_equal(logsumexp(c(1000,100,10000)),
                 10000)
    expect_equal(logsumexp(c(-1000, -999, -1000)),
                 -998.4485552861)
})
