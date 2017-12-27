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

test_that("L2 prox and group lass prox are the same for one group", {
    set.seed(1234)
    d <- 25
    x <- rnorm(d)
    prox1 <- prox_l2(x, 1)
    prox2 <- prox_group(x, list("1"=1), list("1"=1:d))
    expect_equal(prox1, prox2)
}
)

          
