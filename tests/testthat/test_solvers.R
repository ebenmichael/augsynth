context("Synth QP solvers agree (osqp vs frank_wolfe)")

library(tidyverse)

data(basque, package = "Synth")
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~ 0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)


test_that("solvers agree on the QP objective directly", {

    set.seed(42)
    n0 <- 40
    t0 <- 15
    X0 <- matrix(rnorm(n0 * t0), n0, t0)
    X1 <- matrix(rnorm(t0), t0, 1)
    V <- diag(t0)

    w_qp <- augsynth:::synth_qp_osqp(X1, X0, V)
    w_fw <- augsynth:::synth_qp_frank_wolfe(X1, X0, V)

    objective <- function(w) sum((t(X0) %*% w - X1)^2)

    ## frank_wolfe respects the simplex and reaches the osqp objective
    expect_equal(sum(w_fw), 1, tolerance = 1e-10)
    expect_true(all(w_fw >= -1e-10))
    expect_lt(objective(w_fw), objective(w_qp) * (1 + 1e-4) + 1e-10)
})


test_that("frank_wolfe matches osqp without ridge augmentation", {

    syn_qp <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "None", scm = TRUE, t_int = 1975)
    syn_fw <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "None", scm = TRUE, t_int = 1975,
                              solver = "frank_wolfe")

    expect_equal(sum(syn_fw$weights), 1, tolerance = 1e-8)
    expect_true(all(syn_fw$weights >= -1e-8))

    ## same imbalance (the objective) and the same ATT path
    expect_equal(syn_fw$l2_imbalance, syn_qp$l2_imbalance, tolerance = 1e-4)
    expect_lt(max(abs(predict(syn_fw, att = TRUE) -
                        predict(syn_qp, att = TRUE))), 1e-2)
})


test_that("frank_wolfe matches osqp with ridge augmentation (fixed lambda)", {

    syn_qp <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "Ridge", scm = TRUE, t_int = 1975,
                              lambda = 8)
    syn_fw <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "Ridge", scm = TRUE, t_int = 1975,
                              lambda = 8, solver = "frank_wolfe")

    expect_equal(syn_fw$lambda, syn_qp$lambda)
    expect_lt(max(abs(syn_fw$weights - syn_qp$weights)), 1e-2)
    expect_lt(max(abs(predict(syn_fw, att = TRUE) -
                        predict(syn_qp, att = TRUE))), 1e-2)
})


test_that("frank_wolfe matches osqp with ridge augmentation (CV lambda)", {

    syn_qp <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "Ridge", scm = TRUE, t_int = 1975)
    syn_fw <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "Ridge", scm = TRUE, t_int = 1975,
                              solver = "frank_wolfe")

    ## the CV refits run through the requested solver; the selected lambdas
    ## (and so the estimates) should agree on this problem
    expect_equal(syn_fw$lambda, syn_qp$lambda, tolerance = 1e-8)
    expect_lt(max(abs(predict(syn_fw, att = TRUE) -
                        predict(syn_qp, att = TRUE))), 1e-2)
})


test_that("frank_wolfe matches osqp with a fixed effect", {

    syn_qp <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "Ridge", scm = TRUE, t_int = 1975,
                              fixedeff = TRUE, lambda = 8)
    syn_fw <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "Ridge", scm = TRUE, t_int = 1975,
                              fixedeff = TRUE, lambda = 8,
                              solver = "frank_wolfe")

    expect_lt(max(abs(predict(syn_fw, att = TRUE) -
                        predict(syn_qp, att = TRUE))), 1e-2)
})


test_that("a custom solver function can be plugged in", {

    ## passing the osqp solver as a bare function must reproduce the default
    syn_default <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                                   progfunc = "None", scm = TRUE, t_int = 1975)
    syn_custom <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                                  progfunc = "None", scm = TRUE, t_int = 1975,
                                  solver = function(X1, X0, V) {
                                      augsynth:::synth_qp_osqp(X1, X0, V)
                                  })

    expect_equal(syn_custom$weights, syn_default$weights, tolerance = 1e-10)
})


test_that("the solver is remembered by inference refits", {

    syn_fw <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc = "Ridge", scm = TRUE, t_int = 1975,
                              lambda = 8, solver = "frank_wolfe")
    expect_equal(syn_fw$extra_args$solver, "frank_wolfe")

    ## summary (jackknife+ / conformal machinery) replays extra_args
    expect_error(summary(syn_fw, inf_type = "jackknife+"), NA)
})
