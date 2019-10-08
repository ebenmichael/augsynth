context("Testing lambda tuning if ridge is true.")


library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
  filter(regionno != 1)


test_that("Lambda sequence is generated correctly", {
  syn <- augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="Ridge", weightfunc="SCM")
  lambdas <- syn$lambdas
  expect_equivalent(lambdas[length(lambdas)] / lambdas[1], 1/1000)
  expect_equivalent(lambdas[2] / lambdas[1], lambdas[3] / lambdas[2])
})

test_that("Smallest lambda is chosen", {
  syn <- augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="Ridge", weightfunc="SCM")
  expect_equivalent(syn$lambda, min(syn$lambdas))
})
