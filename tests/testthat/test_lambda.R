context("Testing lambda tuning if ridge is true.")


library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
  filter(regionno != 1)


test_that("Lambda sequence is generated correctly", {
  syn <- augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="Ridge", scm=T)
  lambdas <- syn$lambdas
  expect_equivalent(lambdas[length(lambdas)] / lambdas[1], 1e-8)
  expect_equivalent(lambdas[2] / lambdas[1], lambdas[3] / lambdas[2])
})

test_that("Smallest lambda is chosen", {
  syn <- augsynth(gdpcap ~ trt, regionno, year, 1975, basque, 
                  progfunc="Ridge", scm=T, 
                  opts_out=list(min_1se = F))
  expect_equivalent(syn$lambda, syn$lambdas[which.min(syn$lambda_errors)])
})


test_that("Largest lambda within 1 SE of minimum is chosen", {
  syn <- augsynth(gdpcap ~ trt, regionno, year, 1975, basque, 
                  progfunc="Ridge", scm=T, 
                  opts_out=list(min_1se = T))
  min_idx <- which.min(syn$lambda_errors)
  min_1se <- max(syn$lambdas[syn$lambda_errors <= 
                              syn$lambda_errors[min_idx] + 
                              syn$lambda_errors_se[min_idx]])
  expect_equivalent(syn$lambda, min_1se)
})