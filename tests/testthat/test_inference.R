


context("Test inference features of single augsynth objects")

basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
  filter(regionno != 1)

test_that("Default model doesn't contain inference", {
  syn <- augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="Ridge", scm = T)
  expect_equivalent(is.null(syn$results), TRUE)
})

