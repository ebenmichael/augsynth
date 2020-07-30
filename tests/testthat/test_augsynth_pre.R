context("Testing that top level API runs the right functions")

library(Synth)


test_that("augsynth runs single_synth when there is a single treated unit", {

  data(basque)
  basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                              regionno != 17 ~0,
                                              regionno == 17 ~ 1)) %>%
      filter(regionno != 1)

  syn <- augsynth(gdpcap ~ trt, regionno, year, basque,
                  progfunc = "None", scm = T, t_int = 1975)

  syn_single <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                                progfunc = "None", scm = T, t_int = 1975)

  expect_equal(syn$weights, syn_single$weights)
})


test_that("augsynth finds the treated time when is a single treated unit", {

  data(basque)
  basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                              regionno != 17 ~0,
                                              regionno == 17 ~ 1)) %>%
      filter(regionno != 1)

  syn <- augsynth(gdpcap ~ trt, regionno, year, basque,
                  progfunc = "None", scm = T, t_int = 1975)

  syn2 <- augsynth(gdpcap ~ trt, regionno, year, basque,
                  progfunc = "None", scm = T)

  expect_equal(syn$weights, syn2$weights)
})


test_that("augsynth runs single_synth when there is simultaneous adoption", {

  data(basque)
  basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                              !regionno %in% c(16, 17) ~0,
                                              regionno %in% c(16, 17) ~ 1)) %>%
      filter(regionno != 1)

  syn <- augsynth(gdpcap ~ trt, regionno, year, basque,
                  progfunc = "None", scm = T, t_int = 1975)

  syn_single <- single_augsynth(gdpcap ~ trt, regionno, year, basque,
                                progfunc = "None", scm = T, t_int = 1975)

  expect_equal(syn$weights, syn_single$weights)
})


test_that("augsynth runs multisynth when there is staggered adoption", {

  data(basque)
  basque <- basque %>% mutate(trt = case_when((regionno == 17) & (year == 1975) ~ 1,
                                              (regionno == 16) & (year == 1980) ~ 1,
                                              TRUE ~ 0)) %>%
      filter(regionno != 1)

  syn <- augsynth(gdpcap ~ trt, regionno, year, basque,
                  progfunc = "None", scm = T, t_int = 1975)

  syn_multi <- multisynth(gdpcap ~ trt, regionno, year, basque)

  expect_equal(syn$weights, syn_multi$weights)
})