context("Testing multisynth with covariates")
set.seed(1011)

library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when((regionno == 17) & (year == 1975) ~ 1,
                                              (regionno == 16) & (year == 1980) ~ 1,
                                              TRUE ~ 0)) %>%
      filter(regionno != 1)


regions <- basque %>% distinct(regionno) %>% pull(regionno)

test_that("Getting eligible donor units by exact matching works", {

  # binary variable to split on
  fake_bin <- sample(c(0, 1), length(regions), replace = T)
  basque %>%
    inner_join(
      data.frame(regionno = regions, Z = fake_bin) %>%
        mutate(Z = case_when(regionno == 17 ~ 0,
                             regionno == 16 ~ 1,
                             TRUE ~ Z)
              ),
               by = "regionno") -> basque2

  msyn <- multisynth(gdpcap ~ trt | Z, regionno, year, basque2, nu = 0,
                     scm = T)

  # check that there is actually no weight on donors with different Z
  expect_equal(sum(msyn$weights[fake_bin == 1, 1]), 1, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 1]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 1, 2]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 2]), 1, tolerance = 1e-6)


  # again with fixed effect
  msyn <- multisynth(gdpcap ~ trt | Z, regionno, year, basque2, nu = 0,
                     scm = T, fixedeff = T)
  # check that there is actually no weight on donors with different Z
  expect_equal(sum(msyn$weights[fake_bin == 1, 1]), 1, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 1]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 1, 2]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 2]), 1, tolerance = 1e-6)
})