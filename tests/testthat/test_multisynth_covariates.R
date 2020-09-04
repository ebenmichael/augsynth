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
                     scm = T, how_match = "exact")

  # check that there is actually no weight on donors with different Z
  expect_equal(sum(msyn$weights[fake_bin == 1, 1]), 1, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 1]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 1, 2]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 2]), 1, tolerance = 1e-6)


  # again with fixed effect
  msyn <- multisynth(gdpcap ~ trt | Z, regionno, year, basque2, nu = 0,
                     scm = T, fixedeff = T, how_match = "exact")
  # check that there is actually no weight on donors with different Z
  expect_equal(sum(msyn$weights[fake_bin == 1, 1]), 1, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 1]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 1, 2]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 2]), 1, tolerance = 1e-6)
})


test_that("K-NN finds the right number of neighbors", {

  # variables to match on
  Z <- matrix(rnorm(length(regions) * 3), ncol = 3)
  basque %>%
    inner_join(
      data.frame(regionno = regions,
                 Z1 = Z[, 1], Z2 = Z[, 2], Z3 = Z[, 3]),
      by = "regionno") -> basque2
  
  dat <- format_data_stag(quo(gdpcap), quo(trt), quo(regionno),
                          quo(year), basque2)
  k <- 3
  donors <- get_eligible_donors(dat$trt, F, 100)
  knn_donors <- get_knn_donors(dat$trt, Z, donors, k)
  expect_true(all(sapply(knn_donors, sum) == k))

  k <- 20
  expect_warning(get_knn_donors(dat$trt, Z, donors, k))
})

test_that("Getting eligible donor units by knn matching works", {

  # variables to match on
  Z <- matrix(rnorm(length(regions) * 3), ncol = 3)
  basque %>%
    inner_join(
      data.frame(regionno = regions,
                 Z1 = Z[, 1], Z2 = Z[, 2], Z3 = Z[, 3]),
      by = "regionno") -> basque2

  # error if no k is supplied
  expect_error(multisynth(gdpcap ~ trt | Z1 + Z2 + Z3, regionno, 
                          year, basque2,
                          scm = T, how_match = "knn"),
              "Number of neighbors for knn not selected, please choose k.")

  k <- 5
  msyn <- multisynth(gdpcap ~ trt | Z1 + Z2 + Z3, regionno, year, 
                     basque2, scm = T, how_match = "knn", k = k)

  # check that all but k units recieve exactly 0 weight
  expect_equal(sum(msyn$weights[, 1] != 0), k, tolerance = 1e-12)
  expect_equal(sum(msyn$weights[, 2] != 0), k, tolerance = 1e-12) 

  

  # again with fixed effect
    msyn <- multisynth(gdpcap ~ trt | Z1 + Z2 + Z3, regionno, year,
                       basque2, scm = T, fixedeff = T, how_match = "knn", k = k)
  # check that all but k units recieve exactly 0 weight
  expect_equal(sum(msyn$weights[, 1] != 0), k, tolerance = 1e-12)
  expect_equal(sum(msyn$weights[, 2] != 0), k, tolerance = 1e-12) 

  # without synth weights, weights are uniform
  k <- 2
  unimatch <- multisynth(gdpcap ~ trt | Z1 + Z2 + Z3, regionno, year, basque2,
                     scm = T, how_match = "knn", k = k, lambda = 1e10)

  expect_equal(unimatch$weights[unimatch$weights != 0 ], rep(1 / k, 2 * k))

  # matching with more neighbors is worse
  unimatch2 <- multisynth(gdpcap ~ trt | Z1 + Z2 + Z3, regionno, year, basque2,
                     scm = T, how_match = "knn", k = 2 * k, lambda = 1e10)

  trtZ <- Z[regions %in% c(16, 17),]
  imbal1 <- sqrt(sum(sapply(1:2, 
                function(i) sum(unimatch$weights[,i] * (trtZ[i,] - Z) ^ 2 ))))
  imbal2 <- sqrt(sum(sapply(1:2, 
                function(i) sum(unimatch2$weights[,i] * (trtZ[i,] - Z) ^ 2 ))))

  expect_lt(imbal1, imbal2)

})


test_that("Getting eligible donor units by exact and knn matching works", {

  # binary variable to split on
  fake_bin <- sample(c(0, 1), length(regions), replace = T)

  # variables to match on
  Z <- matrix(rnorm(length(regions) * 3), ncol = 3)
  basque %>%
    inner_join(
      data.frame(regionno = regions,
                 Z1 = Z[, 1], Z2 = Z[, 2], Z3 = Z[, 3],
                 Z_bin = fake_bin) %>%
        mutate(Z_bin = case_when(regionno == 17 ~ 0,
                             regionno == 16 ~ 1,
                             TRUE ~ Z_bin)),
      by = "regionno") -> basque2

  # error if no k is supplied
  expect_error(multisynth(gdpcap ~ trt | Z1 + Z2 + Z3 | Z_bin, regionno, 
                          year, basque2,
                          scm = T, how_match = "knn"),
              "Number of neighbors for knn not selected, please choose k.")

  k <- 3
  msyn <- multisynth(gdpcap ~ trt | Z1 + Z2 + Z3 | Z_bin, regionno, year, 
                     basque2, scm = T, how_match = "knn", k = k)
  
  # check that there is actually no weight on donors with different Z
  expect_equal(sum(msyn$weights[fake_bin == 1, 1]), 1, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 1]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 1, 2]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 2]), 1, tolerance = 1e-6)
  
  # check that all but k units recieve exactly 0 weight
  expect_equal(sum(msyn$weights[, 1] != 0), k, tolerance = 1e-12)
  expect_equal(sum(msyn$weights[, 2] != 0), k, tolerance = 1e-12) 
  
  # again with fixed effect
    msyn <- multisynth(gdpcap ~ trt | Z1 + Z2 + Z3 | Z_bin, regionno, year,
                       basque2, scm = T, fixedeff = T, how_match = "knn", k = k)
  # check that all but k units recieve exactly 0 weight
  expect_equal(sum(msyn$weights[, 1] != 0), k, tolerance = 1e-12)
  expect_equal(sum(msyn$weights[, 2] != 0), k, tolerance = 1e-12)

  # check that there is actually no weight on donors with different Z
  expect_equal(sum(msyn$weights[fake_bin == 1, 1]), 1, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 1]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 1, 2]), 0, tolerance = 1e-6)
  expect_equal(sum(msyn$weights[fake_bin == 0, 2]), 1, tolerance = 1e-6) 

  k <- 3
  # without synth weights, weights are uniform
  unimatch <- multisynth(gdpcap ~ trt | Z1 + Z2 + Z3 | Z_bin, regionno, year, basque2,
                     scm = T, how_match = "knn", k = k, lambda = 1e10)

  expect_equal(unimatch$weights[unimatch$weights != 0 ], rep(1 / k, 2 * k))

  # matching without exact gives better matches
  unimatch2 <- multisynth(gdpcap ~ trt | Z1 + Z2 + Z3, regionno, year, basque2,
                     scm = T, how_match = "knn", k = k, lambda = 1e10)

  trtZ <- Z[regions %in% c(16, 17),]
  imbal1 <- sqrt(sum(sapply(1:2, 
                function(i) sum(unimatch$weights[,i] * (trtZ[i,] - Z) ^ 2 ))))
  imbal2 <- sqrt(sum(sapply(1:2, 
                function(i) sum(unimatch2$weights[,i] * (trtZ[i,] - Z) ^ 2 ))))

  expect_lt(imbal2, imbal1)
})