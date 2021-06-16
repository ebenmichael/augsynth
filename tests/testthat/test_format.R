context("Test data formatting")

library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)
                            
test_that("format_data creates matrices with the right dimensions", {
    
    dat <- format_data(quo(gdpcap), quo(trt), quo(regionno), quo(year),1975, basque)

    test_dim <- function(obj, d) {
        expect_equivalent(dim(obj), d)
        }

    test_dim(dat$X, c(17, 20))
    expect_equivalent(length(dat$trt), 17)
    test_dim(dat$y, c(17, 23))
}
)


test_that("format_synth creates matrices with the right dimensions", {
    
    dat <- format_data(quo(gdpcap), quo(trt), quo(regionno), quo(year),1975, basque)
    syn_dat <- format_synth(dat$X, dat$trt, dat$y)
    test_dim <- function(obj, d) {
        expect_equivalent(dim(obj), d)
        }

    test_dim(syn_dat$Z0, c(20, 16))
    test_dim(syn_dat$Z1, c(20, 1))

    test_dim(syn_dat$Y0plot, c(43, 16))
    test_dim(syn_dat$Y1plot, c(43, 1))

    expect_equivalent(syn_dat$Z1, syn_dat$X1)
    expect_equivalent(syn_dat$Z0, syn_dat$X0)
}
)


test_that("multisynth throws errors when there aren't enough pre-treatment times",
  {
    basque2 <- basque %>%
      mutate(trt = case_when(
        regionno == 16 ~ 1,
        year >= 1975 & regionno == 17 ~ 1,
        TRUE ~ 0)
        ) %>%
      filter(regionno != 1)

  # error from always treated unit
  expect_warning(
    expect_error(
      multisynth(gdpcap ~ trt, regionno, year, basque2)
    )
  )

  basque2 <- basque %>%
      mutate(trt = case_when(
        regionno == 16 & year >= 1956 ~ 1,
        year >= 1975 & regionno == 17 ~ 1,
        TRUE ~ 0)
        ) %>%
      filter(regionno != 1)

  # error from one pre-treatment outcome and fixedeff = T
  expect_warning(
    expect_error(multisynth(gdpcap ~ trt, regionno, year, basque2))
  )

  # no error from one pre-treatment outcome and fixedeff = F, just warning
  expect_warning(multisynth(gdpcap ~ trt, regionno, year, basque2, fixedeff = F))

  })


  test_that("formatting for staggered adoption doesn't care about order of time in data",
  {
    basque2 <- basque %>%
      # slice(sample(1:n())) %>%
      mutate(trt = case_when((regionno == 17) & (year >= 1975) ~ 1,
                              (regionno == 16) & (year >= 1980) ~ 1,
                                TRUE ~ 0))

      dat <- format_data_stag(quo(gdpcap), quo(trt), quo(regionno), quo(year), basque2)

      # true treatment times
      true_trt <- c(1975, 1980) - min(basque$year)

      expect_equal(true_trt, sort(dat$trt[is.finite(dat$trt)]))

    basque2 <- basque %>%
        slice(sample(1:n())) %>%
        mutate(trt = case_when((regionno == 17) & (year >= 1975) ~ 1,
                                (regionno == 16) & (year >= 1980) ~ 1,
                                  TRUE ~ 0))

    dat <- format_data_stag(quo(gdpcap), quo(trt), quo(regionno), quo(year), basque2)

    expect_equal(true_trt, sort(dat$trt[is.finite(dat$trt)]))

  })