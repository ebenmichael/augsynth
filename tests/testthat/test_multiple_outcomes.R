context("Generally testing the workflow for synth with multiple outcomes")

library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1),
                            gdpcap_sq = gdpcap ^ 2) %>%
    filter(regionno != 1)

test_that("augsynth and augsynth_multiout are the same without augmentation", {

    syn1 <- augsynth_multiout(gdpcap + gdpcap_sq  ~ trt, regionno, year, 1975, basque,
                    progfunc="None", scm=T)
    syn2 <- augsynth(gdpcap + gdpcap_sq  ~ trt, regionno, year, basque,
                    progfunc="None", scm=T)

    # weights are the same
    expect_equal(c(syn1$weights), c(syn2$weights), tolerance=3e-4)

    # estimates are the same
    expect_equal(c(predict(syn1, att=F)), c(predict(syn2, att = F)), tolerance=5e-5)


    ## level of balance is same
    expect_equal(syn1$l2_imbalance, syn2$l2_imbalance, tolerance=1e-5)        
})

test_that("augsynth and augsynth_multiout are the same with ridge augmentation", {

    syn1 <- augsynth_multiout(gdpcap + gdpcap_sq  ~ trt, regionno, year, 1975, basque,
                    progfunc="Ridge", scm=T, lambda = 10)
    syn2 <- augsynth(gdpcap + gdpcap_sq  ~ trt, regionno, year, basque,
                    progfunc="Ridge", scm=T, lambda = 10)

    # weights are the same
    expect_equal(c(syn1$weights), c(syn2$weights), tolerance=3e-4)

    # estimates are the same
    expect_equal(c(predict(syn1, att=F)), c(predict(syn2, att = F)), tolerance=5e-5)


    ## level of balance is same
    expect_equal(syn1$l2_imbalance, syn2$l2_imbalance, tolerance=1e-5)        
})

test_that("augsynth and augsynth_multiout are the same with fixed effects augmentation", {

    syn1 <- augsynth_multiout(gdpcap + gdpcap_sq  ~ trt, regionno, year, 1975, basque,
                    progfunc="None", scm=T, fixedeff = T)
    syn2 <- augsynth(gdpcap + gdpcap_sq  ~ trt, regionno, year, basque,
                    progfunc="None", scm=T, fixedeff = T)

    # weights are the same
    expect_equal(c(syn1$weights), c(syn2$weights), tolerance=3e-4)

    # estimates are the same
    expect_equal(c(predict(syn1, att=F)), c(predict(syn2, att = F)), tolerance=5e-5)


    ## level of balance is same
    expect_equal(syn1$l2_imbalance, syn2$l2_imbalance, tolerance=1e-5)
})

test_that("single_augsynth and augsynth_multiout are the same for one outcome", {
    syn1 <- augsynth_multiout(gdpcap  ~ trt, regionno, year, 1975, basque,
                    progfunc="None", scm=T, combine_method = "concat")
    syn2 <- augsynth(gdpcap  ~ trt, regionno, year, basque,
                    progfunc="None", scm=T)

    # weights are the same
    expect_equal(c(syn1$weights), c(syn2$weights), tolerance=3e-4)

    # estimates are the same
    expect_equal(c(predict(syn1, att=F)), unname(predict(syn2, att = F)), tolerance=5e-5)


    ## level of balance is same
    expect_equal(syn1$l2_imbalance, syn2$l2_imbalance, tolerance=1e-5)            
})

test_that("single_augsynth and augsynth_multiout are the same for one outcome with ridge augmentation",{
    syn1 <- augsynth_multiout(gdpcap  ~ trt, regionno, year, 1975, basque,
                    progfunc="Ridge", scm=T, combine_method = "concat")
    syn2 <- augsynth(gdpcap  ~ trt, regionno, year, basque,
                    progfunc="Ridge", scm=T)

    # weights are the same
    expect_equal(c(syn1$weights), c(syn2$weights), tolerance=3e-4)

    # estimates are the same
    expect_equal(c(predict(syn1, att=F)), unname(predict(syn2, att = F)),
                 tolerance=5e-5)


    ## level of balance is same
    expect_equal(syn1$l2_imbalance, syn2$l2_imbalance, tolerance=1e-5)            
})

test_that("single_augsynth and augsynth_multiout are the same for one outcome with fixed effect augmentation", {
    syn1 <- augsynth_multiout(gdpcap  ~ trt, regionno, year, 1975, basque,
                    progfunc="None", scm=T, fixedeff = T, combine_method = "concat")
    syn2 <- augsynth(gdpcap  ~ trt, regionno, year, basque,
                    progfunc="None", scm=T, fixedeff = T)

    # weights are the same
    expect_equal(c(syn1$weights), c(syn2$weights), tolerance=3e-4)

    # estimates are the same
    expect_equal(c(predict(syn1, att=F)), unname(predict(syn2, att = F)), tolerance=5e-5)


    ## level of balance is same
    expect_equal(syn1$l2_imbalance, syn2$l2_imbalance, tolerance=1e-5)            
})


test_that("Averaging outcomes with augsynth_multiout gives correct results without fixed effects", {

  sds <- basque %>% filter(trt == 0, year < 1975) %>%
          summarise(across(c(gdpcap, gdpcap_sq), sd)) %>%
          rename(gdpcap_sd = gdpcap, gdpcap_sq_sd = gdpcap_sq)

  basque %>%
    bind_cols(sds) %>%
    mutate(avg = gdpcap / gdpcap_sd + gdpcap_sq / gdpcap_sq_sd,
           avg2 = gdpcap + gdpcap_sq) -> bas_avg
  

  syn1 <- augsynth_multiout(gdpcap + gdpcap_sq ~ trt, regionno, year, 1975, basque,
                            progfunc="None", scm=T, fixedeff = F, combine_method = "avg")
  syn2 <- augsynth(avg  ~ trt, regionno, year, bas_avg,
                    progfunc="None", scm=T, fixedeff = F)

  # weights are the same
  expect_equal(c(syn1$weights), c(syn2$weights), tolerance=3e-4)
})

test_that("Averaging outcomes with augsynth_multiout gives correct results with fixed effects", {

  sds <- basque %>% filter(trt == 0, year < 1975) %>%
          summarise(across(c(gdpcap, gdpcap_sq), sd)) %>%
          rename(gdpcap_sd = gdpcap, gdpcap_sq_sd = gdpcap_sq)

  basque %>%
    bind_cols(sds) %>%
    mutate(avg = gdpcap / gdpcap_sd + gdpcap_sq / gdpcap_sq_sd,
           avg2 = gdpcap + gdpcap_sq) -> bas_avg
  

  syn1 <- augsynth_multiout(gdpcap + gdpcap_sq ~ trt, regionno, year, 1975, basque,
                            progfunc="None", scm=T, fixedeff = T,
                            combine_method = "avg")
  syn2 <- augsynth(avg  ~ trt, regionno, year, bas_avg,
                    progfunc="None", scm=T, fixedeff = T)

  # weights are the same
  expect_equal(c(syn1$weights), c(syn2$weights), tolerance=1e-3)
})



test_that("Concatenating outcomes with augsynth_multiout gives correct results without fixed effects", {


  sds <- basque %>% filter(trt == 0, year < 1975) %>%
          summarise(across(c(gdpcap, gdpcap_sq), sd)) %>%
          rename(gdpcap_sd = gdpcap, gdpcap_sq_sd = gdpcap_sq)

  basque %>%
    bind_cols(sds) %>%
    mutate(gdpcap = gdpcap / gdpcap_sd, gdpcap_sq = gdpcap_sq / gdpcap_sq_sd) %>%
    select(gdpcap, gdpcap_sq, trt, year, regionno) %>%
    pivot_longer(-c(regionno, year, trt)) %>%
    mutate(year = ifelse(name == "gdpcap", year, year - 0.5)) -> bas_cat  

  syn1 <- augsynth_multiout(gdpcap + gdpcap_sq ~ trt, regionno, year, 1975, basque,
                            progfunc="None", scm=T, fixedeff = F,
                            combine_method = "concat")
  syn2 <- augsynth(value  ~ trt, regionno, year, bas_cat,
                    progfunc="None", scm=T, fixedeff = F)

  # weights are the same
  expect_equal(c(syn1$weights), c(syn2$weights), tolerance=5e-4)
})

