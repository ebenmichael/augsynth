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
                    progfunc="None", scm=T)
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
                    progfunc="Ridge", scm=T)
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
                    progfunc="None", scm=T, fixedeff = T)
    syn2 <- augsynth(gdpcap  ~ trt, regionno, year, basque,
                    progfunc="None", scm=T, fixedeff = T)

    # weights are the same
    expect_equal(c(syn1$weights), c(syn2$weights), tolerance=3e-4)

    # estimates are the same
    expect_equal(c(predict(syn1, att=F)), unname(predict(syn2, att = F)), tolerance=5e-5)


    ## level of balance is same
    expect_equal(syn1$l2_imbalance, syn2$l2_imbalance, tolerance=1e-5)            
})
