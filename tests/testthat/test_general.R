context("Generally testing the workflow")


library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)


                            
test_that("SCM gives the right answer", {

    syn <- augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="None", weightfunc="SCM")

    ## average att estimate is as expected
    expect_equal(-.36, mean(summary(syn)$att$Estimate), tolerance=1e-4)

    ## average se estimate is as expected
    expect_equal(0.0945,
                 mean(summary(syn)$att$Std.Error, na.rm=T),
                 tolerance=1e-3)

    ## level of balance is as expected
    expect_equal(.377, syn$l2_imbalance, tolerance=1e-3)

}
)


test_that("Ridge ASCM gives the right answer", {

    asyn <- augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="Ridge",
                     weightfunc="SCM", lambda=8)

    ## average att estimate is as expected
    expect_equal(-.363, mean(summary(asyn)$att$Estimate), tolerance=1e-3)

    ## average se estimate is as expected
    expect_equal(0.162,
                 mean(summary(asyn)$att$Std.Error, na.rm=T),
                 tolerance=1e-3)

    ## level of balance is as expected
    expect_equal(.373, asyn$l2_imbalance, tolerance=1e-3)

}
)




test_that("Ridge ASCM with covariates gives the right answer", {

    covsyn <- augsynth(gdpcap ~ trt| invest + popdens, regionno, year, 
                        1975, basque, progfunc="None", weightfunc="SCM")

    ## average att estimate is as expected
    expect_equal(-.089, mean(summary(covsyn)$att$Estimate), tolerance=1e-3)

    ## average se estimate is as expected
    expect_equal(0.5403,
                 mean(summary(covsyn)$att$Std.Error, na.rm=T),
                 tolerance=1e-3)

    ## level of balance is as expected
    expect_equal(.376, covsyn$l2_imbalance, tolerance=1e-3)

}
)
