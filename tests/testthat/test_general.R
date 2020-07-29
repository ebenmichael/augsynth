context("Generally testing the workflow for augsynth")


library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)


                            
test_that("SCM gives the right answer", {

    syn <- augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="None", scm=T, t_int=1975)
    ## average att estimate is as expected
    expect_equal(-.3686, mean(summary(syn, inf = F)$att$Estimate), tolerance=1e-4)


    ## level of balance is as expected
    expect_equal(.377, syn$l2_imbalance, tolerance=1e-3)

}
)

test_that("SCM finds the correct t_int and gives the right answer", {

    syn1 <- augsynth(gdpcap ~ trt, regionno, year, basque,
                     progfunc="None", scm=T)
    syn2 <- augsynth(gdpcap ~ trt, regionno, year, basque,
                     progfunc = "None", scm = T, t_int = 1975)
    ## average att estimate is as expected
    expect_equal(mean(summary(syn1, inf = F)$att$Estimate), 
                 mean(summary(syn2, inf = F)$att$Estimate), tolerance=1e-4)
    
    ## level of balance is as expected
    expect_equal(syn1$l2_imbalance, syn2$l2_imbalance, tolerance=1e-3)
    
}
)


test_that("Ridge ASCM gives the right answer", {

    asyn <- augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="Ridge",
                     scm=T, lambda=8)

    ## average att estimate is as expected
    expect_equal(-.3696, mean(summary(asyn, inf = F)$att$Estimate), tolerance=1e-3)

    ## level of balance is as expected
    expect_equal(.373, asyn$l2_imbalance, tolerance=1e-3)

}
)




test_that("SCM after residualizing covariates gives the right answer", {

  covsyn_resid <- augsynth(gdpcap ~ trt | invest + popdens,
                      regionno, year, basque,
                      progfunc = "None", scm = T,
                      residualize = T)

  ## average att estimate is as expected
  expect_equal(-.1443,
                mean(summary(covsyn_resid, inf = F)$att$Estimate),
                tolerance = 1e-3)

  ## level of balance is as expected
  expect_equal(.3720, covsyn_resid$l2_imbalance, tolerance=1e-3)

  # perfect auxiliary covariate balance
  expect_equal(0, covsyn_resid$covariate_l2_imbalance, tolerance=1e-3)

}
)

test_that("Ridge ASCM with covariates jointly gives the right answer", {

    covsyn_noresid <- augsynth(gdpcap ~ trt | invest + popdens,
                       regionno, year, basque,
                       progfunc = "None", scm = T,
                       residualize = F)

    ## average att estimate is as expected
    expect_equal(-.3345,
                 mean(summary(covsyn_noresid, inf = F)$att$Estimate),
                 tolerance = 1e-3)

    ## level of balance is as expected
    expect_equal(0.659, covsyn_noresid$l2_imbalance, tolerance=1e-3)

    # covariate balance is as expected
    expect_equal(0.884, covsyn_noresid$covariate_l2_imbalance, tolerance=1e-3)


}
)


test_that("Ridge ASCM after residualizing covariates gives the right answer", {

    covascm_resid <- augsynth(gdpcap ~ trt | invest + popdens,
                       regionno, year, basque,
                       progfunc = "Ridge", scm = T,
                       lambda = 1,
                       residualize = T)

    ## average att estimate is as expected
    expect_equal(-.123,
                 mean(summary(covascm_resid, inf = F)$att$Estimate),
                 tolerance = 1e-3)

    ## level of balance is as expected
    expect_equal(.347, covascm_resid$l2_imbalance, tolerance=1e-3)

    # perfect auxiliary covariate balance
    expect_equal(0, covascm_resid$covariate_l2_imbalance, tolerance=1e-3)

}
)

test_that("Ridge ASCM with covariates jointly gives the right answer", {

    covascm_noresid <- augsynth(gdpcap ~ trt | invest + popdens,
                       regionno, year, basque,
                       progfunc = "Ridge", scm = T,
                       lambda = 1,
                       residualize = F)

    ## average att estimate is as expected
    expect_equal(-.267,
                 mean(summary(covascm_noresid, inf = F)$att$Estimate),
                 tolerance = 1e-3)

    ## level of balance is as expected
    expect_equal(0.419, covascm_noresid$l2_imbalance, tolerance=1e-3)

    # covariate balance is as expected
    expect_equal(0.084, covascm_noresid$covariate_l2_imbalance, tolerance=1e-3)

}
)


test_that("Warning given when inputting an unused argument", {

    expect_warning(
      augsynth(gdpcap ~ trt| invest + popdens, regionno, year, basque, 
               progfunc="Ridge", scm=T, lambda=8, t_int = 1975, 
               bad_param = "Unused input parameter"),
    )
})
