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




test_that("Ridge ASCM with covariates gives the right answer", {

    covsyn <- augsynth(gdpcap ~ trt | invest + popdens,
                       regionno, year, basque,
                       progfunc = "None", scm = T)

    ## average att estimate is as expected
    expect_equal(-.1443,
                 mean(summary(covsyn, inf = F)$att$Estimate),
                 tolerance = 1e-3)

    ## level of balance is as expected
    expect_equal(.3720, covsyn$l2_imbalance, tolerance=1e-3)

}
)

test_that("Warning given when inputting an unused argument", {

    expect_warning(
      augsynth(gdpcap ~ trt| invest + popdens, regionno, year, basque, 
               progfunc="Ridge", scm=T, lambda=8, t_int = 1975, 
               bad_param = "Unused input parameter"),
    )
})
