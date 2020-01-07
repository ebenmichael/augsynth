context("Test time cohort vs unit level analysis")

library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)


                            
test_that("multisynth at the unit level and time cohort level give the same answer for a single treated unit and no augmentation", {

    msyn_unit <- multisynth(gdpcap ~ trt, regionno, year, basque, nu = 0,
                            time_cohort = F, scm = T,
                            eps_rel = 1e-5, eps_abs = 1e-5)
    msyn_time <- multisynth(gdpcap ~ trt, regionno, year, basque, nu = 0,
                            time_cohort = T, scm = T,
                            eps_rel = 1e-5, eps_abs = 1e-5)
    # weights are the same-ish
    expect_equal(c(msyn_unit$weights), c(msyn_time$weights), tolerance=3e-2)

    # estimates are the same-ish
    expect_equal(c(predict(msyn_unit, att=F)),
                 c(predict(msyn_time, att=F)),
                 tolerance=5e-3)


    ## level of balance is same-ish expected
    expect_equal(msyn_unit$ind_l2, msyn_time$ind_l2, tolerance=1e-3)

}
)