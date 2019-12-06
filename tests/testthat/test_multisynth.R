context("Generally testing the workflow for multisynth")

library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)


                            
test_that("augsynth and multisynth give the same answer for a single treated unit and no augmentation", {

    syn <- augsynth(gdpcap ~ trt, regionno, year, 1975, basque,
                    progfunc="None", scm=T)
    msyn <- multisynth(gdpcap ~ trt, regionno, year, basque, nu = 0,
                       scm=T, eps_rel=1e-5, eps_abs=1e-5)
    
    # weights are the same-ish
    expect_equal(c(syn$weights), c(msyn$weights[-16]), tolerance=3e-2)

    # estimates are the same-ish
    pred_msyn <- predict(msyn, att=F)[,1]
    pred_msyn <- pred_msyn[-length(pred_msyn)]
    expect_equal(c(predict(syn, att=F)), pred_msyn, tolerance=5e-3)


    ## level of balance is same-ish expected
    expect_equal(syn$l2_imbalance, msyn$ind_l2, tolerance=1e-3)

}
)


test_that("Pooling doesn't matter for a single treated unit", {

    nopool <- multisynth(gdpcap ~ trt, regionno, year, basque, nu = 0,
                         scm=T, eps_rel=1e-5, eps_abs=1e-5)
    allpool <- multisynth(gdpcap ~ trt, regionno, year, basque, nu = 1,
                          scm=T, eps_rel=1e-5, eps_abs=1e-5)

    # weights are the same
    expect_equal(nopool$weights, allpool$weights)

    # estimates are the same
    expect_equal(predict(nopool), predict(allpool))


    ## level of balance is same-ish expected
    expect_equal(allpool$ind_l2, nopool$ind_l2)

}
)





                            
test_that("Separate synth is the same as fitting separate synths", {


    basque2 <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                                !regionno %in% c(16, 17) ~ 0,
                                                regionno %in% c(16, 17) ~ 1)) %>%
        filter(regionno != 1)


    basque2  %>% filter(regionno != 16) %>% 
        augsynth(gdpcap ~ trt, regionno, year, 1975, .,
                    progfunc="None", scm=T) -> scm17
    basque2  %>% filter(regionno != 17) %>% 
        augsynth(gdpcap ~ trt, regionno, year, 1975, .,
                    progfunc="None", scm=T) -> scm16
    
    msyn <- multisynth(gdpcap ~ trt, regionno, year, basque2, nu = 0,
                       scm=T, eps_rel=1e-5, eps_abs=1e-5)
    
    # weights are the same-ish
    expect_equal(c(scm17$weights), c(msyn$weights[-c(15, 16), 2]), tolerance=3e-2)
    # expect_equal(c(scm16$weights), c(msyn$weights[-c(15, 16), 1]), tolerance=3e-2)
    
    # estimates are the same-ish
    pred_msyn <- predict(msyn, att=F)
    pred_msyn <- pred_msyn[-nrow(pred_msyn), ]
    expect_equal(c(predict(scm17, att=F)), pred_msyn[, 3], tolerance=5e-3)
    expect_equal(c(predict(scm16, att=F)), pred_msyn[, 2], tolerance=5e-3)
}
)

