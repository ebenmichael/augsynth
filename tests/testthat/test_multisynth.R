context("Generally testing the workflow for multisynth")

library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)


                            
test_that("augsynth and multisynth give the same answer for a single treated unit and no augmentation", {

    syn <- single_augsynth(gdpcap ~ trt, regionno, year, 1975, basque,
                    progfunc="None", scm=T, fixedeff = F)
    msyn <- multisynth(gdpcap ~ trt, regionno, year, basque, nu = 0,
                       fixedeff = F, scm=T, eps_rel=1e-5, eps_abs=1e-5)
    
    # weights are the same-ish
    expect_equal(c(syn$weights), c(msyn$weights[-16]), tolerance=3e-4)

    # estimates are the same-ish
    pred_msyn <- predict(msyn, att=F)[,1]
    pred_msyn <- pred_msyn[-length(pred_msyn)]
    expect_equal(unname(predict(syn, att=F)), pred_msyn, tolerance=5e-5)


    ## level of balance is same-ish expected
    expect_equal(syn$l2_imbalance, msyn$avg_l2, tolerance=1e-5)

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
        single_augsynth(gdpcap ~ trt, regionno, year, 1975, .,
                    progfunc="None", scm=T) -> scm17
    basque2  %>% filter(regionno != 17) %>% 
        single_augsynth(gdpcap ~ trt, regionno, year, 1975, .,
                    progfunc="None", scm=T) -> scm16
    
    msyn <- multisynth(gdpcap ~ trt, regionno, year, basque2, nu = 0,
                       scm=T, eps_rel=1e-5, eps_abs=1e-5, fixedeff = F)
    
    # weights are the same-ish
    sscm_weights <- unname(c(scm17$weights))
    mscm_weights <- unname(c(msyn$weights[-c(15, 16), 2]))
    expect_equal(sscm_weights, mscm_weights, tolerance=3e-2)
    expect_equal(rownames(scm17$weights), rownames(as.matrix(msyn$weights[-c(15, 16), 2])))
    # expect_equal(c(scm16$weights), c(msyn$weights[-c(15, 16), 1]), tolerance=3e-2)
    
    # estimates are the same-ish
    pred_msyn <- predict(msyn, att=F)
    pred_msyn <- pred_msyn[-nrow(pred_msyn), ]
    expect_equal(unname(predict(scm17, att=F)), pred_msyn[, 3], tolerance=5e-3)
    expect_equal(unname(predict(scm16, att=F)), pred_msyn[, 2], tolerance=5e-3)
}
)

test_that("Limiting number of lags works", {


    basque2 <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                                !regionno %in% c(16, 17) ~ 0,
                                                regionno %in% c(16, 17) ~ 1)) %>%
        filter(regionno != 1)

    expect_error(
      multisynth(gdpcap ~ trt, regionno, year, basque2, nu = 0,
                 scm=T, eps_rel=1e-5, eps_abs=1e-5, n_lags =3),
      NA
    )
}
)

test_that("L2 imbalance computed correctly", {

  basque2 <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                              !regionno %in% c(16, 17) ~ 0,
                                              regionno %in% c(16, 17) ~ 1)) %>%
      filter(regionno != 1)

  msyn <- multisynth(gdpcap ~ trt, regionno, year, basque2,
                scm=T, eps_rel=1e-5, eps_abs=1e-5)

  glbl <- sqrt(mean(msyn$imbalance[,1]^2))
  ind <- sqrt(mean(
    apply(msyn$imbalance[, -1], 2,
          function(x) sum(x ^ 2) / sum(x != 0))))
  avg_ind <- mean(apply(msyn$imbalance[,-1, drop = F], 2,
              function(x) sqrt(sum(x ^ 2))))
  expect_equal(glbl, msyn$global_l2)
  expect_equal(avg_ind, msyn$avg_l2)
  expect_equal(ind, msyn$ind_l2)
})

test_that("V matrix is equivalent to hard thresholding", {


  basque2 <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                              !regionno %in% c(16, 17) ~ 0,
                                              regionno %in% c(16, 17) ~ 1)) %>%
      filter(regionno != 1)

  V <- c(numeric(10), rep(1,5))
  msyn1 <- multisynth(gdpcap ~ trt, regionno, year, basque2, nu = 0,
                scm=T, eps_rel=1e-8, eps_abs=1e-8, n_lags = 15, V = V)

  msyn2 <- multisynth(gdpcap ~ trt, regionno, year, basque2, nu = 0,
                scm=T, eps_rel=1e-8, eps_abs=1e-8, n_lags = 5)

  expect_equal(msyn1$weights, msyn2$weights, tolerance = 1e-5)
  expect_equal(msyn1$global_l2, msyn2$global_l2, tolerance = 1e-5)
  expect_equal(msyn1$avg_l2, msyn2$avg_l2, tolerance = 1e-5)
}
)

test_that("V matrix is the same for single and multi synth", {

  V <- exp(seq(log(1e-3), log(1), length.out = 20))

  syn <- augsynth(gdpcap ~ trt, regionno, year, basque, progfunc = "none",
                scm=T, V = V)

  msyn <- multisynth(gdpcap ~ trt, regionno, year, basque,
                scm=T, eps_rel=1e-8, eps_abs=1e-8, V = V,
                fixed = F, nu = 0)

  expect_equal(as.numeric(syn$weights), as.numeric(msyn$weights[-16, ]), tolerance = 1e-3)
}
)


                            
test_that("multisynth doesn't depend on unit order", {

    basque2 <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                                !regionno %in% c(16, 17) ~ 0,
                                                regionno %in% c(16, 17) ~ 1)) %>%
        filter(regionno != 1)

    msyn <- multisynth(gdpcap ~ trt, regionno, year, basque2, nu = 0,
                       fixedeff = F, scm=T, eps_rel=1e-5, eps_abs=1e-5)

    msyn2 <- multisynth(gdpcap ~ trt, regionno, year,
                       basque2 %>% arrange(desc(regionno)), nu = 0,
                       fixedeff = F, scm=T, eps_rel=1e-5, eps_abs=1e-5)

    
    # weights are the same
    expect_equal(c(msyn$weights), c(msyn2$weights))

    # estimates are the same
    expect_equal(predict(msyn), predict(msyn2))

}
)


                            
test_that("multisynth doesn't depend on time order", {

    basque2 <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                                !regionno %in% c(16, 17) ~ 0,
                                                regionno %in% c(16, 17) ~ 1)) %>%
        filter(regionno != 1)

    msyn <- multisynth(gdpcap ~ trt, regionno, year, basque2, nu = 0,
                       fixedeff = F, scm=T, eps_rel=1e-5, eps_abs=1e-5)

    msyn2 <- multisynth(gdpcap ~ trt, regionno, year,
                       basque2 %>% arrange(desc(year)), nu = 0,
                       fixedeff = F, scm=T, eps_rel=1e-5, eps_abs=1e-5)

    
    # weights are the same
    expect_equal(c(msyn$weights), c(msyn2$weights))

    # estimates are the same
    expect_equal(predict(msyn), predict(msyn2))

}
)
