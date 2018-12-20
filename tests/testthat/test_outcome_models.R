context("Testing that augmenting synth with different models loads and runs")



library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)


                            
test_that("Augmenting synth with glmnet runs", {

    if(!requireNamespace("glmnet", quietly = TRUE)) {
        ## should fail because glmnet isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="EN", weightfunc="SCM"),
                     "you must install the glmnet package")

        ## install glmnet
        install.packages("glmnet")
    }

    ## should run because glmnet is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="EN", weightfunc="SCM"),
                 NA)    
}
)



test_that("Augmenting synth with random forest runs", {

    if(!requireNamespace("randomForest", quietly = TRUE)) {
        ## should fail because randomForest isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="RF", weightfunc="SCM"),
                     "you must install the randomForest package")

        ## install randomForest
        install.packages("randomForest")
    }

    ## should run because randomForest is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="RF", weightfunc="SCM"),
                 NA)    
}
)




test_that("Augmenting synth with gsynth runs", {

    if(!requireNamespace("gsynth", quietly = TRUE)) {
        ## should fail because gsynth isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="GSYN", weightfunc="SCM"),
                     "you must install the gsynth package")

        ## install gsynth
        devtools::install_github("xuyiqing/gsynth")
    }

    ## should run because gsynth is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="GSYN", weightfunc="SCM"),
                 NA)    
}
)



test_that("Augmenting synth with MCPanel runs", {

    if(!requireNamespace("MCPanel", quietly = TRUE)) {
        ## should fail because MCPanel isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="MCP", weightfunc="SCM"),
                     "you must install the MCPanel package")

        ## install MCPanel
        devtools::install_github("susanathey/MCPanel")
    }

    ## should run because MCPanel is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="MCP", weightfunc="SCM"),
                 NA)    
}
)




test_that("Augmenting synth with CausalImpact runs", {

    if(!requireNamespace("CausalImpact", quietly = TRUE)) {
        ## should fail because CausalImpact isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="CausalImpact", weightfunc="SCM"),
                     "you must install the CausalImpact package")

        ## install CausalImpact
        install.packages("CausalImpact")
    }

    ## should run because CausalImpact is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="CausalImpact", weightfunc="SCM"),
                 NA)    
}
)




test_that("Augmenting synth with keras runs", {

    if(!requireNamespace("keras", quietly = TRUE)) {
        ## should fail because keras isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="seq2seq", weightfunc="SCM"),
                     "you must install the keras package")

        ## install keras
        install.packages("keras")
    }

    ## should run because keras is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, 1975, basque, progfunc="seq2seq", weightfunc="SCM"),
                 NA)    
}
)
