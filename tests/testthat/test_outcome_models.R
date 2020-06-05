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
        expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="EN", scm=T),
                     "you must install the glmnet package")

        ## install glmnet
        install.packages("glmnet")
    }

    ## should run because glmnet is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="EN", scm=T),
                 NA)    
}
)



test_that("Augmenting synth with random forest runs", {

    if(!requireNamespace("randomForest", quietly = TRUE)) {
        ## should fail because randomForest isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="RF", scm=T),
                     "you must install the randomForest package")

        ## install randomForest
        install.packages("randomForest")
    }

    ## should run because randomForest is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="RF", scm=T),
                 NA)    
}
)




test_that("Augmenting synth with gsynth runs and produces the correct result", {

    if(!requireNamespace("gsynth", quietly = TRUE)) {
        ## should fail because gsynth isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="GSYN", scm=T),
                     "you must install the gsynth package")

        ## install gsynth
        install.packages("gsynth")
    }

    ## should run because gsynth is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="GSYN", scm=T),
                 NA) 
    asyn_gsyn <- augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="GSYN", scm=F)
    expect_equal(summary(asyn_gsyn)$average_att$Estimate, -0.1444637, tolerance=1e-4) 
}
)



test_that("Augmenting synth with MCPanel runs", {

    if(!requireNamespace("MCPanel", quietly = TRUE)) {
        ## should fail because MCPanel isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="MCP", scm=T),
                     "you must install the MCPanel package")
    } else {
        ## should run because MCPanel is installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="MCP", scm=T),
                     NA)    
    }

    
}
)




test_that("Augmenting synth with CausalImpact runs", {

    if(!requireNamespace("CausalImpact", quietly = TRUE)) {
        ## should fail because CausalImpact isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="CausalImpact", scm=T),
                     "you must install the CausalImpact package")

        ## install CausalImpact
        install.packages("CausalImpact")
    }

    ## should run because CausalImpact is installed
    expect_error(augsynth(gdpcap ~ trt, regionno, year, basque, progfunc="CausalImpact", scm=T),
                 NA)    
}
)
