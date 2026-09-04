context("Testing that augmenting synth with different models loads and runs")



library(Synth)
data(basque)
basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)



test_that("Augmenting synth with gsynth runs and produces the correct result", {

    if(!requireNamespace("gsynth", quietly = TRUE)) {
        ## should fail because gsynth isn't installed
        expect_error(augsynth(gdpcap ~ trt, regionno, year, basque,
                              progfunc="GSYN", scm=T),
                     "you must install the gsynth package")

        ## install gsynth
        install.packages("gsynth", repos = "http://cran.us.r-project.org")
    }

    ## should run because gsynth is installed
    expect_error(
      augsynth(gdpcap ~ trt, regionno, year, basque,
                                progfunc = "GSYN", scm = T, CV = 0, r = 4),
      NA)
    asyn_gsyn <- augsynth(gdpcap ~ trt, regionno, year, basque,
                          progfunc = "GSYN", scm = F, CV = 0, r = 4)
    expect_equal(summary(asyn_gsyn, inf_type = 'none')$average_att$Estimate,
                 -0.1444637, tolerance=1e-4)
}
)
