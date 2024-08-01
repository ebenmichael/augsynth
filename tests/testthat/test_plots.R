
context("Test plotting features of single augsynth objects")

library(Synth)
data(basque)

basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~ 0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)

syn <- augsynth(gdpcap ~ trt, regionno, year, basque, progfunc = "Ridge", scm = T)

test_that("All plot types for plot.summary.augsynth() are working", {

    test_results <- c()

    for (it in c('conformal', 'jackknife', 'jackknife+', 'permutation', 'permutation_rstat')){
        s_syn <- summary(syn, inf_type = it)
        for (pt in c('estimate only', 'estimate', 'placebo', 'outcomes', 'outcomes raw average')) {
            if((!it %in% c('permutation', 'permutation_rstat')) & (pt == 'placebo')) {
                invisible()
            } else {
                p <- plot(s_syn, plot_type = pt)
                result <- 'ggplot' %in% class(p)
                test_results <- c(test_results, result)
            }
        }
    }

    expect_true(all(test_results))
})

test_that("All plot types for plot.augsynth() are working", {

    test_results <- c()

    for (it in c('conformal', 'jackknife', 'jackknife+', 'permutation', 'permutation_rstat')){
        for (pt in c('estimate only', 'estimate', 'placebo', 'outcomes', 'outcomes raw average')) {
            p <- plot(syn, plot_type = pt, inf_type = it)
            result <- 'ggplot' %in% class(p)
            test_results <- c(test_results, result)
        }
    }

    expect_true(all(test_results))
})

