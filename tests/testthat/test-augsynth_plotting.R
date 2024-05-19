
library( tidyverse )

test_that("generall plotting functions do not crash", {

    data(basque, package="Synth")

    basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                                regionno != 17 ~0,
                                                regionno == 17 ~ 1)) %>%
        filter(regionno != 1)



    syn <- single_augsynth(gdpcap ~ trt, regionno, year, 1975, basque,
                           progfunc="None", scm=T, fixedeff = F)

    syn

    sum <- summary( syn, inf_type = "permutation" )
    sum

    p0 <- plot( sum )
    p0

    p <- plot( sum, plot_type = "estimate" )
    p
    #expect_equal( p0, p )
    pEO <- plot( syn, plot_type = "estimate only" )
    pEO
    p <- plot( sum, plot_type = "outcomes" )
    p
    p <- plot( sum, plot_type = "outcomes raw average" )
    p
    p <- plot( sum, plot_type = "placebo" )
    p

    sum$inf_type

    # These should be the same
    p2 <- plot( syn )
    expect_true( is.ggplot(p2 ) )

    sump <- summary( syn )
    p3 <- plot( sump )

    expect_true( is.ggplot( p3 ) )
})
