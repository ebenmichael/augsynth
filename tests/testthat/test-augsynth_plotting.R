
library( tidyverse )

test_that("generall plotting functions do not crash", {

    data(basque, package="Synth")

    basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                                regionno != 17 ~0,
                                                regionno == 17 ~ 1)) %>%
        filter(regionno != 1, year > 1960 & year < 1990 )



    syn <- single_augsynth(gdpcap ~ trt, regionno, year, 1975, basque,
                           progfunc="None", scm=T, fixedeff = F)

    syn

    sum <- summary( syn, inf_type = "permutation" )
    sum

    # Check that it will add inference on the fly
    auto_add_inf <- augsynth:::augsynth_plot_from_results( syn, inf_type = "permutation_rstat" )
    expect_true( is.ggplot( auto_add_inf ) )



    # use testthat to check print.augsynth.summary writes to console
    expect_output( print( sum ), "Permutation inference" )

    p0 <- plot( sum )
    p0

    p <- plot( sum, plot_type = "estimate" )
    p
    #expect_equal( p0, p )
    pEO <- plot( syn, plot_type = "estimate only" )
    pEO
    p <- plot( syn, plot_type = "outcomes" )
    p
    p <- plot( syn, plot_type = "outcomes raw average" )
    p
    p <- plot( syn, plot_type = "placebo" )
    p

    sum$inf_type

    # These two plots should be the same
    p2 <- plot( syn )
    expect_true( is.ggplot(p2 ) )

    sump <- summary( syn )
    p3 <- plot( sump )
    expect_true( is.ggplot( p3 ) )

    expect_output( print( sump ), "Conformal inference" )

    dt = donor_table( sump )
    expect_equal( sum( dt$weight ), 1, tolerance = 0.000001 )

    # TODO: Check the "cv" flag -- has to rerun conformal inference?
    #plt <- plot( syn, inf_type="jackknife", cv=TRUE )
    #plt
    #expect_true( is.ggplot( plt ) ) )

})




