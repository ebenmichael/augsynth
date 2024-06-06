
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

    # Testing some of the permutation plotting functions
    # (Note ci95 is not really used by core package?)
    ci95plot <- augsynth:::ci95_rstat_plot( syn, inf_type = "permutation_rstat" )
    expect_true( is.ggplot( ci95plot ) )

    tst_plt <- augsynth:::permutation_plot( syn, inf_type = "permutation_rstat" )
    expect_true( is.ggplot( tst_plt ) )


    sum <- summary( syn, inf_type = "permutation" )
    sum
    expect_equal( sum$inf_type, "permutation" )

    # Check that it will add inference on the fly
    auto_add_inf <- plot( syn, inf_type = "permutation_rstat" )
    auto_add_inf
    expect_true( is.ggplot( auto_add_inf ) )

    plt <- augsynth:::augsynth_outcomes_plot( sum )
    expect_true( is.ggplot( plt ) )

    # use testthat to check print.augsynth.summary writes to console
    expect_output( print( sum ), "Permutation inference" )

    p0 <- plot( sum )
    p0

    p <- plot( sum, plot_type = "estimate" )
    p
    expect_equal( p0, p )

    pEO <- plot( sum, plot_type = "estimate only" )
    pEO

    pEO1 <- plot( syn, plot_type = "estimate only", inf_type = "permutation" )
    expect_equal( pEO, pEO1 )

    po <- plot( sum, plot_type = "outcomes" )
    po
    po1 = plot( syn, plot_type = "outcomes", inf_type = "permutation" )
    po1
    expect_equal( po, po1 )

    praw <- plot( sum, plot_type = "outcomes raw average" )
    praw
    praw1 <- plot( syn, plot_type = "outcomes raw average", inf_type="permutation" )
    praw1
    expect_equal( praw, praw1 )

    ppla <- plot( sum, plot_type = "placebo" )
    ppla
    ppla1 <- plot( syn, plot_type = "placebo", inf_type="permutation" )
    ppla1
    expect_equal( ppla, ppla1 )


    # These two plots should be the same
    p2 <- plot( syn )
    p2
    expect_true( is.ggplot(p2 ) )

    sump <- summary( syn )
    p3 <- plot( sump )
    p3
    expect_true( is.ggplot( p3 ) )

    expect_output( print( sump ), "Conformal inference" )

    dt = donor_table( sump )
    expect_equal( sum( dt$weight ), 1, tolerance = 0.000001 )

    # TODO: Check the "cv" flag -- has to rerun conformal inference?
    #plt <- plot( syn, inf_type="jackknife", cv=TRUE )
    #plt
    #expect_true( is.ggplot( plt ) ) )

})




