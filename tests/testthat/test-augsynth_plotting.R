library( tidyverse )
data(basque, package = "Synth")

context("Test plotting features of single augsynth objects")

basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~ 0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)

syn <- augsynth(gdpcap ~ trt, regionno, year, basque, progfunc = "Ridge", scm = T)

test_that("All plot types for plot.summary.augsynth() are working", {

    for (it in c('conformal', 'jackknife', 'jackknife+', 'permutation', 'permutation_rstat')){
        s_syn <- summary(syn, inf_type = it)
        for (pt in c('estimate only', 'estimate', 'placebo', 'outcomes', 'outcomes raw average')) {
            if((!it %in% c('permutation', 'permutation_rstat')) & (pt == 'placebo')) {
                expect_error(plot(s_syn, plot_type = pt))
            } else {
                p <- plot(s_syn, plot_type = pt)
                expect_true('ggplot' %in% class(p))
            }
        }
    }

})

test_that("All plot types for plot.augsynth() are working", {

    for (it in c('conformal', 'jackknife', 'jackknife+', 'permutation', 'permutation_rstat')){
        for (pt in c('estimate only', 'estimate', 'placebo', 'outcomes', 'outcomes raw average')) {
            p <- plot(syn, plot_type = pt, inf_type = it)
            expect_true('ggplot' %in% class(p))
        }
    }

})


test_that("General plotting functions do not crash", {

    syn <- single_augsynth(gdpcap ~ trt, regionno, year, 1975, basque,
                           progfunc="None", scm=T, fixedeff = F)

    # Testing some of the permutation plotting functions
    tst_plt <- augsynth:::permutation_plot( syn, inf_type = "permutation_rstat" )
    expect_true( is.ggplot( tst_plt ) )

    sum <- summary( syn, inf_type = "permutation" )
    expect_equal( sum$inf_type, "permutation" )

    # Check that it will add inference on the fly
    expect_message( auto_add_inf <- plot( syn, inf_type = "permutation_rstat" ) )
    expect_true( is.ggplot( auto_add_inf ) )

    plt <- augsynth:::augsynth_outcomes_plot( sum )
    expect_true( is.ggplot( plt ) )

    # use testthat to check print.augsynth.summary writes to console
    expect_output( print( sum ), "Permutation inference" )

    p0 <- plot( sum )
    p <- plot( sum, plot_type = "estimate" )
    expect_equal( p0, p )

    pEO <- plot( sum, plot_type = "estimate only" )
    expect_message(
        pEO1 <- plot( syn, plot_type = "estimate only", inf_type = "permutation" )
    )
    expect_equal( pEO, pEO1 )

    po <- plot( sum, plot_type = "outcomes" )
    po1 = plot( syn, plot_type = "outcomes", inf_type = "permutation" )
    expect_equal( po, po1 )

    praw <- plot( sum, plot_type = "outcomes raw average" )
    praw1 <- plot( syn, plot_type = "outcomes raw average", inf_type="permutation" )
    expect_equal( praw, praw1 )

    ppla <- plot( sum, plot_type = "placebo" )
    ppla1 <- plot( syn, plot_type = "placebo", inf_type="permutation" )
    expect_equal( ppla, ppla1 )

    # These two plots should be the same
    p2 <- plot( syn )
    expect_true( is.ggplot(p2 ) )

    sump <- summary( syn )
    p3 <- plot( sump )
    expect_true( is.ggplot( p3 ) )

    expect_output( print( sump ), "Conformal inference" )

    dt = sump$donor_table
    expect_equal( sum( dt$weight ), 1, tolerance = 0.000001 )

})
