library( tidyverse )
data(basque, package = "Synth")

context("Test plotting features of single augsynth objects")

basque <- basque %>% mutate(trt = case_when(year < 1975 ~ 0,
                                            regionno != 17 ~ 0,
                                            regionno == 17 ~ 1)) %>%
    filter(regionno != 1)


syn <- single_augsynth(gdpcap ~ trt, regionno, year, 1975, basque,
                       progfunc = "None", scm = T, fixedeff = F)
sum <- summary( syn, inf_type = "permutation" )


test_that("Confirm equivalence of plotting estimate across augsynth and summary objects", {

    # Testing some of the permutation plotting functions
    tst_plt <- augsynth:::permutation_plot( syn, inf_type = "permutation_rstat" )
    expect_true( is_ggplot( tst_plt ) )


    p0 <- plot( sum )
    p <- plot( sum, plot_type = "estimate" )
    expect_equal( p0, p )

    pEO <- plot( sum, plot_type = "estimate only" )
    expect_message(
        pEO1 <- plot( syn, plot_type = "estimate only", inf_type = "permutation" )
    )
    expect_equal( pEO, pEO1 )
})



test_that("Confirm equivalence of plotting outcomes across augsynth and summary objects", {

    po <- plot( sum, plot_type = "outcomes" )
    po1 = plot( syn, plot_type = "outcomes", inf_type = "permutation" )
    expect_equal( po, po1 )

    praw <- plot( sum, plot_type = "outcomes raw average" )
    praw1 <- plot( syn, plot_type = "outcomes raw average", inf_type="permutation" )
    expect_equal( praw$data, praw1$data )

    # These two plots should be the same
    p2 <- plot( syn )
    expect_true( is_ggplot(p2 ) )

})


test_that("Requesting a placebo plot defaults to permutation inference", {

    syn <- single_augsynth(gdpcap ~ trt, regionno, year, 1975, basque,
                           progfunc = "None", scm = T, fixedeff = F)

    sum <- summary( syn, inf_type = "permutation" )

    ppla <- plot( sum, plot_type = "placebo" )
    ppla1 <- plot( syn, plot_type = "placebo", inf_type="permutation" )

    expect_equal( ppla, ppla1 )
    expect_true(all(ppla$data$ATT == ppla1$data$ATT))
})


test_that("Check that plotting functions will add inference on the fly", {

    syn <- single_augsynth(gdpcap ~ trt, regionno, year, 1975, basque,
                           progfunc = "None", scm = T, fixedeff = F)

    expect_message( auto_add_inf <- plot( syn, inf_type = "permutation_rstat" ) )
    expect_true( is_ggplot( auto_add_inf ) )
})


test_that("Check that we can use plotting helper functions from summary", {

    sum <- summary(syn)

    plt <- augsynth:::augsynth_outcomes_plot( sum )
    expect_true( is_ggplot( plt ) )

})
