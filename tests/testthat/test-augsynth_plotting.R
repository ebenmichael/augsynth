

library(tidyverse)
data(basque, package = "Synth")

context("Test plotting features of single augsynth objects")

# Setup small dataset and some augsynth summaries ----

basque <- basque %>%
    mutate(
        trt = case_when(year < 1975 ~ 0,
                        regionno != 17 ~0,
                        regionno == 17 ~ 1)
    ) %>%
    filter(year > 1960, year < 1980) %>%
    filter(regionno > 7 )



syn <- single_augsynth(
    gdpcap ~ trt, regionno, year, 1975, basque,
    progfunc = "None", scm = TRUE, fixedeff = FALSE
)

sum_perm <- summary(syn, inf_type = "permutation")
sum_conf <- summary(syn, inf_type = "conformal")


# The tests ----

test_that("'estimate' plots work for both augsynth and summary objects", {

    p_sum <- plot(sum_perm)
    p_sum_explicit <- plot(sum_perm, plot_type = "estimate")
    p_sum

    expect_s3_class(p_sum, "ggplot")
    expect_equal(p_sum, p_sum_explicit)

    expect_message(
        p_syn <- plot(syn, plot_type = "estimate", inf=TRUE, inf_type="permutation")
    )
    expect_s3_class(p_syn, "ggplot")
    p_sum
    p_sum_explicit
    p_syn

    expect_equal( p_sum, p_syn )
    expect_equal(layer_data(p_sum), layer_data(p_syn), tolerance = 0.000001)

    # Checking that inf=FALSE works, and does not trigger
    # a warning about long run time.
    p_syn2 <- plot( sum_perm, inf=FALSE )
    p_syn2

    p_syn3 <- plot(syn, plot_type = "estimate", inf=FALSE )
    p_syn3
    expect_equal( p_syn2, p_syn3 )
    expect_equal( p_sum, p_syn3)
    expect_equal( layer_data(p_syn2), layer_data(p_syn3) )

})




test_that("'outcomes' plots work for both augsynth and summary objects", {

    p_sum <- plot(sum_conf, plot_type = "outcomes" )
    p_sum_noCI <- plot(sum_perm, plot_type = "outcomes", inf=FALSE )
    p_sum +
        scale_y_continuous( limits=c( 0, 10 ) )

    p_sum_noCI +
        scale_y_continuous( limits=c( 0, 10 ) )

    expect_message(
        p_syn <- plot(syn, plot_type = "outcomes", inf_type = "permutation")
    )

    p_sum
    p_syn

    expect_s3_class(p_sum, "ggplot")
    expect_s3_class(p_syn, "ggplot")

    expect_equal(p_sum$data, p_syn$data)

    p_sum2 <- plot( sum_conf, plot_type = "outcomes" )
    p_sum2
    expect_s3_class(p_sum2, "ggplot")
    expect_equal(p_sum2$data, p_syn$data)
})



test_that("Permutation plotting helper respects requested placebo type", {

    p_att <- augsynth:::permutation_plot(sum_perm, inf_type = "permutation")
    p_rstat <- augsynth:::permutation_plot(sum_perm, inf_type = "permutation_rstat")

    expect_s3_class(p_att, "ggplot")
    expect_s3_class(p_rstat, "ggplot")
    p_att
    p_rstat

    expect_true("ATT" %in% names(p_att$data))
    expect_true("rstat" %in% names(p_rstat$data))

    expect_false(identical(p_att$data$ATT, p_rstat$data$rstat))
})


test_that("'placebo' plots from augsynth objects preserve permutation_rstat when requested", {

    expect_message(
        p_perm <- plot(syn, plot_type = "placebo", inf_type = "permutation")
)
    expect_message(
    p_rstat <- plot(syn, plot_type = "placebo", inf_type = "permutation_rstat")
)

    expect_s3_class(p_perm, "ggplot")
    expect_s3_class(p_rstat, "ggplot")

    expect_equal(rlang::as_label(p_perm$mapping$y), "ATT")
    expect_equal(rlang::as_label(p_rstat$mapping$y), "rstat")

    expect_equal(p_perm$labels$y, "Estimate (ATT)")
    expect_equal(p_rstat$labels$y, "Estimate (RMSPE-adjusted ATT)")
})


test_that("'placebo' plots from augsynth objects default to permutation inference", {

    expect_error(
        p_default <- plot(sum_conf, plot_type = "placebo"),
        "Placebo plots are only available for permutation-based inference"
    )

    expect_message(
        expect_warning(
        p_default <- plot(syn, plot_type = "placebo", inf_type = "conformal"),
        "Placebo plots are only available for permutation-based inference"
    )
    )
    p_default
    expect_s3_class(p_default, "ggplot")

    expect_message(
        p_perm <- plot(syn, plot_type = "placebo", inf_type = "permutation")
    )
    expect_s3_class(p_perm, "ggplot")

    p_perm2 <- plot(sum_perm, plot_type = "placebo" )

    expect_equal(rlang::as_label(p_default$mapping$y), "ATT")
    expect_equal(rlang::as_label(p_perm$mapping$y), "ATT")
    expect_equal(p_default$data, p_perm$data)

    expect_equal(p_default, p_perm2)
    expect_equal(layer_data(p_default), layer_data(p_perm2), tolerance = 0.0001 )
})




test_that("Plotting augsynth objects adds inference on the fly", {

    expect_message(
        p <- plot(syn, inf_type = "permutation_rstat")
    )
    expect_s3_class(p, "ggplot")
})


test_that("Summary plots ignore inf_type passed through dots", {

    expect_warning(
        p <- plot(sum_perm, inf_type = "jackknife"),
        "ignored for summary.augsynth objects"
    )
    expect_s3_class(p, "ggplot")
})


test_that("Outcomes plotting helper works on summary objects", {

    plt <- augsynth:::augsynth_outcomes_plot(sum_perm)
    expect_s3_class(plt, "ggplot")
})


test_that( "'cv' plots work for augsynth objects", {

    expect_error(
        p <- plot(syn, plot_type = "cv"),
        "Cross-validation errors and lambdas are not available in this augsynth object"
    )

    asyn <- augsynth(gdpcap ~ trt, regionno, year, basque, t_int=1975,
                     progfunc = "Ridge", scm = T)
    p2 <- plot( asyn, cv=TRUE )
    p2
    expect_s3_class(p2, "ggplot")

    p3 <- plot( asyn, plot_type="cv" )
    p3
    expect_s3_class(p3, "ggplot")
    expect_equal( p2, p3 )
    expect_equal( layer_data(p2), layer_data(p3) )


    sum_ridge <- summary( asyn, inf_type = "permutation")
    p3 <- plot( sum_ridge, plot_type="cv" )
    p3

    expect_equal( p2, p3 )
    expect_equal( layer_data(p2), layer_data(p3) )

    sum_ridge <- summary( asyn, inf_type = "permutation")
    p3 <- plot( sum_ridge, cv=TRUE )
    p3

    expect_equal( p2, p3 )
    expect_equal( layer_data(p2), layer_data(p3) )
})



