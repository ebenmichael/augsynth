library(tidyverse)
data(basque, package = "Synth")

context("Test plotting features of single augsynth objects")

basque <- basque %>%
    mutate(
        trt = case_when(
            year < 1975 ~ 0,
            regionno != 17 ~ 0,
            regionno == 17 ~ 1
        )
    ) %>%
    filter(regionno != 1)

syn <- single_augsynth(
    gdpcap ~ trt, regionno, year, 1975, basque,
    progfunc = "None", scm = TRUE, fixedeff = FALSE
)

sum_perm <- summary(syn, inf_type = "permutation")
sum_conf <- summary(syn, inf_type = "conformal")


test_that("Estimate plots work for both augsynth and summary objects", {

    p_sum <- plot(sum_perm)
    p_sum_explicit <- plot(sum_perm, plot_type = "estimate")

    expect_s3_class(p_sum, "ggplot")
    expect_equal(p_sum, p_sum_explicit)

    expect_message(
        p_syn <- plot(syn, plot_type = "estimate", inf_type = "permutation")
    )
    expect_s3_class(p_syn, "ggplot")
})


test_that("Outcomes plots work for both augsynth and summary objects", {

    p_sum <- plot(sum_perm, plot_type = "outcomes")
    p_syn <- plot(syn, plot_type = "outcomes", inf_type = "permutation")

    expect_s3_class(p_sum, "ggplot")
    expect_s3_class(p_syn, "ggplot")

    expect_equal(p_sum$data, p_syn$data)
})


test_that("Permutation plotting helper respects requested placebo type", {

    p_att <- augsynth:::permutation_plot(sum_perm, inf_type = "permutation")
    p_rstat <- augsynth:::permutation_plot(sum_perm, inf_type = "permutation_rstat")

    expect_s3_class(p_att, "ggplot")
    expect_s3_class(p_rstat, "ggplot")

    expect_true("ATT" %in% names(p_att$data))
    expect_true("rstat" %in% names(p_rstat$data))

    expect_false(identical(p_att$data$ATT, p_rstat$data$rstat))
})


test_that("Placebo plots from augsynth objects preserve permutation_rstat when requested", {

    p_perm <- plot(syn, plot_type = "placebo", inf_type = "permutation")
    p_rstat <- plot(syn, plot_type = "placebo", inf_type = "permutation_rstat")

    expect_s3_class(p_perm, "ggplot")
    expect_s3_class(p_rstat, "ggplot")

    expect_equal(rlang::as_label(p_perm$mapping$y), "ATT")
    expect_equal(rlang::as_label(p_rstat$mapping$y), "rstat")

    expect_equal(p_perm$labels$y, "Estimate (ATT)")
    expect_equal(p_rstat$labels$y, "Estimate (RMSPE-adjusted ATT)")
})


test_that("Placebo plots from augsynth objects default to permutation inference", {

    expect_warning(
        p_default <- plot(syn, plot_type = "placebo", inf_type = "conformal"),
        "Placebo plots are only available for permutation-based inference"
    )
    expect_s3_class(p_default, "ggplot")

    p_perm <- plot(syn, plot_type = "placebo", inf_type = "permutation")
    expect_s3_class(p_perm, "ggplot")

    expect_equal(rlang::as_label(p_default$mapping$y), "ATT")
    expect_equal(rlang::as_label(p_perm$mapping$y), "ATT")
    expect_equal(p_default$data, p_perm$data)
})


test_that("Placebo plots from augsynth objects preserve permutation_rstat when requested", {

    p_perm <- plot(syn, plot_type = "placebo", inf_type = "permutation")
    p_rstat <- plot(syn, plot_type = "placebo", inf_type = "permutation_rstat")

    expect_s3_class(p_perm, "ggplot")
    expect_s3_class(p_rstat, "ggplot")

    expect_equal(rlang::as_label(p_perm$mapping$y), "ATT")
    expect_equal(rlang::as_label(p_rstat$mapping$y), "rstat")

    expect_equal(p_perm$labels$y, "Estimate (ATT)")
    expect_equal(p_rstat$labels$y, "Estimate (RMSPE-adjusted ATT)")
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
