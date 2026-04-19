

library( tidyverse )

data(basque,  package="Synth")

head( basque)


basque <- basque %>%
    mutate(trt = case_when(year < 1975 ~ 0,
                           regionno != 17 ~ 0,
                           regionno == 17 ~ 1)) %>%
    filter(regionno != 1) %>%
    dplyr::select( -c( sec.agriculture:invest ) )

head( basque )

table( basque$trt, basque$regionno )



test_that( "Placebo distribtion works", {

    syn <- augsynth(gdpcap ~ trt, regionno, year,
                    data = basque, scm = TRUE, progfunc = "none" )
    tt <- treated_table(syn)


    b2 = basque %>%
        mutate( trt = 0 + (regionno == 7 ) * (year>=1975) )
    syn7 <- augsynth(gdpcap ~ trt, regionno, year,
                    data = b2, scm = TRUE, progfunc = "none")

    don1 <- donor_table(syn)
    don2 <- donor_table(syn7)

    # Do the shared donor units have same RMSPE
    expect_equal( dim( don1 ), dim( don2 ) )
    fulldone = don1 %>%
        inner_join( don2, by = "regionno", suffix = c("_syn1", "_syn7") )
    expect_equal( fulldone$RMSPE_syn1, fulldone$RMSPE_syn7, tolerance = 0.000001 )

    gaps1 <- augsynth:::get_placebo_gaps(syn)
    gaps7 <- augsynth:::get_placebo_gaps(syn7)

    # The estimated impacts for the treated unit should be the same as
    # the estimated impact when the treated unit is considered a donor
    # unit.
    g1 <- gaps1 %>%
        pivot_longer( cols = `1955`:`1997`,
                      names_to = "year",
                      values_to = "est" )
    g7 <- gaps7 %>%
        pivot_longer( cols = `1955`:`1997`,
                      names_to = "year",
                      values_to = "est" )
    gs <- g1 %>%
        inner_join( g7, by = c("ID","year"), suffix = c("_syn1", "_syn7") )
    expect_equal( gs$est_syn1, gs$est_syn7, tolerance = 0.000001 )

    expect_equal( tt$ATT, g7$est[g7$ID == 17], tolerance = 0.000001 )

    n_unit = length(unique(basque$regionno))
    n_year = length(unique(basque$year))

    expect_equal( nrow(gaps7), n_unit )
    expect_equal( ncol(gaps7), 1 + n_year )
    expect_equal( dim( syn ), c( n_unit, n_year ) )

    pc <- augsynth:::add_placebo_distribution( syn )
    rs <- pc$results$permutations
    expect_equal( names(rs), c("placebo_dist", "MDES_table") )
    expect_equal( nrow( rs$placebo_dist ), n_year * n_unit )

    expect_equal( syn$unit_var, "regionno" )
    expect_equal( syn$time_var, "year" )

})



test_that( "Placebo distribtion works with ridge", {

    syn <- augsynth(gdpcap ~ trt, regionno, year,
                    data = basque, scm = TRUE, progfunc = "ridge", lambda = 0.1 )
    tt <- treated_table(syn)

    b2 = basque %>%
        mutate( trt = 0 + (regionno == 7 ) * (year>=1975) )
    syn7 <- augsynth(gdpcap ~ trt, regionno, year,
                     data = b2, scm = TRUE, progfunc = "ridge", lambda = 0.1 )

    # The estimated impacts for the treated unit should be the same as
    # the estimated impact when the treated unit is considered a donor
    # unit.
    gaps7 <- augsynth:::get_placebo_gaps(syn7)
    g17 <- gaps7 %>%
        dplyr::filter( ID == 17 ) %>%
        pivot_longer( cols = `1955`:`1997`,
                      names_to = "year",
                      values_to = "est" )
    expect_equal( tt$ATT - g17$est, rep(0,nrow(tt)), tolerance = 0.000001 )


})





test_that("Inference carries through in summary objects", {

    syn <- augsynth(gdpcap ~ trt, regionno, year,
                    data = basque, scm = TRUE, progfunc = "none" )

    sum <- summary( syn, inf_type = "permutation" )
    expect_equal( sum$inf_type, "permutation" )
})


