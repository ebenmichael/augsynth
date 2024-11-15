

library( tidyverse )

data(basque,  package="Synth")

head( basque)

#basque <- basque %>%
#    dplyr::select( regionno, year, gdpcap ) %>%
#    as_tibble()

basque <- basque %>%
    mutate(trt = case_when(year < 1975 ~ 0,
                           regionno != 17 ~ 0,
                           regionno == 17 ~ 1)) %>%
    filter(regionno != 1) %>%
    dplyr::select( -c( sec.agriculture:invest ) )

head( basque )

table( basque$trt, basque$regionno )



test_that( "MDES_table corresponds to default treatment table", {

    syn <- augsynth(gdpcap ~ trt, regionno, year,
                    data=basque, scm = TRUE)

    summ <- summary(syn, inf_type = 'permutation')

    mm <- summ$permutations$MDES_table[1:7] %>%
        select( sort( names(.)))
    mm
    tt <- treated_table(syn) %>%
        select( sort( names(.)))
    tt

    expect_equal( as.data.frame(tt)[c("ATT","raw_average","Yhat", "tx")],
                  as.data.frame(mm)[c("ATT","raw_average","Yhat", "tx")] )
} )

test_that( "Placebo distribtion works", {

    syn <- augsynth(gdpcap ~ trt, regionno, year,
                    data=basque, scm = TRUE, progfunc= "none" )
    tt <- treated_table(syn)
    tt

    donor_table(syn)

    b2 = basque %>%
        mutate( trt = 0 + (regionno == 7 ) * (year>=1975) )
    syn7 <- augsynth(gdpcap ~ trt, regionno, year,
                    data=b2, scm = TRUE, progfunc= "none")

    # These should be the same?  Or close to the same?
    # (They were not under ridge, probably due to the cross-validation procedure.)
    gaps7 <- augsynth:::get_placebo_gaps(syn7)
    g17 <- gaps7 %>%
        dplyr::filter( ID == 17 ) %>%
        pivot_longer( cols = `1955`:`1997`,
                      names_to = "year",
                      values_to = "est" )
    expect_equal( tt$ATT - g17$est, rep(0,nrow(tt)), tolerance = 0.000001 )

    #as.numeric( gaps[ gaps$ID == 10, -1 ] )
    #- as.numeric( gaps7[ gaps7$ID == 10, -1 ] )

    n_unit = length(unique(basque$regionno))
    expect_equal( nrow(gaps7), n_unit )
    n_year = length(unique(basque$year))

    expect_equal( ncol(gaps7), 1 + n_year )
    expect_equal( dim( syn ), c( n_unit, n_year ) )

    pc <- augsynth:::add_placebo_distribution( syn )
    rs <- pc$results$permutations
    expect_equal( names(rs), c("placebo_dist", "MDES_table") )
    expect_equal( nrow( rs$placebo_dist ), n_year * n_unit )

    expect_equal( syn$unit_var, "regionno" )
    expect_equal( syn$time_var, "year" )
})






