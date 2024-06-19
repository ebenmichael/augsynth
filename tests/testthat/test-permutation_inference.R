

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
    filter(regionno != 1)

table( basque$trt, basque$regionno )



test_that( "MDES_table corresponds to default treatment table", {

    syn <- augsynth(gdpcap ~ trt, regionno, year,
                    data=basque, scm = TRUE, inf_type = "permutation" )

    mm <- syn$results$permutations$MDES_table[1:7] %>%
        select( sort( names(.)))
    mm
    tt <- treated_table(syn) %>%
        select( sort( names(.)))
    tt
    expect_equal( as.data.frame(tt), as.data.frame(mm) )


    syn
    expect_equivalent(!is.null(syn$results), TRUE)
    gaps <- augsynth:::get_placebo_gaps(syn)

    n_unit = length(unique(basque$regionno))
    expect_equal( nrow(gaps), n_unit )
    n_year = length(unique(basque$year))

    expect_equal( ncol(gaps), 1 + n_year )
    expect_equal( dim( syn ), c( n_unit, n_year ) )

    pc <- augsynth:::add_placebo_distribution( syn )
    rs <- pc$results$permutations
    expect_equal( names(rs), c("placebo_dist", "MDES_table") )
    expect_equal( nrow( rs$placebo_dist ), n_year * n_unit )

    expect_equal( syn$unit_var, "regionno" )
    expect_equal( syn$time_var, "year" )
})






