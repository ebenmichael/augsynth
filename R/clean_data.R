########################################################
## Collection of scripts to clean different data sources
########################################################


clean_basque <- function(df, sim_num=-1) {
    #' Clean up the basque data that comes from the Synth package
    #' @export
    cleaned <- df %>%
        select(regionno, year, gdpcap) %>%
        rename(unit=regionno, outcome=gdpcap, time=year) %>%
        filter(unit > 1) %>%
        mutate(outcome_id=1,
               treated = unit == 17,
               potential_outcome = ifelse(treated, "Y(1)", "Y(0)"),
               synthetic = "N",
               sim_num=sim_num)

    metadata <- data.frame(t_int=1975, sim_num=-1)
    return(list(outcomes=cleaned,
                metadata=metadata))
}
