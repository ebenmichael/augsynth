############################################
## Scripts to plot aspects of the simulations
############################################


plot_outcomes <- function(outcomes, metadata, trt_unit=NULL) {
    #' Plot the outcomes from a study/simulation
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #' @param trt_unit Unit to count as treated, defaults to using
    #'                 the treated column
    #' @export
    
    if(!is.null(trt_unit)) {
        outcomes <- outcomes %>% 
            mutate(treated=(trt_unit == unit)) %>%
            filter(!(treated==FALSE &
                    synthetic=="N" & 
                    potential_outcome == "Y(1)")
                   )
    }

    # join outcomes with metadata on sim number
    p <- outcomes %>%
        inner_join(metadata) %>%
        ## create interaction variables for grouping and coloring
        mutate(grouping=interaction(treated, synthetic, unit, potential_outcome),
               label=interaction(treated, synthetic, potential_outcome),
               label=plyr::revalue(label,
                                   c("TRUE.N.Y(1)"="Outcome Under Treatment",
                                     "TRUE.N.Y(0)"="Outcome Under Control",
                                     "FALSE.N.Y(0)"="Donor Pool",
                                     "TRUE.Y.Y(0)"="Synthetic Control",
                                     "TRUE.DR.Y(0)"="Double Robust"
                                     )
                                   ),
               synthetic = ifelse(synthetic %in% c("Y", "DR"),
                                  "Y",
                                  "N")
               ) %>%                 
    # plot the outcomes
    ggplot() +
        geom_line(aes(x=time, y=outcome, # plot outcomes vs time
                      group = grouping, color=label,
                      alpha=treated, linetype=synthetic), size=2) +
        geom_vline(aes(xintercept=t_int)) + 
        scale_color_manual(values=c("Donor Pool"="#888888",
                                    "Outcome Under Treatment"="#FDB515",
                                    "Synthetic Control"="#ED4E33",
                                    "Outcome Under Control"="#3B7EA1",
                                    "Double Robust"="#B9D3B6")) +
    scale_alpha_manual(values=c(0.05, 1)) +
        guides(alpha=FALSE, linetype=FALSE)
    if("syn_method" %in% names(outcomes)) {
        p <- p + facet_grid(syn_method ~ outcome_id)
    } else {
        p <- p + facet_wrap( ~outcome_id)
    }
    p <- p + theme_bw()
    return(p)
}

plot_att <- function(outcomes, metadata, trt_unit=NULL) {
    #' Plot the estimate of the att
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #' @param trt_unit Unit to count as treated, defaults to using
    #'                 the treated column
    #' @export
    
    if(!is.null(trt_unit)) {
        outcomes <- outcomes %>% 
            mutate(treated=(trt_unit == unit)) %>%
            filter(!(treated==FALSE &
                    synthetic=="N" & 
                    potential_outcome == "Y(1)")
                   )
    }

    # join outcomes with metadata on sim number
    p <- outcomes %>%
        inner_join(metadata) %>%
        ## create interaction variables for grouping and coloring
        mutate(grouping=interaction(treated, synthetic, unit, potential_outcome),
               label=interaction(treated, synthetic, potential_outcome),
               label=plyr::revalue(label,
                                   c("TRUE.N.Y(1)"="Outcome Under Treatment",
                                     "TRUE.N.Y(0)"="Outcome Under Control",
                                     "FALSE.N.Y(0)"="Donor Pool",
                                     "TRUE.Y.Y(0)"="Synthetic Control",
                                     "TRUE.DR.Y(0)"="Double Robust"
                                     )
                                   ),
               ) %>%
        filter(treated==TRUE) %>%
        select(-label, -grouping, -potential_outcome) %>%
        spread(synthetic, outcome) %>%
        mutate(att=N-Y) %>%
    # plot the outcomes
    ggplot() +
        geom_line(aes(x=time, y=att), size=2, color="#B9D3B6") + 
        geom_vline(aes(xintercept=t_int))
    if("syn_method" %in% names(outcomes)) {
        p <- p + facet_grid(syn_method ~ outcome_id)
    } else {
        p <- p + facet_wrap( ~outcome_id)
    }
    p <- p + ggtitle("ATT Estimate") + theme_bw()
    return(p)
    }
