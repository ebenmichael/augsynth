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
    tmpdf <- outcomes %>%
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
               synthetic = ifelse(synthetic == "N", "N", "Y")
               )
    #return(tmpdf)
    p <- tmpdf %>%                 
        ## plot the outcomes
        ggplot() +
        geom_line(aes(x=time, y=outcome, # plot outcomes vs time
                      group = grouping, color=treated,
                      alpha=treated, linetype=synthetic), size=1.25) +
        geom_vline(aes(xintercept=t_int), linetype=3) + 
        scale_color_manual(values=c("#888888","#2b2b2b"))  +
    scale_alpha_manual(values=c(0.05, 1))
        guides(linetype=FALSE, alpha=FALSE)
    if("syn_method" %in% names(outcomes) & "outcome_id" %in% names(outcomes)) {
        p <- p + facet_grid(syn_method ~ outcome_id, scales="free")
    } else if("outcome_id" %in% names(outcomes)){
        p <- p + facet_wrap( ~outcome_id, scales="free")
    } else if("syn_method" %in% names(outcomes)) {
        p <- p + facet_wrap( ~syn_method, scales="free")
    }
    p <- p + theme_bw()
    return(p)
}


compute_att <- function(outcomes, metadata, trt_unit=NULL) {
    #' Compute the ATT in the post-period
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

    if(is.null(trt_unit)) {
        trt_unit <- (outcomes %>% filter(treated == TRUE) %>% distinct(unit))$unit
    }

    # join outcomes with metadata on sim number
    tmpdf <- outcomes %>%
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
        filter(unit == trt_unit)

    if("syn_method" %in% names(outcomes)) {
        tmpdf  <- tmpdf %>% select(outcome_id, unit, time, t_int,
                                   synthetic, outcome, syn_method)  %>%
            spread(synthetic, outcome) %>%
            gather(syntype, est, -unit, -time, -outcome_id,
                   -t_int, -N, -syn_method)
    } else {
        tmpdf <- tmpdf %>% select(outcome_id, unit, time,
                                  t_int, synthetic, outcome) %>%
            spread(synthetic, outcome) %>%
            gather(syntype, est, -unit, -time, -outcome_id,
                   -t_int, -N)
    }
    tmpdf <- tmpdf %>%
        mutate(syntype= plyr::revalue(syntype,
                                      c("DR"="Adjusted",
                                        "Y"="SC")))
    tmpdf <- tmpdf %>%
        mutate(att=N-est)

    if("syn_method" %in% names(outcomes)) {
        tmpdf <- tmpdf %>% mutate(syntype = syn_method)
    }
    return(tmpdf)
}


plot_att <- function(outcomes, metadata, trt_unit=NULL) {
    #' Plot the estimate of the att
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #' @param trt_unit Unit to count as treated, defaults to using
    #'                 the treated column
    #' @export

    tmpdf <- compute_att(outcomes, metadata, trt_unit)
    
    p <- tmpdf %>%
    # plot the outcomes
    ggplot() +
        geom_line(aes(x=time, y=att, color=syntype), size=2) + 
    geom_vline(aes(xintercept=t_int)) +
        geom_hline(yintercept=0, linetype="dashed") +
        guides(alpha=FALSE, linetype=FALSE)

    p <- p + ggtitle("ATT Estimate") + theme_bw()
    return(p)
    }
