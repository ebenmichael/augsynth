############################################
## Scripts to plot aspects of the simulations
############################################


plot_outcomes <- function(outcomes, metadata, trt_unit=NULL) {
    #' Plot the outcomes from a study/simulation
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #' @param trt_unit Unit to count as treated, defaults to using
    #'                 the treated column
    
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
                                     "TRUE.Y.Y(0)"="Synthetic Control"
                                     )
                                   )
               ) %>%                 
    # plot the outcomes
    ggplot() +
        geom_line(aes(x=time, y=outcome, # plot outcomes vs time
                      group = grouping, color=label,
                      alpha=treated)) +
        geom_vline(aes(xintercept=t_int)) + 
        scale_color_manual(values=c("Donor Pool"="#888888",
                                    "Outcome Under Treatment"="#FDB515",
                                    "Synthetic Control"="#ED4E33",
                                    "Outcome Under Control"="#3B7EA1")) +
    scale_alpha_manual(values=c(0.05, 1)) +
        guides(alpha=FALSE) +
        if("syn_method" %in% names(outcomes)) {
            facet_grid(syn_method ~ outcome_id)
        } else {
            facet_wrap( ~outcome_id)
        } + 
        theme_bw()

    return(p)
}


plot_mse <- function(results, metadata, col) {
    #' Plot the results of a simulation
    #' @param results Dataframe of the results
    #' @param metadata Dataframe with metdata for each run
    #' @param col MSE column to plot
    #'
    #' @return plot and summarized data

    ## rename column
    results$mse <- as.numeric(results[,col])

    plot_df <- results %>%
        inner_join(metadata) %>% # join with metadata on sim_num
        select(corr, outcome_id, syn_method, mse) %>%
        group_by(corr, outcome_id, syn_method) %>%
        summarise(mean=mean(mse),
                  med=median(mse),
                  sd=sd(mse),
                  se=sd / sqrt(n())) %>%
        ungroup()

    p <- plot_df %>%
        ggplot(aes(x=corr, y=mean)) +
        #geom_boxplot(aes(x=corr, y=mse, group=interaction(corr, syn_method),
        #                 color=syn_method), outlier.shape = NA) +
        #ylim(0,1) +
        geom_line(aes(group=syn_method,
                      color=syn_method)) +
        geom_errorbar(aes(ymin=(mean-1.96 * se),
                          ymax=(mean+1.96*se),
                          color=syn_method)) +
        facet_wrap(~ outcome_id) + 
        theme_bw()
    return(list(p=p, summ=plot_df))
}


plot_w_diff <- function(w_diff, metadata, norm_name) {
    #' Plot the difference in weights of a simulation
    #' @w_diff Dataframe with stats on weight differences
    #' @param metadata Dataframe with metdata for each run
    #' @param norm_name Name of the norm to plot
    #'
    #' @return plot and summarized data

    ## rename the norm column
    plot_df <- w_diff
    names(plot_df)[names(plot_df) == norm_name] <- "value"

    ## summarise values
    plot_df <- plot_df %>%
        #gather(variable, value, -pair, -unit, -sim_num) %>%
        inner_join(metadata)%>%
        select(pair, corr, value) %>%
        group_by(pair, corr) %>%
        summarise(mean=mean(value),
                  med=median(value),
                  sd = sd(value),
                  se = sd / sqrt(n())) %>%
        ungroup() %>%
        ## separate out the paring
        separate(pair, c("pair1", "pair2")) %>%
        mutate(pair1 = factor(pair1,
                              levels=c("sep1", "sep2", "index",
                                       "joint", "uniform")),
               pair2 = factor(pair2,
                              levels=c("sep1", "sep2", "index",
                                       "joint", "uniform")))
    ## refactor the variables
    
    ## plot
    p <- plot_df %>%
        ggplot() +
        geom_line(aes(x=corr, y=mean)) +
        #geom_errorbar(aes(ymin=(mean-1.96 * se),
        #                  ymax=(mean+1.96*se))) +
        facet_grid(pair1~pair2) +
        xlab("Correlation between outcomes") +
        ylab("Distance") +
        ggtitle(paste(norm_name,
                      "distance between weights fitted by different methods")) +
        theme_bw()

    return(list(p=p, summ=plot_df))
}
