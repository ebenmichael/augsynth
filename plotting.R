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
                      alpha=treated, linetype=synthetic), size=2) +
        geom_vline(aes(xintercept=t_int)) + 
        scale_color_manual(values=c("Donor Pool"="#888888",
                                    "Outcome Under Treatment"="#FDB515",
                                    "Synthetic Control"="#ED4E33",
                                    "Outcome Under Control"="#3B7EA1")) +
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
        select(corr, n_units, t_total, outcome_id, syn_method, mse) %>%
        group_by(corr, n_units, t_total, syn_method) %>%
        summarise(mean=mean(mse),
                  sd=sd(mse),
                  se=sd / sqrt(n())) %>%
        ungroup()

    p <- plot_df %>% 
        ggplot(aes(x=corr, y=mean)) +
        #geom_boxplot(aes(x=corr, y=mse, group=interaction(corr, syn_method),
        #                 color=syn_method), outlier.shape = NA) +
        #ylim(0,1) +
        geom_line(aes(group=syn_method,
                      color=syn_method), size=1.5) +
        geom_errorbar(aes(ymin=(mean-1.96 * se),
                          ymax=(mean+1.96*se),
                          color=syn_method)) +
        facet_grid(n_units~t_total, labeller = label_both)
        theme_bw() +
        theme(legend.title=element_blank())

    return(list(p=p, summ=plot_df))
}

plot_basque <- function(results, metadata, col) {
    #' Plot the results of simulating from the basque dataset
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
        group_by(outcome_id, corr, syn_method) %>%
        summarise(mean=mean(mse),
                  sd=sd(mse),
                  se=sd / sqrt(n())) %>%
        ungroup()

    p <- plot_df %>% 
        ggplot(aes(x=corr, y=mean)) +
        #geom_boxplot(aes(x=corr, y=mse, group=interaction(corr, syn_method),
        #                 color=syn_method), outlier.shape = NA) +
        #ylim(0,1) +
        geom_line(aes(group=syn_method,
                      color=syn_method), size=1.5) +
        geom_errorbar(aes(ymin=(mean-1.96 * se),
                          ymax=(mean+1.96*se),
                          color=syn_method)) +
        facet_wrap(~outcome_id) + 
        theme_bw() +
        theme(legend.title=element_blank())

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
        select(n_units, t_total, pair, corr, value) %>%
        group_by(n_units, t_total, pair, corr) %>%
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
                                       "joint", "uniform"))) %>%
        filter(pair1 == "sep1", pair2 != "uniform") %>%
        rename(Weight = pair2) %>%
        mutate(Weight = plyr::revalue(Weight, c("sep2" = "Outcome 2",
                                                "joint" = "Joint",
                                                "uniform" = "Uniform",
                                                "index" = "Index")))
    
    ## refactor the variables
    ## plot
    p <- plot_df %>%
        ggplot() +
        geom_line(aes(x=corr, y=mean, group=Weight, color=Weight),
                  size=1.5) +
        #geom_errorbar(aes(ymin=(mean-1.96 * se),
        #                  ymax=(mean+1.96*se))) +
        facet_grid(n_units ~ t_total, labeller = label_both) +
        xlab("Correlation between outcomes") +
        ylab("Distance to Outcome 1") +
        ggtitle(paste(norm_name,
                      "distance between weights fitted by different methods")) +
        theme_bw()

    return(list(p=p, summ=plot_df))
}
