################################################################################
## Scripts to format panel data into matrices
################################################################################

format_synth <- function(outcomes, metadata, trt_unit=1) {
    #' Get the outcomes data into the correct form for Synth::synth
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return List of data to use as an argument for Synth::synth,
    #'         whether the unit was actually treated

    ## check if the unit is *actually* treated
    is_treated <- outcomes %>%
        filter(unit==trt_unit) %>%
        distinct(treated) %>%
        as.logical()
    if(is_treated) {
        ## filter out the counterfactual
        outcomes <- outcomes %>% filter(!((potential_outcome == "Y(0)") &
                                          (unit == trt_unit)))
    } else {
        ## filter out the actually treated unit
        outcomes <- outcomes %>% filter(!treated)
    }
    ## get the number of the treated unit
    distinct_units <- outcomes %>% distinct(unit, .keep_all=TRUE)
    ctrl_units <- distinct_units %>% filter(unit != trt_unit) %>%
        select(unit) %>% as.matrix() %>% as.vector()
    ## get the pre treatment times
    t0 <- metadata$t_int
    times <- outcomes %>% distinct(time)
    pre_trt <- times  %>%
        filter(time < t0) %>% as.matrix() %>% as.vector()

    ## Add fake predictor to make the function go through
    ## no option for no predictors
    synth_data <- outcomes %>%
        mutate(pred1=1)
    ## fit synthetic controls weights
    synth_data <-
        Synth::dataprep(foo=synth_data,
                        predictors=c("pred1"),
                        dependent="outcome",
                        unit.variable="unit",
                        time.variable="time",
                        treatment.identifier=trt_unit,
                        controls.identifier=ctrl_units,
                        time.predictors.prior=pre_trt,
                        time.optimize.ssr=pre_trt,
                        time.plot=as.matrix(times)
                        )
    return(list(synth_data=synth_data,
                is_treated=is_treated))    
    }



format_synth_multi <- function(outcomes, metadata, outcome_col, trt_unit=1) {
    #' Get multiple outcomes data as matrices
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param outcome_col Column name which identifies outcomes
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return List of data to use as an argument for Synth::synth,
    #'         whether the unit was actually treated

    ## create data_out for each outcome
    data_outs <- by(outcomes,
                    outcomes[[outcome_col]],
                    function(subset) format_synth(subset, metadata, trt_unit))

    ## quick and dirty combination of elements of list
    ## TODO: find a better way
    cat_mat <- function(name) {
        do.call(rbind, lapply(data_outs,
                              function(x) x$synth_data[[name]]))
    }


    ## combine all matrices
    data_out <- list()
    data_out$synth_data$Z0 <- cat_mat("Z0")
    data_out$synth_data$Z1 <- cat_mat("Z1")
    data_out$synth_data$Y0plot <- cat_mat("Y0plot")
    data_out$synth_data$Y1plot <- cat_mat("Y1plot")

    ## keep track of if simulation or not
    data_out$is_treated <- data_outs[[1]]$is_treated

    ## keep track of outcome names and locations in vector
    groups <- list()
    curridx <- 1
    for(i in 1:length(data_outs)) {
        len <- dim(data_outs[[i]]$synth_data$Z0)[1]
        if(is.null(len)) len <- 0
        vals <- seq(curridx, curridx + len - 1)
        groups[[names(data_outs)[i]]] <- vals
        curridx <- curridx + len
    }
    data_out$groups <- groups
    return(data_out)
}


format_data <- function(outcomes, metadata, trt_unit=1, outcome_col=NULL) {
    #' Format "long" panel data into "wide" matrices to fit synthetic controls
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param outcome_col Column name which identifies outcomes,
    #'                    if NULL then only one outcome is assumed
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return List of data to use as an argument for Synth::synth,
    #'         whether the unit was actually treated
    #' @export

    ## count number of treated units
    n_t <- outcomes %>% distinct(unit, treated) %>%
        filter(treated) %>%
        count() %>%
        as.numeric()
    print(n_t)
    ## if there is more than one treated unit, average them together
    if(n_t > 1) {
        print("GREATER")
        trtavg <- outcomes %>% filter(treated) %>%
            group_by(time, treated, outcome_id,
                     sim_num, potential_outcome,
                     synthetic) %>%
            summarise(outcome = mean(outcome)) %>%
            mutate(unit=1) %>% 
            data.frame()
        ctrls <- outcomes %>% filter(!treated)
        outcomes <- rbind(trtavg, ctrls)
        trt_unit <- 1
    }

    
    if(is.null(outcome_col)) {
        return(format_synth(outcomes, metadata, trt_unit))
    } else {
        return(format_synth_multi(outcomes, metadata, outcome_col, trt_unit))
    }
}
