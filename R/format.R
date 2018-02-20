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


    ## get the number of the treated unit
    ctrl_units <- outcomes %>%
        distinct(unit, .keep_all=TRUE) %>%
        filter(unit != trt_unit) %>%
        select(unit) %>%
        as.matrix() %>%
        as.vector()
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
                is_treated=TRUE))    
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


format_data <- function(outcomes, metadata, trt_unit=1, outcome_col=NULL,
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated")) {
    #' Format "long" panel data into "wide" matrices to fit synthetic controls
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 1
    #' @param outcome_col Column name which identifies outcomes,
    #'                    if NULL then only one outcome is assumed
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return List of data to use as an argument for Synth::synth,
    #'         whether the unit was actually treated
    #' @export


    ## create new dataframe with renamed columns
    newdf <- outcomes %>%
        rename_(.dots=cols) %>%
        mutate(synthetic="N",
               "potential_outcome"=ifelse(treated,"Y(1)", "Y(0)"))
                                        # add in extra columns

    ## count number of treated units
    n_t <- outcomes %>% distinct(unit, treated) %>%
        filter(treated) %>%
        count() %>%
        as.numeric()

    ## if there is more than one treated unit, average them together
    if(n_t > 1) {
        trtavg <- newdf %>% filter(treated) %>%
            group_by_at(setdiff(names(outcomes), c(outcome, unit)))
        trtavg <- trtavg %>%
            summarise(outcome = mean(outcome)) %>%
            mutate(unit=-1) %>% 
            data.frame()
        ctrls <- newdf %>% filter(!treated)
        newdf <- rbind(trtavg, ctrls)
        trt_unit <- -1
    }

    
    if(is.null(outcome_col)) {
        out <- format_synth(newdf, metadata, trt_unit)
    } else {
        out <- format_synth_multi(newdf, metadata, outcome_col, trt_unit)
    }
    
    ## include averaged outcomes
    out$outcomes <- newdf
    out$trt_unit <- trt_unit

    ## include mapping back to original names
    oldcols <- lapply(1:length(cols), function(i) names(cols)[i])
    names(oldcols) <- sapply(1:length(cols), function(i) cols[[i]])
    out$oldcols <- oldcols
    return(out)
}


format_ipw <- function(outcomes, metadata, outcome_col=NULL,
                       cols=list(unit="unit", time="time",
                                 outcome="outcome", treated="treated")) {
    #' Format "long" panel data into "wide" matrices to fit IPW
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param outcome_col Column name which identifies outcomes,
    #'                    if NULL then only one outcome is assumed
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return List of pre and post period outcomes, and treatment indicator
    #' @export

    newdf <- outcomes %>%
        rename_(.dots=cols) %>%
        mutate(synthetic="N",
               "potential_outcome"=ifelse(treated,"Y(1)", "Y(0)"))
                                        # add in extra columns

    ## get data into wide form
    wide <- newdf %>%
        filter(time < metadata$t_int)
    if(!is.null(outcome_col)) {
        wide$cov <- interaction(wide$time, wide[[outcome_col]])
    } else {
        wide$cov <- wide$time
    }
    wide <- wide %>%
        select(unit, cov, outcome, potential_outcome, treated) %>%
        spread(cov, outcome) %>%
        select(-potential_outcome, -unit) %>%
        as.matrix()

    ## separate out into covariates and treatment vector
    X <- wide[,-1]
    trt <- wide[,1]

    ## get post-period outcomes
    post <- newdf %>%
        filter(time >= metadata$t_int)
    if(!is.null(outcome_col)) {
        post$cov <- interaction(post$time, post[[outcome_col]])
    } else {
        post$cov <- post$time
    }
    post <- post %>%
        select(unit, cov, outcome, potential_outcome, treated) %>%
        spread(cov, outcome) %>%
        select(-potential_outcome, -unit) %>%
        as.matrix()

    y <- post[,-1]

    ## average together treated units
    trtavg <- newdf %>% filter(treated) %>%
        group_by_at(setdiff(names(newdf), c("outcome", "unit"))) %>%
        summarise(outcome = mean(outcome)) %>%
        mutate(unit=-1) %>% 
        data.frame()
    ctrls <- newdf %>% filter(!treated)
    outcomes <- rbind(trtavg, ctrls)
    trt_unit <- -1
    
    return(list(X=X, trt=trt, y=y, outcomes=newdf))
}

