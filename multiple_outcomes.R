###########################################################
## Scripts to run synthetic controls with multiple outcomes
###########################################################
library(tidyverse)
#library(hitandrun)
source("fit_synth.R")
source("entropy.R")

create_index <- function(outcomes, metadata, alpha) {
    #' Create an index of the outcomes with a weighted average
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param alpha Vector of size (n_outcomes - 1) or n_outcomes
    #'              which specifies the weights
    #'
    #' @return index of outcomes

    ## check length of alpha
    n_outcomes <- n_distinct(outcomes$outcome_id)
    if(length(alpha) == (n_outcomes - 1)) {
        alpha <- c(alpha, 1 - sum(alpha))
    } else if(length(alpha) != n_outcomes) {
        stop("alpha must be length n_outcomes or n_outcomes-1")
    }

    ## take the weighted average
    avg <- outcomes %>%
        mutate(grouping = interaction(unit, potential_outcome)) %>%
        group_by(grouping, time) %>%
        summarize(outcome=weighted.mean(outcome, alpha),
                  outcome_id="index",
                  sim_num=first(sim_num),
                  treated=first(treated),
                  potential_outcome=first(potential_outcome),
                  synthetic=first(synthetic),
                  unit=first(unit),
                  potential_outcome=first(potential_outcome)) %>%
        ungroup() %>%
        select(-grouping) %>%
        data.frame()

    return(avg)
}


fit_separate <- function(outcomes, metadata, syn_func, name, trt_unit=1) {
    #' Fit synthetic controls for multiple outcomes separately
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param syn_func Function to get synthetic controls
    #' @param name Name of the synthetic control method
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for both synthetic controls

    ## separately fit synthetic controls
    separate <- by(outcomes, outcomes$outcome_id,
                   function(df) syn_func(df, metadata, trt_unit))

    ## recombine the outcomes
    outcomes <- bind_rows(lapply(separate, function(x) x[[1]]))

    outcomes$syn_method <- paste(name, "separate", sep="_")
    ## get the weights into a list
    weights <- lapply(separate, function(x) x[[2]])
    names(weights) <- sapply(1:length(weights),
                             function(i) paste(paste(name, "sep", sep="_"),
                                               i, sep=""))
    return(list(outcomes=outcomes,
                weights=weights))
}


fit_separate_syn <- function(outcomes, metadata, trt_unit=1) {
    #' Fit synthetic controls for multiple outcomes separately
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for both synthetic controls

    return(fit_separate(outcomes, metadata, get_synth,
                        "synth", trt_unit))
}



fit_separate_ent <- function(outcomes, metadata, trt_unit=1) {
    #' Fit synthetic controls for multiple outcomes separately
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for both synthetic controls

    return(fit_separate(outcomes, metadata, get_entropy,
                        "entropy", trt_unit))
}


concat_synth_out <- function(outcomes, metadata, trt_unit=1) {
    #' Fit synthetic controls for multiple outcomes with random weights
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 1
    #'
    #' @return concatenated output of format_synth
    
    ## get synth formatted data for both outcomes
    synth_out <- by(outcomes, outcomes$outcome_id,
                    function(df) format_synth(df, metadata, trt_unit))

    ## row bind matrices and vectors
    ## Note: Treats both outcomes as the same thing
    synth_data <- list(
        X0 = do.call(rbind, lapply(synth_out, function(x) x$synth_data$X0)),
        X1 = do.call(rbind, lapply(synth_out, function(x) x$synth_data$X1)),
        Z0 = as.matrix(do.call(rbind, lapply(synth_out, function(x) x$synth_data$Z0))),
        Z1 = do.call(rbind, lapply(synth_out, function(x) x$synth_data$Z1)),
        Y0plot = do.call(rbind, lapply(synth_out, function(x) x$synth_data$Y0plot)),
        Y1plot = do.call(rbind, lapply(synth_out, function(x) x$synth_data$Y1plot))
    )

    data_out = list(synth_data=synth_data,
                    is_treated=synth_out[[1]]$is_treated)

    return(data_out)
}


fit_random <- function(outcomes, metadata, trt_unit=1) {
    #' Fit synthetic controls for multiple outcomes with random weights
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 1
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for both synthetic controls

    ## concatenate all the outcomes into one
    data_out <- concat_synth_out(outcomes, metadata, trt_unit)

    ## sample a random weight from the simplex and impute a syntehtic control
    controls <- data_out$synth_data$Y0plot
    n_c <- dim(controls)[2]
    weights <- t(simplex.sample(n_c, 1)$samples)

    ## package into input for impute_controls
    out <- list()
    out$weights <- weights
    out$controls <- controls

    ## impute the controls
    imputed <- impute_controls(outcomes, out, trt_unit)

    ## finalize output
    outcomes <- imputed$outcomes
    outcomes$syn_method = "random"
    weights <- imputed$weights

    return(list(outcomes=outcomes,
                weights=weights))
}


fit_joint <- function(outcomes, metadata, syn_func, name, trt_unit=1) {
    #' Fit synthetic controls for multiple outcomes with the same
    #' weights across different outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param syn_func Function to get synthetic controls
    #' @param name Name of synthetic control method
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for both synthetic controls


    ## make sure that the outcomes are arranged in increading outcome id
    outcomes <- outcomes %>% arrange(outcome_id)
    
    ## concateneate all the outcomes into one
    data_out <- concat_synth_out(outcomes, metadata, trt_unit)
    
    ## fit the weights jointly
    out <- syn_func(data_out)

    ## impute the controls
    imputed <- impute_controls(outcomes, out, trt_unit)

    ## finalize output
    outcomes <- imputed$outcomes
    outcomes$syn_method <- paste(name, "joint", sep="_")
    weights <- list(imputed$weights)
    names(weights) <- paste(name, "joint", sep="_")

    return(list(outcomes=outcomes,
                weights=weights))
}


fit_joint_syn <- function(outcomes, metadata, trt_unit=1) {
    #' Fit synthetic controls for multiple outcomes with the same
    #' weights across different outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for both synthetic controls

    return(fit_joint(outcomes, metadata, fit_synth_formatted,
                     "synth", trt_unit))
}


fit_joint_ent <- function(outcomes, metadata, trt_unit=1) {
    #' Fit entropy regularized synthetic controls for multiple outcomes with
    #' the same weights across different outcomes
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for both synthetic controls

    return(fit_joint(outcomes, metadata, fit_entropy_formatted,
                     "entropy", trt_unit))
}


fit_uniform <- function(outcomes, metadata, trt_unit=1) {
    #' Fit synthetic controls for multiple outcomes with uniform weights
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 1
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for both synthetic controls

    ## concatenate all the outcomes into one
    data_out <- concat_synth_out(outcomes, metadata, trt_unit)

    ## sample a random weight from the simplex and impute a syntehtic control
    controls <- data_out$synth_data$Y0plot
    n_c <- dim(controls)[2]
    weights <- rep(1/n_c, times=n_c)

    ## package into input for impute_controls
    out <- list()
    out$weights <- weights
    out$controls <- controls

    ## impute the controls
    imputed <- impute_controls(outcomes, out, trt_unit)

    ## finalize output
    outcomes <- imputed$outcomes
    outcomes$syn_method = "uniform"
    weights <- imputed$weights

    return(list(outcomes=outcomes,
                weights=weights))
}


fit_index <- function(outcomes, metadata, syn_func, name,
                      trt_unit=1, alpha=NULL) {
    #' Fit synthetic controls for multiple outcomes together with an index
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param syn_func Function to get synthetic controls
    #' @param name Name of the synthetic control method
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha Vector of size (n_outcomes - 1) or n_outcomes
    #'              which specifies the weights,
    #'              if no input then alpha = 1/n_units    
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for synthetic controls

    ## set alpha to uniform if no input
    n_outcomes <- n_distinct(outcomes$outcome_id)
    if(is.null(alpha)) {
        alpha <- rep(1 / n_outcomes, n_outcomes)
    }

    ## combine the multiple outcomes into an index
    outcomes_index <- create_index(outcomes, metadata, alpha)
    
    ## Fit synthetic controls to index, just keep weights
    syn_out <- syn_func(outcomes_index, metadata, trt_unit)
    weights <- syn_out$weights
    is_treated <- syn_out$is_treated

    ## apply weights from index set to outcomes
    syn_outcomes <- data.frame(weights) %>% # convert weights to dataframe
        mutate(unit=as.integer(rownames(weights))) %>% # get unit numbers
        inner_join(outcomes, by="unit") %>% # join with the outcomes table
        group_by(time, outcome_id) %>% # build synthetic outcome for id,time
        summarize(outcome=sum(w.weight * outcome), # compute synthetic control
                  unit=trt_unit,
                  synthetic="Y", # mark as synthetic
                  treated=is_treated, # mark it as treated
                  potential_outcome="Y(0)", # mark it as a control
                  sim_num=first(sim_num)
                  ) %>%
        ungroup()
    ## combine with outcomes
    outcomes <- rbind(outcomes, syn_outcomes)
    outcomes$syn_method <- paste(name, "index", sep="_")
    weights <- list(weights)
    names(weights) <- paste(name, "index", sep="_")
    return(list(outcomes=outcomes,
                weights=weights))
}


fit_index_syn <- function(outcomes, metadata, trt_unit=1, alpha=NULL) {
   #' Fit synthetic controls for multiple outcomes together with an index
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param syn_func Function to get synthetic controls
    #' @param name Name of the synthetic control method
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha Vector of size (n_outcomes - 1) or n_outcomes
    #'              which specifies the weights,
    #'              if no input then alpha = 1/n_units    
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for synthetic controls

    return(fit_index(outcomes, metadata, fit_synth, "synth",
                     trt_unit, alpha))
}


fit_index_ent <- function(outcomes, metadata, trt_unit=1, alpha=NULL) {
    #' Fit entropy regularized synthetic controls for multiple outcomes
    #' together with an index
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param syn_func Function to get synthetic controls
    #' @param name Name of the synthetic control method
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param alpha Vector of size (n_outcomes - 1) or n_outcomes
    #'              which specifies the weights,
    #'              if no input then alpha = 1/n_units    
    #'
    #' @return both outcomes with additional synthetic controls added
    #'         weights for synthetic controls

    return(fit_index(outcomes, metadata, fit_entropy, "entropy",
                     trt_unit, alpha))
}
