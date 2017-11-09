######################################
## Scripts to evaluate simulation runs
######################################
library(tidyverse)
library(parallel)
source("multiple_outcomes.R")
source("simulate_outcomes.R")
source("fit_synth.R")

eval_run <- function(outcomes, metadata) {
    #' Compute evaluation statistics for a simulation run
    #' based on predicting the outcome for the treated unit
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #'
    #' @return Dataframe of evaluation statistics

    n_units <- n_distinct(outcomes$unit)

    ## evaluate the synthetic control
    res_list <- eval_run_unit(outcomes,
                              metadata,
                              1)
    ## extract the different kinds of results
    error <- res_list$error
    w_diff <- res_list$w_diff
    
    ## compute average mse over all time periods
    ## and squared error at the end of the time period
    mse <- error %>%
        group_by(sim_num, syn_method, outcome_id) %>%
        arrange(time) %>% # order by time
        summarise(mse_per_time=mean(diff^2),
                  final_sq_error=last(diff)^2) %>%
        ungroup()
    return(list(mse=mse,
                w_diff=w_diff))
}

loocv <- function(outcomes, metadata, ncores=1) {
    #' Compute leave out one esimates of MSE from a data set
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #' @param ncores Number of cores
    #'
    #' @return Dataframe of evaluation statistics

    ##get rid of the treated unit
    outcomes <- outcomes %>% filter(treated == FALSE)

    units <- (outcomes %>% distinct(unit))[[1]]

    ## evaluate the imputed outcomes of the estimators
    res_list <- mclapply(units,
                         function(unit) eval_run_unit(outcomes,
                                                    metadata,
                                                    unit), mc.cores=ncores)
    ## extract the different kinds of results
    error <- bind_rows(lapply(res_list, function(x) x$error))
    w_diff <- bind_rows(lapply(res_list, function(x) x$w_diff))
    ## compute loocv estimate for MSE, 
    mse <- error %>%
        group_by(sim_num, syn_method, outcome_id, unit) %>%
        summarise(mse_per_time=mean(diff^2),
                  final_sq_error=last(diff^2)) %>%
        group_by(sim_num, syn_method, outcome_id) %>%
        summarise(mse_per_time=mean(mse_per_time),
                  final_sq_error=mean(final_sq_error)) %>%
        ungroup() %>%
        as.data.frame()
    return(list(mse=mse,
                w_diff=w_diff))
}


combine_methods <- function(outcomes, metadata, trt_unit) {
    #' Fit synthetic controls separately and with an index
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #' @param trt_unit Unit that is treated (target for regression)
    #'
    #" @return synthetic controls fit both ways
    ## fit synthetic controls separately
    sep_out <- fit_separate(outcomes, metadata, trt_unit)


    ## fit the synthetic controls jointly
    joint_out <- fit_joint(outcomes, metadata, trt_unit)
    
    ## fit the average of the outcomes
    avg_out <- fit_index(outcomes, metadata, trt_unit)

    ## Use random weights
    #rand_out <- fit_random(outcomes, metadata, trt_unit)

    ## Use uniform weights
    unif_out <- fit_uniform(outcomes, metadata, trt_unit)
    
    ## combine outcomes
    outcomes <- rbind(sep_out$outcomes, joint_out$outcomes,
                      avg_out$outcomes, #rand_out$outcomes,
                      unif_out$outcomes)

    ## combine weights
    weights <- sep_out$weights
    names(weights) <- sapply(1:length(weights),
                             function(i) paste("sep", i, sep=""))
    weights$index <- avg_out$weights
    weights$joint <- joint_out$weights
    #weights$random <- rand_out$weights
    weights$uniform <- unif_out$weights

    return(list(outcomes=outcomes,
                weights=weights))
}


eval_run_unit <- function(outcomes, metadata, trt_unit) {
    #' Compute evaluation statistics for a simulation run
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #' @param trt_unit Unit that is treated (target for regression)
    #'
    #' @return Dataframe of evaluation statistics

    ## get the time of treatment
    t0 <- metadata$t_int
    comb_out <- combine_methods(outcomes, metadata, trt_unit)
    outcomes <- comb_out$outcomes
    weights <- comb_out$weights


    ## Get the error 
    error <- compute_error(outcomes, trt_unit, t0)

    ## compute distance between weights
    w_diff <- compute_weight_stats(weights)
    w_diff$unit <- trt_unit
    w_diff$sim_num <- metadata$sim_num
    return(list(error=error,
                w_diff=w_diff))
}


compute_error <- function(outcomes, trt_unit, t0) {
    #'Compute the error at each time point for one type of synthetic control
    #' @param outcomes Tidy dataframe with the outcomes with synthetic controls
    #' @param trt_unit Unit that is treated (target for regression)
    #' @param t0 Time of treatment
    #'
    #' @return The post treatment MSE of synthetic control and true control

    ## compute the error between synthetic and true control
    error <- outcomes %>%
        filter(unit == trt_unit, potential_outcome == "Y(0)", time>t0) %>%
        spread(synthetic, outcome) %>%
        mutate(diff=(Y-N)) %>%
        select(sim_num, unit, time, outcome_id, syn_method, diff)
    return(error)
}


compute_weight_stats <- function(weights) {
    #' Compute statitsics on the weights of the output
    #'
    #' @param weights List of weights for different methods
    #'
    #' @return Dataframe of statistics


    ## round the weights
    rounded_w <- lapply(weights, function(x) round(x, 3))
    w_names <- names(rounded_w)

    data <- data.frame(pair=character(0), zero_norm=numeric(0),
                       one_norm=numeric(0), two_norm=numeric(0))
    
    data <- bind_rows(
        lapply(1:(length(rounded_w) -1), function(i)
            bind_rows(
                lapply((i+1):length(rounded_w), function(j)
                    list(pair = paste(w_names[i], w_names[j], # name of pair
                                      sep="-"),
                         ## difference in support
                         zero_norm = sum((rounded_w[[i]] -
                                          rounded_w[[j]]) != 0),
                         ## L1 distance = 2 * total varaiation
                         one_norm = sum(abs(rounded_w[[i]] -
                                            rounded_w[[j]])),
                         ## L2 distance
                         two_norm = sum((rounded_w[[i]] -
                                         rounded_w[[j]])^2)
                         )
                    )
            )
            )
        )
    return(data)
}


####################### EVAL FOR FACTOR MODEL ######################

extract_list_results <- function(x) {
    #' Helper function to extract results from a list
    ## parse out lists that aren't length 3 (e.g. something failed)
    x <- x[lapply(x, length) == 3]
    mse  <- do.call(rbind, lapply(x, function(y) y[[1]]))
    w_diff <- do.call(rbind, lapply(x, function(y) y[[2]]))
    metadata <- do.call(rbind, lapply(x, function(y) y[[3]]))
    return(list(mse=mse,
                w_diff=w_diff,
                metadata=metadata))
}


eval_factor_run <- function(n_units, t_total, t_int, d, lambda, corr,
                             sig= 1, effect_type="constant", sim_num=1) {
    #' Generate data from a factor model, add a correlated outcome and evaluate
    #' @param n_units Number of units (first is always treated
    #' @param t_total Total number of time steps
    #' @param t_int Time of intervention
    #' @param d Dimension of predictiors
    #' @param lambda Effect size
    #' @param corr Correlation between outcomes
    #' @param sig Standard deviation of second outcome, defaults to 1    
    #' @param effect_type=c("constant","linear") Type of effect,
    #'                                           default: constant
    #' @param sim_num Simulation number, defaults to 1
    #'
    #' @return simulation results and metadata

    print(sim_num)
    ## simulate data
    #multi <- sim_multi_factor(n_units, t_total, t_int, d, lambda,
    #                          corr, sig, effect_type, sim_num)
    multi <- sim_factor_model2(n_units, t_total, t_int, d, lambda,
                               corr, 10, effect_type, sim_num)
    ## evaluate methods
    results <- eval_run(multi$outcomes, multi$metadata)

    ## return the results and the metadata
    return(list(mse=results$mse,
                w_diff=results$w_diff,
                metadata=multi$metadata))
}

eval_factor <- function(n_units, t_total, t_int, d, lambda, corrs, n_sims,
                             sig=1, effect_type="constant", sim_num=1, n_cores=1) {
    #' Evalaute methods on a factor model with different correlations in parallel
    #' @param n_units Number of units (first is always treated
    #' @param t_total Total number of time steps
    #' @param t_int Time of intervention
    #' @param d Dimension of predictiors
    #' @param lambda Effect size
    #' @param corrs Correlations to evalaute at
    #' @param n_sims Number of simulations per correlation
    #' @param sig Standard deviation of second outcome, defaults to 1    
    #' @param effect_type=c("constant","linear") Type of effect,
    #'                                           default: constant
    #' @param sim_num Simulation number, defaults to 1
    #' @param n_cores Number of cores to use
    #'
    #' @return data.frame for errors, weight differences, and metadata

    ## scoped function to combine every other element of a list
    extract_results <- function(x) {
        mse <- bind_rows(x[seq(1, length(x), 3)])
        w_diff <- bind_rows(x[seq(2, length(x), 3)])
        metadata <- bind_rows(x[seq(3, length(x), 3)])

        return(list(mse=mse,
                    w_diff=w_diff,
                    metadata=metadata))
    }
    
    sims <- 1:(n_sims * length(corrs))

    print(n_sims)
    ## for each value of the correaltion run n_sim simulations and evaluate
    return(extract_list_results(
        lapply(1:length(corrs),
               function(i)
                   extract_list_results(
                       mclapply(
                           1:n_sims,
                           function(j)
                               eval_factor_run(n_units,
                                               t_total,
                                               t_int,
                                               d,
                                               lambda,
                                               corrs[i],
                                               sig,
                                               effect_type,
                                               sim_num=(n_sims *
                                                        (i-1)
                                                   + j)),
                           mc.cores=n_cores)))))
}


################ EVAL FOR REAL DATA #################################

eval_loocv_run <- function(data, corr, sig, sim_num){
    #' Simulate correlated outcomes from panel data and evaluate with loocv
    #' @param data Panel data to simulate from
    #' @param corr Correlation between outcomes
    #' @param sig Standard deviation of second outcome, defaults to 1    
    #' @param sim_num Simulation number

    print(sim_num)
    ## add simulated correalted data
    mo <- add_outcome(data$outcomes, data$metadata, corr, sig, real_data=TRUE,
                      sim_num=sim_num)

    ## evaluate methods
    results <- loocv(mo$outcomes, mo$metadata, ncores=1)

    ## return the results and the metadata
    return(list(mse=results$mse,
                w_diff=results$w_diff,
                metadata=mo$metadata))
    }

eval_loocv <- function(data, corrs, n_sims, sig= 1, n_cores=1) {
    #' Generate correlated outcomes from real panel data
    #' Evaluate weights with LOOCV
    #' @param data Panel data to simulate from (outcomes and metadata)
    #' @param corrs Correlations to evalaute at
    #' @param n_sims Number of simulations per correlation
    #' @param sig Standard deviation of second outcome, defaults to 1
    #' @param n_cores Number of cores to use    
    #'
    #' @return simulation results and metadata
    
    sims <- 1:(n_sims * length(corrs))

    print(n_sims)
    ## for each value of the correaltion run n_sim simulations and evaluate
    return(extract_list_results(
        lapply(1:length(corrs),
               function(i)
                   extract_list_results(
                       mclapply(
                           1:n_sims,
                           function(j)
                               eval_loocv_run(data,
                                              corrs[i],
                                               sig,
                                               sim_num=(n_sims *
                                                        (i-1)
                                                   + j)),
                           mc.cores=n_cores)))))
}
