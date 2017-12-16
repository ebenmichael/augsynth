################################################################################
## Simulations to compare DR estimate and synth estimate
################################################################################
source("evaluate_sim.R")



dr_syn_diff <- function(outcomes, metadata, trt_unit, alpha_w, alpha_o) {
    #' Compute the difference between the DR estimate and the synth estimate
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe with info about the simulation
    #' @param trt_unit Unit that is treated (target for regression)
    #' @param alpha_w regularization parameter for weights
    #' @param alpha_o regularization parameter for outcome model
    #'
    #' @return Dataframe of evaluation statistics

    ## get the time of treatment
    t0 <- metadata$t_int


    ## fit synth weights and DR estimate
    dr_out <- get_dr(outcomes, metadata, trt_unit, alpha_w, alpha_o)


    ## get the distance between pre-period treated and synthetic control
    dists <- compute_dists(dr_out$outcomes, trt_unit, t0)


    ## compute norms of outcome regression parameters
    norms <- compute_norms(dr_out$outparams)

    ## get the distance between the synth estimate and the DR estimate
    drdist <- compute_est_dist(dr_out$outcomes, trt_unit, t0)

    ## combine it all

    alldists <- bind_cols(drdist, norms)
    alldists$l1_w <- dists$l1
    alldists$l2_w <- dists$l2
    alldists$linf_w <- dists$linf
    alldists$alpha_w <- alpha_w
    alldists$alpha_o <- alpha_o

    return(alldists)
}


compute_dists <- function(outcomes, trt_unit, t0) {
    #' compute 1, 2, and infity norms of distance between synth and treated
    #' @param outcomes Tidy dataframe with the outcomes with synthetic controls
    #' @param trt_unit Unit that is treated (target for regression)
    #' @param t0 Time of treatment
    #'
    #' @return l1, l2, and linf norms

    ## compute the error between synthetic control and treated unit
    ## in pre period
    error <- outcomes %>%
        filter(unit == trt_unit, potential_outcome == "Y(0)", time<t0) %>%
        spread(synthetic, outcome) %>%
        mutate(diff=(Y-N))
    l2 <- sqrt(sum(error$diff^2))
    l1 <- sum(abs(error$diff))
    linf <- max(abs(error$diff))
    dists <- list(l1=l1, l2=l2, linf=linf)
    return(dists)
}


compute_est_dist <- function(outcomes, trt_unit, t0) {
    #' Compute the distance between DR and synth
    #' @param outcomes Tidy dataframe with the outcomes with synthetic controls
    #' @param trt_unit Unit that is treated (target for regression)
    #' @param t0 Time of treatment
    #'
    #' @return Distance between DR and synth

    ## compute the error between synthetic control and treated unit
    ## in pre period
    dists <- outcomes %>%
        filter(unit == trt_unit, potential_outcome == "Y(0)", time>=t0) %>%
        spread(synthetic, outcome) %>%
        mutate(diff=(Y-DR)) %>%
        select(diff)
    return(dists)
}


compute_norms <- function(outparams) {
    #' Compute l1, l2 and linf norms of each row in outparams
    norm_fs <- list(l1_out=function(x) sum(abs(x)),
                    l2_out=function(x) sqrt(sum(x^2)),
                    linf_out=function(x) max(abs(x)))

    norms <- bind_rows(lapply(norm_fs,
                    function(f) apply(outparams, 1, f)))

    
    return(norms)
}



diff_sim <- function(n_units, t_total, t_int, d, lambda, alpha_ws, alpha_os,
                     n_sims, sig=1, effect_type="constant",
                     sim_num=1, n_cores=1) {
    #' Simulate from factor model and compute difference between DR and synth
    #' @param n_units Number of units (first is always treated
    #' @param t_total Total number of time steps
    #' @param t_int Time of intervention
    #' @param d Dimension of predictiors
    #' @param lambda Effect size
    #' @param alpha_ws Regularization weights for synth
    #' @param alpha_os Regularization weights for outcome model
    #' @param n_sims Number of simulations per correlation
    #' @param sig Standard deviation of second outcome, defaults to 1    
    #' @param effect_type=c("constant","linear") Type of effect,
    #'                                           default: constant
    #' @param sim_num Simulation number, defaults to 1
    #' @param n_cores Number of cores to use
    #'
    #' @return data.frame for errors, weight differences, and metadata

    sim_num <- 0

    ## for each value of regularization weights, simulate and fit
    return(bind_rows(
        lapply(1:length(alpha_ws),
               function(i)
                   bind_rows(
                       lapply(1:length(alpha_os),
                              function(j)
                                  bind_rows(
                                      mclapply(
                                      1:n_sims,
                                      function(k) {
                                          mo <- sim_factor_model2(n_units,
                                                                  t_total,
                                                                  t_int,
                                                                  d,
                                                                  lambda
                                                                  )

                                          out <-
                                              dr_syn_diff(mo$outcomes %>%
                                                          filter(outcome_id == 1),
                                                          mo$metadata, 1,
                                                          alpha_ws[i],
                                                          alpha_os[j]
                                                          )
                                          out$sim_num <- (i-1) * length(alpha_ws) +
                                              (j-1) * length(alpha_os) + k
                                          return(out)
                                          },
                                      mc.cores=n_cores
                                      )
                                  )
                              )
                   )
               )
    )
    )
}
