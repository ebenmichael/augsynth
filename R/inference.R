################################################################################
## Code for inference
################################################################################

######## RANDOMIZATION INFERENCE

est_att <- function(metadata, fitfunc, trt_unit) {
    #' Fit weights and get att estimates
    #' @param metadata with treatment time
    #' @param fitfunc Partially applied fitting function which takes in
    #'                the number of the treated unit
    #' @param trt_unit Number of the treated unit
    #'
    #' @return Dataframe with ATT estimates

    ## fit the method
    fit <- fitfunc(trt_unit)

    ## compute the att
    suppressMessages(suppressWarnings(
        att <- compute_att(fit$outcomes, metadata, trt_unit=trt_unit)
    ))
    
    return(att %>% select(time, att) %>% mutate(unit = trt_unit))
    
}


compute_stat <- function(att, pretimes, posttimes, statfunc) {
    #' Compute a test statistic statfunc({Y_t | t in times})
    #' @param att ATT estimates in a dataframe
    #' @param pretimes Vector of times in pre-period
    #' @param posttimes Vector of times in post-period
    #' @param statfunc Function to compute test stat, takes in pre and post period
    #'
    #' @return Value of test statistic

    ## limit to times
    post <- (att %>% filter(time %in% posttimes) %>% select(att))$att
    pre <- (att %>% filter(time %in% pretimes) %>% select(att))$att

    ## cat(pre, "\n")
    ## cat(post, "\n")
    ## cat("\n")
    
    ## compute test stat
    tstat <- statfunc(pre, post)
    
    return(tstat)
    
}


firpo_inf <- function(metadata, fitfunc, units, trt_unit,
                      pretimes, posttimes, statfuncs, include_treat) {
    #' Get the permutation distribution of the test stat
    #' assuming equal probability of treatment
    #' @param metadata with treatment time
    #' @param fitfunc Partially applied fitting function which takes in
    #'                the number of the treated unit
    #' @param units Numbers of the units to permute treatment around
    #' @param trt_unit Treated unit
    #' @param pretimes Vector of times in pre-period
    #' @param posttimes Vector of times in post-period
    #' @param statfuncs Function to compute test stats
    #' @param include_treat Whether to include the treated unit
    #'
    #' @return att estimates, test statistics, p-values

    ## compute atts
    atts <- bind_rows(lapply(units, function(u) est_att(metadata, fitfunc, u)))

    ## compute test statistics
    stats <- lapply(statfuncs,
                    function(statfunc)
                        by(atts, atts$unit,
                           function(df)
                               compute_stat(data.frame(df),
                                            pretimes, posttimes,
                                            statfunc)))
    stats <- lapply(stats, function(s)
        data.frame(list(unit=as.numeric(names(s)),
                        stat=sapply(s, function(x) x))))

    ## compute the "p value"
    pvals <- lapply(stats,
                    function(s) {
                        est <- mean(s[s$unit == trt_unit, 2])
                        pval <- sum(s$stat >= est) / dim(s)[1]
                    })

    return(list(atts=atts, stats=stats, pvals=pvals))
}


firpo_inf_synth <- function(outcomes, metadata, trt_unit,
                            pos=FALSE,
                            statfuncs=c(rmse_ratio,mean_abs, abs_tstat),
                            include_treat=FALSE,
                            cols=list(unit="unit", time="time",
                                      outcome="outcome", treated="treated")) {
    #' Get the permutation distribution of SC estimate test statistics
    #' assuming equal probability of treatment
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param pos Whether to only consider units with positive weights, default False
    #' @param statfuncs Function to compute test stats
    #' @param include_treat Whether to include the treated unit in randomization distribution
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return att estimates, test statistics, p-values
    #' @export

    synfunc <- function(u) {
        if(!include_treat) {
            if(u == trt_unit) {
                get_synth(outcomes,
                          metadata, u, cols=cols)
            } else {
                get_synth(outcomes[outcomes[cols$unit] != trt_unit,],
                          metadata, u, cols=cols)
            }
        } else {
            get_synth(outcomes,
                      metadata, u, cols=cols)
        }
    }

    ## ## create a fitting function for synth
    ## synfunc <- function(u) {

    ##     get_synth(outcomes,
    ##               metadata, u, cols=cols)
    ## }

    ## fit synth once
    sc <- get_synth(outcomes, metadata, trt_unit, cols)

    ## use all pre and post periods and units
    times <- unique(sc$outcomes$time)
    pretimes <- times[which(times < metadata$t_int)]
    posttimes <- times[which(times >= metadata$t_int)]
    
    units <- unique(sc$outcomes$unit)
    if(pos) {
        ctrls <- units[which(units != trt_unit)]
        ctrls <- ctrls[which(round(sc$weights, 3) > 0)]
        units <- c(trt_unit, ctrls)
    }
    ## permuate treatment label
    inf <- firpo_inf(metadata, synfunc, units, trt_unit,
                     pretimes, posttimes, statfuncs, include_treat)

    inf$sc <- sc
    return(inf)

}



wpermtest <- function(metadata, fitfunc, units, weights, trt_unit,
                      pretimes, posttimes, statfuncs) {
    #' Get the permutation distribution of the test stat
    #' assuming equal probability of treatment
    #' @param metadata with treatment time
    #' @param fitfunc Partially applied fitting function which takes in
    #'                the number of the treated unit
    #' @param units Numbers of the units to permute treatment around
    #' @param weights Weights for the permutation distribution
    #' @param trt_unit Treated unit
    #' @param pretimes Vector of times in pre-period
    #' @param posttimes Vector of times in post-period
    #' @param statfuncs Function to compute test stats
    #'
    #' @return att estimates, test statistics, p-values
    
    ## compute atts
    atts <- bind_rows(lapply(units, function(u) est_att(metadata, fitfunc, u)))

    ## compute test statistics
    stats <- lapply(statfuncs,
                    function(statfunc)
                        by(atts, atts$unit,
                           function(df)
                               compute_stat(data.frame(df),
                                            pretimes, posttimes,
                                            statfunc)))
    stats <- lapply(stats, function(s)
        data.frame(list(unit=as.numeric(names(s)),
                        stat=sapply(s, function(x) x))))

    ## estimate marginal probabilities
    n_sim <- 10000
    sims <- sapply(1:n_sim,
                   function(x)
                       sample(units, 1, prob=weights))
    probs <-sapply(units, function(x) sum(sims == x)) / n_sim

    ## compute the "p value"
    pvals <- lapply(stats,
                    function(s) {
                        est <- mean(s[s$unit == trt_unit, 2])
                        pval <- sum(probs * (s$stat >= est))
                    })

    return(list(atts=atts, stats=stats, pvals=pvals))
}



wpermtest_sc <- function(outcomes, metadata, trt_unit,
                        statfuncs=c(rmse_ratio, mean_abs, abs_tstat), 
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated")) {
    #' Get the weighted permutation distribution of SC estimate test statistics
    #' estimating p-scores with logit-link synth
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param statfuncs Function to compute test stats
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return att estimates, test statistics, p-values
    #' @export

    ipw <- format_ipw(outcomes, metadata, cols=cols)
    ## create a fitting function for synth
    synfunc <- function(u) {
        if(u == trt_unit) {
            get_synth(outcomes,
                      metadata, u, cols=cols)
        } else {
            get_synth(outcomes[outcomes[cols$unit] != trt_unit,],
                      metadata, u, cols=cols)
        }
    }

    ## fit synth once
    sc <- get_synth(outcomes, metadata, trt_unit, cols)

    ## get estimated propensity scores
    suppressMessages(ent <- get_entropy(outcomes, metadata, trt_unit, eps=sc$primal_obj,
                       cols=cols))

    pscores <- 1 / (1 + exp(-ipw$X %*% ent$dual))
    
    ## use all pre and post periods and units
    times <- unique(sc$outcomes$time)
    pretimes <- times[which(times < metadata$t_int)]
    posttimes <- times[which(times >= metadata$t_int)]
    units <- unique(sc$outcomes$unit)

    odds <- exp(ipw$X %*% ent$dual)
    ## permute treatment label
    inf <- wpermtest(metadata, synfunc, units, pscores, trt_unit,
                     pretimes, posttimes, statfuncs)

    inf$ent <- ent
    inf$pscores <- pscores
    return(inf)

}



wpermtest2 <- function(metadata, fitfunc, units, weights, trt_unit,
                      pretimes, posttimes, statfuncs) {
    #' Get the permutation distribution of the test stat
    #' assuming equal probability of treatment
    #' @param metadata with treatment time
    #' @param fitfunc Partially applied fitting function which takes in
    #'                the number of the treated unit
    #' @param units Numbers of the units to permute treatment around
    #' @param weights Weights for the permutation distribution
    #' @param trt_unit Treated unit
    #' @param pretimes Vector of times in pre-period
    #' @param posttimes Vector of times in post-period
    #' @param statfuncs Function to compute test stats
    #'
    #' @return att estimates, test statistics, p-values
    
    ## compute atts
    atts <- bind_rows(lapply(units, function(u) est_att(metadata, fitfunc, u)))

    ## compute test statistics
    stats <- lapply(statfuncs,
                    function(statfunc)
                        by(atts, atts$unit,
                           function(df)
                               compute_stat(data.frame(df),
                                            pretimes, posttimes,
                                            statfunc)))
    stats <- lapply(stats, function(s)
        data.frame(list(unit=as.numeric(names(s)),
                        stat=sapply(s, function(x) x))))

    ## compute the "p value"
    pvals <- lapply(stats,
                    function(s) {
                        est <- mean(s[s$unit == trt_unit, 2])
                        if(max(s[s$unit != trt_unit,]$stat) < est) {
                            pval <- max(weights)
                        } else {
                            pval <- sum(weights * (s[s$unit != trt_unit,]$stat >= est))
                        }
                    })

    return(list(atts=atts, stats=stats, pvals=pvals))
}



wpermtest2_sc <- function(outcomes, metadata, trt_unit,
                        statfuncs=c(rmse_ratio, mean_abs, abs_tstat), 
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated")) {
    #' Get the weighted permutation distribution of SC estimate test statistics
    #' estimating p-scores with logit-link synth
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param statfuncs Function to compute test stats
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return att estimates, test statistics, p-values
    #' @export

    ipw <- format_ipw(outcomes, metadata, cols=cols)
    ## create a fitting function for synth
    synfunc <- function(u) {
        if(u == trt_unit) {
            get_synth(outcomes,
                      metadata, u, cols=cols)
        } else {
            get_synth(outcomes[outcomes[cols$unit] != trt_unit,],
                      metadata, u, cols=cols)
        }
    }

    ## fit synth once

    sc <- get_synth(outcomes, metadata, trt_unit, cols)

    ## get estimated propensity scores
    suppressMessages(ent <- get_entropy(outcomes, metadata, trt_unit, eps=sc$primal_obj,
                       cols=cols))

    pscores <- 1 / (1 + exp(-ipw$X %*% ent$dual))
    ## use all pre and post periods and units
    times <- unique(sc$outcomes$time)
    pretimes <- times[which(times < metadata$t_int)]
    posttimes <- times[which(times >= metadata$t_int)]
    units <- unique(sc$outcomes$unit)

    ctrls <- units[which(units != trt_unit)]
    posunits <- which(round(sc$weights, 3) > 0)
    ctrls <- ctrls[posunits]
    units <- c(trt_unit, ctrls)


    odds <- exp(ipw$X %*% ent$dual)
    ## permute treatment label
    inf <- wpermtest2(metadata, synfunc, units, ent$weights[posunits], trt_unit,
                     pretimes, posttimes, statfuncs)

    inf$ent <- ent
    inf$pscores <- pscores
    return(inf)

}



importance_test <- function(metadata, fitfunc, units, probs, trt_unit,
                      pretimes, posttimes, statfuncs, n_sim) {
    #' Estimate p-value with non-uniform probabilities of treatment with
    #' importance sampling
    #' @param metadata with treatment time
    #' @param fitfunc Partially applied fitting function which takes in
    #'                the number of the treated unit
    #' @param units Numbers of the units to permute treatment around
    #' @param probs Propensity scores
    #' @param trt_unit Treated unit
    #' @param pretimes Vector of times in pre-period
    #' @param posttimes Vector of times in post-period
    #' @param statfuncs Function to compute test stats
    #' @param n_sim Number of montecarlo samples
    #'
    #' @return att estimates, test statistics, p-values
    ## compute atts
    atts <- bind_rows(lapply(units, function(u) est_att(metadata, fitfunc, u)))

    ## compute test statistics
    stats <- lapply(statfuncs,
                    function(statfunc)
                        by(atts, atts$unit,
                           function(df)
                               compute_stat(data.frame(df),
                                            pretimes, posttimes,
                                            statfunc)))
    stats <- lapply(stats, function(s) sapply(s, function(x) x))
    ## treat one unit, compute test statistic, reweight by propensities

    sims <- sapply(1:n_sim,
                   function(x) {
                       tr <- sample(length(units[-trt_unit]), 1)
                       tvec <- numeric(length(probs))
                       tvec[tr] <- 1
                       prob <- tvec * log(probs) + (1-tvec) * log(1-probs)
                       prob <- exp(sum(prob - max(prob)))
                       c(tr, prob)
                   })

    ## normalize probabilities
    sims <- t(sims)
    sims[,2] <- sims[,2] / sum(sims[,2])
    
    ## compute the p value
    pvals <- lapply(stats,
                   function(stat) {
                       notrt <- stat[-trt_unit]
                       sum((notrt[sims[,1]] > stat[trt_unit]) * sims[,2])
                   })

    return(list(atts=atts, stats=stats, pvals=pvals))
}



importance_test_sc <- function(outcomes, metadata, trt_unit,
                               statfuncs=c(rmse_ratio, mean_abs, abs_tstat),
                               n_sim=1000,
                               cols=list(unit="unit", time="time",
                                         outcome="outcome", treated="treated")) {
    #' Get the weighted permutation distribution of SC estimate test statistics
    #' estimating p-scores with logit-link synth
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param statfuncs Function to compute test stats
    #' @param n_sim Number of montecarlo samples
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return att estimates, test statistics, p-values
    #' @export

    ipw <- format_ipw(outcomes, metadata, cols=cols)
    ## create a fitting function for synth
    synfunc <- function(u) {
        if(u == trt_unit) {
            get_synth(outcomes,
                      metadata, u, cols=cols)
        } else {
            get_synth(outcomes[outcomes[cols$unit] != trt_unit,],
                      metadata, u, cols=cols)
        }
    }

    ## fit synth once

    sc <- get_synth(outcomes, metadata, trt_unit, cols)

    ## get estimated propensity scores
    probs <- sc$weights / (1 + sc$weights)
    
    ## use all pre and post periods and units
    times <- unique(sc$outcomes$time)
    pretimes <- times[which(times < metadata$t_int)]
    posttimes <- times[which(times >= metadata$t_int)]
    units <- unique(sc$outcomes$unit)

    ## permute treatment label
    inf <- importance_test(metadata, synfunc, units, probs, trt_unit,
                           pretimes, posttimes, statfuncs, n_sim)

    return(inf)

}







weighted_test <- function(metadata, fitfunc, units, probs, trt_unit,
                      pretimes, posttimes, statfuncs) {
    #' Estimate p-value with non-uniform probabilities of treatment with
    #' SCM weights
    #' @param metadata with treatment time
    #' @param fitfunc Partially applied fitting function which takes in
    #'                the number of the treated unit
    #' @param units Numbers of the units to permute treatment around
    #' @param probs Propensity scores
    #' @param trt_unit Treated unit
    #' @param pretimes Vector of times in pre-period
    #' @param posttimes Vector of times in post-period    
    #' @param statfuncs Function to compute test stats
    #'
    #' @return att estimates, test statistics, p-values
    ## compute atts
    atts <- bind_rows(lapply(units, function(u) est_att(metadata, fitfunc, u)))
    ## compute test statistics
    stats <- lapply(statfuncs,
                    function(statfunc)
                        by(atts, atts$unit,
                           function(df)
                               compute_stat(data.frame(df),
                                            pretimes, posttimes,
                                            statfunc)))
    stats <- lapply(stats, function(s) sapply(s, function(x) x))
    ## treat one unit, compute test statistic, reweight by propensities

    ## compute and normalize probabilities
    probs <- sapply(1:length(probs),
                    function(i) probs[i] * prod(1-probs[-i]))
    probs <- probs / sum(probs)
    ## compute the p value
    pvals <- lapply(stats,
                    function(stat) {
                        if(length(probs) == length(units) - 1) {
                            notrt <- stat[-trt_unit]
                            sum((notrt >= stat[trt_unit]) * probs)
                        } else {
                            sum((stat >= stat[trt_unit]) * probs)
                        }
                   })

    return(list(atts=atts, stats=stats, pvals=pvals, probs=probs))
}




valid_test <- function(metadata, fitfunc, units, probs, trt_unit,
                      pretimes, posttimes, statfuncs) {
    #' Estimate p-value with non-uniform probabilities of treatment with
    #' SCM weights
    #' @param metadata with treatment time
    #' @param fitfunc Partially applied fitting function which takes in
    #'                the number of the treated unit
    #' @param units Numbers of the units to permute treatment around
    #' @param probs Propensity scores
    #' @param trt_unit Treated unit
    #' @param pretimes Vector of times in pre-period
    #' @param posttimes Vector of times in post-period    
    #' @param statfuncs Function to compute test stats
    #'
    #' @return att estimates, test statistics, p-values
    #' @export 
    ## compute atts
    atts <- bind_rows(lapply(units, function(u) est_att(metadata, fitfunc, u)))
    ## compute test statistics
    stats <- lapply(statfuncs,
                    function(statfunc)
                        by(atts, atts$unit,
                           function(df)
                               compute_stat(data.frame(df),
                                            pretimes, posttimes,
                                            statfunc)))
    stats <- lapply(stats, function(s) sapply(s, function(x) x))
    ## treat one unit, compute test statistic, reweight by propensities

    ## compute and normalize probabilities
    probs <- sapply(1:length(probs),
                    function(i) probs[i] * prod(1-probs[-i]))
    probs <- probs / sum(probs)
    ## compute the p value
    pvals <- lapply(stats,
                   function(stat) {
                       sum((stat >= stat[trt_unit]) * probs)
                   })

    return(list(atts=atts, stats=stats, pvals=pvals, probs=probs))
}


weighted_test_sc <- function(outcomes, metadata, trt_unit,
                             statfuncs=c(rmse_ratio, mean_abs, abs_tstat),
                             cols=list(unit="unit", time="time",
                                       outcome="outcome", treated="treated")) {
    #' Get the weighted permutation distribution of SC estimate test statistics
    #' estimating p-scores with logit-link synth
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param statfuncs Function to compute test stats
    #' @param n_sim Number of montecarlo samples
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return att estimates, test statistics, p-values
    #' @export

    ipw <- format_ipw(outcomes, metadata, cols=cols)
    ## create a fitting function for synth
    synfunc <- function(u) {
        if(u == trt_unit) {
            get_synth(outcomes,
                      metadata, u, cols=cols)
        } else {
            get_synth(outcomes[outcomes[cols$unit] != trt_unit,],
                      metadata, u, cols=cols)
        }
    }

    ## fit synth once

    sc <- get_synth(outcomes, metadata, trt_unit, cols)

    ## get estimated propensity scores
    probs <- sc$weights / (1 + sc$weights)    
    units <- unique(sc$outcomes$unit)
    times <- unique(sc$outcomes$time)
    pretimes <- times[which(times < metadata$t_int)]
    posttimes <- times[which(times >= metadata$t_int)]
    
    ## permute treatment label
    inf <- weighted_test(metadata, synfunc, units, probs, trt_unit,
                         pretimes, posttimes, statfuncs)

    return(inf)

}


weighted_test_bal <- function(outcomes, metadata, trt_unit, hyperparam,
                             link=c("logit", "linear", "pos-linear"),
                             regularizer=c(NULL, "l1", "l2", "ridge", "linf"),
                             normalized=TRUE,
                             statfuncs=c(rmse_ratio, mean_abs, abs_tstat),
                             cols=list(unit="unit", time="time",
                                       outcome="outcome", treated="treated")) {
    #' Get the weighted permutation distribution of SC estimate test statistics
    #' estimating p-scores with logit-link synth
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param statfuncs Function to compute test stats
    #' @param n_sim Number of montecarlo samples
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return att estimates, test statistics, p-values
    #' @export

    ipw <- format_ipw(outcomes, metadata, cols=cols)
    ## create a fitting function for synth
    synfunc <- function(u) {
        if(u == trt_unit) {
            get_synth(outcomes,
                      metadata, u, cols=cols)
        } else {
            get_synth(outcomes[outcomes[cols$unit] != trt_unit,],
                      metadata, u, cols=cols)
        }
    }

    ## fit synth once

    bal <- get_balancer(outcomes, metadata, trt_unit, hyperparam,
                        link, regularizer, normalized, cols=cols)
    ## get back the intercept
    alpha <- -logsumexp(ipw$X[ipw$trt==0,] %*% bal$dual)

    ## get estimated p scores
    probs <- 1 / (1 + exp(-alpha - ipw$X %*% bal$dual))
    
    units <- unique(bal$outcomes$unit)
    times <- unique(bal$outcomes$time)
    pretimes <- times[which(times < metadata$t_int)]
    posttimes <- times[which(times >= metadata$t_int)]
    
    ## permute treatment label
    inf <- weighted_test(metadata, synfunc, units, probs, trt_unit,
                         pretimes, posttimes, statfuncs)

    return(inf)

}



#' test stat is ratio of RMSEs in post and pre
#' @export
rmse_ratio <- function(pre, post) {
    sqrt(mean(post^2)) /
        sqrt(mean(pre^2))
}

#' mean absolute difference
#' @export
mean_abs <- function(pre, post) {
    mean(abs(post))
}

#' absolute value of the t statistic
#' @export
abs_tstat <- function(pre, post) {
    abs( mean(post) / (sd(post) / sqrt(length(post))))
}


## test stat that is diff in diff
diff_n_diff <- function(pre, post) {
    return(mean(post) - mean(pre))
}

## test stat which is mean of absolute differences in post - pre
abs_diff <- function(pre, post) {
    return(mean(post) - mean(abs(pre)))
}

## test stat which is difference in means / rmse in pre period
norm_diff <- function(pre, post) {
    return(mean(post) / sqrt(mean(pre^2)))
}

## test stat which is the difference in mse in post and pre
mse_diff <- function(pre, post) {
    return(mean(post^2) - mean(pre^2))
}

## test stat is last value
last <- function(pre, post) {
    return(post[length(post)])
}




### Chernozhukov method
#' Weighted Least Squares Standard Errors
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param ns Number of samples from the permutation distribution to draw
#' @param trt_unit Treated unit
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' @param order of test statistic
#'
#' @return outcomes with additional synthetic control added and weights
#' @export
cherno_test <- function(outcomes, metadata, ns=1000,  q=c(2,1), trt_unit=1, 
                        cols=list(unit="unit", time="time",
                                  outcome="outcome", treated="treated")) {


    ## format data for synth
    syn_data <- format_data(outcomes, metadata, trt_unit, cols=cols)

    ## pre and post outcomes

    
    trtmat <- syn_data$synth_data$Y0plot
    ctrlmat <- syn_data$synth_data$Y1plot

    t0 <- nrow(syn_data$synth_data$Z0)
    t_final <- nrow(syn_data$synth_data$Y0plot)

    teststats <- matrix(0, nrow=ns, ncol=length(q))

    for(i in 1:ns) {

        ## sample from permutation distribution
        reorder <- sample(1:t_final, t_final)

        ## fit synth with reordered time periods

        new_synth_data <- syn_data

        new_synth_data$synth_data$X0 <- syn_data$synth_data$Y0plot[reorder,,drop=FALSE][1:t0,]
        new_synth_data$synth_data$X1 <- syn_data$synth_data$Y1plot[reorder,,drop=FALSE][1:t0,]

        syn <- fit_synth_formatted(new_synth_data)

        ## get treatment effect estimates

        att <- syn_data$synth_data$Y1plot[reorder,,drop=FALSE][(t0+1):t_final,,drop=FALSE] -
            syn_data$synth_data$Y0plot[reorder,,drop=FALSE][(t0+1):t_final,,drop=FALSE] %*% syn$weights

        teststats[i,] <- sapply(1:length(q),
                                function(j) mean(abs(att)^q[j])^(1/q[j]))
    }

    ## compute test stat for actual data
    syn <- fit_synth_formatted(syn_data)
    real_att <-  syn_data$synth_data$Y1plot[(t0+1):t_final,,drop=FALSE] -
        syn_data$synth_data$Y0plot[(t0+1):t_final,,drop=FALSE] %*% syn$weights
    real_teststat <- sapply(1:length(q),
                            function(j) mean(abs(real_att)^q[j])^(1/q[j]))
    pvals <- sapply(1:length(q), function(i) mean(teststats[,i] >= real_teststat[i]))
    
    return(pvals)
}




##### SAMPLING INFERENCE

standard_error_ <- function(metadata, fitfunc, units, trt_unit,
                           posttimes, weights=NULL) {
    #' Internal function to estimate the variance of the SC estimate
    #' @param metadata with treatment time
    #' @param fitfunc Partially applied fitting function which takes in
    #'                the number of the treated unit
    #' @param units Numbers of the control units to include when estimating the se
    #' @param trt_unit Treated unit
    #' @param posttimes Vector of times in post-period
    #' @param weights Weights to use in SC estimate, default: NULL (Doudchenko-Imbens 2017)
    #'
    #' @return att estimates, test statistics, p-values
    
    ## compute atts
    atts <- lapply(c(trt_unit, units), function(u) est_att(metadata, fitfunc, u))
    trt_att <- est_att(metadata, fitfunc, trt_unit)
    ## compute error for each time period and each unit
    errs <- sapply(atts[-1], function(att) att[att$time %in% posttimes,]$att)

    ## get standard error estimates
    if(is.null(weights)) {
        comb_func <- function(x) sqrt(sum(x^2)) / length(units)
    } else {
        ## 4 different estimators
        
        ## second moment
        comb_func <- function(x) sqrt(mean(x^2)) / sqrt(length(units))
        ses1 <- c(rep(NA, dim(trt_att)[1] - length(posttimes)), apply(errs, 1, comb_func))
        att1 <- trt_att$att
        
        ## weighted second moment
        comb_func <- function(x) sqrt(sum(x^2 * weights)) / sqrt(sum(weights)^2/sum(weights^2))
        ## comb_func <- function(x) sqrt(sum(x^2 * weights^2)) /
        ##                              sqrt(sum(weights^2)) * sqrt(sum(weights^2))
        ses2 <- c(rep(NA, dim(trt_att)[1] - length(posttimes)), apply(errs, 1, comb_func))
        att2 <- trt_att$att
        
        ## unweighted variance
        ## comb_func <- function(x) sqrt(mean((x - mean(x))^2)) / sqrt(length(units))
        comb_func <- function(x) sqrt(mean(x^2) - mean(x)^2) / sqrt(length(units))
        ses3 <- c(rep(NA, dim(trt_att)[1] - length(posttimes)), apply(errs, 1, comb_func))

        ## bias adjustment
        att3 <- trt_att$att - c(rep(0, dim(trt_att)[1] - length(posttimes)),
                                apply(errs, 1, mean))
                                    

        ## weighted variance
        ## comb_func <- function(x) sqrt(sum(x^2 * weights^2) / sum(weights^2) -
        ##                               sum(x * weights)^2 / sum(weights)) * sqrt(sum(weights^2))

        ## comb_func <- function(x) sqrt(sum((x-sum(x*weights))^2 * weights)) / sqrt(sum(weights)^2/sum(weights^2))

        comb_func <- function(x) sqrt(sum(weights * x^2) / sum(weights) -
                                      sum(weights * x)^2 / sum(weights)^2) / sqrt(sum(weights)^2/sum(weights^2))
        ses4 <- c(rep(NA, dim(trt_att)[1] - length(posttimes)), apply(errs, 1, comb_func))

        ## bias adjustment
        att4 <- trt_att$att - c(rep(0, dim(trt_att)[1] - length(posttimes)),
                                apply(errs, 1, function(x) sum(x * weights) / sum(weights)))

        ses <- c(ses1, ses2, ses3, ses4)
        weighted <- rep(c(rep(FALSE, length(ses1)), rep(TRUE, length(ses1))), 2)
        centered = c(rep(FALSE, 2 * length(ses1)), rep(TRUE, 2 * length(ses1)))
        
    }

    ## combine into one dataframe
    trt_att <- rbind(trt_att, trt_att, trt_att, trt_att)
    trt_att$att <- c(att1, att2, att3, att4)
    trt_att$se <- ses
    trt_att$weighted <- weighted
    trt_att$centered <- centered
    return(trt_att)
}



#' Estimate the variance of the SC estimate
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param trt_unit Treated unit
#' @param use_weights Whether to use weights in se estimate, default: FALSE
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' 
#' @return att estimates, test statistics, p-values
#' @export
di_standard_error <- function(outcomes, metadata, trt_unit=1, use_weights=FALSE,
                            cols=list(unit="unit", time="time",
                                      outcome="outcome", treated="treated")) {


    ## format data once
    data_out <- format_data(outcomes, metadata, trt_unit, cols=cols)

    n_c <- dim(data_out$synth_data$Z0)[2]

    t0 <- dim(data_out$synth_data$Z0)[1]
    t_final <- dim(data_out$synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)

    ## iterate over control units
    for(i in 1:n_c) {
        new_data_out <- data_out
        new_data_out$synth_data$Z0 <- data_out$synth_data$Z0[, -i]
        new_data_out$synth_data$Y0plot <- data_out$synth_data$Y0plot[, -i]

        new_data_out$synth_data$Z1 <- data_out$synth_data$Z0[, i, drop=FALSE]
        new_data_out$synth_data$Y1plot <- data_out$synth_data$Y0plot[, i, drop=FALSE]

        ## get synth weights
        syn <- fit_synth_formatted(new_data_out)

        ## estimate satt
        errs[i,] <- new_data_out$synth_data$Y1plot[(t0+1):t_final,] -
            new_data_out$synth_data$Y0plot[(t0+1):t_final,] %*% syn$weights        
    }

    ## att on actual sample
    syn <- fit_synth_formatted(data_out)
    att <- as.numeric(data_out$synth_data$Y1plot -
            data_out$synth_data$Y0plot %*% syn$weights)
    
    ## standard errors
    if(use_weights) {
        se <- sqrt(t(errs^2) %*% syn$weights)
    } else {
        se <- sqrt(apply(errs^2, 2, mean))
    }

    ## combine into dataframe
    out <- outcomes %>% distinct(time)
    out$att <- att

    out$se <- c(rep(NA, t0), se)
    
    return(out)
}







#### BOOTSTRAP


#' Bootstrap standard errors for estimated counterfactual with synth
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param n_boot Number of bootstrap samples
#' @param trt_unit Treated unit
#' @param pred_int Whether to add outcome control variance to make pred interval
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' 
#' @return att estimates, test statistics, p-values
#' @export
bootstrap_sc <- function(outcomes, metadata, n_boot, trt_unit=1, pred_int=F,
                         cols=list(unit="unit", time="time",
                                   outcome="outcome", treated="treated")) {

    ## format data once
    data_out <- format_data(outcomes, metadata, trt_unit, cols=cols)

    n_c <- dim(data_out$synth_data$Z0)[2]

    t0 <- dim(data_out$synth_data$Z0)[1] + 1
    t_final <- dim(data_out$synth_data$Y0plot)[1]
    atts <- matrix(0, n_boot, length(t0:t_final))
    ## fit on bootstrap samples
    for(b in 1:n_boot) {
        ## resample controls
        ctrls <- sample(1:n_c, n_c, replace=TRUE)

        new_data_out <- data_out
        new_data_out$synth_data$Z0 <- data_out$synth_data$Z0[, ctrls]
        new_data_out$synth_data$Y0plot <- data_out$synth_data$Y0plot[, ctrls]

        ## get synth weights
        syn <- fit_synth_formatted(new_data_out)

        ## estimate satt
        atts[b,] <- new_data_out$synth_data$Y1plot[t0:t_final,] -
            new_data_out$synth_data$Y0plot[t0:t_final,] %*% syn$weights
    }
    
    ## att on actual sample
    syn <- fit_synth_formatted(data_out)
    att <- as.numeric(data_out$synth_data$Y1plot -
                      data_out$synth_data$Y0plot %*% syn$weights)
    yhat <- as.numeric(data_out$synth_data$Y0plot %*% syn$weights)
    ## standard errors
    se <- apply(atts, 2, sd)

    if(pred_int) {
        ## estimate of control potential outcome variance
        sig2 <- sapply(t0:t_final, function(i)
            sum(syn$weights * (data_out$synth_data$Y0plot[i,] - yhat[i])^2))
        se <- sqrt(se^2 + sig2)
    }

    ## return(list(att, atts))
    ## bias
    bias <- as.numeric(apply(atts, 2, mean)) - att[t0:t_final]
                       
    ## combine into dataframe
    out <- outcomes %>% distinct(time)
    out$att <- att

    out$se <- c(rep(NA, t0-1), se)
    out$bias <- c(rep(NA, t0-1), bias)

    ## quantiles
    att_post <- att[t0:t_final]
    centered <- atts - t(matrix(att_post, nrow=length(att_post), ncol=dim(atts)[1]))
    out$cilength <- c(rep(NA, t0-1), apply(abs(centered), 2, function(x) quantile(x, .95)))

    out$lower <- c(rep(NA, t0-1), apply(atts, 2, function(x) quantile(x, .025)))
    out$upper <- c(rep(NA, t0-1), apply(atts, 2, function(x) quantile(x, .975)))

    out <- out %>% mutate(lower=2*att-upper, upper=2*att-lower)

    return(out)
}



#' Bootstrap standard errors for estimated counterfactual with synth
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param n_boot Number of bootstrap samples
#' @param trt_unit Treated unit
#' @param pred_int Whether to add outcome control variance to make pred interval
#' @param lambda Ridge hyper-parameter, if NULL then CV
#' @param scm Whether to augment SCM with ridge
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' 
#' @return att estimates, test statistics, p-values
#' @export
bootstrap_ridgeaug <- function(outcomes, metadata, n_boot, trt_unit=1, pred_int=F,
                               lambda=NULL, scm=T, 
                               cols=list(unit="unit", time="time",
                                         outcome="outcome", treated="treated")) {

    ## format data once
    data_out <- format_data(outcomes, metadata, trt_unit, cols=cols)

    ipw_out <- format_ipw(outcomes, metadata, cols=cols)
    n_c <- dim(data_out$synth_data$Z0)[2]

    t0 <- dim(data_out$synth_data$Z0)[1] + 1
    t_final <- dim(data_out$synth_data$Y0plot)[1]

    ## att on actual sample
    aug <- fit_ridgeaug_formatted(ipw_out, data_out, lambda, scm)
    att <- as.numeric(data_out$synth_data$Y1plot -
                      data_out$synth_data$Y0plot %*% aug$weights)
    yhat <- as.numeric(data_out$synth_data$Y0plot %*% aug$weights)    


    ## use one lambda
    lambda <- aug$lambda
    
    atts <- matrix(0, n_boot, length(t0:t_final))
    ## fit on bootstrap samples
    for(b in 1:n_boot) {
        ## resample controls
        ctrls <- sample(1:n_c, n_c, replace=TRUE)

        new_data_out <- data_out
        new_data_out$synth_data$Z0 <- data_out$synth_data$Z0[, ctrls]
        new_data_out$synth_data$Y0plot <- data_out$synth_data$Y0plot[, ctrls]

        new_ipw_out <- ipw_out
        new_ipw_out$X[ipw_out$trt==0,] <- ipw_out$X[ipw_out$trt==0,][ctrls,]
        if(ncol(ipw_out$y) == 1) {
            new_ipw_out$y[ipw_out$trt==0] <- ipw_out$y[ipw_out$trt==0][ctrls]
        } else {
            new_ipw_out$y[ipw_out$trt==0,] <- ipw_out$y[ipw_out$trt==0,][ctrls,]
        }
        
        ## get synth weights
        aug <- fit_ridgeaug_formatted(new_ipw_out, new_data_out, lambda, scm)

        ## estimate satt
        atts[b,] <- new_data_out$synth_data$Y1plot[t0:t_final,] -
            new_data_out$synth_data$Y0plot[t0:t_final,] %*% aug$weights
    }
    

    ## standard errors
    se <- apply(atts, 2, sd)

    if(pred_int) {
        ## estimate of control potential outcome variance
        sig2 <- sapply(t0:t_final, function(i)
            sum(aug$weights * (data_out$synth_data$Y0plot[i,] - yhat[i])^2))
        se <- sqrt(se^2 + sapply(sig2, function(x) max(x, 0)))
    }

    ## return(list(att, atts))
    ## bias
    bias <- as.numeric(apply(atts, 2, mean)) - att[t0:t_final]
                       
    ## combine into dataframe
    out <- outcomes %>% distinct(time)
    out$att <- att

    out$se <- c(rep(NA, t0-1), se)
    out$bias <- c(rep(NA, t0-1), bias)

    ## quantiles
    att_post <- att[t0:t_final]
    centered <- atts - t(matrix(att_post, nrow=length(att_post), ncol=dim(atts)[1]))
    out$cilength <- c(rep(NA, t0-1), apply(abs(centered), 2, function(x) quantile(x, .95)))

    out$lower <- c(rep(NA, t0-1), apply(atts, 2, function(x) quantile(x, .025)))
    out$upper <- c(rep(NA, t0-1), apply(atts, 2, function(x) quantile(x, .975)))

    out <- out %>% mutate(lower=2*att-upper, upper=2*att-lower)

    return(out)
}



#' Bootstrap standard errors for estimated counterfactual with synth
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param n_boot Number of bootstrap samples
#' @param hyperparam Regularization hyperparameter
#' @param trt_unit Treated unit
#' @param link Link function for weights
#' @param regularizer Dual of balance criterion
#' @param normalized Whether to normalize the weights
#' @param outcome_col Column name which identifies outcomes, if NULL then
#'                    assume only one outcome
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' @param opts Optimization options
#'        \itemize{
#'          \item{MAX_ITERS }{Maximum number of iterations to run}
#'          \item{EPS }{Error tolerance}}
#'
#' @return outcomes with additional synthetic control added and weights
#' @export
bootstrap_bal <- function(outcomes, metadata, n_boot, hyperparam, trt_unit=1,
                         link=c("logit", "linear", "pos-linear"),
                         regularizer=c(NULL, "l1", "l2", "ridge", "linf"),
                         normalized=TRUE,
                         cols=list(unit="unit", time="time",
                                   outcome="outcome", treated="treated"),
                         opts=list()) {

    ## format data once
    data_out <- format_ipw(outcomes, metadata, trt_unit, cols=cols)

    n_c <- sum(data_out$trt==0)
    n_t <- sum(data_out$trt==1)

    t0 <- dim(data_out$X)[2]+1
    t_final <- t0 + dim(data_out$y)[2] - 1
    atts <- matrix(0, n_boot, length(t0:t_final))

    ## controls and trt matrix
    Xc <- data_out$X[data_out$trt==0,,drop=FALSE]
    Xt <- data_out$X[data_out$trt==1,,drop=FALSE]

    ## outcomes
    yc <- data_out$y[data_out$trt==0,,drop=FALSE]
    yt <- data_out$y[data_out$trt==1,,drop=FALSE]


    ## fit on bootstrap samples
    for(b in 1:n_boot) {
        ## resample controls
        ctrls <- sample(1:n_c, n_c, replace=TRUE)
        ## resample treated units
        trts <- sample(1:n_t, n_t, replace=TRUE)        
        new_data_out <- data_out
        new_data_out$X <- rbind(Xt[trts,,drop=FALSE], Xc[ctrls,,drop=FALSE])
        new_data_out$trt <- c(rep(1, n_t), rep(0,n_c))
        new_data_out$y <- rbind(yt[trts,,drop=FALSE], yc[ctrls,,drop=FALSE])

        ## get balancer
        capture.output(bal <- fit_balancer_formatted(new_data_out$X, new_data_out$trt,
                                      link=link, regularizer=regularizer,
                                      hyperparam=hyperparam, normalized=normalized,
                                      opts=opts))

        ## estimate satt
        atts[b,] <- colMeans(new_data_out$y[new_data_out$trt==1,,drop=FALSE]) -
            as.numeric(t(new_data_out$y[new_data_out$trt==0,,drop=FALSE]) %*% bal$weights)

    }
    
    ## att on actual sample
    capture.output(bal <- fit_balancer_formatted(data_out$X, data_out$trt,
                                  link=link, regularizer=regularizer,
                                  hyperparam=hyperparam, normalized=normalized,
                                  opts=opts))
    att_pre <- colMeans(new_data_out$X[data_out$trt==1,,drop=FALSE]) -
        as.numeric(t(data_out$X[data_out$trt==0,,drop=FALSE]) %*% bal$weights)
    
    att_post <- colMeans(data_out$y[data_out$trt==1,,drop=FALSE]) -
        as.numeric(t(data_out$y[data_out$trt==0,,drop=FALSE]) %*% bal$weights)
    att <- c(att_pre, att_post)
    ## standard errors
    se <- apply(atts, 2, sd)

    ## return(list(att, atts))
    ## bias
    bias <- as.numeric(apply(atts, 2, mean)) - att[t0:t_final]
                       
    ## combine into dataframe
    out <- outcomes %>% distinct(time)
    out$att <- att

    out$se <- c(rep(NA, t0-1), se)
    out$bias <- c(rep(NA, t0-1), bias)

    ## quantiles
    centered <- atts - t(matrix(att_post, nrow=length(att_post), ncol=dim(atts)[1]))
    ## out$lower <- c(rep(NA, t0-1), apply(abs(centered), 2, function(x) quantile(x, .025)))
    out$cilength <- c(rep(NA, t0-1), apply(abs(centered), 2, function(x) quantile(x, .95)))

    out$lower <- c(rep(NA, t0-1), apply(atts, 2, function(x) quantile(x, .025)))
    out$upper <- c(rep(NA, t0-1), apply(atts, 2, function(x) quantile(x, .975)))

    out <- out %>% mutate(lower=2*att-upper, upper=2*att-lower)
    
    return(out)
}


#' Compute model-based weighted least squares SE estimates of the ATT
#' @param X Matrix of pre-period outcomes
#' @param y Matrix of post-period outcomes
#' @param trt Treatment indicator
#' @param weights Balancing weights
#' @param pred_int Whether to add outcome control variance to make pred interval
#' @param hc Type of HC variance estimate (-1 for homoskedastic)
wls_se_ <- function(X=NULL, y, trt, weights, pred_int, hc) {


    ## ## combine back into a panel structure
    ## n <- nrow(y)
    ## ids <- 1:n
    ## if(is.null(X)) {
    ##     t0 <- 0
    ## } else {
    ##     t0 <- dim(X)[2]
    ## }
    ## t_final <- t0 + dim(y)[2]


    ## new_weights <- numeric(0)
    ## new_weights[trt==0] <- weights
    ## new_weights[trt==1] <- 1

    ## if(!is.null(X)) {
    ##     pnl1 <- data.frame(X)
    ##     colnames(pnl1) <- 1:t0
    ##     pnl1 <- pnl1 %>% mutate(trt=trt, post=1, id=ids, weight=new_weights) %>%
    ##         gather(time, val, -trt, -post, -id, -weight)
        
    ## } else {
    ##     pnl1 <- data.frame(time=NULL, val=NULL, trt=NULL, post=NULL, id=NULL, weight=NULL)
    ## }

    ## pnl2 <- data.frame(y)
    ## colnames(pnl2) <- (t0+1):t_final
    ## pnl2 <- pnl2 %>% mutate(trt=trt, post=1, id=ids, weight=new_weights) %>%
    ##     gather(time, val, -trt, -post, -id, -weight)

    ## pnl <- bind_rows(pnl1, pnl2) %>%
    ##     mutate(time=factor(time,
    ##                        levels=1:t_final)) %>%
    ##     filter(trt==0)

    ## ## get att estimates with WLS
    ## if(t_final > 1) {
    ##     ## fit <- pnl %>% 
    ##     ##     lm(val ~ time + time:trt -1,
    ##     ##        .,
    ##     ##        weights=.$weight)
    ##     fit <- pnl %>% 
    ##         lm(val ~ time -1,
    ##            .,
    ##            weights=.$weight)
    ## } else {
    ##     ## fit <- lm(val ~ trt,
    ##     ##           pnl,
    ##     ##           weights=pnl$weight)
    ##     fit <- lm(val ~ 1,
    ##               pnl,
    ##               weights=pnl$weight)        
    ## }
    
    ## ## att <- as.numeric(coef(fit)[(t_final+1):(2*t_final)])
    ## yhat <- as.numeric(coef(fit)[1:t_final])
    ## ## se <- as.numeric(coef(summary(fit))[,2][(t_final+1):(2*t_final)])

    ## ## se <- sqrt(diag(sandwich::vcovHC(fit, type="HC", sandwich=T)))[(t_final+1):(2*t_final)]
    ## se <- sqrt(diag(sandwich::vcovHC(fit, type="HC", sandwich=T)))
    comb <- cbind(X, y)
    yhat <- as.numeric(t(comb[trt==0,]) %*% weights)
    att <- as.numeric(colMeans(comb[trt==1,,drop=F]) - yhat)
    comb <- comb[trt==0,]

    if(hc == 0) {
        ## heteroskedastic robust variance estimate (huber white)
        se <- sqrt(sapply(1:length(att),
                          function(i) sum(weights^2 * (comb[,i] - yhat[i])^2)))
    } else if(hc == 2) {
        se <- sqrt(sapply(1:length(att),
                          function(i) sum(weights^2/(1-weights) * (comb[,i] - yhat[i])^2)))
    } else if(hc == 3) {
        se <- sqrt(sapply(1:length(att),
                          function(i) sum((weights/(1-weights))^2 * (comb[,i] - yhat[i])^2)))
    }else if(hc == -1){
        ## homoskedastic estimate, different for each time
        vars_t <- sapply(1:length(att), function(i) var(comb[,i]))
        se <- sqrt(vars_t * sum(weights^2))
    }

    if(pred_int) {
        sig2 <- sapply(1:length(att), function(i) sum(weights * (comb[,i] - yhat[i])^2))
        ## sig2 <- summary(fit)$sigma^2
    } else {
        sig2 <- 0
    }
    ## add variance for prediction interval
    ## sig2 <- summary(fit)$sigma^2 / sum(trt==1)
    se <- sqrt(se^2 + sapply(sig2, function(x) max(x,0)))

    return(list(att=att, se=se))    
}



#' Weighted Least Squares Standard Errors
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param trt_unit Treated unit
#' @param pred_int Whether to add outcome control variance to make pred interval
#' @param hc Type of HC variance estimate (-1 for homoskedastic)
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#'
#' @return outcomes with additional synthetic control added and weights
#' @export
wls_se_synth <- function(outcomes, metadata, trt_unit=1, pred_int=F,
                         hc=0, 
                         cols=list(unit="unit", time="time",
                                   outcome="outcome", treated="treated")) {
    
    ## fit synth
    syn <- get_synth(outcomes, metadata, trt_unit, cols)
    
    ## format for WLS
    ipw_dat <- format_ipw(outcomes, metadata, NULL, cols)

    ## get att and standard error estimates
    att_se <- wls_se_(ipw_dat$X, ipw_dat$y, ipw_dat$trt, syn$weights, pred_int, hc)

    ## 
    ## combine into dataframe
    out <- outcomes %>% distinct(time)

    out$att <- att_se$att
    out$se <- att_se$se
    ## out$att <- c(rep(NA, ncol(ipw_dat$X)),
    ##              att_se$att)

    ## out$se <- c(rep(NA, ncol(ipw_dat$X)),
    ##             att_se$se)

    return(out)
}





#' Weighted Least Squares Standard Errors
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param trt_unit Treated unit
#' @param pred_int Whether to add outcome control variance to make pred interval
#' @param hc Type of HC variance estimate (-1 for homoskedastic)
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#'
#' @return outcomes with additional synthetic control added and weights
#' @export
wls_se_ridgeaug <- function(outcomes, metadata, trt_unit=1, pred_int=F, hc=0,
                            lambda=NULL, scm=T, 
                            cols=list(unit="unit", time="time",
                                      outcome="outcome", treated="treated")) {
    
    ## fit synth
    aug <- get_ridgeaug(outcomes, metadata, trt_unit, lambda, scm, cols)

    ## format for WLS
    ipw_dat <- format_ipw(outcomes, metadata, NULL, cols)

    ## get att and standard error estimates
    att_se <- wls_se_(ipw_dat$X, ipw_dat$y, ipw_dat$trt, aug$weights, pred_int, hc)

    ## 
    ## combine into dataframe
    out <- outcomes %>% distinct(time)

    out$att <- att_se$att
    out$se <- att_se$se
    ## out$att <- c(rep(NA, ncol(ipw_dat$X)),
    ##              att_se$att)

    ## out$se <- c(rep(NA, ncol(ipw_dat$X)),
    ##             att_se$se)

    return(out)
}



#' Bootstrap standard errors for estimated counterfactual with synth
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param hyperparam Regularization hyperparameter
#' @param trt_unit Treated unit
#' @param link Link function for weights
#' @param regularizer Dual of balance criterion
#' @param normalized Whether to normalize the weights
#' @param outcome_col Column name which identifies outcomes, if NULL then
#'                    assume only one outcome
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' @param opts Optimization options
#'        \itemize{
#'          \item{MAX_ITERS }{Maximum number of iterations to run}
#'          \item{EPS }{Error tolerance}}
#'
#' @return outcomes with additional synthetic control added and weights
#' @export
wls_se_bal <- function(outcomes, metadata, hyperparam, trt_unit=1,
                         link=c("logit", "linear", "pos-linear"),
                         regularizer=c(NULL, "l1", "l2", "ridge", "linf"),
                         normalized=TRUE,
                         cols=list(unit="unit", time="time",
                                   outcome="outcome", treated="treated"),
                         opts=list()) {


    ## format data once
    ipw_dat <- format_ipw(outcomes, metadata, NULL, cols)

    
    capture.output(bal <- fit_balancer_formatted(ipw_dat$X, ipw_dat$trt,
                                                 link=link, regularizer=regularizer,
                                                 hyperparam=hyperparam, normalized=normalized,
                                                 opts=opts))

    ## get att and standard error estimates
    att_se <- wls_se_(ipw_dat$X, ipw_dat$y, ipw_dat$trt, bal$weights)

    ## 
    ## combine into dataframe
    out <- outcomes %>% distinct(time)
    out$att <- att_se$att

    out$se <- att_se$se

    return(out)
    }




#' Compute standard errors from survey glm
#' @param X Matrix of pre-period outcomes
#' @param y Matrix of post-period outcomes
#' @param trt Treatment indicator
#' @param weights Balancing weights
#' @param pred_int Whether to add outcome control variance to make pred interval
svyglm_se_ <- function(X=NULL, y, trt, weights, pred_int) {


    ## combine back into a panel structure
    n <- nrow(y)
    ids <- 1:n
    if(is.null(X)) {
        t0 <- 0
    } else {
        t0 <- dim(X)[2]
    }
    t_final <- t0 + dim(y)[2]


    new_weights <- numeric(0)
    new_weights[trt==0] <- weights
    new_weights[trt==1] <- 1

    if(!is.null(X)) {
        pnl1 <- data.frame(X)
        colnames(pnl1) <- 1:t0
        pnl1 <- pnl1 %>% mutate(trt=trt, post=1, id=ids, weight=new_weights) %>%
            gather(time, val, -trt, -post, -id, -weight)
        
    } else {
        pnl1 <- data.frame(time=NULL, val=NULL, trt=NULL, post=NULL, id=NULL, weight=NULL)
    }

    
    pnl2 <- data.frame(y)
    colnames(pnl2) <- (t0+1):t_final
    pnl2 <- pnl2 %>% mutate(trt=trt, post=1, id=ids, weight=new_weights) %>%
        gather(time, val, -trt, -post, -id, -weight)

    pnl <- bind_rows(pnl1, pnl2) %>%
        mutate(time=factor(time,
                           levels=1:t_final)) %>%
        filter(trt==0)

    ## make desing
    design <- survey::svydesign(id=~1, weights=~weight, data=pnl)
    
    ## get att estimates with WLS
    if(t_final > 1) {
        ## fit <- pnl %>% 
        ##     lm(val ~ time + time:trt -1,
        ##        .,
        ##        weights=.$weight)
        fit <- survey::svyglm(val ~ time -1, design)
    } else {
        ## fit <- lm(val ~ trt,
        ##           pnl,
        ##           weights=pnl$weight)
        fit <- survey::svyglm(val ~ 1, design)
    }
    ## att <- as.numeric(coef(fit)[(t_final+1):(2*t_final)])
    yhat <- as.numeric(coef(fit))
    ## se <- as.numeric(coef(summary(fit))[,2][(t_final+1):(2*t_final)])

    se <- coef(summary(fit))[,2]

    comb <- cbind(X, y)
    att <- colMeans(comb[trt==1,,drop=F]) - yhat
    comb <- comb[trt==0,]
    if(pred_int) {
        sig2 <- sapply(1:length(att), function(i) sum(weights * (comb[,i] - yhat[i])^2))
    } else {
        sig2 <- 0
    }
    se <- sqrt(se^2 + sig2)

    return(list(att=att, se=se))    
}



#' svyglm standard errors
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param trt_unit Treated unit
#' @param pred_int Whether to add outcome control variance to make pred interval
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#'
#' @return outcomes with additional synthetic control added and weights
#' @export
svyglm_se_synth <- function(outcomes, metadata, trt_unit=1, pred_int=F,
                         cols=list(unit="unit", time="time",
                                   outcome="outcome", treated="treated")) {
    
    ## fit synth
    syn <- get_synth(outcomes, metadata, trt_unit, cols)

    ## format for WLS
    ipw_dat <- format_ipw(outcomes, metadata, NULL, cols)

    ## get att and standard error estimates
    att_se <- svyglm_se_(ipw_dat$X, ipw_dat$y, ipw_dat$trt, syn$weights, pred_int)

    ## 
    ## combine into dataframe
    out <- outcomes %>% distinct(time)

    out$att <- att_se$att
    out$se <- att_se$se
    ## out$att <- c(rep(NA, ncol(ipw_dat$X)),
    ##              att_se$att)

    ## out$se <- c(rep(NA, ncol(ipw_dat$X)),
    ##             att_se$se)

    return(out)
}






#' svyglm standard errors
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param trt_unit Treated unit
#' @param pred_int Whether to add outcome control variance to make pred interval
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#'
#' @return outcomes with additional synthetic control added and weights
#' @export
svyglm_se_synth <- function(outcomes, metadata, trt_unit=1, pred_int=F,
                         cols=list(unit="unit", time="time",
                                   outcome="outcome", treated="treated")) {
    
    ## fit synth
    syn <- get_synth(outcomes, metadata, trt_unit, cols)

    ## format for WLS
    ipw_dat <- format_ipw(outcomes, metadata, NULL, cols)

    ## get att and standard error estimates
    att_se <- svyglm_se_(ipw_dat$X, ipw_dat$y, ipw_dat$trt, syn$weights, pred_int)

    ## 
    ## combine into dataframe
    out <- outcomes %>% distinct(time)

    out$att <- att_se$att
    out$se <- att_se$se
    ## out$att <- c(rep(NA, ncol(ipw_dat$X)),
    ##              att_se$att)

    ## out$se <- c(rep(NA, ncol(ipw_dat$X)),
    ##             att_se$se)

    return(out)
}



#' Use leave out one estimates (placebo gaps) to estimate unit-level variance
#' Do this for ridge-augmented synth
#' @param outcomes Tidy dataframe with the outcomes and meta data
#' @param metadata Dataframe of metadata
#' @param trt_unit Treated unit
#' @param lambda Ridge hyper-parameter, if NULL use CV
#' @param scm Include SCM or not
#' @param use_weights Whether to use weights in se estimate, default: FALSE
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' 
#' @return att estimates, test statistics, p-values
#' @export
loo_se_ridgeaug <- function(outcomes, metadata, trt_unit=1, lambda=NULL,
                            scm=T, use_weights=T,
                            cols=list(unit="unit", time="time",
                                      outcome="outcome", treated="treated")) {


    ## format data once
    data_out <- format_data(outcomes, metadata, trt_unit, cols=cols)
    ipw_dat <- format_ipw(outcomes, metadata, cols=cols)

    n_c <- dim(data_out$synth_data$Z0)[2]

    t0 <- dim(data_out$synth_data$Z0)[1]
    t_final <- dim(data_out$synth_data$Y0plot)[1]
    errs <- matrix(0, n_c, t_final - t0)

    ## att on actual sample
    aug_t <- fit_ridgeaug_formatted(ipw_dat, data_out, lambda, scm)
    att <- as.numeric(data_out$synth_data$Y1plot -
            data_out$synth_data$Y0plot %*% aug_t$weights)

    lam <- aug$lambda
    
    ## iterate over control units
    for(i in 1:n_c) {

        ## reset synth data to make a control a treated
        new_data_out <- data_out
        new_data_out$synth_data$Z0 <- data_out$synth_data$Z0[, -i]
        new_data_out$synth_data$Y0plot <- data_out$synth_data$Y0plot[, -i]

        new_data_out$synth_data$Z1 <- data_out$synth_data$Z0[, i, drop=FALSE]
        new_data_out$synth_data$Y1plot <- data_out$synth_data$Y0plot[, i, drop=FALSE]

        ## reset ipw data to change treatment assignment
        new_ipw_dat <- ipw_dat
        new_ipw_dat$X <- ipw_dat$X[ipw_dat$trt==0,]
        new_ipw_dat$y <- ipw_dat$y[ipw_dat$trt==0,]
        new_ipw_dat$trt <- numeric(nrow(new_ipw_dat$X))
        new_ipw_dat$trt[i] <- 1

        ## get ridge_aug weights
        aug <- fit_ridgeaug_formatted(new_ipw_dat, new_data_out, lam, scm)

        ## estimate satt
        errs[i,] <- new_data_out$synth_data$Y1plot[(t0+1):t_final,] -
            new_data_out$synth_data$Y0plot[(t0+1):t_final,] %*% aug$weights
    }

    ## standard errors
    if(use_weights) {
        se <- sqrt(t(errs^2) %*% aug_t$weights^2)
    } else {
        se <- sqrt(apply(errs^2, 2, mean))
    }

    ## combine into dataframe
    out <- outcomes %>% distinct(time)
    out$att <- att

    out$se <- c(rep(NA, t0), se)
    
    return(out)
}
