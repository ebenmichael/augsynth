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
                      pretimes, posttimes, statfuncs) {
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
                            cols=list(unit="unit", time="time",
                                      outcome="outcome", treated="treated")) {
    #' Get the permutation distribution of SC estimate test statistics
    #' assuming equal probability of treatment
    #' @param outcomes Tidy dataframe with the outcomes and meta data
    #' @param metadata Dataframe of metadata
    #' @param trt_unit Unit that is treated (target for regression), default: 0
    #' @param pos Whether to only consider units with positive weights, default False
    #' @param statfuncs Function to compute test stats
    #' @param cols Column names corresponding to the units,
    #'             time variable, outcome, and treated indicator
    #'
    #' @return att estimates, test statistics, p-values
    #' @export

    synfunc <- function(u) {
        if(u == trt_unit) {
            get_synth(outcomes,
                      metadata, u, cols=cols)
        } else {
            get_synth(outcomes[outcomes[cols$unit] != trt_unit,],
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
                     pretimes, posttimes, statfuncs)

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
#' @param fitfunc Partially applied fitting function which takes in
#'                the number of the treated unit
#' @param units Numbers of the units to include when estimating the se
#' @param trt_unit Treated unit
#' @param use_weights Whether to use weights in se estimate, default: FALSE
#' @param cols Column names corresponding to the units,
#'             time variable, outcome, and treated indicator
#' 
#' @return att estimates, test statistics, p-values
#' @export
standard_error <- function(outcomes, metadata, trt_unit, use_weights=FALSE,
                            cols=list(unit="unit", time="time",
                                      outcome="outcome", treated="treated")) {
        
    ## partially applied synth fitting function
    ## synfunc <- function(u) {
    ##     if(u == trt_unit) {
    ##         get_synth(outcomes,
    ##                   metadata, u, cols=cols)
    ##     } else {
    ##         get_synth(outcomes[outcomes[cols$unit] != trt_unit,],
    ##                   metadata, u, cols=cols)
    ##     }
    ## }

    ## create a fitting function for synth
    synfunc <- function(u) {

        get_synth(outcomes,
                  metadata, u, cols=cols)
    }

    ## fit synth once
    sc <- get_synth(outcomes, metadata, trt_unit, cols)

    ## use all pre and post periods and units
    times <- unique(sc$outcomes$time)
    pretimes <- times[which(times < metadata$t_int)]
    posttimes <- times[which(times >= metadata$t_int)]
    units <- unique(sc$outcomes$unit)
    units <- units[units != trt_unit]

    ## use synth weights
    if(use_weights) {
        weights <- sc$weights
        ## restrict to units with positive weight
        ## units <- units[round(weights,3) > 0]
        ## weights <- weights[round(weights,3) > 0]
    } else {
        weights <- NULL
    }
    
    ## get standard error
    se <- standard_error_(metadata, synfunc, units, trt_unit,
                          posttimes, weights)
    ## se$sc <- sc
    return(se)

}
