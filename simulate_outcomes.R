###############################
## Scripts to simulate outcomes
###############################
library(tidyverse)

sim_factor_model <- function(n_units, t_total, t_int, d, lambda,
                             effect_type="constant", sim_num=1,
                             frac_same=0.15) {
    #' Generate data from a factor model as in ADH 2010 (5) and Firpo 2017 (21)
    #' @param n_units Number of units (first is always treated
    #' @param t_total Total number of time steps
    #' @param t_int Time of intervention
    #' @param d Dimension of predictiors
    #' @param lambda Effect size
    #' @param effect_type=c("constant","linear") Type of effect,
    #'                                           default: constant
    #' @param sim_num Simulation number, defaults to 1
    #' @param frac_same Fraction of units that are similar to the treated unit
    #'
    #' @return \itemize{
    #'           \item{outcomes}{Tidy dataframe of outcomes}
    #'           \item{metadata}{Metadata about the simulation}
    #'         }

    ## create two sets of constants and coefficients, one for
    ## units the treted unit and units that are similar, one
    ## for other units that are different.
    
    ## generate constants and coefficients

    ## random noise for outcomes
    u <- matrix(rnorm(t_total * n_units), n_units, t_total)

    ## random noise for predictors
    v <- array(rnorm(t_total * n_units * d), dim=c(n_units, d, t_total))

    ## scale for outcome auto-correlation (one for each cluster)
    delta <- matrix(runif(t_total * 2, min=-1, max=1), 2, t_total)
    delta[,1] <- delta[,2] + 1
    ## scale for covariate correlation with previous outcome
    kappa <- matrix(runif(t_total * 2, min=-1, max=1), 2, t_total)

    ## coefficients of covariate effect on outcome
    beta <- array(runif(t_total * d * 2, min=-1/d, max=1/d),
                  dim=c(d, t_total, 2))
    beta[,,1] <- 2 * beta[,,1]

    ## covariate autocorrelation
    rho <- array(runif(t_total * d * d * 2, min=-1/d, max=1/d),
                 dim=c(d, d, t_total, 2))
    #### generate data
    ## TODO: Do this not in a for loop

    ## size of treated unit "cluster"
    n_clus <- floor(frac_same * n_units)
    
    ## initialize covariates
    ## Do we need to keep these covariates around?
    Z  <- array(0, dim=c(n_units, d, t_total))

    ## initialize  outcomes
    y <- matrix(0, n_units, t_total)

    ## compute the first values specially
    t <- 1
    Z[,,t] <- v[,,t]
    ## this computes inner products along the first axis of Z
    y[1:n_clus,t] <- Z[1:n_clus,,t] %*% beta[,t, 1] + u[1:n_clus,t]
    y[(n_clus+1):n_units,t] <- Z[(n_clus+1):n_units,,t] %*% beta[,t, 2] +
        u[(n_clus+1):n_units,t]


    ## iterate over time steps
    for(t in 2:t_total) {
        ## update covariates
        Z[1:n_clus,,t] <- kappa[1, t-1] * y[1:n_clus,t-1] +
            Z[1:n_clus,,t-1] %*% rho[,,t-1, 1] +
            v[1:n_clus,,t]

        Z[(n_clus+1):n_units,,t] <- kappa[2, t-1] * y[(n_clus+1):n_units,t-1] +
            Z[(n_clus+1):n_units,,t-1] %*% rho[,,t-1, 2] +
            v[(n_clus+1):n_units,,t]        
        ## update outcomes
        y[1:n_clus,t] <- delta[1, t-1] * y[1:n_clus,t-1] +
            Z[1:n_clus,,t] %*% beta[,t, 1] + u[1:n_clus,t]
        y[(n_clus+1):n_units,t] <- delta[2, t-1] * y[(n_clus+1):n_units,t-1] +
            Z[(n_clus+1):n_units,,t] %*% beta[,t, 2] + u[(n_clus+1):n_units,t]

    }

    ## create a new outcome for the first unit under treatment
    if(effect_type == "linear") {
    ## add in linear treatment effect for first unit
        effect <- lambda * sd(y[1,1:t_int]) *
            (1:t_total - t_int) * (1:t_total > t_int)
    } else if(effect_type == "constant") {
        effect <- lambda * sd(y[1,1:t_int]) * (1:t_total > t_int)
        
    } else {
        stop("effect must be one of c(\"linear\", \"constant\")")
    }

    y <- rbind(y[1,] + effect, y)
    # create a data frame
    outcomes <- data.frame(y)
    names(outcomes) <- 1:t_total
    outcomes <- outcomes %>% mutate(outcome_id=1,
                                    unit=c(1, 1:n_units), # number the units (0=trt)
                                    sim_num=sim_num, # keep simulation number
                                    treated=c(rep(TRUE, 2), # treatment indicator
                                              rep(FALSE, n_units-1)),
                                    # Is this Y_jt(0) or Y_jt(1)
                                    potential_outcome= c("Y(1)",
                                                         rep("Y(0)", n_units)),
                                    # is this a synthetic control
                                    synthetic="N"
                                    )
                                    
    # melt the data frame to use with tidyverse
    outcomes <- outcomes %>% gather(time, outcome, 1:t_total) %>%
        mutate(time=as.integer(time))

    ## create a dataframe with meta data about the simulation
    metadata <- data.frame(list(n_units=n_units,
                                t_total=t_total,
                                t_int=t_int,
                                d=d,
                                lambda=lambda,
                                effect_type=effect_type,
                                sim_num=sim_num,
                                frac_same=frac_same))
    return(list(outcomes=outcomes, metadata=metadata))
}


add_outcome <- function(outcomes, metadata, corr, sig=1, outcome_id=1,
                        real_data=FALSE, sim_num=NULL) {
    #' Simulate a second outcome correlated with the first outcome
    #' @param outcomes Tidy dataframe with the outcomes
    #' @param metadata Dataframe of metadata    
    #' @param corr Correlation between outcomes
    #' @param sig Standard deviation of second outcome, defaults to 1
    #' @param outcome_id Id for outcome to base new outcome on, defaults to 1
    #' @param real_data Is this a real dataset or a simulated one?
    #'                  i.e. do we know the treatment effect?
    #' @param sim_num Simulation number, if null then uses data sim num
    #'
    #' @return \itemize{
    #'           \item{outcomes}{Tidy dataframe of outcomes}
    #'           \item{metadata}{Metadata about the simulation}
    #'         }


    alpha <- sqrt(0.5 * (corr + 1))
    beta <- sqrt(1-corr^2)
    ## if simulated data, then simulated new control variables and add treatment
    if(!real_data) {
        ## simulate the extra outcome under control
        outcomes2 <- outcomes %>%
            filter(!(treated == TRUE & potential_outcome == "Y(1)")) %>%
            mutate(outcome= rnorm(n(),
                                  mean=corr * outcome,
                                  sd=sqrt(1 - corr ^ 2) * sd(outcome)),
                   outcome_id = outcome_id + 1)

        ## add the same treatment effect
        lambda <- metadata$lambda
        effect_type <- metadata$effect_type
        t_int <- metadata$t_int
        t_total <- metadata$t_total
        ## create a new outcome for the first unit under treatment
        control <- outcomes2 %>%
            filter(treated == TRUE & potential_outcome == "Y(0)")
        control <- control$outcome
        if(effect_type == "linear") {
            ## add in linear treatment effect for first unit
            effect <- lambda * sd(control[1:t_int]) *
                (1:t_total - t_int) * (1:t_total > t_int)
        } else if(effect_type == "constant") {
            effect <- lambda * sd(control[1:t_int]) * (1:t_total > t_int)
            
        } else {
            stop("effect must be one of c(\"linear\", \"constant\")")
        }

        outcomes3 <- outcomes2 %>%
            filter(treated)

        outcomes3$outcome <- control + effect
        outcomes3$potential_outcome <- "Y(1)"
        
        outcomes <- rbind(outcomes, outcomes2, outcomes3)
    } else {
        ## if it's real data, just simulated a correlated outcome for all units
        outcomes2 <- outcomes %>%
            mutate(outcome= rnorm(n(),
                                  mean= corr * sig *   
                                      (outcome) * sd(outcome),
                                  sd=sqrt(1 - corr ^ 2) * sig),
            #mutate(outcome = outcome - 2 * alpha * rnorm(n()),
                   outcome_id = outcome_id + 1)
    
        outcomes <- rbind(outcomes, outcomes2)
    }
    ## add correlation and sd info
    ## TODO: Is it necessary to create a new dataframe?
    new_metadata <- data.frame(metadata)
    new_metadata$corr <- corr
    new_metadata$sig <- sig


    if(!is.null(sim_num)) {
        outcomes$sim_num <- sim_num
        new_metadata$sim_num <- sim_num
    }
    return(list(outcomes=outcomes, metadata=new_metadata))    
}


sim_multi_factor <- function(n_units, t_total, t_int, d, lambda, corr,
                             sig= 1, effect_type="constant", sim_num=1) {
    #' Generate data from a factor model and add a correlated outcome
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
    #' @return \itemize{
    #'           \item{outcomes}{Tidy dataframe of outcomes}
    #'           \item{metadata}{Metadata about the simulation}
    #'         }

    out <- sim_factor_model(n_units, t_total, t_int, d, lambda,
                            effect_type=effect_type, sim_num=sim_num)

    out <- add_outcome(out$outcomes, out$metadata, corr, sig)
    return(out)
}


sim_factor_model2 <- function(n_units, t_total, t_int, d, lambda, corr=0, fac_size=10,
                              effect_type="constant", sim_num=1,
                              frac_same=0.3) {
    #' Generate data from a factor model as in ADH 2010 (5) and Firpo 2017 (21)
    #' @param n_units Number of units (first is always treated
    #' @param t_total Total number of time steps
    #' @param t_int Time of intervention
    #' @param d Dimension of predictors
    #' @param lambda Effect size
    #' @param corr Correlation between outcomes
    #' @param fac_size Dimension of hidden factors
    #' @param effect_type=c("constant","linear") Type of effect,
    #'                                           default: constant
    #' @param sim_num Simulation number, defaults to 1
    #' @param frac_same Fraction of units that are similar to the treated unit
    #'
    #' @return \itemize{
    #'           \item{outcomes}{Tidy dataframe of outcomes}
    #'           \item{metadata}{Metadata about the simulation}
    #'         }

    ## size of treated unit "cluster"
    n_clus <- floor(frac_same * n_units)
    
    ## create two sets of constants and coefficients, one for
    ## units the treted unit and units that are similar, one
    ## for other units that are different.
    
    ## generate constants and coefficients

    ## random noise 
    eps <- matrix(rnorm(t_total * n_units), n_units, t_total)

    eps <- array(MASS::mvrnorm(t_total * n_units, mu = rep(0,2),
                               Sigma = matrix(c(1,corr,corr,1), 2, 2)),
                 dim=c(n_units, t_total, 2))


    ## common factor across units
    delta <- runif(t_total, min=-1, max=1)
    delta <- rep(0, t_total)
    ## coefficients of covariate effect on both outcomes
    theta <- array(MASS::mvrnorm(t_total * d, mu = rep(0,2),
                               Sigma = matrix(c(1,corr,corr,1), 2, 2)),
                   dim=c(d, t_total, 2))
    
    ## generate covariates (two sets of two clusters)
    Z <- matrix(0, n_units, d)
    Z[1:n_clus, ] <- matrix(rnorm(n_clus * d), n_clus, d) + 2
    Z[(n_clus + 1):n_units, ] <- matrix(rnorm((n_units - n_clus) * d),
                                         n_units - n_clus, d)
    
    ## factor loadings for both outcomes
    lam <- array(MASS::mvrnorm(t_total * fac_size, mu = rep(0,2),
                               Sigma = matrix(c(1,corr,corr,1), 2, 2)),
                 dim=c(fac_size, t_total, 2))

    ## hidden factors (two sets for two clusters)
    mu <- matrix(0, n_units, fac_size)
    mu[1:n_clus, ] <- matrix(rnorm(n_clus * fac_size), n_clus, fac_size) + 2
    mu[(n_clus + 1):n_units, ] <- matrix(rnorm((n_units - n_clus) * fac_size),
                                         n_units - n_clus, fac_size)
    #### generate data
    ## TODO: Do this not in a for loop

    ## initialize  outcomes
    y <- array(0, dim=c(n_units, t_total, 2))

    ## iterate over time steps
    for(t in 1:t_total) {
        ## update outcomes
        y[,t,1] <- delta[t] + Z %*% theta[,t,1] + mu %*% lam[,t,1] + eps[,t,1]
        y[,t,2] <- delta[t] + Z %*% theta[,t,2] + mu %*% lam[,t,2] + eps[,t,2]
    }


    ## create a new outcome for the first unit under treatment
    if(effect_type == "linear") {
    ## add in linear treatment effect for first unit
        effect <- lambda * sd(y[1,1:t_int,]) *
            (1:t_total - t_int) * (1:t_total > t_int)
    } else if(effect_type == "constant") {
        effect <- lambda * sd(y[1,1:t_int,]) * (1:t_total > t_int)
        
    } else {
        stop("effect must be one of c(\"linear\", \"constant\")")
    }

    #y <- rbind(y[1,,] + effect, y)
    # create a data frame
    outcomes <- data.frame(do.call(rbind, lapply(1:2, function(i) y[,,i])))
    names(outcomes) <- 1:t_total
    outcomes <- outcomes %>% mutate(outcome_id=c(rep(1, n_units), rep(2, n_units)),
                                    unit=rep(1:n_units,2), # number the units (0=trt)
                                    sim_num=sim_num, # keep simulation number
                                    treated=c(rep(TRUE, 1), # treatment indicator
                                              rep(FALSE, n_units-1),
                                              rep(TRUE, 1), # treatment indicator
                                              rep(FALSE, n_units-1)),
                                    # Is this Y_jt(0) or Y_jt(1)
                                    potential_outcome= rep("Y(0)", 2 * n_units),
                                    # is this a synthetic control
                                    synthetic="N"
                                    )

    ## add in the treatment effect
    trt_outcomes <- data.frame(
        do.call(rbind,
                lapply(1:2, function(i) y[1,,i] + effect)))
    names(trt_outcomes) <- 1:t_total
    trt_outcomes <- trt_outcomes %>% mutate(outcome_id=1:2,
                                            unit=rep(1,2), # number the units (0=trt)
                                            sim_num=sim_num, # keep simulation number
                                            treated=rep(TRUE, 2),
                                        # Is this Y_jt(0) or Y_jt(1)
                                            potential_outcome= rep("Y(1)", 2),
                                        # is this a synthetic control
                                            synthetic="N"
                                    )

    outcomes <- rbind(outcomes, trt_outcomes)
    # melt the data frame to use with tidyverse
    outcomes <- outcomes %>% gather(time, outcome, 1:t_total) %>%
        mutate(time=as.integer(time))

    ## create a dataframe with meta data about the simulation
    metadata <- data.frame(list(n_units=n_units,
                                t_total=t_total,
                                t_int=t_int,
                                d=d,
                                lambda=lambda,
                                effect_type=effect_type,
                                sim_num=sim_num,
                                frac_same=frac_same,
                                fac_size=fac_size,
                                corr=corr))
    return(list(outcomes=outcomes, metadata=metadata))
}
