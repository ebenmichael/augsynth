###############################
## Scripts to simulate outcomes
###############################

sim_factor_model<- function(n_units, t_total, t_int, d, lambda, corr=0, fac_size=10,
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
    #' @param effect_type Type of effect
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
