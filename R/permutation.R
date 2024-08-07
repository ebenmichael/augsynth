

#' Run augsynth letting each unit be the treated unit.  This
#' implements the abadie placebo test (e.g., from tobacco paper) but
#' with the augsynth package.
#'
#' @param tx.id Name/ID of treatment group (not stored in ascm)
#'
#' @return Dataframe With each row corresponding to a unit (including
#'   the original treated unit) and each column corresponding to a
#'   time point.  Each entry is the estimated impact for that unit at
#'   that time, after fitting augsynth on that unit as the tx unit and
#'   the other units as possible controls.
#'
#' @noRd
get_placebo_gaps = function( ascm, att = TRUE ) {

    tx_id <- ascm$trt_unit
    wide_data <- ascm$data
    synth_data <- ascm$data$synth_data
    Z <- wide_data$Z
    control_ids = which( wide_data$trt == 0 )
    all_ids = 1:length(wide_data$trt)
    t_final = length( wide_data$time )

    ests <- vapply(all_ids, function(i) {
        new_data <- swap_treat_unit(wide_data, Z, i)
        new_ascm <- do.call(augsynth:::fit_augsynth_internal,
                            c(list(wide = new_data$wide,
                                   synth_data = new_data$synth_data,
                                   Z = new_data$Z,
                                   progfunc = ascm$progfunc, scm = ascm$scm,
                                   fixedeff = ascm$fixedeff),
                              ascm$extra_args))
        est <- predict(new_ascm, att = att )
        est
    }, numeric(t_final) )

    # Verify we recover the original treated unit by seeing if the
    # estimates match our original fit model's estimates
    dim( ests )
    pds = as.numeric( predict( ascm, att = att ) )
    pds
    stopifnot( all( round( ests[ , which( wide_data$trt == 1 ) ] - pds, digits=7 ) == 0 ) )

    ests = as.data.frame( t( ests ) )
    dim( ests )
    ests$ID = tx_id
    ests$ID[ control_ids ] = rownames( ascm$weights )

    ests %>% dplyr::select( ID, everything() )
}

#### Generate donor unit from fit synth model ####

#' Take inner data and change which unit is marked as 'treated'
#'
#' @noRd
swap_treat_unit = function (wide_data, Z, i) {
    wide_data$trt <- rep( 0, length( wide_data$trt ) )
    wide_data$trt[i] = 1

    X0 <- wide_data$X[wide_data$trt == 0, , drop = F]
    x1 <- matrix( wide_data$X[wide_data$trt == 1, , drop = F], ncol = 1 )
    y0 <- wide_data$y[wide_data$trt == 0, , drop = F]
    y1 <- matrix( wide_data$y[wide_data$trt == 1, , drop = F], ncol = 1 )
    new_synth_data <- list()
    new_synth_data$Z0 <- t(X0)
    new_synth_data$X0 <- t(X0)
    new_synth_data$Z1 <- x1
    new_synth_data$X1 <- x1

    wide_data = wide_data[ c( "trt", "X", "y" ) ]

    return(list(wide_data = wide_data,
                synth_data = new_synth_data,
                Z = Z))
}



#' Calculate MDES
#'
#' Use our method for Calculating SEs, p-values, and MDEs by looking
#' at the distribution of impacts across the donor units.
#'
#' @param lest Entity by year estimate estimates for all treatment and control entities
#'   (original tx and all the comparisons as pseudo-treated).  Columns
#'   of entity ID (tx_col), treatment flag (trt), time period (time_col),
#'   and estimated tx-synthetic control difference (ATT).
#' @param treat_year The time period in which the treatment was introduced.
#' @param  tx_col The name to the column containing the treatment variable.
#' @param time_col The name of the column containing the time variable.
#' @return List of three dataframes, one of MDES, one of RMSPEs, and
#'   one of info on the r-statistics (the ATTs divided by the RMSPEs).
#'
#' @noRd
calculate_MDES_table = function( lest, treat_year, tx_col, time_col ) {


    stopifnot( all( c(  tx_col, "trt", time_col, "ATT" ) %in%
                        names( lest ) ) )

    RMSPEs = calculate_RMSPE( lest, treat_year, tx_col, time_col )

    # merge calculated pre-intervention RMSPE for each county to data set
    # of gaps for each year for each county

    # Divide that year's gap by pre-intervention RMSPE.
    rstatcalc <- full_join( lest, RMSPEs, by = tx_col ) %>%
        mutate( rstat=ATT/RMSPE )

    # Calculate the permutation p-values (permuting the r statistics)
    pvalues = rstatcalc %>%
        group_by( !!as.name(time_col) ) %>%
        summarise( p_rstat = mean( abs( rstat[trt == 1] ) <= abs( rstat ) ),
                   SE_rstat = sd( rstat[ trt != 1 ] ) * RMSPE[ trt == 1 ],
                   p_gap = mean( abs( ATT[ trt == 1 ] ) <= abs( ATT ) ),
                   SE_gap = sd( ATT[ trt != 1 ] ),
                   .groups = "drop" )

    # attach R statistic for the treated tract.
    rstattract <- rstatcalc %>%
        filter(trt==1) %>%
        dplyr::select( -!!as.name(tx_col), -trt )

    # Make table of main results with permutation p-values and add
    # indicator of which years are treated and which not.
    main_results <- rstattract %>%
        full_join( pvalues, by=time_col ) %>%
        mutate( tx = ifelse( !!as.name(time_col) >= treat_year, 1, 0 ) )


    # Clean up column ordering and make look nice
    main_results = main_results %>%
        relocate( !!as.name(time_col), tx, Yobs ) %>%
        dplyr::select( -RMSPE )

    # Bundle and return all results
    list( MDES_table = main_results,
          RMSPEs = RMSPEs,
          rstatcalc=rstatcalc )
}


#' Calculate RMSPE for all units on lagged outcomes
#'
#' This calculates the RMSPE of the difference between estimated and
#' observed outcome, averaged across the pre-treatment years, for each
#' county
#'
#' @param lest Entity by year estimate estimates for all treatment and control entities
#'   (original tx and all the comparisons as pseudo-treated).  Columns
#'   of entity ID (tx_col), treatment flag (trt), time period (time_col),
#'   and estimated tx-synthetic control difference (ATT).
#' @param treat_year The time period in which the treatment was introduced.
#' @param  tx_col The name to the column containing the treatment variable.
#' @param time_col The name of the column containing the time variable.
#'
#' @noRd
calculate_RMSPE = function( lest, treat_year, tx_col, time_col ) {
    stopifnot( !is.null( lest$ATT ) )

    RMSPEs = lest %>% filter( !!as.name(time_col) < treat_year ) %>%
        group_by( !!as.name(tx_col) ) %>%
        summarise( RMSPE = sqrt( mean( ATT^2 ) ), .groups = "drop" )

    return(RMSPEs)
}

#' Estimate robust SE
#'
#' Given a list of impact estimates (all assumed to be estimating a
#' target impact of 0), calculate the modified standard deviation by
#' looking an in inner-quartile range rather than the simple standard
#' deviation.  This reduces the effect of small numbers of extreme
#' outliers, and focuses on central variation.  This then gets
#' converted to a standard error (standard deviation of the results,
#' in effect.)
#'
#' @param placebo_estimates The list of impact estimates to use for
#'   estimating the SE.
#' @param k Number of obs to drop from each tail (so 2k obs dropped
#'   total)
#' @param beta Proportion of obs to drop from each tail. E.g., 95\% to
#'   drop 5\% most extreme.  Only one of beta or k can be non-null.
#' @param round_beta TRUE means adjust beta to the nearest integer
#'   number of units to drop (the Howard method).  FALSE means
#'   interpolate quantiles using R's quantile function.
#'
#' @noRd
estimate_robust_SE = function( placebo_estimates, k = NULL, beta=NULL,
                               round_beta = FALSE ) {

    stopifnot( is.null(k) != is.null( beta ) )

    n = length( placebo_estimates )

    if ( is.null( beta ) ) {
        #alpha = (2*k-1) / (2*n)    # Method A
        alpha = k/( n-1 )    # Method B
        beta = 1 - 2*alpha
        q = sort( placebo_estimates )[ c(1+k, n - k) ]
    } else {
        alpha = (1 - beta)/2
        k = alpha * (n-1)  # QUESTION: Why n-1 here????
        if ( round_beta ) {
            k = pmax( 1, round( k ) )
            alpha = k/(n-1)
            beta = 1 - 2*alpha
            q = sort( placebo_estimates )[ c(1+k, n - k) ]
        } else {
            q = quantile( placebo_estimates, probs = c( alpha, 1 - alpha ) )
        }
    }

    del = as.numeric( diff( q ) )
    z = -qnorm( alpha )
    SE_imp = del / (2*z)

    res = data.frame( beta = beta, k = k,
                      q_low = q[[1]], q_high = q[[2]],
                      range = del, z = z, SE_imp = SE_imp )
    res
}


#' Construct an organized dataframe with outcome data in long format
#' from augsynth object
#'
#' @noRd
get_long_data <- function( augsynth ) {

    wide_data <- augsynth$data

    tx_id <- augsynth$trt_unit
    control_ids = which( wide_data$trt == 0 )

    all_ids = 1:nrow(wide_data$X)
    all_ids[ control_ids ] = rownames( augsynth$weights )
    all_ids[ which( wide_data$trt == 1 ) ] = tx_id

    df <- bind_cols(wide_data$X, wide_data$y) %>%
        mutate(!!as.name(augsynth$unit_var) := all_ids) %>%
        pivot_longer(!augsynth$unit_var,
                     names_to = augsynth$time_var,
                     values_to = 'Yobs',
        ) %>%
        mutate(!!as.name(augsynth$time_var) := as.numeric(!!as.name(augsynth$time_var)),
               ever_Tx = ifelse(!!as.name(augsynth$unit_var) == augsynth$trt_unit, 1, 0))

    df
}



add_placebo_distribution <- function(augsynth) {

    # Run permutations
    ests <- get_placebo_gaps(augsynth, att = FALSE)
    time_cols = 2:ncol(ests)

    ests$trt = augsynth$data$trt

    lest = ests %>%
        pivot_longer( cols = all_of( time_cols ),
                      names_to = augsynth$time_var, values_to = "Yhat" ) %>%
        mutate(!!as.name(augsynth$time_var) := as.numeric(!!as.name(augsynth$time_var)))

    df <- get_long_data( augsynth )

    ##### Make dataset of donor units with their weights #####
    units = dplyr::select( df, !!as.name(augsynth$unit_var), ever_Tx ) %>%
        unique()

    tx_col = augsynth$unit_var
    weights = data.frame( tx_col = rownames( augsynth$weights ),
                          weights = augsynth$weights,
                          stringsAsFactors = FALSE)
    colnames(weights)[1] <- augsynth$unit_var

    units = merge( units, weights, by=augsynth$unit_var, all=TRUE )

    # Zero weights for tx units.
    units$weights[ units$ever_Tx == 1 ] = 0

    # confirm that we have placebos for every entity over every observed time period
    stopifnot((ncol(augsynth$data$X) + ncol(augsynth$data$y)) * nrow(augsynth$data$y) == nrow(lest))

    lest = rename( lest, !!as.name(augsynth$unit_var) := ID )

    nn = nrow(lest)
    lest = left_join( lest,
                      df[ c(augsynth$unit_var, augsynth$time_var, "Yobs") ], # issue is that this is calling for original data
                      by = names(lest)[names(lest) %in% names(df[ c(augsynth$unit_var, augsynth$time_var, "Yobs") ])] )
    stopifnot( nrow( lest ) == nn )


    # Impact is difference between observed and imputed control-side outcome
    lest$impact = lest$Yobs - lest$Yhat

    #### Make the actual treatment result information #####

    # Get our observed series (the treated unit)
    T_tract = filter( lest, trt == 1 )

    # The raw donor pool average series
    averages = df %>%
        filter( ever_Tx == 0 ) %>%
        group_by( !!as.name(augsynth$time_var) ) %>%
        summarise( raw_average = mean( .data$Yobs ),
                   .groups = "drop" )
    T_tract = left_join( T_tract, averages, by = augsynth$time_var )

    ##### Calculate ATT and MDES #####

    lest = rename( lest, ATT = impact )
    t0 <- ncol(augsynth$data$X)
    treat_year = augsynth$data$time[t0 + 1]

    res = calculate_MDES_table(lest, treat_year, augsynth$unit_var, augsynth$time_var)

    # Add RMSPE to donor list
    units = left_join( units, res$RMSPEs, by = names(units)[names(units) %in% names(res$RMSPEs)])

    MDES_table = mutate( res$MDES_table,
                         raw_average = T_tract$raw_average,
                         tx = ifelse( !!as.name(augsynth$time_var) >= treat_year, 1, 0 ) ) %>%
        relocate( raw_average, .after = Yhat )

    augsynth$results$permutations <- list( placebo_dist = res$rstatcalc,
                                           MDES_table = MDES_table)

    return(augsynth)
}


#' Generate permutation plots
#'
#' @param results Results from calling augsynth()
#'
#' @param inf_type Inference type (takes a value of 'permutation' or 'permutation_rstat')
#' Type of inference algorithm. Inherits inf_type from `object` or otherwise defaults to "conformal". Options are
#'         \itemize{
#'          \item{"If numeric"}{A multiple of the treated unit's RMSPE above which donor units will be dropped}
#'          \item{"If a character"}{The name or names of donor units to be dropped based on the `unit` parameter
#'           in the augsynth model}
#'         }
#' @export
permutation_plot <- function(augsynth, inf_type = 'permutation') {

    if (!inf_type %in% c('permutation', 'permutation_rstat')) {
        stop("Permutation plots are only available for `permutation` and `permutation_rstat` inference types")
    }

    if(inf_type == 'permutation') {
        measure = "ATT"
        y_lab = "Estimate (ATT)"
    } else {
        measure = 'rstat'
        y_lab = "Estimate (RMPSE-adjusted ATT)"
    }

    if (is.summary.augsynth(augsynth)) {
        placebo_dist <- augsynth$permutations$placebo_dist
    } else if (is.null(augsynth$results$permutations)) {
        augsynth <- add_placebo_distribution(augsynth)
        placebo_dist <- augsynth$results$permutations$placebo_dist
    } else {
        placebo_dist <- augsynth$results$permutations$placebo_dist
    }

    plot_df <- placebo_dist %>%
        mutate(trt_status = factor(trt, levels = c(0, 1), labels = c('Control', 'Treatment')))

    t0 <- ncol(augsynth$data$X)
    treat_year = augsynth$data$time[t0 + 1]

    out_plot <- ggplot2::ggplot(plot_df,
                                aes(x = !!as.name(augsynth$time_var), y = !!as.name(measure),
                                    color = trt_status, linetype = !!as.name(augsynth$unit_var)), size = 0.8) +
        ggplot2::geom_line() +
        ggplot2::geom_vline(lty = 2, xintercept = treat_year) +
        ggplot2::geom_hline(lty = 2, yintercept = 0) +
        ggplot2::scale_color_manual(values = c('Control' = 'gray', 'Treatment' = 'black')) +
        ggplot2::scale_linetype_manual(values = rep('solid', length(unique(plot_df %>% pull(!!as.name(augsynth$unit_var)))))) +
        ggplot2::labs(color = NULL, y = y_lab) +
        ggplot2::guides(linetype = 'none') +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = 'bottom')

    return(out_plot)
}



#' Generate formatted outputs for statistical inference using permutation inference
#'
#' @param augsynth An augsynth object.
#'
#' @noRd
permutation_inf <- function(augsynth, inf_type) {

    t0 <- dim(augsynth$data$synth_data$Z0)[1]
    tpost <- dim(augsynth$data$synth_data$Z0)[1]

    out <- list()
    out$att <- augsynth$results$permutations$MDES_table$ATT

    SEg = NA
    if (inf_type == 'permutation') {
        SEg = augsynth$results$permutations$MDES_table$SE_gap
        pval = augsynth$results$permutations$MDES_table$p_gap
    } else if (inf_type == 'permutation_rstat') {
        SEg = augsynth$results$permutations$MDES_table$SE_rstat
        pval = augsynth$results$permutations$MDES_table$p_rstat
    }
    out$lb <- out$att + (qnorm(0.025) * SEg)
    out$ub <- out$att + (qnorm(0.975) * SEg)
    out$p_val <- pval

    out$lb[c(1:t0)] <- NA
    out$ub[c(1:t0)] <- NA
    out$p_val[c(1:t0)] <- NA
    out$alpha <- 0.05

    return(out)
}


