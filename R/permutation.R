

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
    summarise( p_perm = mean( abs( rstat[trt == 1] ) <= abs( rstat ) ),
               .groups = "drop" )

  # Calculate the estimated SEs of the rstat distribution of future
  # based on the trimmed distribution of placebo estimates (for each
  # time point)
  new_T_tract = rstatcalc %>%
    group_by( !!as.name(time_col) ) %>%
    summarise( ests = estimate_robust_SE( rstat, beta=0.90,
                                          round_beta = TRUE),
               .groups = "drop" ) %>%
    tidyr::unpack( ests )



  # attach R statistic for the treated tract.
  rstattract <- rstatcalc %>%
    filter(trt==1) %>%
    dplyr::select( -!!as.name(tx_col), -trt )

  new_T_tract2 <- full_join(new_T_tract, rstattract, by=time_col)


  # Build the main table of results from the above components.

  main_results = new_T_tract2 %>%
    rename (R_low = q_low,
            R_high = q_high,
            R_range = range,
            SE_rstat = SE_imp ) %>%
    mutate( t       = rstat / SE_rstat,
            p       = 2 * pnorm( -abs(t) ),
            MDE_10   = RMSPE * 2.49 * SE_rstat )


  # Add in the estimated SEs using the ATT, not r-statistics. (This
  # is the un-preferred method, but is basically replicating the
  # above.)
  SE_for_ATT = rstatcalc %>%
    group_by( !!as.name(time_col) ) %>%
    summarise( ests = estimate_robust_SE( ATT, beta=0.90,
                                          round_beta = TRUE),
               .groups = "drop" ) %>%
    tidyr::unpack( ests )
  SE_for_ATT <- full_join(SE_for_ATT, rstattract, by=time_col)

  SE_for_ATT = SE_for_ATT %>%
    rename( SE_gap   = SE_imp ) %>%
    mutate(          t_gap    = ATT / SE_gap,
                     p_gap    = 2 * pnorm( -abs(t_gap) ), # 2-sideds
                     MDE_gap  = 2.49 * SE_gap ) %>%
    dplyr::select( !!as.name(time_col), SE_gap, t_gap, p_gap, MDE_gap )

  main_results = left_join( main_results, SE_for_ATT, by=time_col )


  # Add in permutation p-values and add indicator of which years are
  # treated and which not.

  main_results = left_join( main_results, pvalues, by=time_col ) %>%
    mutate( tx = ifelse( !!as.name(time_col) >= treat_year, 1, 0 ) )


  # Clean up column ordering and make look nice
  main_results = main_results %>%
    relocate( !!as.name(time_col), tx ) %>%
    relocate( beta, k, R_low, R_high, R_range, z, .after = MDE_10 ) %>%
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
#' @param beta Proportion of obs to drop from each tail. E.g., 95% to
#'   drop 5% most extreme.  Only one of beta or k can be non-null.
#' @param round_beta TRUE means adjust beta to the nearest integer
#'   number of units to drop (the Howard method).  FALSE means
#'   interpolate quantiles using R's quantile function.
#'
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



add_placebo_distribution <- function(augsynth){

  # Run permutations
  ests <- get_placebo_gaps(augsynth, att = FALSE)
  time_cols = 2:ncol(ests)

  ests$trt = augsynth$data$trt

  lest = ests %>%
    pivot_longer( cols = all_of( time_cols ),
                  names_to = augsynth$time_var, values_to = "CoLvl" ) %>%
    mutate(!!as.name(augsynth$time_var) := as.numeric(!!as.name(augsynth$time_var)))

  ######## Construct a frame with outcome data in long format - ND added
  df <- bind_cols(augsynth$data$X, augsynth$data$y) %>%
    mutate(!!as.name(augsynth$unit_var) := lest$ID %>% unique()) %>%
    pivot_longer(!augsynth$unit_var,
                 names_to = augsynth$time_var,
                 values_to = 'Yobs',
    ) %>%
    mutate(!!as.name(augsynth$time_var) := as.numeric(!!as.name(augsynth$time_var)),
           ever_Tx = ifelse(!!as.name(augsynth$unit_var) == augsynth$trt_unit, 1, 0))

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
                    by = c( augsynth$unit_var, augsynth$time_var ) )
  stopifnot( nrow( lest ) == nn )


  # Impact is difference between observed and imputed control-side outcome
  lest$impact = lest$Yobs - lest$CoLvl

  #### Make the actual treatment result information #####

  # Get our observed series (the treated unity)
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
  units = left_join( units, res$RMSPEs, by = augsynth$unit_col)

  MDES_table = mutate( res$MDES_table,
                       raw_average = T_tract$raw_average,
                       tx = ifelse( !!as.name(augsynth$time_var) >= treat_year, 1, 0 ) ) %>%
    relocate( raw_average, .after = CoLvl )

  augsynth$results$permutations <- list(placebo_dist = res$rstatcalc,
                                   MDES_table = MDES_table)

  return(augsynth)
}


#' Generate permutation plots
#'
#' @param results Results from calling augsynth()
#' @param inf_type Inference type (takes a value of 'permutation' or 'permutation_rstat')
#' Type of inference algorithm. Inherits inf_type from `object` or otherwise defaults to "conformal". Options are
#'         \itemize{
#'          \item{"If numeric"}{A multiple of the treated unit's RMSPE above which donor units will be dropped}
#'          \item{"If a character"}{The name or names of donor units to be dropped based on the `unit` parameter in the augsynth model}
#'         }
#'
#' with an RMSPE of more than `drop` times the treated units' RMSPE from the donor plot. Defaults to 20x.
#' @export
permutation_plot <- function(augsynth, inf_type = 'permutation') {

  stopifnot(inf_type %in% c('permutation', 'permutation_rstat'))

  if(inf_type == 'permutation') {
    measure = "ATT"
  } else {
    measure = 'rstat'
  }

  if (is.null(augsynth$results$permutations)) {
    augsynth <- add_placebo_distribution(augsynth)
  }

  plot_df <- augsynth$results$permutations$placebo_dist %>%
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
    ggplot2::labs(color = NULL, y = 'Estimate') +
    ggplot2::guides(linetype = 'none') +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = 'bottom')

  return(out_plot)
}


#' Generate augsynth plot with 95% CI based on inference preference
#'
#' @param results Results from calling xxx
#' @param inf_type Inference type (takes a value of 'permutation' or 'permutation_rstat')
#'
ci95_rstat_plot <- function(augsynth, inf_type) {

  stopifnot(inf_type %in% c('permutation', 'permutation_rstat'))

  if(inf_type == 'permutation') {
    measure = "ATT"
  } else {
    measure = 'rstat'
  }

  if (is.null(augsynth$results$permutations)) {
    augsynth <- add_placebo_distribution(augsynth)
  }

  plot_df <- augsynth$results$permutations$MDES_table

  if (measure == 'ATT') {
    plot_df <- plot_df %>%
      mutate(lower = ATT + (qnorm(0.025) * SE_gap),
             upper = ATT + (qnorm(0.975) * SE_gap),
      )
  } else if (measure == 'rstat') {
    plot_df <- plot_df %>%
      mutate(lower = ATT + (qnorm(0.025) * SE_rstat),
             upper = ATT + (qnorm(0.975) * SE_rstat),
      )
  }

  t0 <- ncol(augsynth$data$X)
  treat_year = augsynth$data$time[t0 + 1]

  outplot <- ggplot2::ggplot(data = plot_df, aes(x = !!as.name(augsynth$time_var))) +
    ggplot2::geom_ribbon(data = plot_df %>% filter(!!as.name(augsynth$time_var) >= treat_year),
                         aes(ymin = lower, ymax = upper), alpha = 0.2) +
    ggplot2::geom_line(aes(y = ATT)) +
    ggplot2::geom_vline(lty = 2, xintercept = treat_year) +
    ggplot2::geom_hline(lty = 2, yintercept = 0) +
    ggplot2::scale_color_manual(values = c('Control' = 'gray', 'Treatment' = 'black')) +
    ggplot2::labs(color = NULL, y = 'Estimate') +
    ggplot2::theme_bw()

  return(outplot)
}

#' Generate formatted outputs for statistical inference using permutation inference
#'
#' @param augsynth An augsynth object.
#'
permutation_inf <- function(augsynth, inf_type) {

  t0 <- dim(augsynth$data$synth_data$Z0)[1]
  tpost <- dim(augsynth$data$synth_data$Z0)[1]

  out <- list()
  out$att <- augsynth$results$permutations$MDES_table$ATT

  if (inf_type == 'permutation') {
    out$lb <- augsynth$results$permutations$MDES_table$ATT + (qnorm(0.025) * augsynth$results$permutations$MDES_table$SE_gap)
    out$ub <- augsynth$results$permutations$MDES_table$ATT + (qnorm(0.975) * augsynth$results$permutations$MDES_table$SE_gap)
    out$p_val <- c(rep(0, t0), augsynth$results$permutations$MDES_table$p_gap[c(t0 + 1: tpost)])
  } else if (inf_type == 'permutation_rstat') {
    out$lb <- augsynth$results$permutations$MDES_table$ATT + (qnorm(0.025) * augsynth$results$permutations$MDES_table$SE_rstat)
    out$ub <- augsynth$results$permutations$MDES_table$ATT + (qnorm(0.975) * augsynth$results$permutations$MDES_table$SE_rstat)
    out$p_val <- c(rep(0, t0), augsynth$results$permutations$MDES_table$p[c(t0 + 1: tpost)])
  }

  out$lb[c(1:t0)] <- NA
  out$ub[c(1:t0)] <- NA
  out$p_val[c(1:t0)] <- NA
  out$alpha <- 0.05

  return(out)
}


#' Plot function returning the original level of the outcome variable for the treated unit and its synthetic counterfactual
#' @param augsynth Augsynth object to be plotted
#' @param  measure Whether to plot the synthetic counterfactual or the raw average of donor units
#' @export
augsynth_gap_plot <- function(augsynth, measure = c("synth", "average")) {

  trt_index <- which(augsynth$data$trt == 1)
  df <- bind_cols(augsynth$data$X, augsynth$data$y)
  synth_unit <- t(df[-trt_index, ]) %*% augsynth$weights
  average_unit <- df[-trt_index, ] %>% colMeans()
  treated_unit <- t(df[trt_index, ])

  max_y <- max(c(synth_unit, average_unit, treated_unit))
  min_y <- min(c(synth_unit, average_unit, treated_unit))

  p <- ggplot2::ggplot() +
    ggplot2::geom_line(aes(x = augsynth$data$time, y = treated_unit, linetype = augsynth$trt_unit))

  if ('synth' %in% measure) {
    p <- p +
      ggplot2::geom_line(aes(x = augsynth$data$time, y = synth_unit, linetype = 'Synthetic counterfactual'))
  }

  if ('average' %in% measure) {
    p <- p +
      ggplot2::geom_line(aes(x = augsynth$data$time, y = average_unit, linetype = 'Donor raw average'))
  }

  p <- p +
    ggplot2::scale_y_continuous(limits = c(0, 140)) +
    ggplot2::labs(linetype = NULL,
                  x = augsynth$time_var,
                  y = 'Outcome') +
    ggplot2::ylim(min_y, max_y) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = c(0.75, 0.88))

  return(p)
}


#' Return a summary data frame donor units used in the model, along with their RMSPEs and synthetic weights
#' @param augsynth Augsynth object to be plotted
#' @export
donor_table <- function(augsynth) {

  if (!augsynth$results$inf_type %in% c("permutation", "permutaton_rstat")) {
    augsynth <- add_inference(augsynth, inf_type = 'permutation')
  }

  trt_index <- which(augsynth$data$trt == 1)
  unit_var <- augsynth$unit_var
  RMSPEs <- augsynth$results$permutations$placebo_dist %>%
    select(!!unit_var, RMSPE) %>%
    distinct()
  donor_df <- RMSPEs %>% filter(!!as.name(unit_var) != augsynth$trt_unit) %>%
    mutate(synth_weight = augsynth$weights)

  return(donor_df)

}


#' Return a new augsynth object with specified donor units removed
#' @param augsynth Augsynth object to be plotted
#' @param drop Drop donor units, based on pre-treatment RMSPE or unit name(s)
#' @export
update_augsynth <- function(augsynth, drop = 20){

  inf_type <- augsynth$results$inf_type

  # run placebo tests if necessary
  if (!augsynth$results$inf_type %in% c('permutation', 'permutation_rstat')) {
    augsynth <- add_inference(augsynth, inf_type = 'permutation')
  }

  unit_var <- augsynth$unit_var
  # pre-treatment RMSPE among donors
  donor_RMSPE <- augsynth$results$permutations$placebo_dist %>%
    filter(!!as.name(augsynth$time_var) < augsynth$t_int) %>%
    group_by(!!as.name(augsynth$unit_var)) %>%
    summarise(RMSPE = sqrt(mean(ATT ^ 2)), .groups = "drop")
  # pre-treatment RMSPE for treated unit
  trt_RMSPE <- add_inference(augsynth, inf_type = 'permutation')$results$permutations$placebo_dist %>%
    filter(!!as.name(augsynth$time_var) < augsynth$t_int) %>%
    filter(!!as.name(unit_var) == augsynth$trt_unit) %>%
    pull(RMSPE) %>% unique()

  if (is.numeric(drop)) {
    keep_units <- donor_RMSPE %>% filter(RMSPE / trt_RMSPE <= drop) %>% pull(!!unit_var)
  } else if (is.character(drop)) {
    keep_units <- donor_RMSPE %>% filter((!!as.name(unit_var) %in% drop) == FALSE) %>% pull(!!unit_var) %>% unique()
  }
  keep_units <- c(keep_units, augsynth$trt_unit)

  form <- as.formula(paste(as.character(augsynth$form)[2], as.character(augsynth$form)[1], as.character(augsynth$form)[3]))
  new_data <- as_tibble(augsynth$raw_data, .name_repair = 'unique') %>% filter(!!as.name(unit_var) %in% keep_units)
  new_augsynth <- augsynth(form = form,
                           unit = !!as.name(augsynth$unit_var),
                           time = !!as.name(augsynth$time_var),
                           data = new_data,
                           inf_type = inf_type
  )

  return(new_augsynth)
}
