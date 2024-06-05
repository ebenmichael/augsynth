
# This file contains methods to monitor and modify set of donor units
# in augmented synthetic control method


#' Return a summary data frame donor units used in the model with
#' their synthetic weights.
#'
#' If permutation inference has been conducted, table will also
#' include RMSPEs.  This can be forced with include_RMSPE flag.
#'
#' If the augsynth object does not have permutation-based inference
#' results, the function will call that form of inference, in order to
#' calculate the RMSPEs for each donor unit in turn.
#'
#' @param augsynth Augsynth object to be plotted
#' @param include_RMSPE Include RMSPEs in the table even if
#'   permutation inference has not yet been conducted.
#' @param zap_weights all weights smaller than this value will be set
#'   to zero. Set to NULL to keep all weights.
#' @export
donor_table <- function(augsynth, include_RMSPE = FALSE, zap_weights = 0.0000001 ) {

    trt_index <- which(augsynth$data$trt == 1)
    unit_var <- augsynth$unit_var

    tbl = tibble(
        unit = rownames( augsynth$weights ),
        weight = as.numeric(augsynth$weights) )
    names(tbl)[[1]] = unit_var

    if ( include_RMSPE || (!is.null(augsynth$results) && augsynth$results$inf_type %in% c("permutation", "permutation_rstat")) ) {
        if (is.null(augsynth$results) || (!(augsynth$results$inf_type %in% c("permutation", "permutaton_rstat")) ) ) {
            augsynth <- add_inference(augsynth, inf_type = 'permutation')
        }
        RMSPEs <- augsynth$results$permutations$placebo_dist %>%
            select(!!unit_var, RMSPE) %>%
            distinct()
        tbl <- left_join(tbl, RMSPEs, by = unit_var)
    }

    if ( !is.null(zap_weights) ) {
        tbl <- tbl %>% mutate(weight = ifelse(abs(weight) < zap_weights, 0, weight))
    }

    return(tbl)
}


#' Return a new augsynth object with specified donor units removed
#'
#' @param augsynth Augsynth object to be plotted
#'
#' @param drop Drop donor units, based on pre-treatment RMSPE or unit
#'   name(s). Default of 20 means drop units with an RMSPE 20x higher
#'   than the treated unit.  Can also be a character vector of unit
#'   IDs to drop.
#'
#' @export
update_augsynth <- function(augsynth, drop = 20){

    if (is.null(augsynth$results)){
        inf_type = 'none'
    } else {
        inf_type <- augsynth$results$inf_type
    }

    # run placebo tests if necessary
    if (!inf_type %in% c('permutation', 'permutation_rstat')) {
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
    new_data <- as_tibble(augsynth$raw_data, .name_repair = 'unique') %>%
        filter(!!as.name(unit_var) %in% keep_units)

    new_augsynth <- augsynth(form = form,
                             unit = !!as.name(augsynth$unit_var),
                             time = !!as.name(augsynth$time_var),
                             data = new_data,
                             progfunc = augsynth$progfunc,
                             scm = augsynth$scm,
                             fixedeff = augsynth$fixedeff,
                             cov_agg = augsynth$cov_agg,
                             inf_type = inf_type
    )

    return(new_augsynth)
}
