################################################################################
## Main function for the augmented synthetic controls Method
################################################################################


#' Fit Augmented SCM
#' @param form outcome ~ treatment | auxillary covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param data Panel data as dataframe
#' @param t_int Time of intervention (used for single-period treatment only)
#' @param ... Optional arguments
#' \itemize{
#'   \item Single period augsynth with/without multiple outcomes
#'     \itemize{
#'       \item{"progfunc"}{What function to use to impute control outcomes: Ridge=Ridge regression (allows for standard errors), None=No outcome model, EN=Elastic Net, RF=Random Forest, GSYN=gSynth, MCP=MCPanel, CITS=CITS, CausalImpact=Bayesian structural time series with CausalImpact, seq2seq=Sequence to sequence learning with feedforward nets}
#'       \item{"scm"}{Whether the SCM weighting function is used}
#'       \item{"fixedeff"}{Whether to include a unit fixed effect, default F }
#'       \item{"cov_agg"}{Covariate aggregation functions, if NULL then use mean with NAs omitted}
#'     }
#'   \item Multi period (staggered) augsynth
#'    \itemize{
#'          \item{"relative"}{Whether to compute balance by relative time}
#'          \item{"n_leads"}{How long past treatment effects should be estimated for}
#'          \item{"n_lags"}{Number of pre-treatment periods to balance, default is to balance all periods}
#'          \item{"alpha"}{Fraction of balance for individual balance}
#'          \item{"lambda"}{Regularization hyperparameter, default = 0}
#'          \item{"force"}{Include "none", "unit", "time", "two-way" fixed effects. Default: "two-way"}
#'          \item{"n_factors"}{Number of factors for interactive fixed effects, default does CV}
#'         }
#' }
#' 
#' @return augsynth object that contains:
#'         \itemize{
#'          \item{"weights"}{weights}
#'          \item{"data"}{Panel data as matrices}
#'         }
#' @export
#' 
augsynth <- function(form, unit, time, data, t_int=NULL, ...) {

  call_name <- match.call()

  form <- Formula::Formula(form)
  unit_quosure <- enquo(unit)
  time_quosure <- enquo(time)
  

  ## format data
  outcome <- terms(formula(form, rhs=1))[[2]]
  trt <- terms(formula(form, rhs=1))[[3]]

  # check for multiple outcomes
  multi_outcome <- length(outcome) != 1

  ## get first treatment times
  trt_time <- data %>%
      group_by(!!unit_quosure) %>%
      filter(!all(!!trt == 0)) %>%
      summarise(trt_time = min((!!time_quosure)[(!!trt) == 1])) %>%
      mutate(trt_time = replace_na(as.numeric(trt_time), Inf))

  num_trt_years <- sum(is.finite(unique(trt_time$trt_time)))

  if(multi_outcome & num_trt_years > 1) {
    stop("augsynth is not currently implemented for more than one outcome and more than one treated unit")
  } else if(num_trt_years > 1) {
    message("More than one treatment time found. Running multisynth.")
    if("progfunc" %in% names(list(...))) {
      warning("`progfunc` is not an argument for multisynth, so it is ignored")
    }
    return(multisynth(form, !!enquo(unit), !!enquo(time), data, ...)) 
  } else {
    if (is.null(t_int)) {
      t_int <- trt_time %>% filter(is.finite(trt_time)) %>%
        summarise(t_int = max(trt_time)) %>% pull(t_int)
    }
    if(!multi_outcome) {
      message("One outcome and one treatment time found. Running single_augsynth.")
      return(single_augsynth(form, !!enquo(unit), !!enquo(time), t_int,
                             data = data, ...))
    } else {
      message("Multiple outcomes and one treatment time found. Running augsynth_multiout.")
      return(augsynth_multiout(form, !!enquo(unit), !!enquo(time), t_int,
                               data = data, ...))
    }
  }
}
