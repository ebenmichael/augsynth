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
#'   \item Single period augsynth 
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
  wide <- format_data_stag(outcome, trt, unit_quosure, time_quosure, data=data)

  trt_year = -1
  multi = F
  for (u in wide$trt) {
    if (is.finite(u) && u != trt_year) {
      if (trt_year == -1) {
        trt_year = u
      } else {
        multi = T
        break
      }
    }
  }

  if (multi) {
    return(multisynth(form, !!enquo(unit), !!enquo(time), data, ...)) 
  } else {
    if (is.null(t_int)) {
      t_int = data %>% pull(!!time_quosure) %>% unique() %>% sort() %>% `[`(trt_year+1)
    }
    return(single_augsynth(form, !!enquo(unit), !!enquo(time), t_int, data=data, ...))
  }
}
