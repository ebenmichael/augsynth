#' @export
#' 
augsynth <- function(form, unit, time, data, t_int=NULL, ...) {

# augsynth <- function(form, unit, time, data,
#                      scm = T, progfunc=c("Ridge", "None", "EN", "RF", "GSYN", "MCP", "CITS", "CausalImpact", "seq2seq"), t_int = NULL,
#                      fixedeff = FALSE, cov_agg=NULL, ...,
#                      alpha=NULL, lambda=0,
#                      relative=T, n_leads=NULL, n_lags=NULL,
#                      force="two-way",
#                      n_factors=NULL) {
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
