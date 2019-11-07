#' @export
augsynth3 <- function(form, unit, time, data,
                        relative=T, n_leads=NULL, n_lags=NULL,
                        scm = T, progfunc=c("Ridge", "None", "EN", "RF", "GSYN", "MCP",
                                            "CITS", "CausalImpact", "seq2seq"), t_int = NULL,
                        alpha=NULL, lambda=0,
                        force="two-way",
                        n_factors=NULL,
                        opts_weights=NULL) {
  call_name <- match.call()
  
  form <- Formula::Formula(form)
  unit <- enquo(unit)
  time <- enquo(time)
  
  print("e")
  ## format data
  outcome <- terms(formula(form, rhs=1))[[2]]
  trt <- terms(formula(form, rhs=1))[[3]]
  wide <- format_data_stag(outcome, trt, unit, time, data)
  
  trt_year = -1
  multi = F
  for (unit in wide$trt) {
    if (is.finite(unit)) {
      if (trt_year == -1) { # if 2 units receive treatment in the same year, does it go to multi synth or single synth?
        trt_year = unit
      } else {
        multi = T
        break
      }
    }
  }
  if (multi) {
    return(multisynth(form, unit, time, data,
               force="time", n_factors, n_leads))
    
  } else {
    return(augsynth(form, unit, time, t_int, data, progfunc=progfunc, scm=scm))
  }
}
