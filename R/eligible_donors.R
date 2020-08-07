##############################################################################
## Code to get eligible donor units based on covariates
##############################################################################

get_eligible_donors <- function(Z, trt, how = "exact") {

  J <- sum(is.finite(trt))
  trt_idx <- which(is.finite(trt))
  if(how == "exact") {
    if(is.null(Z)) {
    return(lapply(1:J, function(j) rep(TRUE, length(trt))))
  } else {
    return(lapply(1:J, function(j) apply(Z == Z[trt_idx[j]], 1, all)))
  }
  }
}