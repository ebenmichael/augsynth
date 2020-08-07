##############################################################################
## Code to get eligible donor units based on covariates
##############################################################################

get_donors <- function(trt, Z, time_cohort, n_leads, how = "exact", ...) {

  # first get eligible donors by treatment time
  donors <- get_eligible_donors(trt, time_cohort, n_leads)

  # if Z isn't NULL, futher restrict the donors by matching
  if(!is.null(Z)) {
    donors <- get_matched_donors(trt, Z, donors, how, ...)
  }

  return(donors)
}


get_eligible_donors <- function(trt, time_cohort, n_leads) {

    # get treatment times
    if(time_cohort) {
        grps <- unique(trt[is.finite(trt)])
    } else {
        grps <- trt[is.finite(trt)]
    }

    J <- length(grps)

    # only allow weights on donors treated after n_lags
    donors <- lapply(1:J, function(j) trt > n_leads + grps[j])

    return(donors)
}


get_matched_donors <- function(trt, Z, donors, how = "exact", ...) {

  J <- sum(is.finite(trt))
  trt_idx <- which(is.finite(trt))

  if(how == "exact") {
    return(
      lapply(1:J, 
          function(j) donors[[j]] & apply(Z == Z[trt_idx[j]], 1, all)
      )
    )
  } else {
    stop("Option for exact matching must be in ('exact')")
  }
}