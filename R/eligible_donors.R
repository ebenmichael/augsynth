##############################################################################
## Code to get eligible donor units based on covariates
##############################################################################

get_donors <- function(X, y, trt, Z, time_cohort, n_lags,
                       n_leads, how = "knn", 
                       exact_covariates = NULL, ...) {

  # first get eligible donors by treatment time
  donors <- get_eligible_donors(trt, time_cohort, n_leads)

  # get donors with no NA values
  nona_donors <- get_nona_donors(X, y, trt, n_lags, n_leads, time_cohort)

  donors <- lapply(1:length(donors),
                     function(j) donors[[j]] & nona_donors[[j]])

  # if Z isn't NULL, futher restrict the donors by matching
  if(!is.null(Z)) {
    if(ncol(Z) != 0) {
      donors <- get_matched_donors(trt, Z, donors, how, exact_covariates, ...)
    }
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

    # only allow weights on donors treated after n_leads
    donors <- lapply(1:J, function(j) trt > n_leads + grps[j])

    return(donors)
}

#' Get donors that don't have missing outcomes where treated units have outcomes
get_nona_donors <- function(X, y, trt, n_lags, n_leads, time_cohort) {

  n <- length(trt)
  # find na treatment times
  fulldat <- cbind(X, y)
  is_na <- is.na(fulldat[is.finite(trt), , drop = F])
  # aggregate by time cohort
  if(time_cohort) {
    grps <- unique(trt[is.finite(trt)])
    # if doing a time cohort, convert the boolean mask
    finite_trt <- trt[is.finite(trt)]
    is_na <- t(sapply(grps,
                    function(tj) apply(is_na[finite_trt == tj, , drop = F],
                                       2, all)))
  } else {
      grps <- trt[is.finite(trt)]
  }
  not_na <- !is.na(fulldat)
  J <- length(grps)
  lapply(1:J,
             function(j) {
               idxs <- max(grps[j] - n_lags + 1, 1):min(grps[j] + n_leads,
                                                        ncol(fulldat))
               isna_j <- is_na[j, idxs]
               apply(not_na[, idxs, drop = F][, !isna_j, drop = F], 1, all)
        }) -> donors

  return(donors)
}

get_matched_donors <- function(trt, Z, donors, how, exact_covariates = NULL, k = NULL, ...) {

  J <- sum(is.finite(trt))
  trt_idx <- which(is.finite(trt))
  if(is.null(exact_covariates)) {
    if(how == "exact") {
      return(
        lapply(1:J,
            function(j) donors[[j]] & apply(t(Z) == Z[trt_idx[j], ], 2, all)
        )
      )
    } else if(how == "knn") {
        return(get_knn_donors(trt, Z, donors, k))
    } else {
      stop("Option for exact matching must be in ('exact', 'knn')")
    }
  } else {
        if(how == "exact") {
      return(
        lapply(1:J,
            function(j) donors[[j]] & apply(t(Z) == Z[trt_idx[j], 
                                                   exact_covariates], 2, all)
        )
      )
    } else if(how == "knn") {
        donors <- lapply(1:J,
            function(j) { donors[[j]] &
              apply(t(Z[, exact_covariates, drop = F]) == 
                Z[trt_idx[j],exact_covariates], 2, all)
            }
              )
        approx_covs <- which(!colnames(Z) %in% exact_covariates)
        if(length(approx_covs != 0)) {
          return(get_knn_donors(trt, Z[, approx_covs, drop = F], donors, k))
        } else {
          return(donors)
        }
        
    } else {
      stop("Option for exact matching must be in ('exact', 'knn')")
    }
  }

}

get_knn_donors <- function(trt, Z, donors, k) {

  if(is.null(k)) {
    stop("Number of neighbors for knn not selected, please choose k.")
  }
  # knn matching within time cohort
  trt_idxs <- which(is.finite(trt))
  lapply(1:length(trt_idxs), 
        function(j) {
          idx <- trt_idxs[j]
          # idxs for treated units treated at time tj
          Z_tj <- Z[idx, , drop = F]

          # get donors for treated cohort
          donors_tj <- donors[[j]]
          Z_donors_tj <- Z[donors_tj, , drop = F]
          # check that k is less than the number of donors
          # if not, warn and set k to be the number of donors - 1
          if(k >= nrow(Z_donors_tj)) {
            warning(paste("Number of potential donor units is less than",
                          "the number of required matches,",
                          "returning all donors as matches"))
            return(donors_tj)
          }
          # do knn matching
          nn <- FNN::get.knnx(data = Z_donors_tj, query = Z_tj, k = k)
          # keep track of which indices these are
          donors_j <- logical(length(donors_tj))
          true_idx <- which(donors_tj)[nn$nn.index[1, ]]
          donors_j[true_idx] <- TRUE
          return(donors_j)
         }) -> matches
  names(matches) <- trt_idxs
  return(matches)
}