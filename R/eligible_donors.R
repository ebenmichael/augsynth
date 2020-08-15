##############################################################################
## Code to get eligible donor units based on covariates
##############################################################################

get_donors <- function(trt, Z, time_cohort, n_leads, how = "knn", ...) {

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

#' Get donors that don't have missing outcomes where treated units have outcomes
get_nona_donors <- function(X, trt, mask, time_cohort) {

  n <- length(trt)
  # treatment times
  if(time_cohort) {
      grps <- unique(trt[is.finite(trt)])
      which_t <- lapply(grps, function(tj) (1:n)[trt == tj])
      # if doing a time cohort, convert the boolean mask
      mask <- unique(mask)
  } else {
      grps <- trt[is.finite(trt)]
      which_t <- (1:n)[is.finite(trt)]
  }


  J <- length(grps)
  # handle  X differently if it is a list
  if(typeof(X) == "list") {
      lapply(1:J,
             function(j) {
               isna_j <- is.na(mask[j, ])
               apply(X[[j]], 1, function(x) all(!is.na(x)[!isna_j]))
             }) -> donors
  } else {
    lapply(1:J,
             function(j) {
               isna_j <- is.na(mask[j, ])
               apply(X, 1, function(x) all(!is.na(x)[!isna_j]))
             }) -> donors
  }
  return(donors)
}

get_matched_donors <- function(trt, Z, donors, how, k = NULL, ...) {

  J <- sum(is.finite(trt))
  trt_idx <- which(is.finite(trt))

  if(how == "exact") {
    return(
      lapply(1:J,
          function(j) donors[[j]] & apply(Z == Z[trt_idx[j]], 1, all)
      )
    )
  } else if(how == "knn") {
    if(is.null(k)) {
      stop("Number of neighbors for knn not selected, please choose k.")
    }
    return(get_knn_donors(trt, Z, donors, k))
  } else {
    stop("Option for exact matching must be in ('exact', 'knn')")
  }
}

get_knn_donors <- function(trt, Z, donors, k) {

  # knn matching within time cohort
  grps <- unique(trt[is.finite(trt)])

  # create a list of potential neighbors and queries for each time cohort
  lapply(grps,
         function(tj) {
           # idxs for treated units treated at time tj
           tj_idxs <- trt == tj
           Z_tj <- Z[tj_idxs, , drop = F]

           # get donors for treated cohort
           donors_idxs <- which(trt[is.finite(trt)] == tj)
           # just take the first treated unit since all donors are the same
           # for the same treatment cohort
           donors_tj <- donors[[donors_idxs[1]]]
           Z_donors_tj <- Z[donors_tj, , drop = F]

           # do knn matching
           nn <- FNN::get.knnx(data = Z_donors_tj, query = Z_tj, k = k)
           # keep track of which indices these are
           lapply(1:nrow(nn$nn.index),
                  function(j) {
                    donors_j <- logical(length(donors_tj))
                    true_idx <- which(donors_tj)[nn$nn.index[j, ]]
                    donors_j[true_idx] <- TRUE
                    return(donors_j)
                  }) -> matched_donors_tj
          names(matched_donors_tj) <- donors_idxs
          return(matched_donors_tj)
         }) -> matches
  # flatten list
  matches <- unlist(matches, recursive = F)
  # reorder the donors
  matches[as.character(sort(as.integer(names(matches))))]
}