################################################################################
## Scripts to format panel data into matrices
################################################################################

#' Format "long" panel data into "wide" program evaluation matrices
#' @param outcome Name of outcome column
#' @param trt Name of treatment column
#' @param unit Name of unit column
#' @param time Name of time column
#' @param t_int Time of intervention
#' @param data Panel data as dataframe
#' @noRd
#' @return \itemize{
#'          \item{"X"}{Matrix of pre-treatment outcomes}
#'          \item{"trt"}{Vector of treatment assignments}
#'          \item{"y"}{Matrix of post-treatment outcomes}
#'         }
format_data <- function(outcome, trt, unit, time, t_int, data) {

    ## pre treatment outcomes
    X <- data %>%
        filter(!!time < t_int) %>%
        select(!!unit, !!time, !!outcome) %>%
        spread(!!time, !!outcome) %>%
        select(-!!unit) %>%
        as.matrix()
    


    ## post treatment outcomes
    y <- data %>%
        filter(!!time >= t_int) %>%
        select(!!unit, !!time, !!outcome) %>%
        spread(!!time, !!outcome) %>%
        select(-!!unit) %>%
        as.matrix()


    ## treatment status
    trt <- data %>%
        select(!!unit, !!trt) %>%
        group_by(!!unit) %>%
        summarise(trt = max(!!trt)) %>%
        ungroup() %>%
        pull(trt)

    return(list(X=X, trt=trt, y=y))
}


#' Format "long" panel data into "wide" program evaluation matrices
#' @param outcomes Vectors of names of outcome columns
#' @param trt Name of treatment column
#' @param unit Name of unit column
#' @param time Name of time column
#' @param t_int Time of intervention
#' @param data Panel data as dataframe
#' @noRd
#' @return \itemize{
#'          \item{"X"}{List of matrices of pre-treatment outcomes}
#'          \item{"trt"}{Vector of treatment assignments}
#'          \item{"y"}{List of matrices of post-treatment outcomes}
#'         }
format_data_multi <- function(outcomes, trt, unit, time, t_int, data) {


    lapply(outcomes, 
        function(outcome) format_data(outcome, trt, unit, 
                                     time, t_int, data)
          ) -> formats

    # X <- simplify2array(lapply(formats, function(x) x$X))
    # y <- simplify2array(lapply(formats, function(x) x$y))
    X <- lapply(formats, function(x) t(na.omit(t(x$X))))
    y <- lapply(formats, function(x) t(na.omit(t(x$y))))
    trt <- formats[[1]]$trt
    return(list(X = X, trt = trt, y = y))
}




#' Format "long" panel data into "wide" program evaluation matrices with staggered adoption
#' @param outcome Name of outcome column
#' @param trt Name of treatment column
#' @param unit Name of unit column
#' @param time Name of time column
#' @param data Panel data as dataframe
#' @noRd
#' @return \itemize{
#'          \item{"X"}{Matrix of pre-treatment outcomes}
#'          \item{"trt"}{Vector of treatment assignments}
#'          \item{"y"}{Matrix of post-treatment outcomes}
#'         }
format_data_stag <- function(outcome, trt, unit, time, data) {

    # arrange data by time first
    data <- data %>% arrange(!!time)
      
    ## get first treatment times
    trt_time <- data %>%
        group_by(!!unit) %>%
        summarise(trt_time=(!!time)[(!!trt) == 1][1]) %>%
        mutate(trt_time=replace_na(as.numeric(trt_time), Inf))
    

    t_int <- trt_time %>% filter(is.finite(trt_time)) %>%
        summarise(t_int=max(trt_time)) %>% pull(t_int)

    ## ## boolean mask of available data for treatment groups
    ## mask <- data %>% inner_join(trt_time %>%
    ##                             filter(is.finite(trt_time))) %>%
    ##     filter(!!time < t_int) %>%
    ##     mutate(trt=1-!!trt) %>%
    ##     select(!!unit, !!time, trt_time, trt) %>%
    ##     spread(!!time, trt) %>% 
    ##     group_by(trt_time) %>% 
    ##     summarise_all(list(max)) %>%
    ##     arrange(trt_time) %>% 
    ##     select(-trt_time, -!!unit) %>%
    ##     as.matrix()

    ## boolean mask of available data for treatment groups
    mask <- data %>% inner_join(trt_time %>%
                                filter(is.finite(trt_time)),
                                by = rlang::as_name(unit)) %>%
        filter(!!time < t_int) %>%
        mutate(trt=1-!!trt) %>%
        select(!!unit, !!time, trt_time, trt) %>%
        spread(!!time, trt) %>% 
        ## arrange(!!unit) %>% 
        select(-trt_time, -!!unit) %>%
        as.matrix()
    
    # outcomes as a matrix
    Xy <- data %>%
        select(!!unit, !!time, !!outcome) %>%
        spread(!!time, !!outcome) %>%
        select(-!!unit) %>%
        as.matrix()

    pre_times <- data %>% filter(!!time < t_int) %>%
        distinct(!!time) %>% pull(!!time)
    post_times <- data %>% filter(!!time >= t_int) %>%
        distinct(!!time) %>% pull(!!time)
    X <- Xy[, as.character(pre_times), drop = F]
    y <- Xy[, as.character(post_times), drop = F]

    if(nrow(X) != nrow(y)) {
      stop("There are not the same number of units after the last unit is treated as before the last unit is treated")
    }

    t_vec <- data %>% pull(!!time) %>% unique() %>% sort()
    trt <- sapply(trt_time$trt_time, function(x) which(t_vec == x)-1) %>%
        as.numeric() %>%
        replace_na(Inf)
   

    units <- data %>%
        filter(!!time < t_int) %>%
        select(!!unit, !!time, !!outcome) %>%
        spread(!!time, !!outcome) %>%
        pull(!!unit)

    
    return(list(X=X,
                trt=trt,
                y=y,
                mask=mask,
                time = t_vec,
                units=units))
}


#' Format program eval matrices into synth form
#'
#' @param X Matrix of pre-treatment outcomes
#' @param trt Vector of treatment assignments
#' @param y Matrix of post-treatment outcomes
#' @noRd
#' @return List with data formatted as Synth::dataprep
format_synth <- function(X, trt, y) {


    synth_data <- list()

    ## pre-treatment values as covariates
    synth_data$Z0 <- t(X[trt==0,,drop=F])

    ## average treated units together
    synth_data$Z1 <- as.matrix((colMeans(X[trt==1,,drop=F])), ncol=1)

    ## combine everything together also
    synth_data$Y0plot <- t(cbind(X[trt==0,,drop=F], y[trt==0,,drop=F]))
    synth_data$Y1plot <- as.matrix(colMeans(
        cbind(X[trt==1,,drop=F], y[trt==1,,drop=F])), ncol=1)


    ## predictors are pre-period outcomes
    synth_data$X0 <- synth_data$Z0
    synth_data$X1 <- synth_data$Z1

    return(synth_data)
    
}

#' Remove unit means 
#' @param wide_data X, y, trt
#' @param synth_data List with data formatted as Synth::dataprep
#' @noRd
demean_data <- function(wide_data, synth_data) {

    # pre treatment means
    means <- rowMeans(wide_data$X)

    new_wide_data <- list()
    new_X <- wide_data$X - means
    trt <- wide_data$trt

    new_wide_data$X <- new_X
    new_wide_data$y <- wide_data$y - means
    new_wide_data$trt <- trt

    new_synth_data <- list()
    new_synth_data$X0 <- t(new_X[trt == 0,, drop = FALSE])
    new_synth_data$Z0 <- new_synth_data$X0
    new_synth_data$X1 <- as.matrix((colMeans(new_X[trt==1,,drop = F])), 
                                   ncol = 1)
    new_synth_data$Z1 <- new_synth_data$X1


    # estimate post-treatment as pre-treatment means
    mhat <- replicate(ncol(wide_data$X) + ncol(wide_data$y), means)

    return(list(wide = new_wide_data,
                synth_data = new_synth_data,
                mhat = mhat))
}

#' Helper function to extract covariate matrix from data
#' @param form Formula as outcome ~ treatment | covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param t_int Time of intervention
#' @param data Panel data as dataframe
#' @param cov_agg Covariate aggregation function
#' @noRd
extract_covariates <- function(form, unit, time, t_int, data, cov_agg) {

    ## if no aggregation functions, use the mean (omitting NAs)
    if(is.null(cov_agg)) {
        cov_agg <- c(function(x) mean(x, na.rm=T))
    }

    cov_form <- update(formula(delete.response(terms(form, rhs=2, data=data))),
                        ~. - 1) ## ensure that there is no intercept

    ## pull out relevant covariates and aggregate
    pre_data <- data %>% 
        filter(!! (time) < t_int)

    model.matrix(cov_form,
                    model.frame(cov_form, pre_data,
                                na.action=NULL) ) %>%
        data.frame() %>%
        mutate(unit=pull(pre_data, !!unit)) %>%
        group_by(unit) %>%
        summarise_all(cov_agg) -> Z

    # recombine with any missing units and convert to matrix
    data %>% distinct(!!unit) %>%
      rename(unit = !!unit) %>%
      left_join(Z, by = "unit") %>%
      arrange(unit) %>%
      select(-unit) %>%
      as.matrix() -> Z
    
    if(nrow(distinct(data, !!unit))  != nrow(Z)) {
      stop("Some units missing all covariate data")
    }
    return(Z)
}

#' Check that we can actually run multisynth on the data
#' @param wide Output of format_data_stag
#' @param fixedeff Whether to include a unit fixed effect
#' @param n_leads How long past treatment effects should be estimated for, default is number of post treatment periods for last treated unit
#' @param n_lags Number of pre-treatment periods to balance, default is to balance all periods
check_data_stag <- function(wide, fixedeff, n_leads, n_lags) {

  # If there are less than 5 pre-treatment outcomes, give a warning
  less_5 <- wide$units[wide$trt < 5]
  warn_msg <- ""
  if(length(less_5) != 0) {
    warn_msg <- paste0(
      warn_msg,
      "The following units have less than 5 pre-treatment outcomes: (",
      paste(less_5, collapse = ","),
      "). Be cautious!"
    )
  }

  # check if there are any always treated units
  always_trt <- wide$units[wide$trt == 0]

  # If including a fixed effect, check that there is more than one pretreatment
  # outcome for each unit
  n1 <- wide$units[wide$trt == 1]

  err_msg <- ""
  if(length(always_trt) != 0) {
    err_msg <- paste0(
      err_msg,
      "The following units are always treated and should be removed: (",
      paste(always_trt, collapse = ","),
      ")\n")
  }

  if(length(n1) != 0 & fixedeff) {
    if(nchar(err_msg) > 0) {
      err_msg <- paste0(err_msg, "  Also: ")
    }
    err_msg <- paste0(
      err_msg,
      "You are including a fixed effect with `fixedeff = T`, but the ",
      "following units only have one pre-treatment outcome: (",
      paste(n1, collapse = ","),
      "). Either remove these units or set `fixedeff = F`."
    )
  }

  if(nchar(warn_msg) > 0) {
    warning(warn_msg)
  }
  if(nchar(err_msg) > 0) {
    stop(err_msg)
  }

}