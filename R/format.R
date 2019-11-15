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
#'
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
        pull(trt)

    return(list(X=X, trt=trt, y=y))
}




#' Format "long" panel data into "wide" program evaluation matrices with staggered adoption
#' @param outcome Name of outcome column
#' @param trt Name of treatment column
#' @param unit Name of unit column
#' @param time Name of time column
#' @param data Panel data as dataframe
#'
#' @return \itemize{
#'          \item{"X"}{Matrix of pre-treatment outcomes}
#'          \item{"trt"}{Vector of treatment assignments}
#'          \item{"y"}{Matrix of post-treatment outcomes}
#'         }
format_data_stag <- function(outcome, trt, unit, time, data) {


    ## get first treatment times
    trt_time <- data %>%
        group_by(!!unit) %>%
        summarise(trt_time=(!!time)[(!!trt) == 1][1]) %>%
        mutate(trt_time=replace_na(trt_time, Inf))

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
                                filter(is.finite(trt_time))) %>%
        filter(!!time < t_int) %>%
        mutate(trt=1-!!trt) %>%
        select(!!unit, !!time, trt_time, trt) %>%
        spread(!!time, trt) %>% 
        ## arrange(!!unit) %>% 
        select(-trt_time, -!!unit) %>%
        as.matrix()
    
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
#'
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