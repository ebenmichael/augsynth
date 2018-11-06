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
