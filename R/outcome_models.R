################################################################################
## Code to fit various outcome models
################################################################################

#' Use gsynth to fit factor model for E[Y(0)|X]
#'
#' @param X Matrix of covariates/lagged outcomes
#' @param y Matrix of post-period outcomes
#' @param trt Vector of treatment indicator
#' @param r Number of factors to use (or start with if CV==1)
#' @param r.end Max number of factors to consider if CV==1
#' @param force Fixed effects (0=none, 1=unit, 2=time, 3=two-way)
#' @param CV Whether to do CV (0=no CV, 1=yes CV)
#' @param ... optional arguments for outcome model
#' @noRd
#' @return \itemize{
#'           \item{y0hat }{Predicted outcome under control}
#'           \item{params }{Regression parameters}}
fit_prog_gsynth <- function(X, y, trt, r=0, r.end=5, force=3, CV=1, ...) {
    if(!requireNamespace("gsynth", quietly = TRUE)) {
        stop("In order to fit generalized synthetic controls, you must install the gsynth package.")
    }
    extra_params = list(...)
    if (length(extra_params) > 0) {
        warning("Unused parameters when using gSynth: ", paste(names(extra_params), collapse = ", "))
    }
    
    df_x = data.frame(X, check.names=FALSE)
    df_x$unit = rownames(df_x)
    df_x$trt = rep(0, nrow(df_x))
    df_x <- df_x %>% select(unit, trt, everything())
    long_df_x = gather(df_x, time, obs, -c(unit,trt))

    df_y = data.frame(y, check.names=FALSE)
    df_y$unit = rownames(df_y)
    df_y$trt = trt
    df_y <- df_y %>% select(unit, trt, everything())
    long_df_y = gather(df_y, time, obs, -c(unit,trt))
    long_df = rbind(long_df_x, long_df_y)

    ## time comes from gathered column names as character (#107); recode to
    ## consecutive integers since the names need not be numeric strings
    long_df$time <- match(long_df$time, unique(long_df$time))
    gsyn <- gsynth::gsynth(data = long_df, Y = "obs", D = "trt",
                           index = c("unit", "time"), force = force, CV = CV, r = r)

    t0 <- dim(X)[2]
    t_final <- t0 + dim(y)[2]
    n <- dim(X)[1]
    ## get predicted outcomes
    y0hat <- matrix(0, nrow=n, ncol=(t_final-t0))
    if(!is.null(gsyn$est.co)) {
        ## gsynth < 1.3.0
        y0hat[trt==0,]  <- t(gsyn$Y.co[(t0+1):t_final,,drop=FALSE] -
                                 gsyn$est.co$residuals[(t0+1):t_final,,drop=FALSE])

        y0hat[trt==1,] <- gsyn$Y.ct[(t0+1):t_final,]

        params <- gsyn$est.co
        ## add treated prediction for whole pre-period
        params$Y.ct <- gsyn$Y.ct
    } else {
        ## gsynth >= 1.3.0 (fect): est.co is renamed est and Y.ct holds
        ## imputed Y(0) for all units, in gsyn$id order
        unit_cols <- match(df_x$unit, as.character(gsyn$id))
        if(anyNA(unit_cols)) {
            stop("Failed to map units onto the columns of gsynth's Y.ct ",
                 "(no match in gsyn$id for: ",
                 paste(df_x$unit[is.na(unit_cols)], collapse = ", "), ")")
        }
        y0hat[,] <- t(gsyn$Y.ct[(t0+1):t_final, unit_cols, drop=FALSE])

        params <- gsyn$est
        ## add treated prediction for whole pre-period
        params$Y.ct <- gsyn$Y.ct[, gsyn$tr, drop=FALSE]
        if(is.null(params$factor)) {
            ## 0-column matrix so ncol(params$factor) works downstream
            params$factor <- matrix(0, nrow=t_final, ncol=0)
        }
    }

    ## control and treated residuals
    params$ctrl_resids <- params$residuals
    params$trt_resids <- colMeans(cbind(X[trt==1,,drop=FALSE],
                                        y[trt==1,,drop=FALSE])) -
        rowMeans(params$Y.ct)
    return(list(y0hat=y0hat,
                params=params))
}




