################################################################################
## Fitting, plotting, summarizing staggered synth
################################################################################





#' Fit staggered synth
#' @param form outcome ~ treatment | auxillary covariates
#' @param unit Name of unit column
#' @param time Name of time column
#' @param data Panel data as dataframe
#' @param opts_weights Optional options for fitting synth weights
#'
#' @return augsynth object that contains:
#'         \itemize{
#'          \item{"weights"}{Ridge ASCM weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#'          \item{"mhat"}{Outcome model estimate}
#'          \item{"data"}{Panel data as matrices}
#'         }
#' @export
multisynth <- function(form, unit, time, data,
                       relative=T, gap=NULL,
                       lambda=NULL, alpha=0.5,
                       opts_weights=NULL) {
    
    call_name <- match.call()
    
    form <- Formula::Formula(form)
    unit <- enquo(unit)
    time <- enquo(time)
    
    ## format data
    outcome <- terms(formula(form, rhs=1))[[2]]
    trt <- terms(formula(form, rhs=1))[[3]]
    wide <- format_data_stag(outcome, trt, unit, time, data)

    ## if gap is NULL set it to be the size of X
    if(is.null(gap)) {
        gap <- ncol(wide$X) + 1
    }
    
    ## fit multisynth
    opts_weights <- c(opts_weights,
                      list(link="logit",
                           regularizer="nuc",
                           nlambda=20, lambda.min.ratio=1e-3,
                           opts=list()))

    msynth <- do.call(multisynth_, c(list(X=wide$X, trt=wide$trt, mask=wide$mask,
                                        lambda=lambda, alpha=alpha),
                                   opts_weights))    
    msynth$data <- wide
    msynth$data$time <- data %>% distinct(!!time) %>% pull(!!time)
    msynth$call <- call_name
    msynth$relative <- relative
    msynth$gap <- gap
    ##format output
    class(msynth) <- "multisynth"
    return(msynth)
}

#' Get prediction of average outcome under control
#' @param multisynth Fit multisynth object
#' @param relative Whether to aggregate estimates according to calendar or relative time
#' @param att Whether to estimate the ATT or the missing counterfactual
#'
#' @return Vector of predicted post-treatment control averages for each treatment group
predict.multisynth <- function(multisynth, relative=NULL, att=F) {


    if(is.null(relative)) {
        relative <- multisynth$relative
    }
    gap <- tmp$gap
    d <- ncol(multisynth$data$X)
    fulldat <- cbind(multisynth$data$X, multisynth$data$y)
    ttot <- ncol(fulldat)

    grps <- unique(multisynth$data$trt) %>% sort()
    J <- length(grps) - 1
    n1 <- multisynth$data$trt[is.finite(multisynth$data$trt)] %>%
        table() %>% as.numeric()
    fullmask <- cbind(multisynth$data$mask, matrix(0, nrow=J, ncol=(ttot-d)))
    

    ## get weighted average estimates
    mu0hat <- lapply(multisynth$weights,
                     function(w) t(fulldat) %*% w)

    ## estimate the post-treatment values to gett att estimates
    mu1hat <- vapply(1:J,
                     function(j) colMeans(fulldat[multisynth$data$trt ==grps[j],
                                                , drop=FALSE]),
                     numeric(ttot))

    tauhat <- lapply(mu0hat, function(x) mu1hat - x)
    
    ## re-index time if relative to treatment
    if(relative) {
        total_len <- ttot + d - grps[1] ## total length of predictions
        mu0hat <- lapply(mu0hat,
                         function(x) {
                             vapply(1:J,
                                    function(j) {
                                        vec <- c(rep(NA, d-grps[j]),
                                          x[1:grps[j],j],
                                          x[(grps[j]+1):(min(grps[j] + gap, ttot)), j])
                                        c(vec, rep(NA, total_len - length(vec)))
                                    },
                                    numeric(total_len))
                         })

        tauhat <- lapply(tauhat,
                         function(x) {
                             vapply(1:J,
                                    function(j) {
                                        vec <- c(rep(NA, d-grps[j]),
                                                 x[1:grps[j],j],
                                                 x[(grps[j]+1):(min(grps[j] + gap, ttot)), j])
                                        c(vec, rep(NA, total_len - length(vec)))
                                    },
                                    numeric(total_len))
                         })
    ## get the overall average estimate
    lapply(mu0hat,
           function(x) {
               avg <- apply(x, 1, function(z) sum(n1 * z, na.rm=T) / sum(n1 * !is.na(z)))
               cbind(avg, x)
           }) -> mu0hat

    lapply(tauhat,
           function(x) {
               avg <- apply(x, 1, function(z) sum(n1 * z, na.rm=T) / sum(n1 * !is.na(z)))
               cbind(avg, x)
           }) -> tauhat
        
    } else {

        ## only average currently treated units
        lapply(mu0hat, function(x) {
            avg1 <- rowSums(t(fullmask) *  x * n1) /
                rowSums(t(fullmask) *  n1)
            avg2 <- rowSums(t(1-fullmask) *  x * n1) /
                rowSums(t(1-fullmask) *  n1)
            avg <- replace_na(avg1, 0) * apply(fullmask, 2, min) +
                replace_na(avg2,0) * apply(1-fullmask, 2, max)
            cbind(avg, x)
        }) -> mu0hat

        ## only average currently treated units
        lapply(tauhat, function(x) {
            avg1 <- rowSums(t(fullmask) *  x * n1) /
                rowSums(t(fullmask) *  n1)
            avg2 <- rowSums(t(1-fullmask) *  x * n1) /
                rowSums(t(1-fullmask) *  n1)
            avg <- replace_na(avg1, 0) * apply(fullmask, 2, min) +
                replace_na(avg2,0) * apply(1-fullmask, 2, max)
            cbind(avg, x)
        }) -> tauhat        
    }
    

    if(att) {
        return(tauhat)
    } else {
        return(mu0hat)
    }
}


#' Print function for multisynth
#' @export
print.multisynth <- function(multisynth) {
    ## straight from lm
    cat("\nCall:\n", paste(deparse(multisynth$call), sep="\n", collapse="\n"), "\n\n", sep="")

    ## ## print att estimates
    ## att_post <- colMeans(augsynth$data$y[augsynth$data$trt == 1,,drop=F]) -
    ##     predict(augsynth)

    ## cat(paste("Average ATT Estimate: ",
    ##           format(round(mean(att_post),3), nsmall = 3), "\n\n", sep=""))
}

