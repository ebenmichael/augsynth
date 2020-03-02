##############################################################################
## Outcome regression with multiple treated units
##############################################################################

#' Fit a time regression
#' @param X Matrix of outcomes
#' @param trt Vector of treatment status for each unit
#' @param n_leads How long past treatment effects should be estimated for
#' @param reg_param Regularization hyperparameter
#' @param lowlim Lower bound for coefs
#' @param uplim upper bound for coefs
#' @param ... Extra optimization hyperparameters
#' @noRd
#' @return \itemize{
#'           \item{y0hat }{List of predicted outcome under control}
#'           \item{residuals }{List of residuals}
#'           \item{params }{Regression parameters}}
fit_time_reg <- function(X, trt, n_leads, reg_param, lowlim = 0, uplim = 1, ...) {

    grps <- trt[is.finite(trt)]
    J <- length(grps)
    tmax <- max(trt[is.finite(trt)])

    # fit QP
    reg_weights <- fit_time_reg_qp_(X, trt, n_leads, lowlim, uplim, reg_param, ...)

    # get predicted outcomes (repeated as a matrix) and residuals
    y0hat <- lapply(1:J, 
        function(j) {
            # compute time fixed effects from pure controls
            time_eff <- matrix(colMeans(X[!is.finite(trt),]),
                              nrow=nrow(X), ncol=ncol(X),
                              byrow=T)
            Xj <- X - time_eff
            zero_mat <- matrix(0, nrow = nrow(X), ncol = (tmax - grps[j]))
            Xj <- cbind(zero_mat, Xj[, 1:grps[j], drop = F])
            # take out pure control means
            y0hatj <- Xj %*% reg_weights[,j, drop = F]
            matrix(y0hatj, nrow=nrow(X), ncol=ncol(X)) + time_eff
        })

    residuals <- lapply(1:J, function(j) X - y0hat[[j]])

    return(list(y0hat = y0hat,
                residuals = residuals,
                time_weights = reg_weights))
}


#' Fit a time regression
#' @param X Matrix of outcomes
#' @param trt Vector of treatment status for each unit
#' @param n_leads How long past treatment effects should be estimated for
#' @param reg_param Regularization hyperparameter
#' @param lowlim Lower bound for coefs
#' @param uplim upper bound for coefs
#' @param ... Extra optimization hyperparameters
#' @noRd
#' @return reg_weights Fitted regression weights
fit_time_reg_qp_ <- function(X, trt, n_leads, lowlim, uplim, reg_param, ...) {

    grps <- trt[is.finite(trt)]
    J <- length(grps)
    ttot <- ncol(X)
    max_trt <- max(grps)

    # get data in the right form
    data_mats <- collect_data(X, trt, n_leads)

    # create constraint matrices
    constraints <- make_constraints(J, grps, lowlim, uplim)

    # get components of QP
    Qmat <- get_Qmat(data_mats$pre_mats)
    pvec <- get_pvec(data_mats$pre_mats, data_mats$post_vecs)

    I0 <- get_regularization_matrix(J, max_trt, reg_param)
    # add in regularization
    # I0 <- Matrix::bdiag(reg_param1 * Matrix::Diagonal(max_trt), 
    #                     reg_param2 * Matrix::Diagonal(J * max_trt))
    

    
    Qmat <- Qmat + I0

    # fit QP
    settings <- osqp::osqpSettings(verbose = FALSE, ...)
    out <- osqp::solve_osqp(Qmat, pvec, constraints$Amat, 
                            constraints$lvec, constraints$uvec, 
                            pars=settings)

    # collect as matrix
    # reg_weights <- matrix(out$x, ncol = J + 1)
    reg_weights <- matrix(out$x, ncol = J)
    # pooled <- reg_weights[,1]
    # add in common component
    # reg_weights <- reg_weights[, 1] + reg_weights[, -1]
    # reverse to calendar time
    # reg_weights <- reg_weights[nrow(reg_weights):1, ]
    return(reg_weights)
}


#' Organize data for the QP
#' @param X Matrix of outcomes
#' @param trt Vector of treatment status for each unit
#' @param n_leads How long past treatment effects should be estimated for
#' @noRd
collect_data <- function(X, trt, n_leads) {

    grps <- trt[is.finite(trt)]
    J <- length(grps)
    ttot <- ncol(X)
    max_trt <- max(grps)


    # sapply(1:ncol(X), 
    #        function(tj) {
    #            mean(X[trt >= tj])
    #        }) -> ctrl_means

    # X <- t(t(X) - ctrl_means)

    # get pre-treatment matrices
    lapply(grps, function(tj) {
        # donor unit pre tj outcomes
        idxs <- trt > tj + n_leads
        pre_mat <- cbind(#1, 
                         matrix(0, nrow = nrow(X), ncol = (max_trt - tj)),
                         X[, 1:tj, drop = F])
        # subtract out pure control means
        pre_mat <- t(t(pre_mat) - colMeans(pre_mat[!is.finite(trt),,drop = F]))
        # restrict to units that won't be treated w/in n_leads
        pre_mat[idxs,,drop = F]
    }) -> pre_mats


    # get post treatment averages
    lapply(grps, 
           function(tj) {
               # avg of donor units post tj outcomes
               idxs <- trt > tj + n_leads
               donors <- rowMeans(X[, (tj + 1):(tj + n_leads), drop = F])
               # subtract out pure control means
               donors <- donors - mean(donors[!is.finite(trt)])
               # restrict to units that won't be treated w/in n_leads
               donors[idxs]
           }) -> post_vecs
    return(list(pre_mats = pre_mats, post_vecs = post_vecs))
}


get_Qmat <- function(pre_mats) {

    #### matrix in QP
    cov_mats <- lapply(pre_mats, function(x) t(x) %*% x)
    # unit specific covariance matrices
    Qmat <- Matrix::bdiag(cov_mats)
    return(Qmat)
}


get_Qmat_pool <- function(pre_mats) {

    #### matrix in QP
    cov_mats <- lapply(pre_mats, function(x) t(x) %*% x)
    pooled_cov <- Reduce(`+`, cov_mats)
    cov_mats_bind <-do.call(rbind, cov_mats)
    # unit specific covariance matrices
    Qmat <- Matrix::bdiag(cov_mats)
    # pooling terms
    Qmat <- rbind(t(cov_mats_bind), Qmat)
    Qmat <- cbind(rbind(pooled_cov, cov_mats_bind), Qmat)

    return(Qmat)
}


get_pvec <- function(pre_mats, post_vecs) {
    # vector in QP
    lapply(1:length(pre_mats), function(j) {
        t(pre_mats[[j]]) %*% post_vecs[[j]]
    }) -> pvec_list

    pvec <- do.call(c, pvec_list)
    
    return(-2 * pvec)
}


get_pvec_pool <- function(pre_mats, post_vecs) {
    # vector in QP
    lapply(1:length(pre_mats), function(j) {
        t(pre_mats[[j]]) %*% post_vecs[[j]]
    }) -> pvec_list

    pvec_pool <- Reduce(`+`, pvec_list)
    pvec <- do.call(c, pvec_list)
    pvec <- c(pvec_pool, pvec)

    return(-2 * pvec)
}


make_constraints <- function(J, grps, lowlim, uplim) {

    tmax <- max(grps)
    # sum to 1 constraints
    A1 <- Matrix::t(Matrix::bdiag(lapply(1:J, function(j) c(0, rep(1, tmax)))))
    A1 <- Matrix::t(Matrix::bdiag(lapply(1:J, function(j) rep(1, tmax))))
    l1 <- rep(1, J)
    # l1 <- rep(-Inf, J)
    u1 <- rep(1, J)
    # u1 <- rep(Inf, J)

    # upper lower limits
    diag_w_intercept <- Matrix::bdiag(list(0, Matrix::Diagonal(tmax)))[-1, ]
    A2 <- Matrix::bdiag(lapply(1:J, function(j) diag_w_intercept))
    A2 <- Matrix::Diagonal(J * tmax)
    # make sure that only weighting times that exist
    l2 <- sapply(1:J, function(j) {
        c(rep(0, tmax - grps[j]), rep(lowlim, grps[j]))
    })
    
    u2 <- sapply(1:J, function(j) {
        c(rep(0, tmax - grps[j]), rep(uplim, grps[j]))
    })

    # combine

    Amat <- rbind(A1, A2)
    lvec <- c(l1, l2)
    uvec <- c(u1, u2)
    return(list(Amat = Amat, lvec = lvec, uvec = uvec))
}

make_constraints_pool <- function(J, grps, lowlim, uplim) {

    tmax <- max(grps)
    # sum to 1 constraints
    A1 <- cbind(0,
                matrix(1, ncol = tmax, nrow = J),
                Matrix::t(Matrix::bdiag(lapply(1:J, function(j) c(0, rep(1, tmax))))))
    
    l1 <- rep(1, J)
    # l1 <- rep(-Inf, J)
    u1 <- rep(1, J)
    # u1 <- rep(Inf, J)

    # upper lower limits
    diag_w_intercept <- Matrix::bdiag(list(0, Matrix::Diagonal(tmax)))[-1, ]
    pool_A2 <- do.call(rbind, lapply(1:J, function(j) diag_w_intercept))
    A2 <- Matrix::bdiag(lapply(1:J, function(j) diag_w_intercept))
    A2 <- cbind(pool_A2, A2)

    # restrict global intercept to 0
    A2 <- rbind(c(1, numeric(ncol(A2) - 1)), A2)
    
    # make sure that only weighting times that exist
    l2 <- sapply(1:J, function(j) {
        c(rep(0, tmax - grps[j]), rep(lowlim, grps[j]))
    })
    l2 <- c(0, l2)
    
    u2 <- sapply(1:J, function(j) {
        c(rep(0, tmax - grps[j]), rep(uplim, grps[j]))
    })
    u2 <- c(0, u2)
    # l2 <- rep(lowlim, J * tmax)
    # u2 <- rep(uplim, J  * tmax)

    # combine

    Amat <- rbind(A1, A2)
    lvec <- c(l1, l2)
    uvec <- c(u1, u2)
    return(list(Amat = Amat, lvec = lvec, uvec = uvec))
}

get_regularization_matrix <- function(J, max_trt, reg_param) {

    single_reg_mat <- Matrix::bdiag(list(0, Matrix::Diagonal(max_trt)))
    I0 <- Matrix::bdiag(lapply(1:J,function(j) single_reg_mat))
    I0 <- reg_param * Matrix::Diagonal(J * max_trt)
    return(reg_param * I0)
}

get_regularization_matrix_pool <- function(J, max_trt, reg_param1, reg_param2) {

    single_reg_mat <- Matrix::bdiag(list(0, Matrix::Diagonal(max_trt)))
        grp_reg_mats <- Matrix::bdiag(lapply(1:J,function(j) single_reg_mat))
        I0 <- Matrix::bdiag(reg_param1 * single_reg_mat,
                            reg_param2 * grp_reg_mats)

    return(I0)
}