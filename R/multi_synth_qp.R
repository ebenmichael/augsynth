################################################################################
## Solve the multisynth problem as a QP
################################################################################


#' Internal function to fit synth with staggered adoption with a QP solver
#' @param X Matrix of pre-final intervention outcomes, or list of such matrices after transformations
#' @param trt Vector of treatment levels/times
#' @param mask Matrix with indicators for observed pre-intervention times for each treatment group
#' @param n_leads Number of time periods after treatment to impute control values.
#'            For units treated at time T_j, all units treated after T_j + n_leads
#'            will be used as control values. If larger than the number of periods,
#'            only never never treated units (pure controls) will be used as comparison units
#' @param n_lags Number of pre-treatment periods to balance, default is to balance all periods
#' @param relative Whether to re-index time according to treatment date, default T
#' @param nu Hyper-parameter that controls trade-off between overall and individual balance.
#'              Larger values of nu place more emphasis on individual balance.
#'              Balance measure is
#'                nu ||global|| + (1-nu) ||individual||
#'              Default: 0
#' @param lambda Regularization hyper-parameter. Default, 0
#' @param time_cohort Whether to average synthetic controls into time cohorts
#' @param verbose Whether to print logs for osqp
#' @param eps_rel Relative error tolerance for osqp
#' @param eps_abs Absolute error tolerance for osqp
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Matrix of unit weights}
#'          \item{"imbalance"}{Matrix of overall and group specific imbalance}
#'          \item{"global_l2"}{Imbalance overall}
#'          \item{"ind_l2"}{Matrix of imbalance for each group}
#'         }
multisynth_qp <- function(X, trt, mask, n_leads=NULL, n_lags=NULL,
                          relative=T, nu=0, lambda=0, time_cohort = FALSE,
                          donors = NULL,
                          verbose = FALSE, 
                          eps_rel=1e-4, eps_abs=1e-4) {

    n <- if(typeof(X) == "list") dim(X[[1]])[1] else dim(X)[1]
    d <- if(typeof(X) == "list") dim(X[[1]])[2] else dim(X)[2]

    if(is.null(n_leads)) {
        n_leads <- d+1
    } else if(n_leads > d) {
        n_leads <- d+1
    }
    if(is.null(n_lags)) {
        n_lags <- d
    } else if(n_lags > d) {
        n_lags <- d
    }

    ## treatment times
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
    n1 <- sapply(1:J, function(j) length(which_t[[j]]))

    # only allow weights on eligible donors
    # if null, then all donors treated after n_lags are eligible 
    if(is.null(donors)) {
      donors <- lapply(1:J, function(j) trt > n_leads + grps[j])
    } else {
      donors <- lapply(1:J,
                function(j) {
                  (trt > n_leads + grps[j]) & donors[[j]]
                }
            )
    }

    ## handle X differently if it is a list
    if(typeof(X) == "list") {
        x_t <- lapply(1:J, function(j) colSums(X[[j]][which_t[[j]], mask[j,]==1, drop=F]))
        
        # Xc contains pre-treatment data for valid donor units
        Xc <- lapply(1:nrow(mask),
                 function(j) X[[j]][donors[[j]], mask[j,]==1, drop=F])
    } else {
        x_t <- lapply(1:J, function(j) colSums(X[which_t[[j]], mask[j,]==1, drop=F]))        
        
        # Xc contains pre-treatment data for valid donor units
        Xc <- lapply(1:nrow(mask),
                 function(j) X[donors[[j]], mask[j,]==1, drop=F])
    }
    ## make matrices for QP
    n0s <- sapply(Xc, nrow)
    if(any(n0s == 0)) {
      stop("Some treated units have no possible donor units!")
    }
    n0 <- sum(n0s)

    const_mats <- make_constraint_mats(trt, grps, n_leads, n_lags, Xc, d, n1)
    Amat <- const_mats$Amat
    lvec <- const_mats$lvec
    uvec <- const_mats$uvec

    ## quadratic balance measures

    qvec <- make_qvec(Xc, x_t, nu, n_lags, d)

    Pmat <- make_Pmat(Xc, x_t, nu, n_lags, lambda, d)

    ## Optimize
    settings <- do.call(osqp::osqpSettings, 
                        c(list(verbose = verbose, 
                               eps_rel = eps_rel, 
                               eps_abs = eps_abs)))

    out <- osqp::solve_osqp(Pmat, qvec, Amat, lvec, uvec, pars = settings)

    ## get weights
    total_ctrls <- n0 * J
    weights <- matrix(out$x[1:total_ctrls], nrow = n0)
    nj0 <- as.numeric(lapply(Xc, nrow))
    nj0cumsum <- c(0, cumsum(nj0))
    ## compute imbalance
    # if(relative) {
        imbalance <- vapply(1:J,
                            function(j) {
                                dj <- length(x_t[[j]])
                                ndim <- min(dj, n_lags)
                                c(numeric(d-ndim),
                                x_t[[j]][(dj-ndim+1):dj] -
                                    t(Xc[[j]][,(dj-ndim+1):dj]) %*% 
                                    out$x[(nj0cumsum[j] + 1):nj0cumsum[j + 1]])
                            },
                            numeric(d))
    # } else {
    #     imbalance <- vapply(1:J,
    #                         function(j) c(x_t[[j]] -  t(Xc[[j]]) %*% weights[,j],
    #                                       numeric(d-length(x_t[[j]]))),
    #                         numeric(d))
    # }

    avg_imbal <- rowSums(t(t(imbalance)))


    global_l2 <- sqrt(sum(avg_imbal^2)) #/ sqrt(length(avg_imbal))
    ind_l2 <- sum(apply(imbalance, 2, function(x) sqrt(sum(x^2)))) #/ sqrt(prod(dim(imbalance)))
    ## pad weights with zeros for treated units and divide by number of treated units

    vapply(1:J,
           function(j) {
             weightj <-  numeric(n)
             weightj[donors[[j]]] <- out$x[(nj0cumsum[j] + 1):nj0cumsum[j + 1]]
             weightj
           },
           numeric(n)) -> weights

    weights <- t(t(weights) / n1)

    output <- list(weights = weights,
                   imbalance = cbind(avg_imbal, imbalance),
                   global_l2 = global_l2,
                   ind_l2 = ind_l2)

    return(output)
}


#' Create constraint matrices for multisynth QP
#' @param trt Vector of treatment levels/times
#' @param grps Treatment times
#' @param n_leads Number of time periods after treatment to impute control values.
#' @param n_lags Number of pre-treatment periods to balance
#' @param Xc List of outcomes for possible comparison units
#' @param d Max number of lagged outcomes
#' @param n1 Vector of number of treated units per cohort
#' @noRd
#' @return 
#'         \itemize{
#'          \item{"Amat"}{Linear constraint matrix}
#'          \item{"lvec"}{Lower bounds for linear constraints}
#'          \item{"uvec"}{Upper bounds for linear constraints}
#'         }
make_constraint_mats <- function(trt, grps, n_leads, n_lags, Xc, d, n1) {

    J <- length(grps)
    idxs0  <- trt  > n_leads + min(grps)

    n0 <- sum(idxs0)

    ## sum to n1 constraint
    A1 <- do.call(Matrix::bdiag, lapply(1:(J), function(x) rep(1, n0)))
    A1 <- Matrix::bdiag(lapply(1:J, function(j) rep(1, nrow(Xc[[j]]))))
    
    Amat <- as.matrix(Matrix::t(A1))
    Amat <- Matrix::rbind2(Matrix::t(A1), Matrix::Diagonal(nrow(A1)))


    # constraints for transformed weights
    A_trans1 <- do.call(Matrix::bdiag,
                       lapply(1:J,
                        function(j)  {
                            dj <- ncol(Xc[[j]])
                            ndim <- min(dj, n_lags)
                            max_dim <- min(d, n_lags)
                            mat <- Xc[[j]][, (dj - ndim + 1):dj, drop = F]
                            n0 <- nrow(mat)
                            zero_mat <- Matrix::Matrix(0, n0, max_dim - ndim)
                            Matrix::t(cbind(zero_mat, mat))
                       }))
    # sum of total number of pre-periods
    sum_tj <- min(d, n_lags) * J
    A_trans2 <- - Matrix::Diagonal(sum_tj)
    A_trans <- Matrix::cbind2(A_trans1, A_trans2)

    # add in zero columns for transformated weights
    Amat <- Matrix::cbind2(Amat, 
                           Matrix::Matrix(0,
                                          nrow = nrow(Amat),
                                          ncol = sum_tj))
    Amat <- Matrix::rbind2(Amat, A_trans)

    lvec <- c(n1, # sum to n1 constraint
              numeric(nrow(A1)), # lower bound by zero
              numeric(sum_tj) # constrain transformed weights
             ) 
    
    uvec <- c(n1, #sum to n1 constraint
              rep(Inf, nrow(A1)),
              numeric(sum_tj) # constrain transformed weights
              )


    return(list(Amat = Amat, lvec = lvec, uvec = uvec))
}

#' Make the vector in the QP
#' @param Xc List of outcomes for possible comparison units
#' @param x_t List of outcomes for treated units
#' @param nu Hyperparameter between global and individual balance
#' @param n_lags Number of lags to balance
#' @param d Largest number of pre-intervention time periods
#' @noRd
make_qvec <- function(Xc, x_t, nu, n_lags, d) {

    J <- length(x_t)

    qvec <- lapply(1:J,
                   function(j) {
                       dj <- length(x_t[[j]])
                       ndim <- min(dj, n_lags)
                       max_dim <- min(d, n_lags)
                       vec <- x_t[[j]][(dj - ndim + 1):dj] / ndim
                       c(numeric(max_dim - ndim), vec)
                   })

    avg_target_vec <- lapply(x_t,
                            function(xtk) {
                                dk <- length(xtk)
                                ndim <- min(dk, n_lags)
                                max_dim <- min(d, n_lags)
                                c(numeric(max_dim - ndim), 
                                    xtk[(dk - ndim + 1):dk])
                            }) %>% reduce(`+`)
    qvec_avg <- rep(avg_target_vec, J)

    qvec <- - (nu * qvec_avg / n_lags + (1 - nu) * reduce(qvec, c))
    total_ctrls <- lapply(Xc, nrow) %>% reduce(`+`)
    return(c(numeric(total_ctrls), qvec))
}


#' Make the matrix in the QP
#' @param Xc List of outcomes for possible comparison units
#' @param x_t List of outcomes for treated units
#' @param nu Hyperparameter between global and individual balance
#' @param n_lags Number of lags to balance
#' @param lambda Regularization hyperparameter
#' @param d Largest number of pre-intervention time periods
#' @noRd
make_Pmat <- function(Xc, x_t, nu, n_lags, lambda, d) {

    J <- length(x_t)


    ndims <- vapply(1:J,
                    function(j) min(length(x_t[[j]]), n_lags),
                    numeric(1))
    max_dim <- min(d, n_lags)
    total_dim <- sum(ndims)
    total_dim <- max_dim * J
    V1 <- Matrix::bdiag(lapply(ndims, 
                        function(ndim) Matrix::Diagonal(max_dim, 1 / ndim)))

    tile_sparse <- function(d,j) {
        kronecker(Matrix::Matrix(1, nrow = j, ncol = j), Matrix::Diagonal(d))
    }
    V2 <- tile_sparse(max_dim, J) / n_lags
    Pmat <- nu * V2 + (1 - nu) * V1
    # combine
    total_ctrls <- lapply(Xc, nrow) %>% reduce(`+`)
    Pmat <- Matrix::bdiag(Matrix::Matrix(0, nrow = total_ctrls,
                                         ncol = total_ctrls),
                          Pmat)
    I0 <- Matrix::bdiag(Matrix::Diagonal(total_ctrls),
                        Matrix::Matrix(0, nrow = total_dim, ncol = total_dim))
    return(Pmat + lambda * I0)

}