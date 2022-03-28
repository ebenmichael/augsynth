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
#' @param norm_pool Normalizing value for pooled objective, default: number of treated units squared
#' @param norm_sep Normalizing value for separate objective, default: number of treated units
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
multisynth_qp <- function(X, trt, mask, Z = NULL, n_leads=NULL, n_lags=NULL,
                          relative=T, nu=0, lambda=0, V = NULL, time_cohort = FALSE,
                          donors = NULL, norm_pool = NULL, norm_sep = NULL,
                          verbose = FALSE, 
                          eps_rel=1e-4, eps_abs=1e-4) {

    # if Z has no columns then set it to NULL
    if(!is.null(Z)) {
      if(ncol(Z) == 0) {
        Z <- NULL
      }
    }

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
    V <- make_V_matrix(n_lags, V)
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
    if(is.null(norm_sep)) {
      norm_sep <- 1#J
    }
    if(is.null(norm_pool)) {
      norm_pool <- 1#J ^ 2
    }
    n1 <- sapply(1:J, function(j) length(which_t[[j]]))

    # if no specific donors passed in,
    # then all donors treated after n_lags are eligible
    if(is.null(donors)) {
      donors <- get_eligible_donors(trt, time_cohort, n_leads)
    }

    ## handle X differently if it is a list
    if(typeof(X) == "list") {
        x_t <- lapply(1:J, function(j) colSums(X[[j]][which_t[[j]], mask[j,]==1, drop=F]))
        
        # Xc contains pre-treatment data for valid donor units
        Xc <- lapply(1:nrow(mask),
                 function(j) X[[j]][donors[[j]], mask[j,]==1, drop=F])
        
        # std dev of outcomes for first treatment time
        sdx <- sd(X[[1]][is.finite(trt)])
    } else {
        x_t <- lapply(1:J, function(j) colSums(X[which_t[[j]], mask[j,]==1, drop=F]))        
        
        # Xc contains pre-treatment data for valid donor units
        Xc <- lapply(1:nrow(mask),
                 function(j) X[donors[[j]], mask[j,]==1, drop=F])

        # std dev of outcomes
        sdx <- sd(X[is.finite(trt)])
    }

    # get covariates for donors
    if(!is.null(Z)) {
      # scale covariates to have same variance as pure control outcomes
      Z_scale <- sdx * apply(Z, 2, 
        function(z) (z - mean(z[!is.finite(trt)])) / sd(z[!is.finite(trt)]))

      z_t <- lapply(1:J, function(j) colSums(Z_scale[which_t[[j]], , drop = F]))
      Zc <- lapply(1:J, function(j) Z_scale[donors[[j]], , drop = F])
    } else {
      z_t <- lapply(1:J, function(j) c(0))
      Zc <- lapply(1:J, function(j) Matrix::Matrix(0,
                                                   nrow = sum(donors[[j]]),
                                                   ncol = 1))
    }
    dz <- ncol(Zc[[1]])

    # replace NA values with zero
    x_t <- lapply(x_t, function(xtk) tidyr::replace_na(xtk, 0))
    Xc <- lapply(Xc, function(xck) apply(xck, 2, tidyr::replace_na, 0))
    ## make matrices for QP
    n0s <- sapply(Xc, nrow)
    if(any(n0s == 0)) {
      stop("Some treated units have no possible donor units!")
    }
    n0 <- sum(n0s)

    const_mats <- make_constraint_mats(trt, grps, n_leads, n_lags, Xc, Zc, d, n1)
    Amat <- const_mats$Amat
    lvec <- const_mats$lvec
    uvec <- const_mats$uvec

    ## quadratic balance measures

    qvec <- make_qvec(Xc, x_t, z_t, nu, n_lags, d, V, norm_pool, norm_sep)

    Pmat <- make_Pmat(Xc, x_t, dz, nu, n_lags, lambda, d, V, norm_pool, norm_sep)

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
    imbalance <- vapply(1:J,
                        function(j) {
                            dj <- length(x_t[[j]])
                            ndim <- min(dj, n_lags)
                            c(numeric(d-ndim),
                            x_t[[j]][(dj-ndim+1):dj] -
                                t(Xc[[j]][,(dj-ndim+1):dj, drop = F]) %*% 
                                out$x[(nj0cumsum[j] + 1):nj0cumsum[j + 1]])
                        },
                        numeric(d))
    avg_imbal <- rowMeans(t(t(imbalance)))

    Vsq <- t(V) %*% V
    global_l2 <- c(sqrt(t(avg_imbal[(d - n_lags + 1):d]) %*% Vsq %*%
                          avg_imbal[(d - n_lags + 1):d])) / sqrt(d)
    avg_l2 <- mean(apply(imbalance, 2,
                  function(x) c(sqrt(t(x[(d - n_lags + 1):d]) %*% Vsq %*%
                                x[(d - n_lags + 1):d]))))
    ind_l2 <- sqrt(mean(
      apply(imbalance, 2,
      function(x) c(x[(d - n_lags + 1):d] %*% Vsq %*% x[(d - n_lags + 1):d]) /
          sum(x[(d - n_lags + 1):d] != 0))))
    # pad weights with zeros for treated units and divide by number of treated units
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
                   ind_l2 = ind_l2,
                   avg_l2 = avg_l2,
                   V = V)

    if(!is.null(Z)) {
      # imbalance in auxiliary covariates
      z_t <- sapply(1:J, function(j) colMeans(Z[which_t[[j]], , drop = F]))
      imbal_z <- z_t - t(Z) %*% weights
      avg_imbal_z <- rowSums(t(t(imbal_z) * n1)) / sum(n1)
      global_l2_z <- sqrt(sum(avg_imbal_z ^ 2))
      ind_l2_z <- sum(apply(imbal_z, 2, function(x) sqrt(sum(x ^ 2))))
      imbal_z <- cbind(avg_imbal_z, imbal_z)
      rownames(imbal_z) <- colnames(Z)

      output$imbalance_aux <- imbal_z
      output$global_l2_aux <- global_l2_z
      output$ind_l2_aux <- ind_l2_z
    }
    

    
    

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
make_constraint_mats <- function(trt, grps, n_leads, n_lags, Xc, Zc, d, n1) {

    J <- length(grps)
    idxs0  <- trt  > n_leads + min(grps)

    n0 <- sum(idxs0)

    ## sum to n1 constraint
    A1 <- do.call(Matrix::bdiag, lapply(1:(J), function(x) rep(1, n0)))
    A1 <- Matrix::bdiag(lapply(1:J, function(j) rep(1, nrow(Xc[[j]]))))
    
    Amat <- as.matrix(Matrix::t(A1))
    Amat <- Matrix::rbind2(Matrix::t(A1), Matrix::Diagonal(nrow(A1)))

    dz <- ncol(Zc[[1]])
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
    A_trans <- Matrix::cbind2(
      Matrix::cbind2(A_trans1, A_trans2),
      Matrix::Matrix(0, nrow = nrow(A_trans1), ncol = dz * J))

    # constraints for transformed weights on auxiliary covariates
    A_transz <- Matrix::t(Matrix::bdiag(Zc))
    A_transz <- Matrix::cbind2(
      Matrix::cbind2(A_transz, 
                     Matrix::Matrix(0, nrow = nrow(A_transz), ncol = sum_tj)),
      -Matrix::Diagonal(dz * J))

    # add in zero columns for transformated weights
    Amat <- Matrix::cbind2(Amat, 
                           Matrix::Matrix(0,
                                          nrow = nrow(Amat),
                                          ncol = sum_tj + dz * J))
    Amat <- Matrix::rbind2(Matrix::rbind2(Amat, A_trans), A_transz)

    lvec <- c(n1, # sum to n1 constraint
              numeric(nrow(A1)), # lower bound by zero
              numeric(sum_tj), # constrain transformed weights
              numeric(dz * J) # constrain transformed weights
             ) 
    
    uvec <- c(n1, #sum to n1 constraint
              rep(Inf, nrow(A1)),
              numeric(sum_tj), # constrain transformed weights
              numeric(dz * J) # constrain transformed weights
              )


    return(list(Amat = Amat, lvec = lvec, uvec = uvec))
}

#' Make the vector in the QP
#' @param Xc List of outcomes for possible comparison units
#' @param x_t List of outcomes for treated units
#' @param nu Hyperparameter between global and individual balance
#' @param n_lags Number of lags to balance
#' @param d Largest number of pre-intervention time periods
#' @param V Scaling matrix
#' @param norm_pool Normalizing value for pooled objective
#' @param norm_sep Normalizing value for separate objective
#' @noRd
make_qvec <- function(Xc, x_t, z_t, nu, n_lags, d, V, norm_pool, norm_sep) {

    J <- length(x_t)
    Vsq <- t(V) %*% V
    qvec <- lapply(1:J,
                   function(j) {
                       dj <- length(x_t[[j]])
                       ndim <- min(dj, n_lags)
                       max_dim <- min(d, n_lags)
                       vec <- x_t[[j]][(dj - ndim + 1):dj] / ndim
                       Vsq %*% c(numeric(max_dim - ndim), vec)
                   })

    avg_target_vec <- lapply(x_t,
                            function(xtk) {
                                dk <- length(xtk)
                                ndim <- min(dk, n_lags)
                                max_dim <- min(d, n_lags)
                                c(numeric(max_dim - ndim), 
                                    xtk[(dk - ndim + 1):dk])
                            }) %>% reduce(`+`) %*% Vsq
    qvec_avg <- rep(avg_target_vec, J)
    # qvec <- - (nu * qvec_avg / n_lags + (1 - nu) * reduce(qvec, c))
    # qvec <- - (nu * qvec_avg / (J ^ 2 * n_lags) +
    #           (1 - nu) * reduce(qvec, c) / J)
    qvec <- - (nu * qvec_avg / (norm_pool * n_lags * J ^ 2) +
               (1 - nu) * reduce(qvec, c) / (norm_sep * J))

    qvec_avg_z <- z_t %>% reduce(`+`)
    qvec_avg_z <- rep(qvec_avg_z, J)
    # qvec_z <- - (nu * qvec_avg_z + (1 - nu) * reduce(z_t, c)) / length(z_t[[1]])
    # qvec_z <- - (nu * qvec_avg_z / J ^2 +
    #              (1 - nu) * reduce(z_t, c) / J) / length(z_t[[1]])
    qvec_z <- - (nu * qvec_avg_z / (norm_pool * J ^ 2) +
                 (1 - nu) * reduce(z_t, c) / (norm_sep * J)) / length(z_t[[1]])

    total_ctrls <- lapply(Xc, nrow) %>% reduce(`+`)
    return(c(numeric(total_ctrls), qvec, qvec_z))
}


#' Make the matrix in the QP
#' @param Xc List of outcomes for possible comparison units
#' @param x_t List of outcomes for treated units
#' @param nu Hyperparameter between global and individual balance
#' @param n_lags Number of lags to balance
#' @param lambda Regularization hyperparameter
#' @param d Largest number of pre-intervention time periods
#' @param V Scaling matrix
#' @param norm_pool Normalizing value for pooled objective
#' @param norm_sep Normalizing value for separate objective
#' @noRd
make_Pmat <- function(Xc, x_t, dz, nu, n_lags, lambda, d, V,
                      norm_pool, norm_sep) {

    J <- length(x_t)

    Vsq <- t(V) %*% V
    ndims <- vapply(1:J,
                    function(j) min(length(x_t[[j]]), n_lags),
                    numeric(1))
    max_dim <- min(d, n_lags)
    total_dim <- sum(ndims)
    total_dim <- max_dim * J
    V1 <- Matrix::bdiag(lapply(ndims, 
                        function(ndim) Matrix::Diagonal(max_dim, 1 / ndim)))
    V1 <- Matrix::bdiag(lapply(ndims, function(ndim) Vsq / ndim))
    tile_sparse <- function(j) {
        kronecker(Matrix::Matrix(1, nrow = j, ncol = j), Vsq)
    }
    tile_sparse_cov <- function(d, j) {
        kronecker(Matrix::Matrix(1, nrow = j, ncol = j),
                  Matrix::Diagonal(d))
    }
    V2 <- tile_sparse(J) / n_lags
    # Pmat <- nu * V2 + (1 - nu) * V1
    # Pmat <- nu * V2 / J ^ 2 + (1 - nu) * V1 / J
    Pmat <- nu * V2 / (norm_pool * J ^ 2) + (1 - nu) * V1 / (norm_sep * J)
    V1_z <- Matrix::Diagonal(dz * J, 1 / dz)
    V2_z <- tile_sparse_cov(dz, J) / dz
    # Pmat_z <- nu * V2_z + (1 - nu) * V1_z
    # Pmat_z <- nu * V2_z / J ^ 2 + (1 - nu) * V1_z / J
    Pmat_z <- nu * V2_z / (norm_pool * J ^ 2) + (1 - nu) * V1_z / (norm_sep * J)
    # combine
    total_ctrls <- lapply(Xc, nrow) %>% reduce(`+`)
    Pmat <- Matrix::bdiag(Matrix::Matrix(0, nrow = total_ctrls,
                                         ncol = total_ctrls),
                          Pmat, Pmat_z)
    I0 <- Matrix::bdiag(Matrix::Diagonal(total_ctrls),
                        Matrix::Matrix(0, nrow = total_dim + dz * J,
                                          ncol = total_dim + dz * J))
    return(Pmat + lambda * I0)

}