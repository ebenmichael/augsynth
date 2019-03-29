################################################################################
## Solve the multisynth problem as a QP
################################################################################


#' Internal function to fit synth with staggered adoption with a QP solver
#' @param X Matrix of pre-final intervention outcomes
#' @param trt Vector of treatment levels/times
#' @param mask Matrix with indicators for observed pre-intervention times for each treatment group
#' @param gap Number of time periods after treatment to impute control values.
#'            For units treated at time T_j, all units treated after T_j + gap
#'            will be used as control values. If larger than the number of periods,
#'            only never never treated units (pure controls) will be used as comparison units
#' @param relative Whether to re-index time according to treatment date, default T
#' @param alpha Hyper-parameter that controls trade-off between overall and individual balance.
#'              Larger values of alpha place more emphasis on individual balance.
#'              Balance measure is
#'                alpha ||global|| + (1-alpha) ||individual||
#'              Default: 0
#' @param lambda Regularization hyper-parameter. Default, 0
#' 
#'
#' @importMethodsFrom Matrix %*%
multisynth_qp <- function(X, trt, mask, gap=NULL, relative=T, alpha=0, lambda=0) {

    n <- dim(X)[1]

    d <- dim(X)[2]
    if(is.null(gap)) {
        gap <- d+1
    } else if(gap > d) {
        gap <- d+1
    }

    ## average over treatment times/groups
    grps <- sort(unique(trt[is.finite(trt)]))

    J <- length(grps)


    x_t <- lapply(1:J, function(j) colSums(X[trt ==grps[j], mask[j,]==1, drop=F]))    

    ## All possible donor units for all treatment groups
    Xc <- lapply(1:nrow(mask),
                 function(j) X[trt > gap + min(grps), mask[j,]==1, drop=F])

    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))
    
    ## weight the global parameters by the number of treated units still untreated in calander time
    nw <- sapply(1:max(trt[is.finite(trt)]), function(t) sum(trt[is.finite(trt)] >= t))

    ## reverse to get relative absolute time
    if(relative) {
        nw <- rev(nw)
    }


    
    ## make matrices for QP
    n0 <- sum(!is.finite(trt))
    idxs0  <- trt  > gap + min(grps)
    n0 <- sum(idxs0)
    
    ## sum to n1 constraint
    A1 <- do.call(Matrix::bdiag, lapply(1:(J+1), function(x) rep(1, n0)))
    Amat <- as.matrix(Matrix::t(A1))
    Amat <- Matrix::rbind2(Matrix::t(A1), Matrix::Diagonal(nrow(A1)))

    ## sum to global weights
    A2 <- Matrix::kronecker(t(c(-1, rep(1, J))), Matrix::Diagonal(n0))
    Amat <- Matrix::rbind2(Amat, A2)

    ## zero out weights for inelligible units
    A3 <- do.call(Matrix::bdiag,
                  c(lapply(1:J,
                           function(j) {
                               z <- numeric(n0)
                               z[trt[idxs0] <= grps[j] + gap] <- 1
                               z
                           })))

    
    Amat <- Matrix::rbind2(Amat,
                           Matrix::cbind2(Matrix::Matrix(0, ncol=n0, nrow=J),
                                          Matrix::t(A3)))
    
    lvec <- c(sum(n1), n1, # sum to n1 constraint
              numeric(nrow(A1)), # lower bound by zero
              numeric(n0), # sum to global weights
              numeric(J)) # zero out weights for impossible donors
    
    uvec <- c(sum(n1), n1, #sum to n1 constraint
              rep(sum(n1), n0), sapply(n1, function(n1j) rep(n1j, n0)), # upper bound by n1j
              numeric(n0), # sum to global weights
              numeric(J)) # zero out weights for impossible donors

    ## quadratic balance measures

    dvec <- lapply(1:J, function(j) c((1-alpha) / length(x_t[[j]]) * Xc[[j]] %*% x_t[[j]]))

    dvec_avg <- lapply(x_t,
                       function(xtk) {
                           lapply(1:J,
                                  function(j)  {
                                      dk <- length(xtk)
                                      dj <- ncol(Xc[[j]])
                                      ndim <- min(dk, dj)
                                      ## relative time: inner product from the end
                                      if(relative) {
                                          Xc[[j]][,(dj-ndim+1):dj] %*%
                                              xtk[(dk-ndim+1):dk] *
                                              n1[j] / sum(n1)
                                      } else {
                                          ## absolute time: inner product from start
                                          Xc[[j]][,1:ndim] %*%
                                              xtk[1:ndim] *
                                              n1[j] / sum(n1)
                                      }
                                  }) %>% reduce(`+`)
                       }) %>% reduce(`+`)
    
    dvec <- -c(alpha * dvec_avg, reduce(dvec,c))


    V1 <- do.call(Matrix::bdiag,
                  lapply(Xc, function(x) x / sqrt(ncol(x))))

    lapply(1:J,
           function(k) {
               lapply(1:J,
                      function(j) {
                          dk <- ncol(Xc[[k]])
                          dj <- ncol(Xc[[j]])
                          ndim <- min(dk, dj)
                          if(relative) {
                              ## inner product from end
                              Xc[[k]][,(dk-ndim+1):dk] %*%
                                  t(Xc[[j]][,(dj-ndim+1):dj]) *
                                  n1[k] * n1[j] / sum(n1)^2
                          } else {
                              ## inner product from start
                              Xc[[k]][,1:ndim] %*%
                                  t(Xc[[j]][,1:ndim]) *
                                  n1[k] * n1[j] / sum(n1)^2
                          }
                      }) %>% reduce(`+`)
           }) %>% reduce(`+`) -> V2


    Hmat <- Matrix::bdiag(list(alpha * (V2 + lambda * Matrix::Diagonal(nrow(V2))),
    (1-alpha) * (V1 %*% Matrix::t(V1) + lambda * Matrix::Diagonal(nrow(V1)))))



    ## Optimize
    settings <- osqp::osqpSettings(verbose = FALSE, eps_abs=1e-8, eps_rel = 1e-8)
    out <- osqp::solve_osqp(Hmat, dvec, Amat, lvec, uvec, pars=settings)

    ## get weights
    weights <- matrix(out$x[-(1:n0)], nrow=n0)

    ## compute imbalance
    if(relative) {
        imbalance <- vapply(1:J,
                            function(j) c(numeric(d-length(x_t[[j]])),
                                          x_t[[j]] -  t(Xc[[j]]) %*% weights[,j]),
                            numeric(d))
    } else {
        imbalance <- vapply(1:J,
                            function(j) c(x_t[[j]] -  t(Xc[[j]]) %*% weights[,j],
                                          numeric(d-length(x_t[[j]]))),
                            numeric(d))
    }

    avg_imbal <- rowSums(t(t(imbalance)))


    ## pad weights with zeros for treated units and divide by number of treated units
    weights <- matrix(0, nrow=n, ncol=J)
    weights[idxs0,] <- matrix(out$x[-(1:n0)], nrow=n0)
    weights <- t(t(weights) / n1)
    
    output <- list(weights=weights,
                   imbalance=cbind(avg_imbal, imbalance))
    
    return(output)
    
}
