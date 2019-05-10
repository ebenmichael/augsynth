################################################################################
## Solve the multisynth problem as a QP
################################################################################


#' Internal function to fit synth with staggered adoption with a QP solver
#' @param X Matrix of pre-final intervention outcomes, or list of such matrices after transformations
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

    n <- if(typeof(X) == "list") dim(X[[1]])[1] else dim(X)[1]
    d <- if(typeof(X) == "list") dim(X[[1]])[2] else dim(X)[2]

    if(is.null(gap)) {
        gap <- d+1
    } else if(gap > d) {
        gap <- d+1
    }

    ## average over treatment times/groups
    grps <- sort(unique(trt[is.finite(trt)]))

    J <- length(grps)


    ## handle X differently if it is a list
    if(typeof(X) == "list") {
        x_t <- lapply(1:J, function(j) colSums(X[[j]][trt ==grps[j], mask[j,]==1, drop=F]))    
        
        ## All possible donor units for all treatment groups
        Xc <- lapply(1:nrow(mask),
                 function(j) X[[j]][trt > gap + min(grps), mask[j,]==1, drop=F])
    } else {
        x_t <- lapply(1:J, function(j) colSums(X[trt ==grps[j], mask[j,]==1, drop=F]))    
        
        ## All possible donor units for all treatment groups
        Xc <- lapply(1:nrow(mask),
                 function(j) X[trt > gap + min(grps), mask[j,]==1, drop=F])
    }
    
    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))
    
    ## weight the global parameters by the number of treated units still untreated in calander time
    nw <- sapply(1:max(trt[is.finite(trt)]), function(t) sum(trt[is.finite(trt)] >= t))

    ## reverse to get relative absolute time
    if(relative) {
        nw <- rev(nw)
    }


    n1tot <- sum(n1)
    
    ## make matrices for QP
    idxs0  <- trt  > gap + min(grps)
    n0 <- sum(idxs0)    
    const_mats <- make_constraint_mats(trt, grps, gap)
    Amat <- const_mats$Amat
    lvec <- const_mats$lvec
    uvec <- const_mats$uvec

    ## quadratic balance measures
    dvec <- lapply(1:J, function(j) c(Xc[[j]] %*%
                                      ## diag(nw[(d-ncol(Xc[[j]])+1):d] / n1tot) %*%
                                      x_t[[j]]) / ## n1[j]
                                    length(x_t[[j]])
                   )

    dvec_avg <- lapply(1:J,
                       function(j)  {
                           lapply(x_t,
                                  function(xtk) {
                                      
                                      dk <- length(xtk)
                                      dj <- ncol(Xc[[j]])
                                      ndim <- min(dk, dj)
                                      ## relative time: inner product from the end
                                      if(relative) {
                                          Xc[[j]][,(dj-ndim+1):dj] %*%
                                              xtk[(dk-ndim+1):dk]
                                           
                                      } else {
                                          ## absolute time: inner product from start
                                          Xc[[j]][,1:ndim] %*%
                                              xtk[1:ndim]
                                      }
                                  }) %>% reduce(`+`)
                       }) %>% reduce(c)
    
    dvec <- - (alpha * dvec_avg / d + (1 - alpha) * reduce(dvec,c))

    V1 <- do.call(Matrix::bdiag,
                  lapply(1:J, function(j) Xc[[j]] / sqrt(ncol(Xc[[j]]))))

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
                                  t(Xc[[j]][,(dj-ndim+1):dj]) / d
                          } else {
                              ## inner product from start
                              Xc[[k]][,1:ndim] %*%
                                  t(Xc[[j]][,1:ndim]) / d
                          }
                      }) %>% do.call(cbind, .)
           }) %>% do.call(rbind, .) -> V2


    Hmat <- alpha * V2 + (1 - alpha) * V1 %*% Matrix::t(V1) + lambda * Matrix::Diagonal(nrow(V1))

    ## Optimize
    settings <- osqp::osqpSettings(verbose = FALSE, eps_abs=1e-7, eps_rel = 1e-7,
                                   max_iter=5000)
    out <- osqp::solve_osqp(Hmat, dvec, Amat, lvec, uvec, pars=settings)

    ## get weights
    weights <- matrix(out$x, nrow=n0)

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
    weights[idxs0,] <- matrix(out$x, nrow=n0)
    weights <- t(t(weights) / n1
                 )

    output <- list(weights=weights,
                   imbalance=cbind(avg_imbal, imbalance))
    
    return(output)
    
}


#' Create constraint matrices for multisynth QP
#' @param trt Vector of treatment levels/times
#' @param grps Unique treatment time groups
#' @param gap Number of time periods after treatment to impute control values.
#'            For units treated at time T_j, all units treated after T_j + gap
#'            will be used as control values. If larger than the number of periods,
#'            only never never treated units (pure controls) will be used as comparison units#'
#' @return 
#'         \itemize{
#'          \item{"Amat"}{Linear constraint matrix}
#'          \item{"lvec"}{Lower bounds for linear constraints}
#'          \item{"uvec"}{Upper bounds for linear constraints}
#'         }
make_constraint_mats <- function(trt, grps, gap) {

    J <- length(grps)
    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))
    idxs0  <- trt  > gap + min(grps)

    n0 <- sum(idxs0)

    ## sum to n1 constraint
    A1 <- do.call(Matrix::bdiag, lapply(1:(J), function(x) rep(1, n0)))

    Amat <- as.matrix(Matrix::t(A1))
    Amat <- Matrix::rbind2(Matrix::t(A1), Matrix::Diagonal(nrow(A1)))

    ## zero out weights for inelligible units
    A3 <- do.call(Matrix::bdiag,
                  c(lapply(1:J,
                           function(j) {
                               z <- numeric(n0)
                               z[trt[idxs0] <= grps[j] + gap] <- 1
                               z
                           })))
    

    Amat <- Matrix::rbind2(Amat, Matrix::t(A3))

    
    lvec <- c(n1, # sum to n1 constraint
              numeric(nrow(A1)), # lower bound by zero
              numeric(J)) # zero out weights for impossible donors
    
    uvec <- c(n1, #sum to n1 constraint
              rep(Inf, nrow(A1)),
              numeric(J)) # zero out weights for impossible donors


    return(list(Amat=Amat, lvec=lvec, uvec=uvec))
}
