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
#' 
#' @return \itemize{
#'          \item{"weights"}{Matrix of unit weights}
#'          \item{"imbalance"}{Matrix of overall and group specific imbalance}
#'          \item{"global_l2"}{Imbalance overall}
#'          \item{"ind_l2"}{Matrix of imbalance for each group}
#'         }
multisynth_qp <- function(X, trt, mask, n_leads=NULL, n_lags=NULL,
                          relative=T, nu=0, lambda=0, ...) {

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
    grps <- trt[is.finite(trt)]

    J <- length(grps)


    n_t <- sum(is.finite(trt))
    which_t <- (1:n)[is.finite(trt)]

    ## handle X differently if it is a list
    if(typeof(X) == "list") {
        ## x_t <- lapply(1:J, function(j) colSums(X[[j]][trt ==grps[j], mask[j,]==1, drop=F]))    
        x_t <- lapply(1:J, function(j) colSums(X[[j]][which_t[j], mask[j,]==1, drop=F]))
        ## All possible donor units for all treatment groups
        Xc <- lapply(1:nrow(mask),
                 function(j) X[[j]][trt > n_leads + min(grps), mask[j,]==1, drop=F])
    } else {
        ## x_t <- lapply(1:J, function(j) colSums(X[trt ==grps[j], mask[j,]==1, drop=F]))    
        x_t <- lapply(1:J, function(j) colSums(X[which_t[j], mask[j,]==1, drop=F]))        
        ## All possible donor units for all treatment groups
        ## Xc <- lapply(1:nrow(mask),
        ##              function(j) X[trt > n_leads + min(grps), mask[j,]==1, drop=F])
        Xc <- lapply(1:nrow(mask),
                 function(j) X[trt > n_leads + min(grps), mask[j,]==1, drop=F])        
    }
    
    n1 <- sapply(1:J, function(j) 1)
    
    n1tot <- sum(n1)
    
    ## make matrices for QP
    idxs0  <- trt  > n_leads + min(grps)
    n0 <- sum(idxs0)    
    const_mats <- make_constraint_mats(trt, grps, n_leads)
    Amat <- const_mats$Amat
    lvec <- const_mats$lvec
    uvec <- const_mats$uvec

    ## quadratic balance measures
    dvec <- lapply(1:J,
                   function(j) {
                       dj <- length(x_t[[j]])
                       ndim <- min(dj, n_lags)
                       c(Xc[[j]][,(dj-ndim+1):dj, drop=F] %*%
                         x_t[[j]][(dj-ndim+1):dj]) /
                           ndim
                   }
                   )

    dvec_avg <- lapply(1:J,
                       function(j)  {
                           lapply(x_t,
                                  function(xtk) {
                                      
                                      dk <- length(xtk)
                                      dj <- ncol(Xc[[j]])
                                      ndim <- min(dk, dj, n_lags)
                                      ## relative time: inner product from the end
                                      if(relative) {
                                          Xc[[j]][,(dj-ndim+1):dj, drop=F] %*%
                                              xtk[(dk-ndim+1):dk]
                                           
                                      } else {
                                          ## absolute time: inner product from start
                                          Xc[[j]][,1:ndim] %*%
                                              xtk[1:ndim]
                                      }
                                  }) %>% reduce(`+`)
                       }) %>% reduce(c)
    
    dvec <- - (nu * dvec_avg / n_lags + (1 - nu) * reduce(dvec,c))

    V1 <- do.call(Matrix::bdiag,
                  lapply(1:J,
                         function(j) {
                             dj <- length(x_t[[j]])
                             ndim <- min(dj, n_lags)
                             Xc[[j]][,(dj-ndim+1):dj, drop=F] / sqrt(ndim)
                         })
                  )
                         

    lapply(1:J,
           function(k) {
               lapply(1:J,
                      function(j) {
                          dk <- ncol(Xc[[k]])
                          dj <- ncol(Xc[[j]])
                          ndim <- min(dk, dj, n_lags)
                          if(relative) {
                              ## inner product from end
                              Xc[[k]][,(dk-ndim+1):dk, drop=F] %*%
                                  t(Xc[[j]][,(dj-ndim+1):dj, drop=F]) / n_lags
                          } else {
                              ## inner product from start
                              Xc[[k]][,1:ndim] %*%
                                  t(Xc[[j]][,1:ndim]) / n_lags
                          }
                      }) %>% do.call(cbind, .)
           }) %>% do.call(rbind, .) -> V2


    Hmat <- nu * V2 + (1 - nu) * V1 %*% Matrix::t(V1) + lambda * Matrix::Diagonal(nrow(V1))

    ## Optimize
    settings <- osqp::osqpSettings(verbose = FALSE, ...)
    out <- osqp::solve_osqp(Hmat, dvec, Amat, lvec, uvec, pars=settings)

    ## get weights
    weights <- matrix(out$x, nrow=n0)

    ## compute imbalance
    if(relative) {
        imbalance <- vapply(1:J,
                            function(j) {
                                dj <- length(x_t[[j]])
                                ndim <- min(dj, n_lags)
                                c(numeric(d-ndim),
                                x_t[[j]][(dj-ndim+1):dj] -
                                    t(Xc[[j]][,(dj-ndim+1):dj]) %*% weights[,j])
                            },
                            numeric(d))
    } else {
        imbalance <- vapply(1:J,
                            function(j) c(x_t[[j]] -  t(Xc[[j]]) %*% weights[,j],
                                          numeric(d-length(x_t[[j]]))),
                            numeric(d))
    }

    avg_imbal <- rowSums(t(t(imbalance)))


    global_l2 <- sqrt(sum(avg_imbal^2))
    ind_l2 <- sum(apply(imbalance, 2, function(x) sqrt(sum(x^2))))
    ## pad weights with zeros for treated units and divide by number of treated units
    weights <- matrix(0, nrow=n, ncol=J)
    weights[idxs0,] <- matrix(out$x, nrow=n0)
    weights <- t(t(weights) / n1
                 )

    
    output <- list(weights=weights,
                   imbalance=cbind(avg_imbal, imbalance),
                   global_l2=global_l2,
                   ind_l2=ind_l2)
    
    return(output)
    
}


#' Create constraint matrices for multisynth QP
#' @param trt Vector of treatment levels/times
#' @param grps Treatment times
#' @param n_leads Number of time periods after treatment to impute control values.
#'            For units treated at time T_j, all units treated after T_j + n_leads
#'            will be used as control values. If larger than the number of periods,
#'            only never never treated units (pure controls) will be used as comparison units#'
#' @return 
#'         \itemize{
#'          \item{"Amat"}{Linear constraint matrix}
#'          \item{"lvec"}{Lower bounds for linear constraints}
#'          \item{"uvec"}{Upper bounds for linear constraints}
#'         }
make_constraint_mats <- function(trt, grps, n_leads) {

    J <- length(grps)
    n1 <- sapply(1:J, function(j) 1)
    idxs0  <- trt  > n_leads + min(grps)

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
                               z[trt[idxs0] <= grps[j] + n_leads] <- 1
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
