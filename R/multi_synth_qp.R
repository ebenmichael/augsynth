################################################################################
## Solve the multisynth problem as a QP
################################################################################


#' Internal function to fit synth with staggered adoption according to time relative to treatment
#' #' @importMethodsFrom Matrix %*%
multisynth_relative_qp <- function(X, trt, mask, gap, alpha, lambda, ...) {

    n <- dim(X)[1]

    d <- dim(X)[2]
    
    ## average over treatment times/groups
    grps <- sort(unique(trt[is.finite(trt)]))

    J <- length(grps)


    x_t <- lapply(1:J, function(j) colMeans(X[trt ==grps[j], mask[j,]==1, drop=F]))
    x_t <- lapply(1:J, function(j) colSums(X[trt ==grps[j], mask[j,]==1, drop=F]))    
    
    ## pure controls, indexed by time
    ## Xc <- X[!is.finite(trt),,drop=FALSE]
    ## All possible donor units for all treatment groups
    Xc <- lapply(1:nrow(mask),
                 function(j) X[trt > gap + min(grps), mask[j,]==1, drop=F])

    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))


    ## compute global average according to time relative to treatment
    sum_t <- reduce(x_t, function(x, y) c(numeric(length(y)-length(x)),x) + y)
    
    ## weight the global parameters by the number of treated units still untreated in calander time
    nw <- sapply(1:max(trt[is.finite(trt)]), function(t) sum(trt[is.finite(trt)] >= t))
    ## reverse to get relative absolute time
    nw <- rev(nw)    


    
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
    glblzero <- numeric(n0)
    glblzero[is.finite(trt)] <- 1
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
    dvec_avg <- lapply(1:J,
                       function(j) {
                           lapply(1:J,
                                  function(k) {
                                      c(c(numeric(d-length(x_t[[k]])),x_t[[k]]) %*%
                                        t(cbind(matrix(0, nrow=nrow(Xc[[j]]), ncol=(d-ncol(Xc[[j]]))), Xc[[j]])))
                                  }) %>% reduce(`+`)
                       }) %>% reduce(`+`)
    
    dvec_avg <- Xc[[1]] %*% reduce(x_t, function(x, y) x[1:length(x_t[[1]])] + y[1:length(x_t[[1]])])


    dvec_avg <- lapply(x_t,
                       function(xtk) {
                           lapply(1:J,
                                  function(j)  {
                                      ## Xc[[j]][,1:length(x_t[[1]])] %*% xtk[1:length(x_t[[1]])] * n1[j] / sum(n1)
                                      dk <- length(xtk)
                                      dj <- ncol(Xc[[j]])
                                      ndim <- min(dk, dj)
                                      Xc[[j]][,(dj-ndim+1):dj] %*% xtk[(dk-ndim+1):dk] * n1[j] / sum(n1)
                                  }) %>% reduce(`+`)
                       }) %>% reduce(`+`)
    
    dvec <- -c(alpha * dvec_avg, reduce(dvec,c))


    V1 <- do.call(Matrix::bdiag,
                  lapply(Xc, function(x) x / sqrt(ncol(x))))

    lapply(Xc,
           function(x0k) {
               lapply(Xc,
                      function(x0j) {
                          cbind(matrix(0, nrow=nrow(x0k), ncol=(d-ncol(x0k))), x0k) %*%
                              t(cbind(matrix(0, nrow=nrow(x0j), ncol=(d-ncol(x0j))), x0j))
                      }) %>% reduce(`+`)
           }) %>% reduce(`+`) -> V2
    
    V2 <- Xc[[1]] %*% t(Xc[[1]])


    lapply(1:J,
           function(k) {
               lapply(1:J,
                      function(j) {
                          ## Xc[[k]][,1:length(x_t[[1]])] %*%
                          ##     t(Xc[[j]][,1:length(x_t[[1]])]) * n1[k] * n1[j] / sum(n1)^2
                          dk <- ncol(Xc[[k]])
                          dj <- ncol(Xc[[j]])
                          ndim <- min(dk, dj)
                          Xc[[k]][,(dk-ndim+1):dk] %*%
                              t(Xc[[j]][,(dj-ndim+1):dj]) * n1[k] * n1[j] / sum(n1)^2
                      }) %>% reduce(`+`)
           }) %>% reduce(`+`) -> V2
    
    Hmat <- Matrix::bdiag(list(alpha * V2, (1-alpha) * V1 %*% Matrix::t(V1))) +
        lambda * Matrix::Diagonal(nrow(V1) + nrow(V2))


    print(length(dvec))
    print(dim(Hmat))
    print(dim(Amat))
    print(length(lvec))
    print(length(uvec))

    settings <- osqp::osqpSettings(verbose = TRUE, eps_abs=1e-8, eps_rel = 1e-8)
    out <- osqp::solve_osqp(Hmat, dvec, Amat, lvec, uvec, pars=settings)

    ## get weights
    weights <- matrix(out$x[-(1:n0)], nrow=n0)

    ## compute imbalance
    imbalance <- vapply(1:J,
                         function(j) c(numeric(d-length(x_t[[j]])),
                                       x_t[[j]] -  t(Xc[[j]]) %*% weights[,j]),
                        numeric(d))

    avg_imbal <- rowSums(t(t(imbalance)))
    
    output <- list(weights=weights,
                   imbalance=cbind(avg_imbal, imbalance))
    
    return(output)
    
}
