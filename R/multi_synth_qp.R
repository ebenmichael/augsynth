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
    Xc <- X[!is.finite(trt),,drop=FALSE]
    Xc <- lapply(1:nrow(mask),
                 function(j) Xc[,mask[j,]==1, drop=F])

    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))


    ## compute global average according to time relative to treatment
    sum_t <- reduce(x_t, function(x, y) c(numeric(length(y)-length(x)),x) + y)
    
    ## weight the global parameters by the number of treated units still untreated in calander time
    nw <- sapply(1:max(trt[is.finite(trt)]), function(t) sum(trt[is.finite(trt)] >= t))
    ## reverse to get relative absolute time
    nw <- rev(nw)    


    
    ## make matrices for QP
    n0 <- sum(!is.finite(trt))

    ## sum to one constraint
    A1 <- do.call(Matrix::bdiag, lapply(1:(J+1), function(x) rep(1, n0)))
    Amat <- as.matrix(Matrix::t(A1))
    Amat <- Matrix::rbind2(Matrix::t(A1), Matrix::Diagonal(nrow(A1)))
    
    ## sum to global weights
    A2 <- Matrix::kronecker(t(c(-1, rep(1, J))), Matrix::Diagonal(n0))
    Amat <- Matrix::rbind2(Amat, A2)

    
    lvec <- c(sum(n1), n1, numeric(nrow(A1)), numeric(n0))
    uvec <- c(sum(n1), n1, rep(sum(n1), n0), sapply(n1, function(n1j) rep(n1j, n0)), numeric(n0))

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

    print(dvec_avg)

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
    print(dvec_avg)    
    
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


    ## if(alpha == 1) {
    ##     dvec <- -dvec_avg
    ##     Hmat <- V2
    ##     ## sum to one constraint
    ##     A1 <- t(rep(1, n0))
    ##     Amat <- Matrix::rbind2(A1, Matrix::Diagonal(ncol(A1)))
    
    ##     lvec <- c(sum(n1), numeric(n0))
    ##     uvec <- c(sum(n1), rep(sum(n1), n0))
    ## }
    print(length(dvec))
    print(dim(Hmat))
    print(dim(Amat))
    print(length(lvec))
    print(length(uvec))

    ## print(dvec)
    ## print(Hmat)
    ## dvec <- lapply(1:J, function(j) c((1-alpha) / length(x_t[[j]]) * Xc[[j]] %*% x_t[[j]]) +
    ##                                 c(alpha * n1[j] / sum(n1) * Xc[[j]] %*%
    ##                                   ## diag(nw[rev(mask[j,]==1)] / sum(n1)) %*%
    ##                                   sum_t[mask[j,]==1]))
    ## dvec <- - reduce(dvec, c)



    ## V1 <- do.call(Matrix::bdiag,
    ##               lapply(Xc, function(x) x / sqrt(ncol(x))))

    ## Xc %>%
    ##     lapply(function(x) cbind(matrix(0, nrow=nrow(x), ncol=(d-ncol(x))), x)) %>%
    ##     do.call(rbind, .) -> V2

    ## Xc %>%
    ##     lapply(function(x) x[,1:length(x_t[[1]])]) %>%
    ##     do.call(rbind, .) -> V2

    
    ## just do a sparse matrix if no global constraint
    ## if(alpha == 0) {
    ##     Hmat <- V1 %*% Matrix::t(V1) + lambda * Matrix::Diagonal(nrow(V1))
    ## } else {
    ##     Hmat <- V1 %*% Matrix::t(V1) + alpha * V2 %*% ## (diag(nw) / sum(n1)) %*% 
    ##         t(V2) / sum(n1)^2 + lambda * diag(nrow(V1))
    ## }

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





#' Internal function to fit synth with staggered adoption according to absolute time
#' #' @importMethodsFrom Matrix %*%
multisynth_absolute_qp <- function(X, trt, mask, gap, alpha, lambda, ...) {

    n <- dim(X)[1]

    d <- dim(X)[2]
    
    ## average over treatment times/groups
    grps <- sort(unique(trt[is.finite(trt)]))

    J <- length(grps)


    x_t <- lapply(1:J, function(j) colSums(X[trt ==grps[j], mask[j,]==1, drop=F]))    
    
    ## pure controls, indexed by time
    Xc <- X[!is.finite(trt),,drop=FALSE]
    Xc <- lapply(1:nrow(mask),
                 function(j) Xc[,mask[j,]==1, drop=F])

    n1 <- sapply(1:J, function(j) sum(trt==grps[j]))


    ## compute global average according to time relative to treatment
    sum_t <- reduce(x_t, function(x, y) c(x, numeric(length(y)-length(x))) + y)
    
    ## weight the global parameters by the number of treated units still untreated in calander time
    nw <- sapply(1:max(trt[is.finite(trt)]), function(t) sum(trt[is.finite(trt)] >= t))
    ## reverse to get relative absolute time
    nw <- rev(nw)    


    
    ## make matrices for QP
    n0 <- sum(!is.finite(trt))

    ## sum to one constraint
    A1 <- do.call(Matrix::bdiag, lapply(1:(J), function(x) rep(1, n0)))
    Amat <- as.matrix(Matrix::t(A1))
    Amat <- Matrix::rbind2(Matrix::t(A1), Matrix::Diagonal(nrow(A1)))


    lvec <- c(n1, numeric(nrow(A1)))
    uvec <- c(n1, sapply(n1, function(n1j) rep(n1j, n0)))

    ## quadratic balance measures
    dvec <- lapply(1:J, function(j) c((1-alpha) / length(x_t[[j]]) * Xc[[j]] %*% x_t[[j]]) +
                                    c(alpha * n1[j] / sum(n1) * Xc[[j]] %*%
                                      ## diag(nw[rev(mask[j,]==1)] / sum(n1)) %*%
                                      sum_t[mask[j,]==1]))
    dvec <- - reduce(dvec, c)
    


    V1 <- do.call(Matrix::bdiag,
                  lapply(Xc, function(x) x / sqrt(ncol(x))))

    Xc %>%
        lapply(function(x) cbind(x, matrix(0, nrow=nrow(x), ncol=(d-ncol(x))))) %>%
        do.call(rbind, .) -> V2

    
    ## just do a sparse matrix if no global constraint
    if(alpha == 0) {
        Hmat <- V1 %*% Matrix::t(V1) + lambda * Matrix::Diagonal(nrow(V1))
    } else {
        Hmat <- V1 %*% Matrix::t(V1) + alpha * V2 %*% ## (diag(nw) / sum(n1)) %*% 
            t(V2) / sum(n1)^2 + lambda * diag(nrow(V1))
    }

    settings <- osqp::osqpSettings(verbose = FALSE, eps_abs=1e-8, eps_rel = 1e-8)
    out <- osqp::solve_osqp(Hmat, dvec, Amat, lvec, uvec, pars=settings)


    ## get weights
    weights <- matrix(out$x, nrow=n0)

    ## compute imbalance
    imbalance <- vapply(1:J,
                        function(j) c(x_t[[j]] -  t(Xc[[j]]) %*% weights[,j],
                                      numeric(d-length(x_t[[j]]))),
                        numeric(d))

    avg_imbal <- rowSums(t(t(imbalance)))

    obj <- sapply(1:J,
                  function(j) sum((x_t[[j]] -  t(Xc[[j]]) %*% weights[,j])^2) / length(x_t[[j]]))
    obj <- (1-alpha) * reduce(obj, `+`) + alpha * avg_imbal %*% diag(nw/sum(nw)) %*% avg_imbal
    
    output <- list(weights=weights,
                   imbalance=cbind(avg_imbal, imbalance),
                   objective=obj)
    
    return(output)
    
}
