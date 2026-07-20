################################################################################
## Solvers for the synthetic control QP
##   min_w (X0' w - X1)' V (X0' w - X1)  s.t.  w >= 0, sum(w) = 1
## Each solver is a function with signature (X1, X0, V) -> numeric weights,
## where X1 is the t0 x 1 target, X0 is the n0 x t0 matrix of control
## outcomes, and V is a t0 x t0 scaling matrix. To add a new solver, add it
## to the switch in get_synth_solver() or pass the function itself as the
## `solver` argument.
################################################################################

#' Resolve a synth solver from a name or a function
#' @param solver "osqp" (default), "frank_wolfe", or a function with
#'        signature (X1, X0, V) returning a weight vector
#' @noRd
get_synth_solver <- function(solver = "osqp") {
    if (is.function(solver)) {
        return(solver)
    }
    switch(match.arg(solver, c("osqp", "frank_wolfe")),
           osqp = synth_qp_osqp,
           frank_wolfe = synth_qp_frank_wolfe)
}

#' Solve the synth QP with OSQP (forms the n0 x n0 Gram matrix)
#' @param X1 Target vector
#' @param X0 Matrix of control outcomes (n0 x t0)
#' @param V Scaling matrix
#' @noRd
synth_qp_osqp <- function(X1, X0, V) {

    Pmat <- X0 %*% V %*% t(X0)
    qvec <- - t(X1) %*% V %*% t(X0)

    n0 <- nrow(X0)
    A <- rbind(rep(1, n0), diag(n0))
    l <- c(1, numeric(n0))
    u <- c(1, rep(1, n0))

    settings = osqp::osqpSettings(verbose = FALSE,
                                  eps_rel = 1e-8,
                                  eps_abs = 1e-8)
    sol <- osqp::solve_osqp(P = Pmat, q = qvec,
                            A = A, l = l, u = u,
                            pars = settings)

    return(sol$x)
}

#' Solve the synth QP via Frank-Wolfe support identification
#'
#' Runs Frank-Wolfe with exact line search (the algorithm of the synthdid
#' package: Arkhangelsky, Athey, Hirshberg, Imbens & Wager 2021, dual
#' BSD-3/GPL >= 2) to identify the sparse active donor set, then solves the
#' QP exactly on that support (OSQP on a support-sized problem), and
#' verifies KKT optimality with a full gradient screen, re-admitting any
#' donors the support restriction missed. The n0 x n0 donor Gram matrix of
#' the full QP is never formed: Frank-Wolfe costs O(n0 t0) per iteration
#' and the polished QP is at most (5 t0) x (5 t0), so memory stays
#' O(n0 t0) -- the constraint that binds with very large donor pools.
#'
#' @param X1 Target vector
#' @param X0 Matrix of control outcomes (n0 x t0)
#' @param V Scaling matrix
#' @param max_iter Frank-Wolfe iteration cap
#' @param support_size Donors kept for the QP polish (top Frank-Wolfe
#'        weights); default min(n0, 5 t0)
#' @noRd
synth_qp_frank_wolfe <- function(X1, X0, V, max_iter = 2000L,
                                 support_size = NULL) {

    ## fold V into the design: A = V^{1/2} X0' (t0 x n0), b = V^{1/2} X1
    if (all(V == diag(diag(V), nrow(V)))) {
        s <- sqrt(diag(V))
        A <- t(X0) * s
        b <- as.vector(X1) * s
    } else {
        E <- eigen(V, symmetric = TRUE)
        R <- diag(sqrt(pmax(E$values, 0)), nrow(V)) %*% t(E$vectors)
        A <- R %*% t(X0)
        b <- as.vector(R %*% X1)
    }
    n0 <- ncol(A)
    t0 <- nrow(A)

    ## Frank-Wolfe with exact line search from uniform weights
    w <- rep(1 / n0, n0)
    Aw <- as.vector(A %*% w)
    val_old <- Inf
    for (it in seq_len(max_iter)) {
        grad <- crossprod(A, Aw - b)               # half-gradient, O(n0 t0)
        i <- which.min(grad)
        dA <- A[, i] - Aw
        denom <- sum(dA^2)
        if (denom <= 0) break
        step <- max(0, min(1, -sum((Aw - b) * dA) / denom))
        if (step <= 1e-12) break
        w <- (1 - step) * w
        w[i] <- w[i] + step
        Aw <- (1 - step) * Aw + step * A[, i]
        val <- sum((Aw - b)^2)
        if (val_old - val < 1e-9 * val) break      # relative-decrease stop
        val_old <- val
    }

    ## exact QP restricted to a donor support (OSQP on a small problem)
    qp_support <- function(As, bs) {
        ns <- ncol(As)
        sol <- osqp::solve_osqp(
            P = crossprod(As), q = -crossprod(As, bs),
            A = rbind(rep(1, ns), diag(ns)),
            l = c(1, numeric(ns)), u = c(1, rep(1, ns)),
            pars = osqp::osqpSettings(verbose = FALSE,
                                      eps_rel = 1e-8, eps_abs = 1e-8))
        ws <- pmax(sol$x, 0)
        ws / sum(ws)
    }

    k <- if (is.null(support_size)) min(n0, 5L * t0)
         else min(n0, support_size)
    ## FW starts uniform and never zeroes untouched donors: its support is
    ## the set pushed above the shrunken-uniform baseline
    support <- which(w > min(w) * (1 + 1e-9))
    if (!length(support)) support <- seq_len(min(n0, t0))
    if (length(support) > k)
        support <- support[order(w[support], decreasing = TRUE)[seq_len(k)]]

    ## polish, then KKT gradient screening: at the optimum every donor
    ## outside the active set has (half-)gradient >= the common active-set
    ## level mu; one O(n0 t0) pass finds donors FW missed, add them and
    ## re-polish until no violations (typically 1-2 rounds)
    for (round in seq_len(10L)) {
        w <- numeric(n0)
        w[support] <- qp_support(A[, support, drop = FALSE], b)
        g <- as.vector(crossprod(A, A %*% w - b))
        active <- support[w[support] > 1e-12]
        mu <- stats::median(g[active])
        tol_g <- 1e-8 * max(abs(g))
        violated <- setdiff(which(g < mu - tol_g), support)
        if (!length(violated)) break
        violated <- violated[order(g[violated])[seq_len(min(length(violated),
                                                            k))]]
        support <- c(support, violated)
    }

    return(w)
}
