#######################################################
# Helper scripts to fit synthetic controls to simulations
#######################################################

#' Fit synthetic controls on outcomes with data in 'wide' format
#' @param wide_data List of outcomes in wide format
#' @param V Matrix to scale the objective by
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Synth weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#' }
fit_synth_wide <- function(wide_data, V = NULL) {

  # get pre treatment outcomes into the right form
  X1 <- colMeans(wide_data$X[wide_data$trt == 1, , drop = F])
  X0 <- wide_data$X[wide_data$trt == 0, , drop = F]

  return(fit_synth_formatted(X1, X0, V))
}

#' Fit synthetic controls on outcomes after formatting data
#' @param X1 Vector of pre-treatment outcomes for treated unit
#' @param X0 Matrix of pre-treatment outcomes for control units
#' @param V Matrix to scale the obejctive by
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Synth weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#' }
fit_synth_formatted <- function(X1, X0, V = NULL) {


    t0 <- dim(X0)[2]
    # if no V, set equal to 1
    if(is.null(V)) {
        V <- diag(rep(1, t0))
    } else if(is.vector(V)) {
        V <- diag(V)
    } else if(ncol(V) == 1 & nrow(V) == t0) {
        V <- diag(c(V))
    } else if(ncol(V) == t0 & nrow(V) == 1) {
        V <- diag(c(V))
    } else if(nrow(V) == t0) {
    } else {
        stop("`V` must be a vector with t0 elements or a t0xt0 matrix")
    }

    # fit synth
    weights <- synth_qp(X1, X0, V)
    # compute the l2 imbalance
    l2_imbalance <- sqrt(sum((t(X0) %*% weights - X1) ^ 2))

    # L2 imbalance scaled by L2 imbalance for uniform weights
    uni_w <- matrix(1 / nrow(X0), nrow = nrow(X0), ncol = 1)
    unif_l2_imbalance <- sqrt(sum((t(X0) %*% uni_w - X1) ^ 2))
    scaled_l2_imbalance <- l2_imbalance / unif_l2_imbalance

    return(list(weights = weights,
                l2_imbalance = l2_imbalance,
                scaled_l2_imbalance = scaled_l2_imbalance))
}

#' Solve the synth QP directly
#' @param X1 Target vector
#' @param X0 Matrix of control outcomes
#' @param V Scaling matrix
#' @noRd
synth_qp <- function(X1, X0, V) {

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
