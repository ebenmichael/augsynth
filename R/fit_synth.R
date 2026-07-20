#######################################################
# Helper scripts to fit synthetic controls to simulations
#######################################################

#' Make a V matrix from a vector (or null)
#'
#' @noRd
make_V_matrix <- function(t0, V) {
  if(is.null(V)) {
        V <- diag(rep(1, t0))
    } else if(is.vector(V)) {
        if(length(V) != t0) {
          stop(paste("`V` must be a vector with", t0, "elements or a", t0,
                     "x", t0, "matrix"))
        }
        V <- diag(V)
    } else if(ncol(V) == 1 & nrow(V) == t0) {
        V <- diag(c(V))
    } else if(ncol(V) == t0 & nrow(V) == 1) {
        V <- diag(c(V))
    } else if(nrow(V) == t0) {
    } else {
        stop(paste("`V` must be a vector with", t0, "elements or a", t0,
                     "x", t0, "matrix"))
    }

  return(V)
}

#' Fit synthetic controls on outcomes after formatting data
#' @param synth_data Panel data in format of Synth::dataprep
#' @param V Matrix to scale the obejctive by
#' @param solver Solver for the synth QP: "osqp" (default), "frank_wolfe",
#'        or a function (X1, X0, V) -> weights; see R/solvers.R
#' @noRd
#' @return \itemize{
#'          \item{"weights"}{Synth weights}
#'          \item{"l2_imbalance"}{Imbalance in pre-period outcomes, measured by the L2 norm}
#'          \item{"scaled_l2_imbalance"}{L2 imbalance scaled by L2 imbalance of uniform weights}
#' }
fit_synth_formatted <- function(synth_data, V = NULL, solver = "osqp") {


    t0 <- dim(synth_data$Z0)[1]
    ## if no  is supplied, set equal to 1

    V <- make_V_matrix(t0, V)

    weights <- synth_qp(synth_data$X1, t(synth_data$X0), V, solver = solver)
    l2_imbalance <- sqrt(sum((synth_data$Z0 %*% weights - synth_data$Z1)^2))

    ## primal objective value scaled by least squares difference for mean
    uni_w <- matrix(1/ncol(synth_data$Z0), nrow=ncol(synth_data$Z0), ncol=1)
    unif_l2_imbalance <- sqrt(sum((synth_data$Z0 %*% uni_w - synth_data$Z1)^2))
    scaled_l2_imbalance <- l2_imbalance / unif_l2_imbalance

    return(list(weights=weights,
                l2_imbalance=l2_imbalance,
                scaled_l2_imbalance=scaled_l2_imbalance))
}

#' Solve the synth QP with the requested solver
#' @param X1 Target vector
#' @param X0 Matrix of control outcomes
#' @param V Scaling matrix
#' @param solver Solver name or function; resolved by get_synth_solver()
#' @noRd
synth_qp <- function(X1, X0, V, solver = "osqp") {

    solve_fn <- get_synth_solver(solver)
    return(solve_fn(X1, X0, V))
}
