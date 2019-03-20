
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "augsynth_types.h"
#include <numeric>
#include <string>

using namespace arma;
using namespace Rcpp;
using namespace std;


//' Gradient for Multiple synthetic controls and absolute time
// [[Rcpp::export]]
arma::mat grad_multisynth_absolute(arma::mat theta, List opts) {


  // treated averages, re-indexed by time relative to treatment
  arma::mat Xt = as<arma::mat>(opts["Xt"]);

  // all data in absolute time
  arma::mat X = as<arma::mat>(opts["X"]);

  // treatment time vector
  arma::vec trt = as<arma::vec>(opts["trt"]);

  // unqiue treatment time vector
  arma::vec unique_trt =  as<arma::vec>(opts["unique_trt"]);

  // number of units at each level
  arma::vec n1 = as<arma::vec>(opts["n1"]);
  int ntot = accu(n1);
  // mask for time periods for treated units that are pre-treatment
  arma::mat bool_mask = as<arma::mat>(opts["mask"]);

  // length of gap for estimating treatment effects
  int gap = as<int>(opts["gap"]);
  
  weightPtr weight_func = *as<wptr>(opts["weight_func"]);
  // arma::mat ipw_weights = as<arma::mat>(opts["ipw_weights"]);

  // initialize gradient as zero
  arma::mat grad = zeros(bool_mask.n_cols, bool_mask.n_rows+1);

  int d = X.n_cols;
  
  // iterate over treated units
  arma::mat restrict_Xc;
  arma::mat weights;
  for(int j=0; j < bool_mask.n_rows; j++) {

    int tj = unique_trt(j);
    // restrict to units that are treated at time > T_j + gap, apply mask
    restrict_Xc = X.rows(find(trt > tj + gap));
    restrict_Xc = restrict_Xc.each_row() % bool_mask.row(j);

    weights = weight_func(restrict_Xc, theta.col(0) + theta.col(j+1)
                          // ipw_weights
                          );
    grad.col(j+1) = restrict_Xc.t() * weights;
    grad.col(0) += restrict_Xc.t() * weights * n1[j]/ntot;
  }

  //combine to get gradient
  return grad - Xt;    
}


// [[Rcpp::export]]
gptr make_grad_multisynth_absolute() {
  return gptr(new gradPtr(grad_multisynth_absolute));
}



//' Gradient for Multiple synthetic controls and relative time
// [[Rcpp::export]]
arma::mat grad_multisynth_relative(arma::mat theta, List opts) {

  // treated averages, re-indexed by time relative to treatment
  arma::mat Xt = as<arma::mat>(opts["Xt"]);

  // all data in absolute time
  arma::mat X = as<arma::mat>(opts["X"]);

  // treatment time vector
  arma::vec trt = as<arma::vec>(opts["trt"]);

  // unqiue treatment time vector
  arma::vec unique_trt =  as<arma::vec>(opts["unique_trt"]);

  // number of units at each level
  arma::vec n1 = as<arma::vec>(opts["n1"]);
  int ntot = accu(n1);
  // mask for time periods for treated units that are pre-treatment
  arma::mat bool_mask = as<arma::mat>(opts["mask"]);

  // length of gap for estimating treatment effects
  int gap = as<int>(opts["gap"]);
  
  weightPtr weight_func = *as<wptr>(opts["weight_func"]);
  // arma::mat ipw_weights = as<arma::mat>(opts["ipw_weights"]);

  // initialize gradient as zero
  arma::mat grad = zeros(bool_mask.n_cols, bool_mask.n_rows+1);

  int d = X.n_cols;
  
  // iterate over treated units
  arma::mat restrict_Xc;
  arma::mat weights;
  for(int j=0; j < bool_mask.n_rows; j++) {

    int tj = unique_trt(j);
    // restrict to units that are treated at time > T_j + gap, apply mask
    restrict_Xc = X.rows(find(trt > tj + gap));
    restrict_Xc = shift(restrict_Xc.each_row() %
                        bool_mask.row(j),
                        d - tj,
                        1);

    weights = weight_func(restrict_Xc, theta.col(0) + theta.col(j+1)
                          // ipw_weights
                          );
    grad.col(j+1) = restrict_Xc.t() * weights;
    grad.col(0) += restrict_Xc.t() * weights * n1[j]/ntot;
  }

  //combine to get gradient
  return grad - Xt;
    
}


// [[Rcpp::export]]
gptr make_grad_multisynth_relative() {
  return gptr(new gradPtr(grad_multisynth_relative));
}





//' Nuclear norm prox
//'
//' @param x Input matrix
//' @param lam Prox scaling factor
//' @param opts List of options (opts["lam"] holds the other scaling
//'
//' @return Singular value soft thresholded X
// [[Rcpp::export]]
mat prox_nuc(mat x, double lam, List opts) {
  mat U; vec s; mat V;
  // SVD then threshold singular values
  bool fail = svd_econ(U, s, V, x);
  int d = x.n_rows;
  int m = x.n_cols;

  // threshold singular values
  lam = lam * as<double>(opts["lam"]);
  s = (s - lam) % (s > lam) + (s + lam) % (s < -lam);

  // If the SVD fails, return a matrix of NA
  // TODO: Figure out if this has weird side effects
  if(1-fail) {
    mat out(x.n_rows, x.n_cols);
    out.fill(datum::nan);
    return out;
  } else {
    return U * diagmat(s) * V.t();
  }
}



//' Weighted squared L2 prox
//'
//' @param x Input matrix
//' @param lam Prox scaling factor
//' @param opts List of options. opts["lam"] holds the other scaling, opts["w"] holds unit specific weights
//'
//' @return Column soft thresholded X
// [[Rcpp::export]]
mat prox_weighted_l2_sq(mat x, double lam, List opts) {
  lam = lam * as<double>(opts["lam"]);
  vec ridge_w = as<vec>(opts["w"]) * lam;
  return x / (1 + ridge_w);
}



//' Prox for weighted ridge on global params and nuclear norm on individual params
//'
//' @param x Input matrix (contains global and local parameters
//' @param lam Prox scaling factor
//' @param opts List of options (opts["alpha"] holds the ratio between global and local balance
//'
//' @return L2 squared prox values
// [[Rcpp::export]]
mat prox_multilevel_weighted_ridge_nuc(mat x, double lam, List opts) {

  double alpha = as<double>(opts["alpha"]);  
  
  // separate out global and local parameters
  mat xglobal = prox_weighted_l2_sq(x.col(0), lam * (1-alpha), opts);
  mat xlocal = prox_nuc(x.cols(1, x.n_cols-1), lam * alpha, opts);

  
  return join_rows(xglobal, xlocal);
  
}

// [[Rcpp::export]]
pptr make_prox_multilevel_weighted_ridge_nuc() {
  return pptr(new proxPtr(prox_multilevel_weighted_ridge_nuc));
}


//' Prox for weighted ridge on global params and nuclear norm on individual params
//'
//' @param x Input matrix (contains global and local parameters
//' @param lam Prox scaling factor
//' @param opts List of options (opts["alpha"] holds the ratio between global and local balance
//'
//' @return L2 squared prox values
// [[Rcpp::export]]
mat prox_multilevel_weighted_ridge_nuc_normalized(mat x, double lam, List opts) {

  // separate out the intercept
  int d = x.n_rows;

  mat alpha = x.row(0);
  mat beta = x.rows(1,d-1);

  // prox on beta
  beta = prox_multilevel_weighted_ridge_nuc(beta, lam, opts);

  return join_vert(alpha, beta);  
}

// [[Rcpp::export]]
pptr make_prox_multilevel_weighted_ridge_nuc_normalized() {
  return pptr(new proxPtr(prox_multilevel_weighted_ridge_nuc_normalized));
}



//' Nuclear norm projection ||x||_* <= lam
//'
//' @param x Input matrix
//' @param lam Constraint on nuclear norm
//' @param opts List of options (opts["lam"] holds the other scaling
//'
//' @return Singular value soft thresholded X
// [[Rcpp::export]]
mat proj_nuc(mat x, double lam, List opts) {
  mat U; vec s; mat V;
  // SVD then threshold singular values
  bool fail = svd_econ(U, s, V, x);
  int d = x.n_rows;
  int m = x.n_cols;

  lam = lam * as<double>(opts["lam"]);
  // if the constraint is already satisfied, just return the matrix
  if(accu(s) <= lam) {
    return x;
  }

  
  // sort singular values to find optimal threshold
  // Rcout << s << "\n";
  vec cumsum_s = cumsum(s);
    // Rcout << cumsum_s << "\n";
  int rho = 0;
  while(s(rho + 1) - (cumsum_s(rho + 1) - lam) / (rho + 2) > 0 & rho < s.n_elem - 2) {
    // Rcout << s(rho + 1) - (cumsum_s(rho + 1) - lam) / (rho + 2)  << "\n";
    rho +=  1;
  }
  // Rcout << rho << "\n";
  double theta = (cumsum_s(rho) - lam) / (rho + 1);

  // threshold singular values
  s = (s - theta) % (s > theta) + (s + theta) % (s < -theta);
  // Rcout << "Constraint: " << lam << "Actual: " << accu(s)  << "\n\n---------\n";
  // If the SVD fails, return a matrix of NA
  // TODO: Figure out if this has weird side effects
  if(1-fail) {
    mat out(x.n_rows, x.n_cols);
    out.fill(datum::nan);
    return out;
  } else {
    return U * diagmat(s) * V.t();
  }
}

// [[Rcpp::export]]
pptr make_proj_nuc() {
  return pptr(new proxPtr(proj_nuc));
}



//' Squared L2 Prox for global parameters, nuclear norm projection for local parameters
//'
//' @param x Input matrix (contains global and local parameters
//' @param lam Prox scaling factor
//' @param opts List of options (opts["alpha"] holds the ratio between global and local balance
//'
//' @return L2 squared prox values
// [[Rcpp::export]]
mat proj_multilevel_weighted_ridge_nuc(mat x, double lam, List opts) {

  double alpha = as<double>(opts["alpha"]);  
  // separate out global and local parameters
  //NOTE: hyperparameter is inverted
  List inv_opts = List::create(_["lam"] = 1/as<double>(opts["lam"]),
                               _["w"]=as<vec>(opts["w"]));
  mat xglobal = prox_weighted_l2_sq(x.col(0), lam * 1/(1-alpha), inv_opts); 
  mat xlocal = proj_nuc(x.cols(1, x.n_cols-1), lam * alpha, opts);

  
  return join_rows(xglobal, xlocal);
  
}

// [[Rcpp::export]]
pptr make_proj_multilevel_weighted_ridge_nuc() {
  return pptr(new proxPtr(proj_multilevel_weighted_ridge_nuc));
}

