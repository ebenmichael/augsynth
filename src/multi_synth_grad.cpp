
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

  // control and treated data
  arma::mat Xc = as<arma::mat>(opts["Xc"]);
  arma::mat Xt = as<arma::mat>(opts["Xt"]);

  // number of units at each level
  arma::vec n1 = as<arma::vec>(opts["n1"]);
  int ntot = accu(n1);
  // mask for time periods for treated units that are pre-treatment
  arma::mat bool_mask = as<arma::mat>(opts["mask"]);
 
  
  weightPtr weight_func = *as<wptr>(opts["weight_func"]);
  // arma::mat ipw_weights = as<arma::mat>(opts["ipw_weights"]);

  // initialize gradient as zero
  arma::mat grad = zeros(bool_mask.n_cols, bool_mask.n_rows+1);


  // iterate over treated units
  arma::mat restrict_Xc;
  arma::mat weights;
  for(int j=0; j < bool_mask.n_rows; j++) {
    // apply mask
    restrict_Xc = Xc.each_row() % bool_mask.row(j);
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
