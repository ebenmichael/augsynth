
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "augsynth_types.h"
#include <numeric>
#include <string>

using namespace arma;
using namespace Rcpp;
using namespace std;


//' Generic balancing loss gradient
// [[Rcpp::export]]
arma::mat balancing_grad_multisynth(arma::mat theta, List opts) {

  // control and treated data
  arma::mat Xc = as<arma::mat>(opts["Xc"]);
  arma::mat Xt = as<arma::mat>(opts["Xt"]);

  // number of units at each level
  arma::vec n1 = as<arma::vec>(opts["n1"]);
  int ntot = accu(n1);
  // mask for time periods for treated units that are pre-treatment
  arma::mat bool_mask = as<arma::mat>(opts["mask"]);
 
  
  weightPtrIPW weight_func = *as<wptripw>(opts["weight_func"]);
  arma::mat ipw_weights = as<arma::mat>(opts["ipw_weights"]);

  // initialize gradient as zero
  arma::mat grad = zeros(bool_mask.n_cols, bool_mask.n_rows+1);


  // iterate over treated units
  arma::mat restrict_Xc;
  arma::mat weights;
  for(int j=0; j < bool_mask.n_rows; j++) {
    // apply mask
    restrict_Xc = Xc.each_row() % bool_mask.row(j);
    weights = weight_func(restrict_Xc, theta.col(0) + theta.col(j+1),
                          ipw_weights);
    grad.col(j+1) = restrict_Xc.t() * weights;

    grad.col(0) += restrict_Xc.t() * weights * n1[j]/ntot;
  }

  //combine to get gradient
  return grad - Xt;
    
}


// [[Rcpp::export]]
gptr make_balancing_grad_multisynth() {
  return gptr(new gradPtr(balancing_grad_multisynth));
}
