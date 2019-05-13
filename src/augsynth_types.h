#include <RcppArmadillo.h>
#include <Rcpp.h>


using namespace arma;
using namespace Rcpp;


typedef double (*lossPtr)(arma::vec x, List opts);
typedef mat (*gradPtr)(mat x, List opts);
typedef mat (*proxPtr)(mat x, double t, List opts);
typedef mat (*weightPtr)(mat X, mat theta);
typedef mat (*weightPtr2)(mat eta);
typedef mat (*weightPtrIPW)(mat X, mat theta, mat q);

typedef XPtr<gradPtr> gptr;
typedef XPtr<proxPtr> pptr;
typedef XPtr<weightPtr> wptr;
typedef XPtr<weightPtr2> wptr2;
typedef XPtr<weightPtrIPW> wptripw;


typedef mat (*fullWeightPtr)(mat Xc, mat theta, wptr weight_func, List opts);
typedef XPtr<fullWeightPtr> fwptr;



typedef double (*kernelPtr)(mat x, mat y, double p);
typedef XPtr<kernelPtr> kptr;
