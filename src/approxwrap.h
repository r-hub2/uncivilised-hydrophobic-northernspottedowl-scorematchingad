#ifndef approxwrap
#define approxwrap

#include "scorematchingad.h"
#include "manifoldtransforms/manifolds.hpp"
#include "utils/approx.hpp"
#include <Rcpp.h>

//' @noRd
//' @describeIn evaltape_internal The value of a recorded function approximated by Taylor expansion.
//' Returns the approximate value of `pfun` at `x`.
//' @details
//' # pTaylorApprox 
//' Approximates the value of a `CppAD` tape at `x` using a Taylor approximation at `centre`. The dynamic parameters of the tape are set by `dynparam`.
//' @param centre For pTaylorApprox. A vector in the domain of the taped function to approximate the value at `x` from.
//' @param order For pTaylorApprox. The order of Taylor expansion to use.
// [[Rcpp::export]]
vecd pTaylorApprox(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                     vecd x, vecd centre,
                     vecd dynparam, size_t order);


#endif
