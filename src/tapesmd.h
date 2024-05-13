#ifndef tapesmd_h
#define tapesmd_h

# include "scorematchingad.h"
# include <cppad/cppad.hpp>
# include "tapebdryweight.h"
# include "manifoldtransforms/bdryweights.hpp"
# include "manifoldtransforms/manifolds.hpp"
# include "utils/PrintFor.hpp"

CppAD::ADFun<double> tapesmd(veca1 u, //a vector. The composition measurement for taping
                             veca1 theta, //a vector of parameters for taping
                             CppAD::ADFun<double> & lltape,
                             transform<a1type> &tran,
                             manifold<a1type> &M,
                             a1type (*h2fun)(const veca1 &, const double &), // the weight function h^2
                             const double & acut, //the acut constraint for the weight functions
                             bool verbose
                             );

//in R store a pointer to the ADFun object
//' @noRd
//' @title The score matching objective calculator.
//' @param xbetain a concatenated vector of sqrt(x) and beta
//' @param n The dimension of x.
//' @param manifoldname The name of the manifold to transform to
//' @param weightname The name of the weight function to use
//' @param acut The constraint a_c in the weight function
//' @return An RCpp::XPtr object pointing to the ADFun
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > ptapesmd(veca1 u_ad,
                                      veca1 theta_ad,
                                      Rcpp::XPtr< CppAD::ADFun<double> > pll,
                                      transform_a1type & tran,
                                      manifold_a1type & man,
                                      std::string weightname,
                                      const double acut,
                                      bool verbose);

#endif
