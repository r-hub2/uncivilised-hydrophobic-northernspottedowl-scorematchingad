#ifndef tapell_h
#define tapell_h

# include <RcppEigen.h>
# include <Rcpp.h>
# include <cppad/cppad.hpp>
# include "scorematchingad.h"
# include "likelihoods/likelihoods.hpp"
# include "utils/PrintFor.hpp"

// declare a function that tapes a log likelihood
CppAD::ADFun<double> tapellcpp(veca1 z, //data measurement tranformed to M manifold
                            veca1 theta, //theta parameter
                               a1type (*llf)(const veca1 &, const veca1 &), //the log likelihood function
                               transform<a1type> & tran, //it seems points must be passed for abstract classes (note error when compiling without the *, and Stefan's demo)
                               Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta, //TRUE (1) values indicate that the corresponding value of theta is not a variable (dynamic or independent)
                               bool verbose
                               );


//' @noRd
//' @title Tape of a log-likelihood calculation 2
//' @param p dimension of measurements
//' @param bd dimension of the parameter vector
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to the ADFun
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > ptapell2(veca1 z_ad, //data measurement on the M manifold
                                     veca1 theta_ad,
                                     Rcpp::XPtr<llPtr> llfXPtr, //the log likelihood function
                                     transform_a1type & tran,
                                     Eigen::Matrix<int, Eigen::Dynamic, 1> fixedtheta,
                                     bool verbose
                                     );

//' @noRd
//' @title Get an XPtr to a named log-likelihood function in source code of package
//' @param llname name of the likelihood function
//' @return An RCpp::XPtr object pointing to a `llPtr` object of the log-likelihood function. Since `llPtr` is itself a pointer object, we have an XPtr pointing to a pointer that points to a function.
// [[Rcpp::export]]
Rcpp::XPtr<llPtr> getllptr(std::string llname);

//' @name evalll
//' @title Evaluate a custom log-likelihood function
//' @description Evaluates a custom log-likelihood function from [`customll()`] without taping. This is useful to check that the custom log-likelihood is behaving.
//' To check a tape of the custom log-likelihood use [`buildsmdtape()`] then [`evaltape()`].
//' @param ll A compiled log-likelihood function created by [`customll()`].
// ( ll is an XPtr to a llPtr object that points to a log-likelihood function )
//' @param u A vector of measurements for an individual
//' @param theta A vector of parameters
//' @return The value of the log-likelihood at `u` with parameters `theta`.
//' @export
// [[Rcpp::export]]
a1type evalll(Rcpp::XPtr<llPtr> ll, const veca1& u, const veca1& theta);

#endif
