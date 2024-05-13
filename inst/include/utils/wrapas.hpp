# ifndef MYCPP_WRAPAS
# define MYCPP_WRAPAS

# include "../scorematchingad.h"
# include <Rcpp.h>
# include <cppad/cppad.hpp>

// definitions for Rcpp::wrap and Rcpp::as for various data types
namespace Rcpp {
  // from an SEXP to a a1type
  template <> a1type as( SEXP inval) {
    double cppval;
    cppval = Rcpp::as<double>(inval);
    a1type out(cppval);
    return(out);
  }

  // a1type to SEXP
  template <> SEXP wrap(const a1type &inval){
    double cppval(CppAD::Value(inval));
    return(Rcpp::wrap(cppval)); //returns SEXP
  } 

  // from an SEXP to an eigen vector if a1type
  template <> veca1 as( SEXP invec) {
    Rcpp::NumericVector insvec(invec);
    veca1 out(insvec.size());
    for (long int i = 0; i<out.size(); i++){
      out[i] = insvec[i];
    }
    return(out);
  }

  // eigen vector of a1type to SEXP
  template <> SEXP wrap(const veca1 &invec){
    Rcpp::NumericVector out(invec.size());
    for (long int i=0; i<out.size(); i++){
      out[i] = CppAD::Value(invec[i]);
    }
    return(Rcpp::wrap(out)); //returns SEXP
  } 



  // from an SEXP to an eigen matrix of a1type
  template <> mata1 as( SEXP inmat) {
    Rcpp::NumericMatrix inmat2;
    inmat2 = Rcpp::as<Rcpp::NumericMatrix>(inmat);
    mata1 out(inmat2.rows(), inmat2.cols());
    for (long int i = 0; i<out.cols(); i++){
      for (long int j = 0; j<out.rows(); j++){
        out(j, i) = inmat2(j, i);
      }
    }
    return(out);
  }

  // eigen matrix of a1type to SEXP
  template <> SEXP wrap(const mata1 &inmat){
    Rcpp::NumericMatrix out(inmat.rows(), inmat.cols());
    for (long int i = 0; i<out.cols(); i++){
      for (long int j = 0; j<out.rows(); j++){
        out[j + i*out.rows()] = CppAD::Value(inmat(j, i));
      }
    }
    return(Rcpp::wrap(out)); //returns SEXP
  } 

}

# endif
