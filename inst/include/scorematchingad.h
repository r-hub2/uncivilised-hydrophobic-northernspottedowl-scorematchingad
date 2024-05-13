# ifndef scorematchingad_TYPES
# define scorematchingad_TYPES

# include <RcppCommon.h> //don't use Rcpp.h until after the wraps and as declaration at the end of this file
# include <iostream>        // standard input/output
# include <vector>          // standard vector
# include <cppad/example/cppad_eigen.hpp>  //load eigen
# include <cppad/cppad.hpp> // the CppAD package

// [[Rcpp::depends(RcppEigen)]] //include RcppEigen here so that Eigen:: below makes sense
#include <RcppEigen.h>

typedef std::vector<double> svecd;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> vecd; //a vector of double values
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> matd;//a matrix of double
typedef CppAD::AD<double> a1type;   // for first (outer) level of taping
typedef Eigen::Matrix<a1type, Eigen::Dynamic, 1> veca1; //a vector of a1type values
typedef Eigen::Matrix<a1type, Eigen::Dynamic, Eigen::Dynamic> mata1;//a matrix of a1types
typedef CppAD::AD<a1type> a2type;  // for second (inner) level of taping
typedef Eigen::Matrix<a2type, Eigen::Dynamic, 1> veca2;

template <typename T>
struct manifold { //exactly like a class, but with default public members https://stackoverflow.com/questions/24196885/can-c-struct-have-member-functions
  //make these members 'pure virtual' means this is an 'abstract class'
  virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &) = 0; //projection matrix for manifold, given vector on the manifold.
  virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dPmatfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &, const int &) = 0; //elementwise derivative of projection matrix for manifold
  //for taping will need to pass copies - so that the coefficients of the tape are not updated by other calls part way through
  virtual ~manifold(){}; //destructor
  std::string name() const { //function giving the manifold-transform name
    std::string out = "NA";
    return(out);
  };
  manifold(){};
};

template <typename T>
struct transform { //exactly like a class, but with default public members https://stackoverflow.com/questions/24196885/can-c-struct-have-member-functions
  //make these members 'pure virtual' means this is an 'abstract class'
  virtual Eigen::Matrix<T, Eigen::Dynamic, 1> toM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &) = 0; //map from simplex to manifold
  virtual Eigen::Matrix<T, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &) = 0; //transformation from manifold to simplex
  virtual T logdetJfromM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &) = 0; //determinant of Jacobian of the tranformation - for correcting the likelihood function as it is a density
  //for taping will need to pass copies - so that the coefficients of the tape are not updated by other calls part way through
  virtual ~transform(){}; //destructor
  std::string name() const { //function giving the transform name
    std::string out = "NA";
    return(out);
  };
  transform(){};
};


// the following is for passing around likelihood functions
typedef a1type (*llPtr)(const veca1&, const veca1&);

typedef manifold<a1type> manifold_a1type;
RCPP_EXPOSED_CLASS_NODECL(manifold_a1type)
typedef transform<a1type> transform_a1type;
RCPP_EXPOSED_CLASS_NODECL(transform_a1type)

namespace Rcpp {
  template <> a1type as( SEXP );
  template <> SEXP wrap(const a1type&);
  template <> veca1 as( SEXP );
  template <> SEXP wrap(const veca1&);
  template <> mata1 as( SEXP );
  template <> SEXP wrap(const mata1&);
}

# endif

