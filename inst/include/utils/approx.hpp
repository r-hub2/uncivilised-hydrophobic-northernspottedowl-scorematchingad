# ifndef utils_approx
# define utils_approx
#include <RcppEigen.h>
#include <cppad/cppad.hpp> // the CppAD package
  template <class Type>
  Eigen::Matrix<Type, Eigen::Dynamic, 1> taylorapprox(
		  CppAD::ADFun<Type> &f,  //a tape with independent values that are points on the manifold (not the parameters)
		  const Eigen::Matrix<Type, Eigen::Dynamic, 1> centre,
		  const size_t order,
		  const Eigen::Matrix<Type, Eigen::Dynamic, 1> x){
    //In CppAD speak consider the input function X(t) to be
    //X(t) = centre + t*(x - centre). So X(0) = centre, X(1) = x.
    //First derivative of X at 0, is x - centre
    //Higher order derivative are all 0.
    Eigen::Matrix<Type, Eigen::Dynamic, 1> out(f.Range());
    Eigen::Matrix<Type, Eigen::Dynamic, 1> diff(x.size());
    out.setZero();
    out += f.Forward(0, centre); //zeroth order - constant component of taylor expansion
    if (order >= 1){
      diff = x - centre; // for some reason forward can't take the lazy evaluation version of x - centre direclty. (x - centre).eval() also works
      out += f.Forward(1, diff); //now out[0] is evaluation of a linear approximation of f
    }
    if (order >= 2){
      for (size_t i=2; i<=order; i++){
	      diff.setZero();
        out += f.Forward(i, diff); //now out[0] is evaluation of a quadratic approximation of f
      }
    }
    return(out);
  }

# endif
