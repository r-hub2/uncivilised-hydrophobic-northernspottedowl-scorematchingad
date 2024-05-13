# ifndef CPPAD_INTERFACE
# define CPPAD_INTERFACE

// things for manipulating and evaluating CppAD::ADFun objects from R

//for content that is Rcpp specific
#include "scorematchingad.h"
#include "utils/wrapas.hpp"  //needed because converting veca1 from R

//' @noRd
//' @name evaltape_internal
//' @title Advanced: Evaluate `CppAD` Tapes via their Pointer
//' @description The recommended method for evaluating tapes is [`evaltape()`].
//' Internally, `evaltape()` and other methods are using the methods documented here.
//' There methods access the tapes using `Rcpp::XPtr` objects and perform evaluations a single point at a time. 
//' @return A vector of numeric values, except `pParameter()`, which returns logical values.
//' @family tape evaluators
//' @param pfun Rcpp::XPtr to an ADFun. Can be obtained as the `ptr` field of an [`ADFun`] object.
//' @param x A vector in the domain of the taped function
//' @param dynparam a vector of the dynamic parameters, if `pfun` has no dynamic parameter than pass `vector("numeric")`.

//' @noRd
//' @describeIn evaltape_internal Evaluates a tape without any differentiation at the given values of `x` and dynparam. 
//' The name `pForward0` is a reference to the zero order `CppAD` method [`forward`](https://cppad.readthedocs.io/forward_zero.html), and the prefix 'p' is because the tape is passed as a pointer.
//' @param pfun Rcpp::XPtr to an ADFun. Can be obtained as the `ptr` field of an [`ADFun`] object.
//' @param x A vector in the domain of the taped function
//' @param dynparam a vector of the dynamic parameters, if `pfun` has no dynamic parameter than pass `vector("numeric")`.
// [[Rcpp::export]]
vecd pForward0(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam);

//' @noRd
//' @describeIn evaltape_internal Evaluates a the Jacobian of a tape using the `CppAD` `Jacobian` method <https://cppad.readthedocs.io/latest/Jacobian.html>. 
// [[Rcpp::export]]
vecd pJacobian(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam);

//' @noRd
//' @describeIn evaltape_internal Evaluates a the Hessian of a tape using the `CppAD` `Hessian` method <https://cppad.readthedocs.io/latest/Hessian.html>, assuming that range space of the taped function has dimension of `1`. 
// [[Rcpp::export]]
vecd pHessian(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam);

//' @noRd
//' @describeIn evaltape_internal Test whether the returned values are constant with respect to the independent values using 
//' `CppAD`'s `Parameter` method <https://cppad.readthedocs.io/latest/fun_property.html>.
//' Returns A vector of logical values. `TRUE` indicates that element of the tape result is constant.
//' @details 
//' # pParameter
//' The `CppAD` function [`Parameter(i)`](https://cppad.readthedocs.io/latest/fun_property.html#parameter) returns `TRUE` when the `i`th component of the range does not depend on the independent value
//' (the `i`th component may still depend on the value of the dynamic parameters - see <https://cppad.readthedocs.io/latest/glossary.html#dynamic> ).
// [[Rcpp::export]]
std::vector<bool> pParameter(Rcpp::XPtr< CppAD::ADFun<double> > pfun);
// According to the help, applying Variable(u) to each return value would be false if u depends on the dynamic parameters and does not depend on the independent variable vector.




//' @noRd
//' @title Tape the Jacobian of CppAD Tape
//' @family tape builders
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with dynamic parameters and independent parameters
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters
//' @description Creates a tape of the Jacobian of function taped by CppAD.
//' When the function returns a real value (as is the case for densities and the score matching objective) the Jacobian is equivalent to the gradient.
//' The `x` vector is used as the value to conduct the taping.
//' @details
//' When the returned tape is evaluated (via say [`pForward0()`], the resultant vector contains the Jacobian in long format (see <https://cppad.readthedocs.io/latest/Jacobian.html>).
//' Suppose the function represented by `pfun` maps from \eqn{n}-dimensional space to \eqn{m}-dimensional space, then
//' the first \eqn{n} elements of vector is the gradient of the first component of function output.
//' The next \eqn{n} elements of the vector is the gradient of the second component of the function output.
//' The Jacobian as a matrix, could then be obtained by [`as.matrix()`] with `byrow = TRUE` and `ncol = n`.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object.
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> >  pTapeJacobian(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam);

//' @noRd
//' @title Tape the Hessian of a CppAD Tape
//' @family tape builders
//' @inheritParams pTapeJacobian
//' @description Creates a tape of the Hessian of a function taped by CppAD.
//' The taped function represented by `pfun` must be scalar-valued (i.e. a vector of length 1).
//' The `x` vector and `dynparam` are used as the values to conduct the taping.
//' @details
//' When the returned tape is evaluated (via say [`pForward0()`], the resultant vector contains the Hessian in long format (see <https://cppad.readthedocs.io/latest/Hessian.html>).
//' Suppose the function represented by `pfun` maps from \eqn{n}-dimensional space to \eqn{1}-dimensional space, then
//' the first \eqn{n} elements of the vector is the gradient of the partial derivative with respect to the first dimension of the function's domain.
//' The next \eqn{n} elements of the vector is the gradient of the partial derivative of the second dimension of the function's domain.
//' The Hessian as a matrix, can be obtained by using [`as.matrix()`] with `ncol = n`.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object.
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> >  pTapeHessian(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam);


//' @noRd
//' @title Tape the Gradient Offset of a Quadratic CppAD Tape
//' @family tape builders
//' @inheritParams pTapeJacobian
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object. The independent argument to the function are the dynamic parameters of `pfun`.
//' @details A quadratic function can be written as
//' \deqn{f(x;\theta) = \frac{1}{2} x^T W(\theta) x + b(\theta)^Tx + c.}
//' The gradient of \eqn{f(x; \theta)} with respect to \eqn{x} is
//' \deqn{\Delta f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T)x + b(\theta).}
//' The Hessian is 
//' \deqn{H f(x; \theta) = \frac{1}{2}(W(\theta) + W(\theta)^T),}
//' which does not depend on \eqn{x},
//' so the gradient of the function can be rewritten as
//' \deqn{\Delta f(x;\theta) = H f(x; \theta) x + b(\theta)^T.}
//' The tape calculates \eqn{b(\theta)} as
//'  \deqn{b(\theta) = \Delta f(x;\theta) - H f(x; \theta) x,}
//' which does not depend on \eqn{x}.
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> >  pTapeGradOffset(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam);


//' @noRd
//' @title Tape the log of Jacobian determinant of a CppAD Tape
//' @family tape builders
//' @param pfun Rcpp::XPtr to an ADFun tape a tape with dynamic parameters and independent parameters
//' @param x A vector in the domain of the taped function.
//' @param dynparam a vector of the dynamic parameters
//' @description Creates a tape of the log of the Jacobian determinant of a function taped by CppAD.
//' The `x` vector is used as the value to conduct the taping.
//' @return A `Rcpp::XPtr` to a CppAD::ADFun object.
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> >  ptapelogdetJ(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam);

//' @noRd
//' @title Switch Dynamic and Independent Values of a Tape
//' @family tape builders
//' @description Convert an ADFun so that the independent values become dynamic parameters
//' and the dynamic parameters become independent values
//' @param pfun An Rcpp::XPtr to an ADFun object (i.e. a tape of a function)
//' @param newvalue The independent value (in the sense after the switch has occurred) at which to tape the ADFun
//' @param newdynparam The value of the dynamic parameters (after the switch) at which to tape the ADFun
//' @return A pointer to an ADFun
// [[Rcpp::export]]
Rcpp::XPtr< CppAD::ADFun<double> > swapDynamic(Rcpp::XPtr< CppAD::ADFun<double> > pfun, veca1 newvalue, veca1 newdynparam);

# endif
