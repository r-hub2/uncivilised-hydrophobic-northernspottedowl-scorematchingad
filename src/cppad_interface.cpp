#include "cppad_interface.h"


Rcpp::XPtr< CppAD::ADFun<double> > swapDynamic(Rcpp::XPtr< CppAD::ADFun<double> > pfun, veca1 newvalue, veca1 newdynparam){
  //check inputs and tape match
  if (pfun->Domain() != newdynparam.size()){Rcpp::stop("Size of newdynparam must match domain size of taped function.");}
  if (pfun->size_dyn_ind() != newvalue.size()){Rcpp::stop("Size of newvalue must match the parameter size of the taped function.");}



  //convert taped object to higher order
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  veca1 y(1);

  //START TAPING
  CppAD::Independent(newvalue, newdynparam);

  pfunhigher.new_dynamic(newvalue); //before switch the newvalue is the dynamic parameter vector
  y = pfunhigher.Forward(0, newdynparam); //before the switch the newdynparam is the independent value

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(newvalue, y);
  //out->optimize(); //meant to remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore. But asserts errors.
  out->check_for_nan(false);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}


vecd pJacobian(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam){
  //check inputs and tape match
  if (pfun->Domain() != x.size()){Rcpp::stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){Rcpp::stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}

  vecd grad(x.size());
  pfun->new_dynamic(dynparam);
  grad = pfun->Jacobian(x);  //treat the Rcpp::XPtr as a regular pointer

  return(grad);
}

vecd pForward0(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam){
  //check inputs and tape match
  if (pfun->Domain() != x.size()){Rcpp::stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){Rcpp::stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}

  vecd out(1);
  pfun->new_dynamic(dynparam);
  out = pfun->Forward(0, x);  //treat the Rcpp::XPtr as a regular pointer

  return(out);
}

vecd pHessian(Rcpp::XPtr< CppAD::ADFun<double> > pfun, vecd x, vecd dynparam){
  //check inputs and tape match
  if (pfun->Domain() != x.size()){Rcpp::stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){Rcpp::stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}

  vecd hess(x.size() * x.size(), 1);
  pfun->new_dynamic(dynparam);
  hess = pfun->Hessian(x, 0);  //treat the Rcpp::XPtr as a regular pointer
  return(hess);
}


Rcpp::XPtr< CppAD::ADFun<double> >  pTapeJacobian(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam){
  // x and dynparam must have elements of a1type so that taping can proceed
  //check inputs and tape match
  if (pfun->Domain() != x.size()){Rcpp::stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){Rcpp::stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}



  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Jacobian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  veca1 jac(pfunhigher.Domain() * pfunhigher.Range());
  jac = pfunhigher.Jacobian(x);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(x, jac);
  //out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.  But asserts errors.
  out->check_for_nan(false);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

Rcpp::XPtr< CppAD::ADFun<double> >  pTapeHessian(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam){
  // x and dynparam must have elements of a1type so that taping can proceed
  //check inputs and tape match
  if (pfun->Domain() != x.size()){Rcpp::stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){Rcpp::stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}


  if (pfun->Range()>1){
    Rcpp::stop("Taped function 'pfun' must return a vector of length 1. Currently 'pfun' returns a vector of length %i.", pfun->Range());
  }

  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Hessian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  veca1 hess(pfunhigher.Domain() * pfunhigher.Domain());
  hess = pfunhigher.Hessian(x, 0);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(x, hess);
  //out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore. But asserts errors.
  out->check_for_nan(false);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

std::vector<bool> pParameter(Rcpp::XPtr< CppAD::ADFun<double> > pfun){
  std::vector<bool> isparameter(pfun->Range());
  for (size_t i = 0; i < pfun->Range(); i++){
    isparameter[i] = pfun->Parameter(i);
  }
  return(isparameter);
}
// According to the help, applying Variable(u) to each return value would be false if u depends on the dynamic parameters and does not depend on the independent variable vector.

Rcpp::XPtr< CppAD::ADFun<double> >  pTapeGradOffset(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam){
  // x and dynparam must have elements of a1type so that taping can proceed
  //check inputs and tape match
  if (pfun->Domain() != x.size()){Rcpp::stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){Rcpp::stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}


  if (pfun->Range()>1){
    Rcpp::stop("Taped function 'pfun' must return a vector of length 1. Currently 'pfun' returns a vector of length %i.", pfun->Range());
  }

  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Hessian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  CppAD::Independent(dynparam);  
  pfunhigher.new_dynamic(dynparam);
  veca1 jac(pfunhigher.Domain());
  jac = pfunhigher.Jacobian(x);
  mata1 hess(pfunhigher.Domain() * pfunhigher.Domain(), 1);
  hess = pfunhigher.Hessian(x, 0);
  //arrange hess into a matrix
  hess.resize(pfunhigher.Domain(),pfunhigher.Domain());

  veca1 gradoffset(pfunhigher.Domain());
  gradoffset = jac - (hess * x);

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(dynparam, gradoffset);
  //out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.  But asserts errors.
  out->check_for_nan(false);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

Rcpp::XPtr< CppAD::ADFun<double> >  ptapelogdetJ(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                    veca1 x, veca1 dynparam){
  // domain and range must have equal size for the determinant of the Jacobian to make sense
  if (pfun->Domain() != pfun->Range()){Rcpp::stop("Domain (size %i) and range (size %i) need to be equal for determinant of Jacobian.", pfun->Domain(), pfun->Range());}
  // x and dynparam must have elements of a1type so that taping can proceed
  //check inputs and tape match
  if (pfun->Domain() != x.size()){Rcpp::stop("Size of input vector %i does not match domain size %i of taped function.", x.size(), pfun->Domain());}
  if (pfun->size_dyn_ind() != dynparam.size()){Rcpp::stop("Size of parameter vector %i does not match parameter size %i of the taped function.", dynparam.size(), pfun->size_dyn_ind());}
  
  //convert taped object to higher order, so that the 'base' type of the tape is a1type, so x and dynparam can be passed into Jacobian()
  CppAD::ADFun<a1type, double> pfunhigher;
  pfunhigher = pfun->base2ad();

  CppAD::Independent(x, dynparam);  //start taping with x as the usual independent parameter and dynparam as the dynamic parameter
  pfunhigher.new_dynamic(dynparam);
  mata1 jacmat(pfunhigher.Domain() * pfunhigher.Range(), 1);
  jacmat = pfunhigher.Jacobian(x);
  jacmat.resize(pfunhigher.Domain(), pfunhigher.Range()); //first row is: du1/dz1, du2/dz1, du3/dz1. Second row is du1/dz2, du2/dz2, du3/dz2
  veca1 logdet(1);
  logdet[0] = CppAD::log(CppAD::abs(jacmat.determinant()));

  //end taping
  CppAD::ADFun<double>* out = new CppAD::ADFun<double>; //returning a pointer
  out->Dependent(x, logdet);
  //out->optimize(); //remove some of the extra variables that were used for recording the ADFun f above, but aren't needed anymore.  But asserts errors.
  out->check_for_nan(false);

  Rcpp::XPtr< CppAD::ADFun<double> > pout(out, true);
  return(pout);
}

