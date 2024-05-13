# include "approxwrap.h"

vecd pTaylorApprox(Rcpp::XPtr< CppAD::ADFun<double> > pfun,
                     vecd x, vecd centre,
                     vecd dynparam, size_t order){
  vecd out(pfun->Range());
  pfun->new_dynamic(dynparam);
  out = taylorapprox(*pfun,
                     centre,
                     order,
                     x);

  return(out);
}

