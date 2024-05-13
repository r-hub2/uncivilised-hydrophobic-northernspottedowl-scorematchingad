# include "tapebdryweight.h"

CppAD::ADFun<double> tapeh2(veca1 z,
                            a1type (*h2fun)(const veca1 &, const double &),
                            const double & acut){
  //tape relationship between z and h2
  CppAD::Independent(z);
  // range space vector
  size_t m = 1;               // number of ranges space variables
  veca1 y(m); // vector of ranges space variables
  y[0] = h2fun(z, acut);
  CppAD::ADFun<double> tape;  //copying the change_parameter example, a1type is used in constructing f, even though the input and outputs to f are both a2type.
  tape.Dependent(z, y);
  return(tape);
}

