#ifndef tapebdryweight
#define tapebdryweight

# include "scorematchingad.h"

CppAD::ADFun<double> tapeh2(veca1 z,
                            a1type (*h2fun)(const veca1 &, const double &),
                            const double & acut);

#endif
