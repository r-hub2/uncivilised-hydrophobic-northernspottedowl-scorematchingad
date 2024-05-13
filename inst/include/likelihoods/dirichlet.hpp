#ifndef likelihoods_dirichlet
#define likelihoods_dirichlet
#include <RcppEigen.h>
// # include <cppad/example/atomic_three/mat_mul.hpp> // for matrix multiplication tapesmo's matrix multiplication don't seem to need this though :S

namespace ll { // namespace for log-likelihood functions

    template <class T>
    T ll_dirichlet(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
	       const Eigen::Matrix<T, Eigen::Dynamic, 1> &beta)
    {   size_t d  = u.size();
        T y(0.);  // initialize summation
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }
   

} // namespace ll

#endif
