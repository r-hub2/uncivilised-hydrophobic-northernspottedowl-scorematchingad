#ifndef likelihoods_vMF
#define likelihoods_vMF
# include <RcppEigen.h>

namespace ll {

  template <class T>
  T ll_vMF(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
                 const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta)
    //lklhood is log(k*mu.x), where mu is a unit vector
    //for this software theta is the vector k*mu, which is unrestricted (except mu shouldn't be zero)
  {
    T y(0.);  // initialize summation
    for(size_t i = 0; i < u.size(); i++)
    {   y   += theta[i] * u[i];
    }
    return y;
  }

}

#endif
