#ifndef likelihoods_FB
#define likelihoods_FB

#include <RcppEigen.h>

namespace ll {

  template <class T>
  T ll_FB(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
               const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta)
    //lklhood is log(u*A*u + km*u) - from Mardia et al 2016
    // A is symmetric, and the sum of the diagonals is zero
    //m*k is any vector
  {
    size_t d  = u.size();
    size_t Binghamthetasize = d - 1 + (d-1)*d/2;
    Eigen::Matrix<T, Eigen::Dynamic, 1> Btheta;
    Btheta = theta.block(0,0, Binghamthetasize, 1);
    Eigen::Matrix<T, Eigen::Dynamic, 1> Ftheta;
    Ftheta = theta.block(Binghamthetasize,0, d, 1);
    T out;
    out = ll_Bingham(u, Btheta);
    out += ll_vMF(u, Ftheta);
    return(out);
  }

}

#endif
