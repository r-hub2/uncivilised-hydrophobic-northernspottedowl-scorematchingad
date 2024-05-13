#ifndef likelihoods_Bingham
#define likelihoods_Bingham

# include <RcppEigen.h>

namespace ll {
  template <class T>
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
  BinghamMatrix(const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta){
    // A is symmetric, and the sum of the diagonals is zero
    //assume the parameter vector theta is encoded as:
    //c(diag(A)[1:(p-1)], A[upper.tri(ALs)]
    size_t d = (-1 + std::sqrt(1 + 4 * 2 * (1+theta.size()))) / 2 + 0.5;//the +0.5 makes sure truncation gets to the correct integer
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Amat(d, d);
    Amat.setZero();
    //populate the diagonal
    size_t vecidx = 0;
    for (size_t row=0; row < d-1; row++){
      Amat(row,row) = theta[vecidx];
      vecidx +=1;
    }
    Amat(d-1, d-1) = -theta.block(0, 0, d - 1, 1).sum(); //so trace of A is zero

    //populate the upper and lower triangles
    //the upper triangle has d-1 rows, the rows have 1 to (d-1) elements. Arithmetic series:
    //(d-1)/2 [2 + (d-1−1)] = (d-1) (2 + d - 2)/2 = (d - 1) d/2
    Eigen::Matrix<T, Eigen::Dynamic, 1> upptriblock((d - 1) * d/2);
    upptriblock = theta.block(d-1, 0, upptriblock.size(), 1);
    vecidx = 0;
    for (size_t col=1; col < d; col++){
      for (size_t row=0; row < col; row++){
        Amat(row, col) = upptriblock[vecidx]; //bug fix - column index is 0!!
        Amat(col, row) = upptriblock[vecidx];
        vecidx +=1;
      }
    }
    return(Amat);
  }

  template <class T>
  T ll_Bingham(const Eigen::Matrix<T, Eigen::Dynamic, 1> &u,
           const Eigen::Matrix<T, Eigen::Dynamic, 1> &theta)
    //lklhood is log(u*A*u) - from Mardia et al 2016
    // A is symmetric, and the sum of the diagonals is zero
  {
    size_t d  = u.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Amat(d, d);
    Amat = BinghamMatrix(theta);

    Eigen::Matrix<T, 1, 1> out_e;
    out_e = u.transpose() * Amat * u;
    T out(out_e[0]);
    return out;
  }
}
#endif
