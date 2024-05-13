#ifndef mantrans_manifold_Hn111
#define mantrans_manifold_Hn111
// code for various tools for the additive log ratio transform
#include <RcppEigen.h>
#include <cppad/cppad.hpp> //for CppAD::log
namespace mantran {//names space for manifold-transformation pair (triplets}
template <typename Type>
struct Hn111 : public manifold<Type> {
  ~Hn111(){};
  Hn111(){};
  
  std::string name() const {
    std::string out = "Hn111";
    return(out);
  }

  //Pmat is the same as for simplex (both planes with normal of (1,1,1,....1)
  // manifold tangent-plane projection matrix P (for isometric(?) embeddings this is closely related to the manifold metric
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);
    Eigen::Matrix<Type, Eigen::Dynamic, 1> ones(n);
    ones.setOnes();
    double nd = n;
    Pmat = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n) - (ones*ones.transpose()/nd);
    return(Pmat);
  }

  //partial derivative of the tangent-plane projection matrix
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> dPmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x, const int &d) override {
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> bvx(n, n);
    bvx.setZero();
    return(bvx);
  }
};
}//namesspace mantran for manifold-transformation pair (triplets}

#endif

