#ifndef mantrans_manifold_Euc
#define mantrans_manifold_Euc
// code for various tools for the additive log ratio transform
#include <RcppEigen.h>
namespace mantran {//names space for manifold-transformation pair (triplets}
template <typename Type>
struct Euc : public manifold<Type> {
  ~Euc(){};
  Euc(){};

  std::string name() const {
    std::string out = "Euc";
    return(out);
  }


  // manifold tangent-plane projection matrix P (for isometric(?) embeddings this is closely related to the manifold metric
  Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
    int n = x.size();
    Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);
    Pmat = Eigen::Matrix<Type, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n);
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
