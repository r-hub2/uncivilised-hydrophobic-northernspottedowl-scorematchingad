#ifndef mantrans_manifold_sim
#define mantrans_manifold_sim
#include <RcppEigen.h>
namespace mantran {
template <typename T>
struct sim : public manifold<T> {
  ~sim(){};
  sim(){};

  std::string name() const {
    std::string out = "sim";
    return(out);
  }

  // manifold tangent-plane projection matrix P (for isometric(?) embeddings this is closely related to the manifold metric
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Pmatfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x) override {
    int n = x.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> Pmat(n, n);
    Eigen::Matrix<T, Eigen::Dynamic, 1> ones(n);
    ones.setOnes();
    double nd = n;
    Pmat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Identity(n,n) - (ones*ones.transpose()/nd);
    return(Pmat);
  }

  //partial derivative of the tangent-plane projection matrix
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> dPmatfun(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x, const int &d) override {
    int n = x.size();
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> bvx(n, n);
    bvx.setZero();
    return(bvx);
  }
};
}
#endif
