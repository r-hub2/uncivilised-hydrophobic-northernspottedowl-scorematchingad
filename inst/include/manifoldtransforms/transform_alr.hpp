#ifndef mantrans_transform_alr
#define mantrans_transform_alr
// code for various tools for the additive log ratio transform
#include <RcppEigen.h>
namespace mantran {//names space for manifold-transformation pair (triplets}
template <typename Type>
struct alr : public transform<Type> {
  ~alr(){};
  alr(){};

  std::string name() const {
    std::string out = "alr";
    return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size() - 1);
     out = x.block(0,0,x.size() - 1, 1) / x[x.size() - 1];
     out = out.array().log();
     return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size() + 1);
     Type one_on_u_p;
     one_on_u_p = x.array().exp().sum() + 1.;
     out << x.array().exp(), 1.;
     out /= one_on_u_p;
     return(out);
  }

  Type logdetJfromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &z) override {
    Eigen::Matrix<Type, Eigen::Dynamic, 1> u(z.size() + 1);
    u = fromM(z);
    Type out;
    out = u.array().log().sum();
    return(out);
  }
};
}//namesspace mantran for manifold-transformation pair (triplets}

#endif
