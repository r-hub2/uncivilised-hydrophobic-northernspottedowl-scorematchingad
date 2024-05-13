#ifndef mantrans_transform_sqrt
#define mantrans_transform_sqrt
// code for various tools for the positive quadrant of the sphere
#include <RcppEigen.h>
namespace mantran {
template <typename Type>
struct sqrt : public transform<Type> {
  ~sqrt(){};
  sqrt(){};

  std::string name() const {
    std::string out = "sqrt";
    return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.cwiseSqrt();
     return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.cwiseProduct(x);
     return(out);
  }

  Type logdetJfromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &z) override {
     Type out;
     out = z.array().log().sum() + 0.6931472 * z.size(); //final number here is log(2)
     return(out);
  }


};
}
#endif
