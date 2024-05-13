#ifndef mantrans_transform_identity
#define mantrans_transform_identity
#include <RcppEigen.h>
namespace mantran {
template <typename T>
struct identity : public transform<T> {
  ~identity(){};
  identity(){};

  std::string name() const {
    std::string out = "";
    return(out);
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> toM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x) override {
    Eigen::Matrix<T, Eigen::Dynamic, 1> out(x.size());
    out = x;
    return(out);
  }

  Eigen::Matrix<T, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &x) override {
    Eigen::Matrix<T, Eigen::Dynamic, 1> out(x.size());
    out = x;
    return(out);
  }

  T logdetJfromM(const Eigen::Matrix<T, Eigen::Dynamic, 1> &z) override {
    T out;
    out = 0.;
    return(out);
  }
};
}
#endif
