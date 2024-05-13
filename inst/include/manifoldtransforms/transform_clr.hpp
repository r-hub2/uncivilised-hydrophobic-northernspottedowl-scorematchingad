#ifndef mantrans_transform_clr
#define mantrans_transform_clr
// code for various tools for the additive log ratio transform
#include <RcppEigen.h>
namespace mantran {//names space for manifold-transformation pair (triplets}
template <typename Type>
struct clr : public transform<Type> {
  ~clr(){};
  clr(){};

  std::string name() const {
    std::string out = "clr";
    return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> toM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.array().log(); //log all elements of x
     Eigen::Matrix<Type, Eigen::Dynamic, 1> sumlog(x.size());//use a matrix so that -= is a known operation
     sumlog.setConstant(out.mean()); //sum logged values - mean would work just as well, but sum has fewer operations (except maybe when dimensions are very large?)
     out -= sumlog; //take the sumlog away from each element
     return(out);
  }

  Eigen::Matrix<Type, Eigen::Dynamic, 1> fromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &x) override {
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out(x.size());
     out = x.array().exp(); //exp all elements of x
     Type sumexp = out.sum(); 
     Eigen::Matrix<Type, Eigen::Dynamic, 1> out2(x.size());
     out2 = out / sumexp; //normalise by sum
     return(out2);
  }
  
  Type logdetJfromM(const Eigen::Matrix<Type, Eigen::Dynamic, 1> &z) override {
    Eigen::Matrix<Type, Eigen::Dynamic, 1> u(z.size());
    u = fromM(z);
    Type out; //use a matrix so that -= is a known operation
    out = u.array().log().sum() + log(u.size());
    return(out);
  }
  //could us Sylvester's determinant theorem for direct value
  // or matrix determinant lemma according to Wikipedia

};
}//namesspace mantran for manifold-transformation pair (triplets}

#endif

