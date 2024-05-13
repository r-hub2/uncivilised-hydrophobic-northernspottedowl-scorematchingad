#ifndef PrintForIncluded_cpp
#define PrintForIncluded_cpp

#include <RcppEigen.h>
#include <cppad/cppad.hpp>

template <class T> //T must be something CppAD::PrintFor can print
void PrintForVec(const char* before, const T & printvec){
  CppAD::PrintFor(before, printvec[0]);
  if (printvec.size() > 1){
    for(size_t i=1; i<printvec.size(); i++){
      CppAD::PrintFor(" ", printvec[i]);
    }
  }
}

template <class T> //T must be something CppAD::PrintFor can print
void PrintForMatrix(const char* before, const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> & printmat){
  CppAD::PrintFor(before, printmat(0,0));
  for(size_t i=1; i<printmat.cols(); i++){
    CppAD::PrintFor(" ", printmat(0, i));
  }
  for(size_t row=1; row<printmat.rows(); row++){
    CppAD::PrintFor("\n", printmat(row,0));
    for(size_t col=1; col<printmat.cols(); col++){
      CppAD::PrintFor(" ", printmat(row, col));
    }
  }
}

#endif
