#ifndef GUARD_vec_h
#define GUARD_vec_h

#include <eigen3/Eigen/Core>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <vector>

using vec3          = Eigen::Vector3d;
using cvec3         = Eigen::Vector3cd;
using Array         = Eigen::Array<double, Eigen::Dynamic, 1>;
using ComplexArray  = Eigen::Array<std::complex<double>, Eigen::Dynamic, 1>;
using ComplexVector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
using Matrix        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ComplexMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Tensor3       = Eigen::Tensor<double,3,Eigen::RowMajor>;
using Tensor4       = Eigen::Tensor<double,4,Eigen::RowMajor>;
using ComplexTensor3 = Eigen::Tensor<std::complex<double>,3,Eigen::RowMajor>;
using ComplexTensor4 = Eigen::Tensor<std::complex<double>,4,Eigen::RowMajor>;
using AB_type       = Eigen::Matrix<std::complex<double>, 2, Eigen::Dynamic, Eigen::RowMajor>;
using Eigen::Ref;

using position_t = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using tmatrix_t = Eigen::TensorMap<ComplexTensor3>;

const double PI = 3.14159265358979323846;

#endif
