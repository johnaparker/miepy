#ifndef GUARD_vec_h
#define GUARD_vec_h

#include <eigen3/Eigen/Core>
#include <vector>

using vec3          = Eigen::Vector3d;
using cvec3         = Eigen::Vector3cd;
using Array         = Eigen::Array<double, Eigen::Dynamic, 1>;
using ComplexArray  = Eigen::Array<std::complex<double>, Eigen::Dynamic, 1>;
using ComplexVector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
using Matrix        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ComplexMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using AB_type       = Eigen::Matrix<std::complex<double>, 2, Eigen::Dynamic, Eigen::RowMajor>;
using Eigen::Ref;

using position_t = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using tmatrix_t  = std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

#endif
