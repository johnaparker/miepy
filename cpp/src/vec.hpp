#ifndef GUARD_vec_h
#define GUARD_vec_h

#include <eigen3/Eigen/Core>

using vec3          = Eigen::Vector3d;
using Array         = Eigen::Array<double, Eigen::Dynamic, 1>;
using ComplexArray  = Eigen::Array<std::complex<double>, Eigen::Dynamic, 1>;
using ComplexVector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
using Matrix        = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using ComplexMatrix = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using AB_type       = Eigen::Matrix<std::complex<double>, 2, Eigen::Dynamic, Eigen::RowMajor>;
using Eigen::Ref;

#endif
