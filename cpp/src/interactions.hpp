#ifndef GUARD_vsh_interactions_h
#define GUARD_vsh_interactions_h

#include <complex>
#include <eigen3/Eigen/Core>
#include "vsh_translation.hpp"
#include <vector>

using position_t = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
using tmatrix_t  = std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;
using ComplexVector = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;

enum class solver {
    bicgstab,
    exact
};

ComplexVector bicgstab(const Ref<const ComplexMatrix>& A, const Ref<const ComplexVector>& b,
        int maxiter = 1000, double tolerance = 1e-5);

ComplexVector solve_linear_system(const Ref<const ComplexMatrix>& agg_tmatrix,
        const Ref<const ComplexVector>& p_src, solver method = solver::bicgstab);

ComplexMatrix sphere_aggregate_tmatrix(const Ref<const position_t>& positions,
        const Ref<const ComplexMatrix>& mie, double k);

ComplexMatrix particle_aggregate_tmatrix(const Ref<const position_t>& positions,
        tmatrix_t tmatrix, double k);

#endif
