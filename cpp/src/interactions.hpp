#ifndef GUARD_vsh_interactions_h
#define GUARD_vsh_interactions_h

#include <complex>
#include <eigen3/Eigen/Core>
#include "vsh_translation.h"
#include <vector>

using position_t = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;
using tmatrix_t  = std::vector<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;

enum class solver {
    bicgstab,
    exact
}

ComplexArray solve_linear_system(const Ref<const ComplexMatrix>& agg_tmatrix,
        const Ref<const ComplexArray>& p_src, solver method = solver::bicgstab);

ComplexMatrix sphere_aggregate_tmatrix(const Ref<const position_t>& positions,
        const Ref<const ComplexMatrix>& mie, double k);

ComplexMatrix particle_aggregate_tmatrix(const Ref<const position_t>& positions,
        tmatrix_t tmatrix, double k);

#endif
