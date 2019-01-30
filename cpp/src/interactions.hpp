#ifndef GUARD_vsh_interactions_h
#define GUARD_vsh_interactions_h

#include <complex>
#include <eigen3/Eigen/Core>
#include <vector>
#include "vec.hpp"

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
        const tmatrix_t& tmatrix, double k);

ComplexMatrix reflection_matrix_nia(const Ref<const position_t>& positions,
        const Ref<const ComplexMatrix>& mie, double k, std::complex<double> reflection, double z);

#endif
