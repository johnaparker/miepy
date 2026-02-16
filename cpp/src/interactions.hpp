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

struct bicgstab_result {
    ComplexVector solution;
    int iterations;
    double residual;
};

bicgstab_result bicgstab_profiled(const Ref<const ComplexMatrix>& A,
        const Ref<const ComplexVector>& b, int maxiter = 1000, double tolerance = 1e-5);

ComplexVector solve_linear_system(const Ref<const ComplexMatrix>& agg_tmatrix,
        const Ref<const ComplexVector>& p_src, solver method = solver::bicgstab);

ComplexVector apply_block_preconditioner(const Ref<const ComplexMatrix>& M_inv_blocks,
        const Ref<const ComplexVector>& x, int block_size);

ComplexVector bicgstab_preconditioned(const Ref<const ComplexMatrix>& A,
        const Ref<const ComplexVector>& b,
        const Ref<const ComplexMatrix>& M_inv_blocks, int block_size,
        int maxiter = 1000, double tolerance = 1e-5);

bicgstab_result bicgstab_preconditioned_profiled(const Ref<const ComplexMatrix>& A,
        const Ref<const ComplexVector>& b,
        const Ref<const ComplexMatrix>& M_inv_blocks, int block_size,
        int maxiter = 1000, double tolerance = 1e-5);

ComplexVector solve_linear_system_preconditioned(const Ref<const ComplexMatrix>& agg_tmatrix,
        const Ref<const ComplexVector>& p_src,
        const Ref<const ComplexMatrix>& M_inv_blocks, int block_size,
        solver method = solver::bicgstab);

ComplexMatrix sphere_aggregate_tmatrix(const Ref<const position_t>& positions,
        const Ref<const ComplexMatrix>& mie, double k);

ComplexMatrix particle_aggregate_tmatrix(const Ref<const position_t>& positions,
        const tmatrix_t& tmatrix, double k);

ComplexMatrix reflection_matrix_nia(const Ref<const position_t>& positions,
        const Ref<const ComplexMatrix>& mie, double k, std::complex<double> reflection, double z);

#endif
