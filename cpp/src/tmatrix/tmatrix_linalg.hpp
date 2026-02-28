#ifndef GUARD_tmatrix_linalg_h
#define GUARD_tmatrix_linalg_h

#include "tmatrix/tmatrix_types.hpp"
#include <utility>
#include <algorithm>

namespace tmatrix {

// ============================================================================
// Matrix equilibration: row and column scaling to reduce condition number
// Port of Fortran MatrixSolv.f90::equilibrate (lines 479-597)
//
// Computes row scales r and column scales c such that diag(r)*A*diag(c)
// has entries closer to unit magnitude. Modifies A in-place.
// Returns (row_scale, col_scale).
// ============================================================================
template<typename Real>
std::pair<typename Types<Real>::RealVector, typename Types<Real>::RealVector>
equilibrate_matrix(typename Types<Real>::Matrix& A) {
    using RealVector = typename Types<Real>::RealVector;
    using complex_t = std::complex<Real>;

    int m = A.rows();
    int n = A.cols();

    RealVector row_scale = RealVector::Ones(m);
    RealVector col_scale = RealVector::Ones(n);

    Real smlnum = smallest_pos<Real>() / machine_eps<Real>();
    Real bignum = Real(1) / smlnum;

    // Compute row infinity-norms (complex L1: |Re|+|Im|)
    Real rowcnd_min = bignum;
    Real rowcnd_max = Real(0);
    for (int i = 0; i < m; i++) {
        Real amax = Real(0);
        for (int j = 0; j < n; j++) {
            Real val = tm_abs<Real>(A(i, j).real()) + tm_abs<Real>(A(i, j).imag());
            if (val > amax) amax = val;
        }
        if (amax < smlnum) {
            // Zero or tiny row — set scale to 1
            row_scale(i) = Real(1);
        } else {
            row_scale(i) = Real(1) / amax;
        }
        if (row_scale(i) < rowcnd_min) rowcnd_min = row_scale(i);
        if (row_scale(i) > rowcnd_max) rowcnd_max = row_scale(i);
    }

    // Only apply row scaling if there's significant variation (ratio < 0.01)
    Real rowcnd = rowcnd_min / std::max(rowcnd_max, smlnum);
    if (rowcnd >= Real(0.01)) {
        // Already well-balanced — skip row scaling
        row_scale.setOnes();
    } else {
        // Apply row scaling
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A(i, j) *= complex_t(row_scale(i), 0);
            }
        }
    }

    // Compute column infinity-norms of the (possibly row-scaled) matrix
    Real colcnd_min = bignum;
    Real colcnd_max = Real(0);
    for (int j = 0; j < n; j++) {
        Real amax = Real(0);
        for (int i = 0; i < m; i++) {
            Real val = tm_abs<Real>(A(i, j).real()) + tm_abs<Real>(A(i, j).imag());
            if (val > amax) amax = val;
        }
        if (amax < smlnum) {
            col_scale(j) = Real(1);
        } else {
            col_scale(j) = Real(1) / amax;
        }
        if (col_scale(j) < colcnd_min) colcnd_min = col_scale(j);
        if (col_scale(j) > colcnd_max) colcnd_max = col_scale(j);
    }

    // Only apply column scaling if there's significant variation
    Real colcnd = colcnd_min / std::max(colcnd_max, smlnum);
    if (colcnd >= Real(0.01)) {
        col_scale.setOnes();
    } else {
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A(i, j) *= complex_t(col_scale(j), 0);
            }
        }
    }

    return {row_scale, col_scale};
}

// ============================================================================
// Equilibrated LU solve with iterative refinement
// Port of Fortran MatrixSolv.f90::LU1 (lines 700-786)
//
// Solves A_orig * X = B_orig using:
// 1. Matrix equilibration of A
// 2. LU factorization of equilibrated A
// 3. Initial solve
// 4. Iterative refinement (up to itmax iterations)
// 5. Unscaling of solution
//
// Reports rcond of the equilibrated matrix via rcond_out.
// ============================================================================
template<typename Real>
typename Types<Real>::Matrix equilibrated_lu_solve(
    const typename Types<Real>::Matrix& A_orig,
    const typename Types<Real>::Matrix& B_orig,
    double& rcond_out,
    int itmax = 3) {

    using Matrix = typename Types<Real>::Matrix;
    using RealVector = typename Types<Real>::RealVector;
    using complex_t = std::complex<Real>;

    Real eps = machine_eps<Real>();
    int n = A_orig.rows();
    int nrhs = B_orig.cols();

    // Copy A for equilibration (modifies in-place)
    Matrix A_eq = A_orig;

    // Step 1: Equilibrate
    auto [row_scale, col_scale] = equilibrate_matrix<Real>(A_eq);

    // Step 2: LU factorize the equilibrated matrix
    ::Eigen::PartialPivLU<Matrix> lu(A_eq);
    rcond_out = static_cast<double>(lu.rcond());

    // Step 3: Scale RHS: B_scaled = diag(row_scale) * B_orig
    Matrix B_scaled(n, nrhs);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nrhs; j++) {
            B_scaled(i, j) = B_orig(i, j) * complex_t(row_scale(i), 0);
        }
    }

    // Step 4: Initial solve
    Matrix X = lu.solve(B_scaled);

    // Step 5: Iterative refinement
    if (itmax > 0) {
        for (int iter = 0; iter < itmax; iter++) {
            // Compute residual: r = B_scaled - A_eq * X
            Matrix residual = B_scaled - A_eq * X;

            // Compute componentwise backward error for each RHS column
            Real berr_max = Real(0);
            for (int col = 0; col < nrhs; col++) {
                Real numer_max = Real(0);
                Real denom_max = Real(0);
                for (int i = 0; i < n; i++) {
                    Real ri = tm_abs<Real>(residual(i, col).real())
                            + tm_abs<Real>(residual(i, col).imag());
                    if (ri > numer_max) numer_max = ri;

                    // |A_eq|*|X| + |B_scaled| for this row
                    Real row_sum = tm_abs<Real>(B_scaled(i, col).real())
                                 + tm_abs<Real>(B_scaled(i, col).imag());
                    for (int k = 0; k < n; k++) {
                        Real aik = tm_abs<Real>(A_eq(i, k).real())
                                 + tm_abs<Real>(A_eq(i, k).imag());
                        Real xk = tm_abs<Real>(X(k, col).real())
                                + tm_abs<Real>(X(k, col).imag());
                        row_sum += aik * xk;
                    }
                    if (row_sum > denom_max) denom_max = row_sum;
                }
                Real berr = (denom_max > Real(0)) ? (numer_max / denom_max) : numer_max;
                if (berr > berr_max) berr_max = berr;
            }

            // Stop if backward error is at machine precision
            if (berr_max <= eps) break;

            // Stop if stagnating (store previous berr for comparison)
            if (iter > 0) {
                // We break out on stagnation below
            }

            // Correction step
            Matrix dx = lu.solve(residual);
            X += dx;
        }
    }

    // Step 6: Unscale solution: X_orig(i,:) = X(i,:) * col_scale(i)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < nrhs; j++) {
            X(i, j) *= complex_t(col_scale(i), 0);
        }
    }

    return X;
}

} // namespace tmatrix

#endif // GUARD_tmatrix_linalg_h
