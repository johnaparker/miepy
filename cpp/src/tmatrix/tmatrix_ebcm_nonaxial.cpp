#include "tmatrix/tmatrix_ebcm_nonaxial.hpp"
#include "tmatrix/tmatrix_special.hpp"
#include "tmatrix/tmatrix_svwf.hpp"
#include <vector>
#include <cmath>

namespace tmatrix {

// ============================================================================
// 3D mixed product: (n x A) . B with n = (n_r, n_theta, n_phi)
// Generalizes the axisymmetric mixt_product which assumes n_phi = 0.
// ============================================================================
template<typename Real>
static std::complex<Real> mixt_product_3d(Real n_r, Real n_theta, Real n_phi,
                                          const std::complex<Real>* A,
                                          const std::complex<Real>* B) {
    using complex_t = std::complex<Real>;
    complex_t nr(n_r, 0), nt(n_theta, 0), np(n_phi, 0);
    // (n x A) = (nt*A[2] - np*A[1], np*A[0] - nr*A[2], nr*A[1] - nt*A[0])
    // (n x A).B = sum of component products
    return (nt * A[2] - np * A[1]) * B[0]
         + (np * A[0] - nr * A[2]) * B[1]
         + (nr * A[1] - nt * A[0]) * B[2];
}

// ============================================================================
// Determine the m-value for a given matrix index in the NFMDS mode ordering
//
// For plus=false: m=0 (Nrank), m=-1 (Nrank), m=+1 (Nrank), m=-2 (Nrank-1), ...
// For plus=true:  m=0 (Nrank), m=+1 (Nrank), m=-1 (Nrank), m=+2 (Nrank-1), ...
//
// Returns the actual m-value (signed) and the n-value for the given index.
// ============================================================================
template<typename Real>
static void index_to_mn(int idx, int Nrank, int Mrank, bool plus,
                        int& m_out, int& n_out) {
    if (idx < Nrank) {
        // m=0 block
        m_out = 0;
        n_out = idx + 1;
        return;
    }

    int remaining = idx - Nrank;
    for (int m = 1; m <= Mrank; m++) {
        int block_size = Nrank - m + 1;
        if (remaining < block_size) {
            // First sub-block
            m_out = plus ? m : -m;
            n_out = m + remaining;
            return;
        }
        remaining -= block_size;
        if (remaining < block_size) {
            // Second sub-block
            m_out = plus ? -m : m;
            n_out = m + remaining;
            return;
        }
        remaining -= block_size;
    }
    // Should not reach here
    m_out = 0;
    n_out = 1;
}

// ============================================================================
// Apply mirror and azimuthal symmetry post-processing
// Port of Proces1.f90::AzimMirSym (lines 453-614)
//
// Rows use plus=false ordering: m = 0, -1, +1, -2, +2, ...
// Cols use plus=true ordering:  m = 0, +1, -1, +2, -2, ...
//
// After the integration, for each matrix entry:
// - Mirror symmetry: multiply diagonal blocks by sp, off-diagonal by sm
// - Azimuthal N-fold symmetry: multiply all blocks by factor
// ============================================================================
template<typename Real>
static void apply_symmetry(typename Types<Real>::Matrix& A,
                           int Nmax, int Nrank, int Mrank,
                           bool mirror, int Nazimutsym) {
    if (Nazimutsym < 2 && !mirror) return;

    using complex_t = std::complex<Real>;
    const complex_t im_unit(0, 1);
    Real dalfa = (Nazimutsym >= 2) ? Real(2) * pi<Real>() / Real(Nazimutsym) : Real(0);

    for (int i = 0; i < Nmax; i++) {
        int ml, n;
        index_to_mn<Real>(i, Nrank, Mrank, false, ml, n);  // rows: plus=false

        for (int j = 0; j < Nmax; j++) {
            int m1l, n1;
            index_to_mn<Real>(j, Nrank, Mrank, true, m1l, n1);  // cols: plus=true

            int m_abs = std::abs(ml);
            int m1_abs = std::abs(m1l);

            if (mirror) {
                Real sp = ((n + n1 + m_abs + m1_abs) % 2 == 0) ? Real(2) : Real(0);
                Real sm = Real(2) - sp;

                A(i, j)             *= complex_t(sp, 0);
                A(i, j + Nmax)      *= complex_t(sm, 0);
                A(i + Nmax, j)      *= complex_t(sm, 0);
                A(i + Nmax, j + Nmax) *= complex_t(sp, 0);
            }

            if (Nazimutsym >= 2) {
                complex_t fact(Real(1), 0);
                for (int p = 2; p <= Nazimutsym; p++) {
                    Real arg = Real((ml + m1l) * (p - 1)) * dalfa;
                    fact += tm_exp<Real>(im_unit * complex_t(arg, 0));
                }

                A(i, j)             *= fact;
                A(i, j + Nmax)      *= fact;
                A(i + Nmax, j)      *= fact;
                A(i + Nmax, j + Nmax) *= fact;
            }
        }
    }
}

// ============================================================================
// Assemble Q-matrix for a non-axisymmetric particle
// Port of Proces1.f90::matrix_Q (lines 9-92) + MatrixQ.f90::matQ (lines 10-91)
//
// index1, index2: Bessel function type (1=regular, 3=outgoing)
//   Q31: index1=3, index2=1 -> sign = -1
//   Q11: index1=1, index2=1 -> sign = +1
//
// Returns [2*Nmax, 2*Nmax] matrix
// ============================================================================
template<typename Real>
static typename Types<Real>::Matrix assemble_Q_nonaxial(
    const NonAxialGeometry<Real>& geom,
    int index1, int index2,
    Real k_real, std::complex<Real> k_int,
    std::complex<Real> n_rel,
    int Mrank, int Nrank, int Nmax,
    int Nint1, int Nint2, bool mirror,
    bool conducting = false) {

    using complex_t = std::complex<Real>;
    using Matrix = typename Types<Real>::Matrix;

    Matrix A = Matrix::Zero(2 * Nmax, 2 * Nmax);

    // Pre-factor: sign * i * k^2 / pi
    Real sign_val = Real(1);
    if (index1 == 3 && index2 == 1) sign_val = Real(-1);
    complex_t f = complex_t(sign_val, 0) * complex_t(0, 1)
                * complex_t(k_real * k_real, 0)
                / complex_t(pi<Real>(), 0);

    // Get quadrature points
    std::vector<std::vector<Real>> param1G, param2G, weightsG;
    geom.quadrature_points(Nint1, Nint2, mirror, param1G, param2G, weightsG);

    int Nparam = static_cast<int>(param1G.size());

    // Count total quadrature points for OpenMP scheduling
    int total_points = 0;
    std::vector<int> patch_offsets(Nparam + 1, 0);
    for (int ip = 0; ip < Nparam; ip++) {
        patch_offsets[ip] = total_points;
        total_points += static_cast<int>(param1G[ip].size());
    }
    patch_offsets[Nparam] = total_points;

    // Parallel accumulation
    #pragma omp parallel
    {
        Matrix A_local = Matrix::Zero(2 * Nmax, 2 * Nmax);
        std::vector<std::vector<complex_t>> mv1, nv1, mv3, nv3;

        #pragma omp for schedule(dynamic)
        for (int flat = 0; flat < total_points; flat++) {
            // Find patch and point index from flat index
            int iparam = 0;
            for (int ip = 0; ip < Nparam; ip++) {
                if (flat < patch_offsets[ip + 1]) {
                    iparam = ip;
                    break;
                }
            }
            int pint = flat - patch_offsets[iparam];

            Real p1 = param1G[iparam][pint];
            Real p2 = param2G[iparam][pint];
            Real weight = weightsG[iparam][pint];

            GeometryPoint3D<Real> pt = geom.evaluate(p1, p2, iparam);
            Real r = pt.r;
            Real theta = pt.theta;
            Real phi = pt.phi;
            Real dA = pt.dA;
            Real nuv_r = pt.n_r;
            Real nuv_theta = pt.n_theta;
            Real nuv_phi = pt.n_phi;

            complex_t zl(k_real * r, 0);       // k_medium * r (real)
            complex_t zc = k_int * complex_t(r, 0);  // k_interior * r

            // Outer SVWFs: plus=false (rows: m = 0,-1,+1,-2,+2,...), sym=false
            svwf_complete<Real>(index1, zl, theta, phi,
                               Mrank, Nrank, Nmax, false, false, mv3, nv3);

            // Inner SVWFs: plus=true (cols: m = 0,+1,-1,+2,-2,...), sym=false
            svwf_complete<Real>(index2, zc, theta, phi,
                               Mrank, Nrank, Nmax, true, false, mv1, nv1);

            complex_t fact = f * complex_t(dA * weight, 0);

            // Accumulate Q-matrix via mixed products (n x A) . B
            for (int i = 0; i < Nmax; i++) {
                complex_t mvl[3] = {mv3[0][i], mv3[1][i], mv3[2][i]};
                complex_t nvl[3] = {nv3[0][i], nv3[1][i], nv3[2][i]};

                for (int j = 0; j < Nmax; j++) {
                    complex_t mvc[3] = {mv1[0][j], mv1[1][j], mv1[2][j]};
                    complex_t nvc[3] = {nv1[0][j], nv1[1][j], nv1[2][j]};

                    // Block (0,0): magnetic-magnetic
                    complex_t v1 = mixt_product_3d<Real>(nuv_r, nuv_theta, nuv_phi, mvc, nvl);
                    complex_t v2 = mixt_product_3d<Real>(nuv_r, nuv_theta, nuv_phi, nvc, mvl);
                    if (!conducting) {
                        A_local(i, j) += (v1 + n_rel * v2) * fact;
                    } else {
                        A_local(i, j) += v2 * fact;
                    }

                    // Block (0,1): magnetic-electric
                    v1 = mixt_product_3d<Real>(nuv_r, nuv_theta, nuv_phi, nvc, nvl);
                    v2 = mixt_product_3d<Real>(nuv_r, nuv_theta, nuv_phi, mvc, mvl);
                    if (!conducting) {
                        A_local(i, j + Nmax) += (v1 + n_rel * v2) * fact;
                    } else {
                        A_local(i, j + Nmax) += v2 * fact;
                    }

                    // Block (1,0): electric-magnetic
                    v1 = mixt_product_3d<Real>(nuv_r, nuv_theta, nuv_phi, mvc, mvl);
                    v2 = mixt_product_3d<Real>(nuv_r, nuv_theta, nuv_phi, nvc, nvl);
                    if (!conducting) {
                        A_local(i + Nmax, j) += (v1 + n_rel * v2) * fact;
                    } else {
                        A_local(i + Nmax, j) += v2 * fact;
                    }

                    // Block (1,1): electric-electric
                    v1 = mixt_product_3d<Real>(nuv_r, nuv_theta, nuv_phi, nvc, mvl);
                    v2 = mixt_product_3d<Real>(nuv_r, nuv_theta, nuv_phi, mvc, nvl);
                    if (!conducting) {
                        A_local(i + Nmax, j + Nmax) += (v1 + n_rel * v2) * fact;
                    } else {
                        A_local(i + Nmax, j + Nmax) += v2 * fact;
                    }
                }
            }
        }

        #pragma omp critical
        A += A_local;
    }

    // Apply mirror and azimuthal symmetry post-processing
    apply_symmetry<Real>(A, Nmax, Nrank, Mrank,
                         mirror && geom.is_mirror_symmetric(),
                         geom.azimuthal_symmetry());

    return A;
}

// ============================================================================
// Main entry point: compute full [2, rmax, 2, rmax] T-matrix
// for a non-axisymmetric particle.
//
// Algorithm (port of TNONAXSYM.f90::TMatrix_Nrank_MrankNONAXSYM):
// 1. Assemble Q31, apply symmetry
// 2. Assemble Q11, apply symmetry
// 3. T_nfmds = Q11 * Q31^{-1}
// 4. Map to MiePy convention
// ============================================================================
template<typename Real>
TmatrixResult compute_nonaxisymmetric_tmatrix(
    const NonAxialGeometry<Real>& geom,
    double k_double, std::complex<double> n_rel_double,
    int lmax, int Nint1, int Nint2,
    bool conducting) {

    using complex_t = std::complex<Real>;
    using Matrix = typename Types<Real>::Matrix;

    int rmax = lmax * (lmax + 2);
    int Nrank = lmax;
    int Mrank = lmax;
    int Nmax = Nrank + Mrank * (2 * Nrank - Mrank + 1);
    // Nmax should equal rmax when Mrank = Nrank = lmax

    Real k_real = Real(k_double);
    complex_t n_rel(Real(n_rel_double.real()), Real(n_rel_double.imag()));
    complex_t k_int = complex_t(k_real, 0) * n_rel;

    bool mirror = geom.is_mirror_symmetric();

    // Assemble Q31 (index1=3, index2=1)
    Matrix Q31 = assemble_Q_nonaxial<Real>(geom, 3, 1, k_real, k_int, n_rel,
                                            Mrank, Nrank, Nmax, Nint1, Nint2,
                                            mirror, conducting);

    // Assemble Q11 (index1=1, index2=1)
    Matrix Q11 = assemble_Q_nonaxial<Real>(geom, 1, 1, k_real, k_int, n_rel,
                                            Mrank, Nrank, Nmax, Nint1, Nint2,
                                            mirror, conducting);

    // T_nfmds = Q11 * Q31^{-1} via LU decomposition
    // Same technique as axisymmetric non-DS path:
    // lu(Q31^T).solve(Q11^T).transpose() = Q11 * Q31^{-1}
    ::Eigen::PartialPivLU<Matrix> lu(Q31.transpose());
    double rcond = static_cast<double>(lu.rcond());
    Matrix T_nfmds = lu.solve(Q11.transpose()).transpose();

    // ========================================================================
    // Convention mapping: NFMDS -> MiePy [2, rmax, 2, rmax]
    //
    // Both rows and columns of T_nfmds are in plus=false ordering
    // (0, -1, +1, -2, +2, ...) because inner indices cancel in Q11 * Q31^{-1}
    // and both outer bases use plus=false.
    //
    // Port of get_tmatrix.py lines 132-174
    // ========================================================================
    std::vector<std::complex<double>> T(2 * rmax * 2 * rmax, std::complex<double>(0, 0));

    auto Tidx = [rmax](int a1, int r1, int a2, int r2) -> int {
        return ((a1 * rmax + r1) * 2 + a2) * rmax + r2;
    };

    // MiePy mode index: r = n^2 + n + m - 1 (0-based)
    auto mode_r = [](int n, int m) -> int {
        return n * n + n + m - 1;
    };

    // Generate the plus=false mode sequence for both rows and columns
    // Sequence: (m=0, n=1..Nrank), (m=-1, n=1..Nrank), (m=+1, n=1..Nrank),
    //           (m=-2, n=2..Nrank), (m=+2, n=2..Nrank), ...
    struct ModeInfo { int m; int n; };
    std::vector<ModeInfo> modes(Nmax);
    {
        int mp = 0, np = 1;
        for (int i = 0; i < Nmax; i++) {
            if (np > Nrank) {
                if (mp >= 0) {
                    mp = -(mp + 1);
                } else {
                    mp = -mp;
                }
                np = std::abs(mp);
                // For m=0, n starts at 1; for |m|>0, n starts at |m|
                if (mp == 0) np = 1;
            }
            modes[i] = {mp, np};
            np++;
        }
    }

    for (int i = 0; i < Nmax; i++) {
        int m1p = modes[i].m;
        int n1p = modes[i].n;
        int r1 = mode_r(n1p, m1p);

        for (int j = 0; j < Nmax; j++) {
            int m2p = modes[j].m;
            int n2p = modes[j].n;
            int r2 = mode_r(n2p, m2p);

            // Phase factor: -(i^(n2-n1))
            int dn = n2p - n1p;
            std::complex<double> factor;
            int dn_mod = ((dn % 4) + 4) % 4;
            switch (dn_mod) {
                case 0: factor = std::complex<double>(-1, 0); break;
                case 1: factor = std::complex<double>(0, -1); break;
                case 2: factor = std::complex<double>(1, 0); break;
                case 3: factor = std::complex<double>(0, 1); break;
            }

            // Sign corrections for odd m pairs
            double f1 = 1.0, f2 = 1.0;
            if (m1p % 2 != 0 && m2p % 2 != 0) {
                // Both m values are odd
                double sign_prod = (m1p * m2p > 0) ? 1.0 : -1.0;
                f1 = std::pow(sign_prod, double(n1p + n2p + 1));
                f2 = std::pow(sign_prod, double(n1p + n2p));
            }

            // Map from NFMDS blocks to MiePy convention
            std::complex<double> t_mm(double(T_nfmds(i, j).real()),
                                      double(T_nfmds(i, j).imag()));
            std::complex<double> t_ee(double(T_nfmds(i + Nmax, j + Nmax).real()),
                                      double(T_nfmds(i + Nmax, j + Nmax).imag()));
            std::complex<double> t_em(double(T_nfmds(i + Nmax, j).real()),
                                      double(T_nfmds(i + Nmax, j).imag()));
            std::complex<double> t_me(double(T_nfmds(i, j + Nmax).real()),
                                      double(T_nfmds(i, j + Nmax).imag()));

            T[Tidx(0, r1, 0, r2)] = t_ee * factor * f1;    // elec-elec
            T[Tidx(1, r1, 1, r2)] = t_mm * factor * f1;    // mag-mag
            T[Tidx(0, r1, 1, r2)] = -t_em * factor * f2;   // elec-mag
            T[Tidx(1, r1, 0, r2)] = -t_me * factor * f2;   // mag-elec
        }
    }

    return {std::move(T), rcond, 0};
}

// Explicit instantiations
template TmatrixResult compute_nonaxisymmetric_tmatrix<double>(
    const NonAxialGeometry<double>&, double, std::complex<double>, int, int, int, bool);

#if MIEPY_HAS_QUAD
template TmatrixResult compute_nonaxisymmetric_tmatrix<__float128>(
    const NonAxialGeometry<__float128>&, double, std::complex<double>, int, int, int, bool);
#endif

} // namespace tmatrix
