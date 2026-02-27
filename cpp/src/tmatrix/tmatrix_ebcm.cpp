#include "tmatrix/tmatrix_ebcm.hpp"
#include "tmatrix/tmatrix_special.hpp"
#include "tmatrix/tmatrix_svwf.hpp"
#include <vector>
#include <cmath>

namespace tmatrix {

// Mixed product: (n x A) . B = (n2*A3 - n3*A2)*B1 + (n3*A1 - n1*A3)*B2 + (n1*A2 - n2*A1)*B3
// In our case n3 (n_phi) = 0 for axisymmetric particles
template<typename Real>
static std::complex<Real> mixt_product(Real n_r, Real n_theta,
                                       const std::complex<Real>* A,
                                       const std::complex<Real>* B) {
    // n = (n_r, n_theta, 0), A = (A_r, A_theta, A_phi), B = (B_r, B_theta, B_phi)
    // (n x A) = (n_theta*A_phi, -n_r*A_phi, n_r*A_theta - n_theta*A_r)
    // (n x A).B = n_theta*A_phi*B_r - n_r*A_phi*B_theta + (n_r*A_theta - n_theta*A_r)*B_phi
    using complex_t = std::complex<Real>;
    complex_t nr(n_r, 0), nt(n_theta, 0);
    return nt * A[2] * B[0] - nr * A[2] * B[1] + (nr * A[1] - nt * A[0]) * B[2];
}

// ============================================================================
// Assemble Q-matrix for a single azimuthal mode m
// Port of Proces1.f90::matrix_Q_m + MatrixQ.f90::matQ_m
//
// index1, index2: Bessel function type (1=regular, 3=outgoing)
//   Q31: index1=3 (outer=outgoing), index2=1 (inner=regular) -> sign = -1
//   Q11: index1=1 (outer=regular), index2=1 (inner=regular) -> sign = +1
//
// Returns 2D matrix of size [2*NmaxL, 2*NmaxC]
// Layout: [[Q_mm, Q_mn], [Q_nm, Q_nn]] where m=magnetic, n=electric
// ============================================================================
template<typename Real>
static typename Types<Real>::Matrix assemble_Q_m(
    const AxialGeometry<Real>& geom,
    int index1, int index2,
    Real k_real, std::complex<Real> k_int,
    std::complex<Real> n_rel,
    int m, int Nrank, int Nmax, int Nint,
    bool mirror, bool use_ds,
    const std::vector<Real>& zRe, const std::vector<Real>& zIm,
    bool conducting = false) {

    using complex_t = std::complex<Real>;
    using Matrix = typename Types<Real>::Matrix;

    // Determine dimensions based on DS mode
    // Both Q31 and Q11 use localized outer (NmaxL = Nmax).
    // DS only changes the inner (column) basis to distributed sources.
    // We truncate to the first Nmax DS sources so both matrices are square
    // [2*Nmax, 2*Nmax]. The DS basis cancels out in T = Q11 * Q31^{-1}.
    int NmaxL = Nmax;
    int NmaxC = Nmax;

    Matrix A = Matrix::Zero(2 * NmaxL, 2 * NmaxC);

    // Pre-factor: sign * i * 2 * k^2
    Real sign_val = Real(1);
    if (index1 == 3 && index2 == 1) sign_val = Real(-1);
    complex_t f = complex_t(sign_val, 0) * complex_t(0, 1) * complex_t(Real(2), 0)
                  * complex_t(k_real * k_real, 0);

    // Get quadrature points
    std::vector<std::vector<Real>> paramG, weightsG;
    geom.quadrature_points(Nint, mirror, paramG, weightsG);

    // Temporary SVWF storage
    std::vector<std::vector<complex_t>> mv1, nv1, mv3, nv3;

    int Nparam = geom.num_curves();
    for (int iparam = 0; iparam < Nparam; iparam++) {
        int Nintl = static_cast<int>(paramG[iparam].size());
        for (int pint = 0; pint < Nintl; pint++) {
            Real param = paramG[iparam][pint];
            Real weight = weightsG[iparam][pint];

            GeometryPoint<Real> pt = geom.evaluate(param, iparam);
            Real r = pt.r;
            Real theta = pt.theta;
            Real dA = pt.dA;
            Real nuv_r = pt.n_r;
            Real nuv_theta = pt.n_theta;

            complex_t zl(k_real * r, 0);    // k_medium * r (real)
            complex_t zc = k_int * complex_t(r, 0);  // k_interior * r

            int m_minus = -m;

            // Compute outer VSWFs (index1, m_minus, medium wavenumber)
            // Always localized — DS only affects the inner (column) basis
            svwf_localized<Real>(index1, zl, theta, m_minus, Nrank, NmaxL, mv3, nv3);

            // Compute inner VSWFs (index2, m, interior wavenumber)
            if (!use_ds) {
                svwf_localized<Real>(index2, zc, theta, m, Nrank, NmaxC, mv1, nv1);
            } else {
                // DS mode: inner uses distributed sources with interior wavenumber
                if (index2 == 1 && index1 == 3) {
                    svwf_distributed<Real>(1, k_int, r, theta, zRe, zIm, m, Nrank, mv1, nv1);
                } else if (index2 == 1 && index1 == 1) {
                    // Q11 with DS: outer=localized regular, inner=distributed regular
                    svwf_distributed<Real>(1, k_int, r, theta, zRe, zIm, m, Nrank, mv1, nv1);
                } else {
                    svwf_localized<Real>(index2, zc, theta, m, Nrank, NmaxC, mv1, nv1);
                }
            }

            complex_t fact = f * complex_t(dA * weight, 0);

            // Assemble Q-matrix via mixed products (n x A) . B
            // Port of MatrixQ.f90::matQ_m
            // For PEC (conducting): only v2 terms, no ind_ref multiplication
            for (int i = 0; i < NmaxL; i++) {
                int ni = (std::abs(m) == 0) ? (i + 1) : (std::abs(m) + i);

                // Extract outer SVWF vectors for row i
                complex_t mvl[3] = {mv3[0][i], mv3[1][i], mv3[2][i]};
                complex_t nvl[3] = {nv3[0][i], nv3[1][i], nv3[2][i]};

                for (int j = 0; j < NmaxC; j++) {
                    int nj = (std::abs(m) == 0) ? (j + 1) : (std::abs(m) + j);

                    complex_t mvc[3] = {mv1[0][j], mv1[1][j], mv1[2][j]};
                    complex_t nvc[3] = {nv1[0][j], nv1[1][j], nv1[2][j]};

                    Real s = ((ni + nj) % 2 == 0) ? Real(1) : Real(-1);

                    // Block (0,0): magnetic-magnetic
                    complex_t v1 = mixt_product<Real>(nuv_r, nuv_theta, mvc, nvl);
                    complex_t v2 = mixt_product<Real>(nuv_r, nuv_theta, nvc, mvl);
                    if (!conducting) {
                        A(i, j) += (v1 + n_rel * v2) * fact;
                        if (mirror) A(i, j) += complex_t(s, 0) * (v1 + n_rel * v2) * fact;
                    } else {
                        A(i, j) += v2 * fact;
                        if (mirror) A(i, j) += complex_t(s, 0) * v2 * fact;
                    }

                    if (m != 0) {
                        // Block (0,1): magnetic-electric
                        v1 = mixt_product<Real>(nuv_r, nuv_theta, nvc, nvl);
                        v2 = mixt_product<Real>(nuv_r, nuv_theta, mvc, mvl);
                        if (!conducting) {
                            A(i, j + NmaxC) += (v1 + n_rel * v2) * fact;
                            if (mirror) A(i, j + NmaxC) -= complex_t(s, 0) * (v1 + n_rel * v2) * fact;
                        } else {
                            A(i, j + NmaxC) += v2 * fact;
                            if (mirror) A(i, j + NmaxC) -= complex_t(s, 0) * v2 * fact;
                        }

                        // Block (1,0): electric-magnetic
                        v1 = mixt_product<Real>(nuv_r, nuv_theta, mvc, mvl);
                        v2 = mixt_product<Real>(nuv_r, nuv_theta, nvc, nvl);
                        if (!conducting) {
                            A(i + NmaxL, j) += (v1 + n_rel * v2) * fact;
                            if (mirror) A(i + NmaxL, j) -= complex_t(s, 0) * (v1 + n_rel * v2) * fact;
                        } else {
                            A(i + NmaxL, j) += v2 * fact;
                            if (mirror) A(i + NmaxL, j) -= complex_t(s, 0) * v2 * fact;
                        }
                    }

                    // Block (1,1): electric-electric
                    v1 = mixt_product<Real>(nuv_r, nuv_theta, nvc, mvl);
                    v2 = mixt_product<Real>(nuv_r, nuv_theta, mvc, nvl);
                    if (!conducting) {
                        A(i + NmaxL, j + NmaxC) += (v1 + n_rel * v2) * fact;
                        if (mirror) A(i + NmaxL, j + NmaxC) += complex_t(s, 0) * (v1 + n_rel * v2) * fact;
                    } else {
                        A(i + NmaxL, j + NmaxC) += v2 * fact;
                        if (mirror) A(i + NmaxL, j + NmaxC) += complex_t(s, 0) * v2 * fact;
                    }
                }
            }
        }
    }

    return A;
}

// ============================================================================
// Compute T-matrix for a single azimuthal mode m
// T_m = -Q11_m * Q31_m^{-1}
// ============================================================================
template<typename Real>
static typename Types<Real>::Matrix tmatrix_m(
    const AxialGeometry<Real>& geom,
    Real k_real, std::complex<Real> k_int,
    std::complex<Real> n_rel,
    int m, int Nrank, int Nmax, int Nint,
    bool mirror, bool use_ds,
    const std::vector<Real>& zRe, const std::vector<Real>& zIm,
    bool conducting = false) {

    using Matrix = typename Types<Real>::Matrix;

    // Q31: outgoing(3) x regular(1) in interior medium
    Matrix Q31 = assemble_Q_m<Real>(geom, 3, 1, k_real, k_int, n_rel,
                                     m, Nrank, Nmax, Nint, mirror, use_ds, zRe, zIm, conducting);

    // Q11: regular(1) x regular(1) in interior medium
    Matrix Q11 = assemble_Q_m<Real>(geom, 1, 1, k_real, k_int, n_rel,
                                     m, Nrank, Nmax, Nint, mirror, use_ds, zRe, zIm, conducting);

    // T_m = Q11 * Q31^{-1}
    // The negation is incorporated into the phase factor during mapping.
    // Both matrices are square [2*Nmax, 2*Nmax] (DS columns truncated to Nmax).
    ::Eigen::PartialPivLU<Matrix> lu(Q31);
    Matrix T_m = Q11 * lu.inverse();

    return T_m;
}

// ============================================================================
// Main entry point: compute full [2, rmax, 2, rmax] T-matrix
// ============================================================================
template<typename Real>
std::vector<std::complex<double>> compute_axisymmetric_tmatrix(
    const AxialGeometry<Real>& geom,
    double k_double, std::complex<double> n_rel_double,
    int lmax, int Nint,
    bool use_ds, bool complex_plane, double eps_z_double,
    bool conducting) {

    using complex_t = std::complex<Real>;

    int rmax = lmax * (lmax + 2);
    int Nrank = lmax;
    bool mirror = geom.is_mirror_symmetric();

    Real k_real = Real(k_double);
    complex_t n_rel(Real(n_rel_double.real()), Real(n_rel_double.imag()));
    complex_t k_int = complex_t(k_real, 0) * n_rel;
    Real eps_z = Real(eps_z_double);

    bool use_ds_internal = use_ds;
    std::vector<Real> zRe, zIm;
    if (use_ds_internal) {
        geom.distributed_sources(Nrank, complex_plane, eps_z, zRe, zIm);
        // DS breaks mirror symmetry exploitation: the sign factors for
        // localized SVWFs under theta -> pi-theta don't apply to DS SVWFs
        // (each source has a distinct position, not related by simple parity).
        mirror = false;
    }

    // Output T-matrix: [2, rmax, 2, rmax] flattened in row-major
    std::vector<std::complex<double>> T(2 * rmax * 2 * rmax, std::complex<double>(0, 0));

    // Helper to index into [2, rmax, 2, rmax]
    auto Tidx = [rmax](int a1, int r1, int a2, int r2) -> int {
        return ((a1 * rmax + r1) * 2 + a2) * rmax + r2;
    };

    // Mode index: r = n*(n+2) + m (0-based index for mode (n,m))
    // where n = 1..lmax, m = -n..n
    auto rindex = [](int n, int m) -> int {
        return n * (n + 2) + m - 1;  // 0-based: (n-1)*(n+2-1) + m + something
        // Actually: r = n*(n+2) - n + m - 1 = n^2 + n + m - 1
        // Wait, let's match miepy convention: r = (n-1)*(n+1) + m + n = n^2 + m - 1 + n
        // miepy mode_indices: for n,m -> r = n*(n+1) + m - 1  (0-based)
        // No: miepy uses r = n*(n+2) + m for 1-based... let me check
    };
    (void)rindex;

    // Use miepy convention: for (n, m), 0-based index r = n*n + n + m - 1
    // This is n*(n+1) + m - 1 (since sum of 2k+1 for k=1..n-1 = n^2-1, plus m offset)
    // Actually from the code: rmax = lmax*(lmax+2), and mode_indices gives
    // r values 0..rmax-1 for n=1..lmax, m=-n..n
    // r = n*n + n + m - 1 = (n-1)*(n+1) + m + n - 1... let me just use n^2 - 1 + n + m
    // r = n^2 + n + m - 1, but that gives r(1,-1)=0, r(1,0)=1, r(1,1)=2 -- yes!
    auto mode_r = [](int n, int m) -> int {
        return n * n + n + m - 1;
    };

    // Loop over azimuthal modes m = 0, 1, ..., lmax
    for (int m = 0; m <= lmax; m++) {
        int Nmax = (m == 0) ? Nrank : (Nrank - m + 1);
        if (Nmax <= 0) continue;

        // Compute T-matrix for this m
        auto T_m = tmatrix_m<Real>(geom, k_real, k_int, n_rel,
                                    m, Nrank, Nmax, Nint, mirror, use_ds_internal, zRe, zIm,
                                    conducting);

        // Map per-m block into full [2, rmax, 2, rmax] T-matrix
        // T_m is always [2*Nmax, 2*Nmax] (pseudo-inverse collapses DS columns)

        for (int i = 0; i < Nmax; i++) {
            int n1 = (m == 0) ? (i + 1) : (m + i);

            for (int j = 0; j < Nmax; j++) {
                int n2 = (m == 0) ? (j + 1) : (m + j);

                // Phase factor from NFMDS to MiePy convention
                // factor = -(i^(n2-n1))
                int dn = n2 - n1;
                std::complex<double> factor;
                int dn_mod = ((dn % 4) + 4) % 4;
                switch (dn_mod) {
                    case 0: factor = std::complex<double>(-1, 0); break;
                    case 1: factor = std::complex<double>(0, -1); break;
                    case 2: factor = std::complex<double>(1, 0); break;
                    case 3: factor = std::complex<double>(0, 1); break;
                }

                // For m >= 0
                if (m == 0) {
                    int r1 = mode_r(n1, 0);
                    int r2 = mode_r(n2, 0);

                    // Magnetic-magnetic: T_m[i, j]
                    T[Tidx(1, r1, 1, r2)] = double(T_m(i, j).real()) * factor
                                           + std::complex<double>(0, 1) * double(T_m(i, j).imag()) * factor;
                    // Electric-electric: T_m[i+Nmax, j+Nmax]
                    T[Tidx(0, r1, 0, r2)] = double(T_m(i + Nmax, j + Nmax).real()) * factor
                                           + std::complex<double>(0, 1) * double(T_m(i + Nmax, j + Nmax).imag()) * factor;
                } else {
                    // Positive m
                    {
                        int r1 = mode_r(n1, m);
                        int r2 = mode_r(n2, m);

                        std::complex<double> t_mm(double(T_m(i, j).real()), double(T_m(i, j).imag()));
                        std::complex<double> t_ee(double(T_m(i + Nmax, j + Nmax).real()), double(T_m(i + Nmax, j + Nmax).imag()));
                        std::complex<double> t_em(double(T_m(i + Nmax, j).real()), double(T_m(i + Nmax, j).imag()));
                        std::complex<double> t_me(double(T_m(i, j + Nmax).real()), double(T_m(i, j + Nmax).imag()));

                        T[Tidx(1, r1, 1, r2)] = t_mm * factor;
                        T[Tidx(0, r1, 0, r2)] = t_ee * factor;
                        T[Tidx(0, r1, 1, r2)] = t_em * factor * double(m > 0 ? 1 : -1);
                        T[Tidx(1, r1, 0, r2)] = t_me * factor * double(m > 0 ? 1 : -1);
                    }

                    // Negative m: use symmetry relation
                    // For m -> -m: sign factors from parity
                    {
                        int r1 = mode_r(n1, -m);
                        int r2 = mode_r(n2, -m);

                        std::complex<double> t_mm(double(T_m(i, j).real()), double(T_m(i, j).imag()));
                        std::complex<double> t_ee(double(T_m(i + Nmax, j + Nmax).real()), double(T_m(i + Nmax, j + Nmax).imag()));
                        std::complex<double> t_em(double(T_m(i + Nmax, j).real()), double(T_m(i + Nmax, j).imag()));
                        std::complex<double> t_me(double(T_m(i, j + Nmax).real()), double(T_m(i, j + Nmax).imag()));

                        // For negative m: diagonal blocks same, off-diagonal flip sign
                        T[Tidx(1, r1, 1, r2)] = t_mm * factor;
                        T[Tidx(0, r1, 0, r2)] = t_ee * factor;
                        T[Tidx(0, r1, 1, r2)] = -t_em * factor * double(m > 0 ? 1 : -1);
                        T[Tidx(1, r1, 0, r2)] = -t_me * factor * double(m > 0 ? 1 : -1);
                    }
                }
            }
        }
    }

    return T;
}

// Explicit instantiations
template std::vector<std::complex<double>> compute_axisymmetric_tmatrix<double>(
    const AxialGeometry<double>&, double, std::complex<double>, int, int, bool, bool, double, bool);

#if MIEPY_HAS_QUAD
template std::vector<std::complex<double>> compute_axisymmetric_tmatrix<__float128>(
    const AxialGeometry<__float128>&, double, std::complex<double>, int, int, bool, bool, double, bool);
#endif

} // namespace tmatrix
