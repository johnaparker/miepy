#include "tmatrix/tmatrix_ebcm.hpp"
#include "tmatrix/tmatrix_special.hpp"
#include "tmatrix/tmatrix_svwf.hpp"
#include "tmatrix/tmatrix_linalg.hpp"
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
    // In DS mode:
    //   NmaxC = Nrank (inner always uses Nrank DS sources)
    //   Q31 (index1=3): NmaxL = Nrank (outer also DS, Nrank sources)
    //   Q11 (index1=1): NmaxL = Nmax  (outer localized, Nmax modes)
    // In non-DS mode: NmaxL = NmaxC = Nmax
    int NmaxL, NmaxC;
    if (use_ds) {
        NmaxC = Nrank;
        NmaxL = (index1 == 3) ? Nrank : Nmax;
    } else {
        NmaxL = Nmax;
        NmaxC = Nmax;
    }

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
            if (use_ds && index1 == 3) {
                // Q31 DS: distributed outgoing for outer, medium wavenumber
                complex_t kc(k_real, 0);
                svwf_distributed<Real>(3, kc, r, theta, zRe, zIm, m_minus, Nrank, mv3, nv3);
            } else {
                // Q11 or non-DS: localized for outer
                svwf_localized<Real>(index1, zl, theta, m_minus, Nrank, NmaxL, mv3, nv3);
            }

            // Compute inner VSWFs (index2, m, interior wavenumber)
            if (use_ds) {
                // DS mode: inner basis uses distributed sources with interior wavenumber
                svwf_distributed<Real>(1, k_int, r, theta, zRe, zIm, m, Nrank, mv1, nv1);
            } else {
                svwf_localized<Real>(index2, zc, theta, m, Nrank, NmaxC, mv1, nv1);
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
// Assemble incident matrix for a single azimuthal mode m
// Maps localized basis -> DS outer basis. Used in DS T-matrix extraction.
// Port of Proces1.f90::incident_matrix_m + MatrixQ.f90::matQinc_m
//
// Dimensions: [2*Nrank rows, 2*Nmax columns]
// Rows: DS outgoing (Hankel), medium wavenumber, -m
// Columns: localized regular (Bessel), medium wavenumber, +m
//
// Key differences from assemble_Q_m:
//   - No n_rel multiplication (just v1 + v2)
//   - No mirror symmetry
//   - Inner uses medium wavenumber (not interior)
//   - Outer always uses distributed outgoing (medium wavenumber)
// ============================================================================
template<typename Real>
static typename Types<Real>::Matrix assemble_incident_m(
    const AxialGeometry<Real>& geom,
    Real k_real,
    int m, int Nrank, int Nmax, int Nint,
    const std::vector<Real>& zRe, const std::vector<Real>& zIm) {

    using complex_t = std::complex<Real>;
    using Matrix = typename Types<Real>::Matrix;

    Matrix A = Matrix::Zero(2 * Nrank, 2 * Nmax);

    // Pre-factor: -i * 2 * k^2 (same sign as Q31: sign = -1)
    complex_t f = complex_t(Real(-1), 0) * complex_t(0, 1) * complex_t(Real(2), 0)
                  * complex_t(k_real * k_real, 0);

    // Get quadrature points (no mirror symmetry)
    std::vector<std::vector<Real>> paramG, weightsG;
    geom.quadrature_points(Nint, false, paramG, weightsG);

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

            complex_t zl(k_real * r, 0);  // k_medium * r

            int m_minus = -m;

            // Outer SVWFs (rows): DS outgoing Hankel, medium wavenumber, -m
            complex_t kc(k_real, 0);
            svwf_distributed<Real>(3, kc, r, theta, zRe, zIm, m_minus, Nrank, mv3, nv3);

            // Inner SVWFs (cols): localized regular Bessel, medium wavenumber, +m
            svwf_localized<Real>(1, zl, theta, m, Nrank, Nmax, mv1, nv1);

            complex_t fact = f * complex_t(dA * weight, 0);

            // Assemble via mixed products — NO n_rel, NO mirror
            for (int i = 0; i < Nrank; i++) {
                complex_t mvl[3] = {mv3[0][i], mv3[1][i], mv3[2][i]};
                complex_t nvl[3] = {nv3[0][i], nv3[1][i], nv3[2][i]};

                for (int j = 0; j < Nmax; j++) {
                    complex_t mvc[3] = {mv1[0][j], mv1[1][j], mv1[2][j]};
                    complex_t nvc[3] = {nv1[0][j], nv1[1][j], nv1[2][j]};

                    // Block (0,0): magnetic-magnetic
                    complex_t v1 = mixt_product<Real>(nuv_r, nuv_theta, mvc, nvl);
                    complex_t v2 = mixt_product<Real>(nuv_r, nuv_theta, nvc, mvl);
                    A(i, j) += (v1 + v2) * fact;

                    if (m != 0) {
                        // Block (0,1): magnetic-electric
                        v1 = mixt_product<Real>(nuv_r, nuv_theta, nvc, nvl);
                        v2 = mixt_product<Real>(nuv_r, nuv_theta, mvc, mvl);
                        A(i, j + Nmax) += (v1 + v2) * fact;

                        // Block (1,0): electric-magnetic
                        v1 = mixt_product<Real>(nuv_r, nuv_theta, mvc, mvl);
                        v2 = mixt_product<Real>(nuv_r, nuv_theta, nvc, nvl);
                        A(i + Nrank, j) += (v1 + v2) * fact;
                    }

                    // Block (1,1): electric-electric
                    v1 = mixt_product<Real>(nuv_r, nuv_theta, nvc, mvl);
                    v2 = mixt_product<Real>(nuv_r, nuv_theta, mvc, nvl);
                    A(i + Nrank, j + Nmax) += (v1 + v2) * fact;
                }
            }
        }
    }

    return A;
}

// ============================================================================
// Compute T-matrix for a single azimuthal mode m
// Non-DS: T_m = Q11 * Q31^{-1}
// DS:     T_m = Q11 * Q31^{-1} * Inc
//
// Returns pair of (T-matrix block, rcond of Q31)
// ============================================================================
template<typename Real>
static std::pair<typename Types<Real>::Matrix, double> tmatrix_m(
    const AxialGeometry<Real>& geom,
    Real k_real, std::complex<Real> k_int,
    std::complex<Real> n_rel,
    int m, int Nrank, int Nmax, int Nint,
    bool mirror, bool use_ds,
    const std::vector<Real>& zRe, const std::vector<Real>& zIm,
    bool conducting = false) {

    using Matrix = typename Types<Real>::Matrix;

    if (use_ds) {
        // DS mode: T_m = Q11 * Q31^{-1} * Inc
        //
        // Q31: [2*Nrank, 2*Nrank] (DS outer × DS inner)
        Matrix Q31 = assemble_Q_m<Real>(geom, 3, 1, k_real, k_int, n_rel,
                                         m, Nrank, Nmax, Nint, mirror, true, zRe, zIm, conducting);

        // Incident matrix: [2*Nrank, 2*Nmax] (DS outer × localized)
        Matrix Inc = assemble_incident_m<Real>(geom, k_real, m, Nrank, Nmax, Nint, zRe, zIm);

        // Solve Q31 * x = Inc → x = Q31^{-1} * Inc  [2*Nrank, 2*Nmax]
        double rcond;
        Matrix x = equilibrated_lu_solve<Real>(Q31, Inc, rcond);

        // Q11: [2*Nmax, 2*Nrank] (localized outer × DS inner)
        Matrix Q11 = assemble_Q_m<Real>(geom, 1, 1, k_real, k_int, n_rel,
                                         m, Nrank, Nmax, Nint, mirror, true, zRe, zIm, conducting);

        // T_m = Q11 * x = Q11 * Q31^{-1} * Inc  [2*Nmax, 2*Nmax]
        Matrix T_m = Q11 * x;
        return {T_m, rcond};
    } else {
        // Non-DS mode: T_m = Q11 * Q31^{-1}
        // Both matrices are square [2*Nmax, 2*Nmax]
        Matrix Q31 = assemble_Q_m<Real>(geom, 3, 1, k_real, k_int, n_rel,
                                         m, Nrank, Nmax, Nint, mirror, false, zRe, zIm, conducting);

        Matrix Q11 = assemble_Q_m<Real>(geom, 1, 1, k_real, k_int, n_rel,
                                         m, Nrank, Nmax, Nint, mirror, false, zRe, zIm, conducting);

        // The negation is incorporated into the phase factor during mapping.
        // T_m = Q11 * Q31^{-1} via: lu(Q31^T).solve(Q11^T)^T
        double rcond;
        Matrix T_m = equilibrated_lu_solve<Real>(
            Matrix(Q31.transpose()), Matrix(Q11.transpose()), rcond).transpose();
        return {T_m, rcond};
    }
}

// ============================================================================
// Main entry point: compute full [2, rmax, 2, rmax] T-matrix
// ============================================================================
template<typename Real>
TmatrixResult compute_axisymmetric_tmatrix(
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

    // Per-m reciprocal condition numbers (each m writes to its own index)
    std::vector<double> rconds(lmax + 1, 1.0);

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
    #pragma omp parallel for schedule(dynamic, 1)
    for (int m = 0; m <= lmax; m++) {
        int Nmax = (m == 0) ? Nrank : (Nrank - m + 1);
        if (Nmax <= 0) continue;

        // Compute T-matrix for this m
        auto [T_m, rcond] = tmatrix_m<Real>(geom, k_real, k_int, n_rel,
                                    m, Nrank, Nmax, Nint, mirror, use_ds_internal, zRe, zIm,
                                    conducting);
        rconds[m] = rcond;

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

    // Find worst (minimum) rcond across all m
    double worst_rcond = 1.0;
    int worst_m = 0;
    for (int m = 0; m <= lmax; m++) {
        if (rconds[m] < worst_rcond) {
            worst_rcond = rconds[m];
            worst_m = m;
        }
    }

    return {std::move(T), worst_rcond, worst_m};
}

// ============================================================================
// Diagnostic: return Q31 and Q11 matrices for a single azimuthal mode m
// Uses the same code path as compute_axisymmetric_tmatrix.
// Returns pair of flattened row-major [2*Nmax, 2*Nmax] complex double arrays.
// ============================================================================
template<typename Real>
std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
diagnostic_Q_matrices_m(
    const AxialGeometry<Real>& geom,
    double k_double, std::complex<double> n_rel_double,
    int m, int lmax, int Nint,
    bool use_ds, bool complex_plane, double eps_z_double,
    bool conducting) {

    using complex_t = std::complex<Real>;
    using Matrix = typename Types<Real>::Matrix;

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
        mirror = false;
    }

    int Nmax = (m == 0) ? Nrank : (Nrank - m + 1);
    if (Nmax <= 0) {
        return {{}, {}};
    }

    // Compute Q31 and Q11 using the same assemble_Q_m as the production code
    Matrix Q31 = assemble_Q_m<Real>(geom, 3, 1, k_real, k_int, n_rel,
                                     m, Nrank, Nmax, Nint, mirror, use_ds_internal, zRe, zIm, conducting);
    Matrix Q11 = assemble_Q_m<Real>(geom, 1, 1, k_real, k_int, n_rel,
                                     m, Nrank, Nmax, Nint, mirror, use_ds_internal, zRe, zIm, conducting);

    // Flatten to std::vector<complex<double>> in row-major order
    // Q31 and Q11 may have different dimensions in DS mode
    int rows31 = Q31.rows(), cols31 = Q31.cols();
    int rows11 = Q11.rows(), cols11 = Q11.cols();
    std::vector<std::complex<double>> Q31_flat(rows31 * cols31);
    std::vector<std::complex<double>> Q11_flat(rows11 * cols11);
    for (int i = 0; i < rows31; i++)
        for (int j = 0; j < cols31; j++)
            Q31_flat[i * cols31 + j] = std::complex<double>(
                double(Q31(i, j).real()), double(Q31(i, j).imag()));
    for (int i = 0; i < rows11; i++)
        for (int j = 0; j < cols11; j++)
            Q11_flat[i * cols11 + j] = std::complex<double>(
                double(Q11(i, j).real()), double(Q11(i, j).imag()));

    return {Q31_flat, Q11_flat};
}

// Explicit instantiations
template TmatrixResult compute_axisymmetric_tmatrix<double>(
    const AxialGeometry<double>&, double, std::complex<double>, int, int, bool, bool, double, bool);

template std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
diagnostic_Q_matrices_m<double>(
    const AxialGeometry<double>&, double, std::complex<double>, int, int, int, bool, bool, double, bool);

#if MIEPY_HAS_QUAD
template TmatrixResult compute_axisymmetric_tmatrix<__float128>(
    const AxialGeometry<__float128>&, double, std::complex<double>, int, int, bool, bool, double, bool);

template std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
diagnostic_Q_matrices_m<__float128>(
    const AxialGeometry<__float128>&, double, std::complex<double>, int, int, int, bool, bool, double, bool);
#endif

} // namespace tmatrix
