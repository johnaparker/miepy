#include "tmatrix/tmatrix_special.hpp"
#include <algorithm>
#include <stdexcept>

namespace tmatrix {

// ============================================================================
// Bessel function constants (matching BesLeg.f90 parameters)
// ============================================================================
template<typename Real>
struct BesselParams {
    static constexpr Real ZeroSinXX = Real(1e-10);
    static constexpr Real FactNBes = Real(400.0);
    static constexpr Real InitBesVal = Real(1e-35);
    static constexpr Real MaxArgBes = Real(1e6);
    static constexpr Real UpperBoundSeq = Real(1e10);
    static constexpr Real LowerBoundSeq = Real(1e-10);
};

// ============================================================================
// bessel_j: Spherical Bessel j_n(z) via backward recurrence
// Port of BesLeg.f90::besel_j (lines 17-122)
// ============================================================================
template<typename Real>
void bessel_j(std::complex<Real> z, int N,
              std::vector<std::complex<Real>>& j,
              std::vector<std::complex<Real>>& jd) {
    using complex_t = std::complex<Real>;
    using P = BesselParams<Real>;

    j.resize(N + 1);
    jd.resize(N + 1);

    const complex_t zero(0, 0);
    const complex_t one(1, 0);

    for (int k = 0; k <= N; k++) {
        j[k] = zero;
        jd[k] = zero;
    }

    Real a0 = tm_abs<Real>(z);

    if (a0 < P::ZeroSinXX) {
        j[0] = one;
    } else if (a0 < P::MaxArgBes) {
        j[0] = tm_sin(z) / z;
        j[1] = j[0] / z - tm_cos(z) / z;

        if (N >= 2) {
            complex_t j0 = j[0];
            complex_t j1 = j[1];

            // Compute starting order for backward recurrence
            int Ma = int(a0 + Real(4) * std::pow(double(a0), 0.33) + Real(2)
                        + std::sqrt(double(Real(101) + a0))) + 20;
            int Mn = N + int(std::sqrt(double(P::FactNBes * N)));
            int M = std::max(Ma, Mn);

            complex_t f2 = zero;
            complex_t f1(P::InitBesVal, 0);

            for (int k = M; k >= 0; k--) {
                complex_t kc(Real(2 * k + 3), 0);
                complex_t f = kc * f1 / z - f2;
                if (k <= N) j[k] = f;
                f2 = f1;
                f1 = f;

                // Rescale to prevent overflow/underflow
                Real af1 = tm_abs<Real>(f1);
                if (af1 > P::UpperBoundSeq) {
                    f2 *= P::LowerBoundSeq;
                    f1 *= P::LowerBoundSeq;
                    for (int l = k; l <= N; l++)
                        j[l] *= P::LowerBoundSeq;
                } else if (af1 < P::LowerBoundSeq) {
                    f2 *= P::UpperBoundSeq;
                    f1 *= P::UpperBoundSeq;
                    for (int l = k; l <= N; l++)
                        j[l] *= P::UpperBoundSeq;
                }
            }

            // Normalize using known j(0) or j(1)
            complex_t scale;
            if (tm_abs<Real>(f1) > tm_abs<Real>(f2)) {
                scale = j0 / f1;
            } else {
                scale = j1 / f2;
            }
            for (int k = 0; k <= N; k++)
                j[k] *= scale;
        }

        // Compute derivatives: jd(n) = z*j(n-1) - n*j(n)
        for (int k = 1; k <= N; k++) {
            complex_t kc(Real(k), 0);
            jd[k] = z * j[k - 1] - kc * j[k];
        }
    } else {
        // Large argument asymptotic
        for (int k = 0; k <= N; k++) {
            Real arg = Real(0.5) * Real(k + 1) * pi<Real>();
            complex_t argc(arg, 0);
            j[k] = tm_cos(z - argc) / z;
        }
        for (int k = 1; k <= N; k++) {
            complex_t kc(Real(k), 0);
            jd[k] = z * j[k - 1] - kc * j[k];
        }
    }
}

// ============================================================================
// bessel_h: Spherical Hankel h^(1)_n(z) via forward recurrence for y_n
// Port of BesLeg.f90::besel_h (lines 124-235)
// ============================================================================
template<typename Real>
void bessel_h(std::complex<Real> z, int N,
              std::vector<std::complex<Real>>& h,
              std::vector<std::complex<Real>>& hd) {
    using complex_t = std::complex<Real>;
    using P = BesselParams<Real>;

    h.resize(N + 1);
    hd.resize(N + 1);

    const complex_t zero(0, 0);
    const complex_t one(1, 0);
    const complex_t im(0, 1);

    complex_t z2 = z * z;
    complex_t z3 = z2 * z;
    Real a0 = tm_abs<Real>(z);

    std::vector<complex_t> j(N + 1), y(N + 1);

    for (int k = 0; k <= N; k++) {
        j[k] = zero;
        y[k] = zero;
        h[k] = zero;
        hd[k] = zero;
    }

    if (a0 < P::ZeroSinXX) {
        j[0] = one;
    } else if (a0 < P::MaxArgBes) {
        j[0] = tm_sin(z) / z;
        j[1] = j[0] / z - tm_cos(z) / z;

        if (N >= 2) {
            complex_t j0 = j[0];
            complex_t j1 = j[1];

            int Ma = int(a0 + Real(4) * std::pow(double(a0), 0.33) + Real(2)
                        + std::sqrt(double(Real(101) + a0))) + 20;
            int Mn = N + int(std::sqrt(double(P::FactNBes * N)));
            int M = std::max(Ma, Mn);

            complex_t f2 = zero;
            complex_t f1(P::InitBesVal, 0);

            for (int k = M; k >= 0; k--) {
                complex_t kc(Real(2 * k + 3), 0);
                complex_t f = kc * f1 / z - f2;
                if (k <= N) j[k] = f;
                f2 = f1;
                f1 = f;

                Real af1 = tm_abs<Real>(f1);
                if (af1 > P::UpperBoundSeq) {
                    f2 *= P::LowerBoundSeq;
                    f1 *= P::LowerBoundSeq;
                    for (int l = k; l <= N; l++)
                        j[l] *= P::LowerBoundSeq;
                } else if (af1 < P::LowerBoundSeq) {
                    f2 *= P::UpperBoundSeq;
                    f1 *= P::UpperBoundSeq;
                    for (int l = k; l <= N; l++)
                        j[l] *= P::UpperBoundSeq;
                }
            }

            complex_t scale;
            if (tm_abs<Real>(f1) > tm_abs<Real>(f2)) {
                scale = j0 / f1;
            } else {
                scale = j1 / f2;
            }
            for (int k = 0; k <= N; k++)
                j[k] *= scale;
        }

        // Spherical Neumann functions via forward recurrence
        y[0] = -tm_cos(z) / z;
        y[1] = (y[0] - tm_sin(z)) / z;
        if (N >= 2) {
            for (int k = 2; k <= N; k++) {
                if (tm_abs<Real>(j[k - 1]) > tm_abs<Real>(j[k - 2])) {
                    y[k] = (j[k] * y[k - 1] - one / z2) / j[k - 1];
                } else {
                    complex_t kc(Real(2 * k - 1), 0);
                    y[k] = (j[k] * y[k - 2] - kc / z3) / j[k - 2];
                }
            }
        }

        // Hankel = j + i*y
        for (int k = 0; k <= N; k++)
            h[k] = j[k] + im * y[k];

        // Derivatives
        for (int k = 1; k <= N; k++) {
            complex_t kc(Real(k), 0);
            hd[k] = z * h[k - 1] - kc * h[k];
        }
    } else {
        // Large argument asymptotic
        for (int k = 0; k <= N; k++) {
            Real arg = Real(0.5) * Real(k + 1) * pi<Real>();
            complex_t argc(arg, 0);
            h[k] = tm_exp(im * (z - argc)) / z;
        }
        for (int k = 1; k <= N; k++) {
            complex_t kc(Real(k), 0);
            hd[k] = z * h[k - 1] - kc * h[k];
        }
    }
}

// ============================================================================
// legendre_normalized: Normalized associated Legendre + pi_nm, tau_nm
// Port of BesLeg.f90::Leg_normalized (lines 744-850)
// ============================================================================
template<typename Real>
void legendre_normalized(Real theta, int m, int Nmax,
                         std::vector<Real>& Pnm,
                         std::vector<Real>& dPnm,
                         std::vector<Real>& pinm,
                         std::vector<Real>& taunm) {
    Pnm.resize(Nmax + 1);
    dPnm.resize(Nmax + 1);
    pinm.resize(Nmax + 1);
    taunm.resize(Nmax + 1);

    Real cth = tm_cos(theta);
    Real sth = tm_sin(theta);

    // Compute Pnm via recurrence
    if (m == 0) {
        Pnm[0] = tm_sqrt(Real(2)) / Real(2);
        if (Nmax >= 1)
            Pnm[1] = tm_sqrt(Real(3) / Real(2)) * cth;
        for (int n = 1; n < Nmax; n++) {
            Real f1 = tm_sqrt(Real(2 * n + 1) / Real(n + 1));
            f1 *= tm_sqrt(Real(2 * n + 3) / Real(n + 1));
            Real f2 = tm_sqrt(Real(2 * n + 3) / Real(2 * n - 1));
            f2 *= Real(n) / Real(n + 1);
            Pnm[n + 1] = f1 * cth * Pnm[n] - f2 * Pnm[n - 1];
        }
    } else {
        for (int n = 0; n < m; n++)
            Pnm[n] = Real(0);

        // P_mm_real: compute Pm = sqrt((2m+1)/2) * prod_{k=1..m} sqrt((m+k)/(4k)) * sin^m(theta)
        Real tamp = Real(1);
        for (int k = 1; k <= m; k++) {
            Real f = tm_sqrt(Real(m + k) / (Real(4) * Real(k)));
            tamp *= f * sth;
        }
        Real Pm = tm_sqrt(Real(2 * m + 1) / Real(2)) * tamp;
        Pnm[m] = Pm;

        for (int n = m; n < Nmax; n++) {
            Real f1 = tm_sqrt(Real(2 * n + 1) / Real(n + 1 - m));
            f1 *= tm_sqrt(Real(2 * n + 3) / Real(n + 1 + m));
            Real f2 = tm_sqrt(Real(2 * n + 3) / Real(2 * n - 1));
            f2 *= tm_sqrt(Real(n - m) / Real(n - m + 1));
            f2 *= tm_sqrt(Real(n + m) / Real(n + m + 1));
            Pnm[n + 1] = f1 * cth * Pnm[n] - f2 * Pnm[n - 1];
        }
    }

    // Compute dPnm (only for m == 0)
    if (m == 0) {
        dPnm[0] = Real(0);
        for (int n = 1; n <= Nmax; n++) {
            Real f2 = tm_sqrt(Real(2 * n + 1) / Real(2 * n - 1));
            Real f1 = Real(n) * f2;
            dPnm[n] = f1 * Pnm[n - 1] + f2 * cth * dPnm[n - 1];
        }
    } else {
        for (int n = 0; n <= Nmax; n++)
            dPnm[n] = Real(0);
    }

    // Compute pinm and taunm
    if (m == 0) {
        for (int n = 0; n <= Nmax; n++) {
            pinm[n] = Real(0);
            taunm[n] = -sth * dPnm[n];
        }
    } else {
        if (m == 1) {
            pinm[0] = Real(0);
            taunm[0] = Real(0);
            pinm[1] = tm_sqrt(Real(3)) / Real(2);
        } else {
            for (int n = 0; n < m; n++) {
                pinm[n] = Real(0);
                taunm[n] = Real(0);
            }
            // pi_mm_real
            Real tamp2 = Real(1);
            for (int k = 1; k <= m - 1; k++) {
                Real f = tm_sqrt(Real(m + k) / (Real(4) * Real(k)));
                tamp2 *= f * sth;
            }
            pinm[m] = tm_sqrt(Real(2 * m + 1)) / Real(2) * tamp2;
        }

        // pinm recurrence
        for (int n = m; n < Nmax; n++) {
            Real f1 = tm_sqrt(Real(2 * n + 1) / Real(n + 1 - m));
            f1 *= tm_sqrt(Real(2 * n + 3) / Real(n + 1 + m));
            Real f2 = tm_sqrt(Real(2 * n + 3) / Real(2 * n - 1));
            f2 *= tm_sqrt(Real(n - m) / Real(n - m + 1));
            f2 *= tm_sqrt(Real(n + m) / Real(n + m + 1));
            pinm[n + 1] = f1 * cth * pinm[n] - f2 * pinm[n - 1];
        }

        // taunm
        for (int n = m; n <= Nmax; n++) {
            Real f1 = Real(n);
            Real f2 = tm_sqrt(Real(2 * n + 1) / Real(2 * n - 1));
            f2 *= tm_sqrt(Real(n - m) / Real(n + m));
            f2 *= Real(n + m);
            taunm[n] = f1 * cth * pinm[n] - f2 * pinm[n - 1];
        }
    }
}

// ============================================================================
// Complex-argument Legendre functions for distributed sources
// Port of BesLeg.f90 lines 852-924
// ============================================================================

template<typename Real>
std::complex<Real> P_mm_complex(std::complex<Real> sint, std::complex<Real> cost, int m) {
    using complex_t = std::complex<Real>;
    if (m == 0) {
        return tm_sqrt(Real(3) / Real(2)) * cost;
    }
    complex_t tamp(1, 0);
    for (int k = 1; k <= m; k++) {
        Real f = tm_sqrt(Real(m + k) / (Real(4) * Real(k)));
        tamp *= complex_t(f, 0) * sint;
    }
    Real f = tm_sqrt(Real(2 * m + 1) / Real(2));
    return complex_t(f, 0) * tamp;
}

template<typename Real>
std::complex<Real> pi_mm_complex(std::complex<Real> sint, int m) {
    using complex_t = std::complex<Real>;
    if (m == 0) {
        return complex_t(0, 0);
    }
    complex_t tamp(1, 0);
    for (int k = 1; k <= m - 1; k++) {
        Real f = tm_sqrt(Real(m + k) / (Real(4) * Real(k)));
        tamp *= complex_t(f, 0) * sint;
    }
    Real f = tm_sqrt(Real(2 * m + 1)) / Real(2);
    return complex_t(f, 0) * tamp;
}

template<typename Real>
std::complex<Real> tau_mm_complex(std::complex<Real> sint, std::complex<Real> cost, int m) {
    using complex_t = std::complex<Real>;
    if (m == 0) {
        return -tm_sqrt(Real(3) / Real(2)) * sint;
    }
    complex_t tamp(1, 0);
    for (int k = 1; k <= m - 1; k++) {
        Real f = tm_sqrt(Real(m + k) / (Real(4) * Real(k)));
        tamp *= complex_t(f, 0) * sint;
    }
    Real f = Real(m) * tm_sqrt(Real(2 * m + 1)) / Real(2);
    return complex_t(f, 0) * cost * tamp;
}

// ============================================================================
// Gauss-Legendre quadrature via Newton-Raphson
// Port of Integr.f90::Gauss_Legendre1 (lines 225-288)
// ============================================================================
template<typename Real>
void gauss_legendre(Real x1, Real x2, int Np,
                    std::vector<Real>& weights,
                    std::vector<Real>& nodes) {
    weights.resize(Np);
    nodes.resize(Np);

    const Real eps = std::is_same<Real, double>::value ? Real(1e-15) : Real(1e-30);
    const int maxIter = 100000;

    Real xm = Real(0.5) * (x1 + x2);
    Real xl = Real(0.5) * (x2 - x1);
    int M = (Np + 1) / 2;

    for (int i = 0; i < M; i++) {
        Real z = tm_cos(pi<Real>() * (Real(i + 1) - Real(0.25)) / (Real(Np) + Real(0.5)));

        Real z1;
        int iter = 0;
        do {
            Real p1 = Real(1);
            Real p2 = Real(0);
            for (int j = 1; j <= Np; j++) {
                Real p3 = p2;
                p2 = p1;
                p1 = (Real(2 * j - 1) * z * p2 - Real(j - 1) * p3) / Real(j);
            }
            Real pp = Real(Np) * (z * p1 - p2) / (z * z - Real(1));
            z1 = z;
            z = z1 - p1 / pp;
            iter++;
        } while (tm_abs(z - z1) > eps && iter < maxIter);

        nodes[i] = xm - xl * z;
        nodes[Np - 1 - i] = xm + xl * z;
        // Recompute Legendre polynomial derivative at converged root
        {
            Real p1 = Real(1), p2 = Real(0);
            for (int j = 1; j <= Np; j++) {
                Real p3 = p2; p2 = p1;
                p1 = (Real(2*j-1)*z*p2 - Real(j-1)*p3)/Real(j);
            }
            Real pp = Real(Np)*(z*p1 - p2)/(z*z - Real(1));
            weights[i] = Real(2) * xl / ((Real(1) - z * z) * pp * pp);
        }
        weights[Np - 1 - i] = weights[i];
    }
}

// ============================================================================
// Explicit template instantiations
// ============================================================================
template void bessel_j<double>(std::complex<double>, int,
    std::vector<std::complex<double>>&, std::vector<std::complex<double>>&);
template void bessel_h<double>(std::complex<double>, int,
    std::vector<std::complex<double>>&, std::vector<std::complex<double>>&);
template void legendre_normalized<double>(double, int, int,
    std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&);
template std::complex<double> P_mm_complex<double>(std::complex<double>, std::complex<double>, int);
template std::complex<double> pi_mm_complex<double>(std::complex<double>, int);
template std::complex<double> tau_mm_complex<double>(std::complex<double>, std::complex<double>, int);
template void gauss_legendre<double>(double, double, int, std::vector<double>&, std::vector<double>&);

#if MIEPY_HAS_QUAD
template void bessel_j<__float128>(std::complex<__float128>, int,
    std::vector<std::complex<__float128>>&, std::vector<std::complex<__float128>>&);
template void bessel_h<__float128>(std::complex<__float128>, int,
    std::vector<std::complex<__float128>>&, std::vector<std::complex<__float128>>&);
template void legendre_normalized<__float128>(__float128, int, int,
    std::vector<__float128>&, std::vector<__float128>&, std::vector<__float128>&, std::vector<__float128>&);
template std::complex<__float128> P_mm_complex<__float128>(std::complex<__float128>, std::complex<__float128>, int);
template std::complex<__float128> pi_mm_complex<__float128>(std::complex<__float128>, int);
template std::complex<__float128> tau_mm_complex<__float128>(std::complex<__float128>, std::complex<__float128>, int);
template void gauss_legendre<__float128>(__float128, __float128, int, std::vector<__float128>&, std::vector<__float128>&);
#endif

} // namespace tmatrix
