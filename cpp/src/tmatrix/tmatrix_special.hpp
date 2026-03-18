#ifndef GUARD_tmatrix_special_h
#define GUARD_tmatrix_special_h

#include "tmatrix_types.hpp"
#include <vector>

namespace tmatrix {

// Spherical Bessel function j_n(z) and derivatives jd_n(z) = d[z*j_n(z)]/dz
// Port of BesLeg.f90 besel_j (lines 17-122)
// n ranges from 0..N, results stored in j[0..N], jd[0..N]
template<typename Real>
void bessel_j(std::complex<Real> z, int N,
              std::vector<std::complex<Real>>& j,
              std::vector<std::complex<Real>>& jd);

// Spherical Hankel function h^(1)_n(z) and derivatives hd_n(z) = d[z*h_n(z)]/dz
// Port of BesLeg.f90 besel_h (lines 124-235)
template<typename Real>
void bessel_h(std::complex<Real> z, int N,
              std::vector<std::complex<Real>>& h,
              std::vector<std::complex<Real>>& hd);

// Normalized associated Legendre functions and angular functions
// Port of BesLeg.f90 Leg_normalized (lines 744-850)
// Pnm[0..Nmax], dPnm[0..Nmax], pinm[0..Nmax], taunm[0..Nmax]
template<typename Real>
void legendre_normalized(Real theta, int m, int Nmax,
                         std::vector<Real>& Pnm,
                         std::vector<Real>& dPnm,
                         std::vector<Real>& pinm,
                         std::vector<Real>& taunm);

// P_mm for complex argument (distributed sources)
// Port of BesLeg.f90 P_mm (lines 852-873)
template<typename Real>
std::complex<Real> P_mm_complex(std::complex<Real> sint, std::complex<Real> cost, int m);

// pi_mm for complex argument
// Port of BesLeg.f90 pi_mm (lines 876-898)
template<typename Real>
std::complex<Real> pi_mm_complex(std::complex<Real> sint, int m);

// tau_mm for complex argument
// Port of BesLeg.f90 tau_mm (lines 901-923)
template<typename Real>
std::complex<Real> tau_mm_complex(std::complex<Real> sint, std::complex<Real> cost, int m);

// Gauss-Legendre quadrature nodes and weights on [x1, x2]
// Port of Integr.f90 Gauss_Legendre1 (lines 225-288)
template<typename Real>
void gauss_legendre(Real x1, Real x2, int N,
                    std::vector<Real>& weights,
                    std::vector<Real>& nodes);

} // namespace tmatrix

#endif
