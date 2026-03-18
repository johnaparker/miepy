#ifndef GUARD_tmatrix_svwf_h
#define GUARD_tmatrix_svwf_h

#include "tmatrix_types.hpp"
#include <vector>

namespace tmatrix {

// Localized vector spherical wave functions (VSWF) for azimuthal mode m
// Port of SVWF.f90::MN (lines 11-74)
//
// index: 1 = regular (Bessel j), 3 = radiating (Hankel h)
// z = k*r (complex wavenumber * radial distance)
// theta: zenith angle
// m: azimuthal mode number
// Nrank: maximum expansion order
// Nmax: dimension of output vectors
//   For m=0: Nmax = Nrank
//   For |m|>0: Nmax = Nrank - |m| + 1
//
// Output: MV[3][Nmax], NV[3][Nmax] (3 = r,theta,phi components)
// Excludes the exp(i*m*phi) factor
template<typename Real>
void svwf_localized(int index, std::complex<Real> z, Real theta,
                    int m, int Nrank, int Nmax,
                    std::vector<std::vector<std::complex<Real>>>& MV,
                    std::vector<std::vector<std::complex<Real>>>& NV);

// Distributed source vector spherical wave functions
// Port of SVWF.f90::MN_DS (lines 128-208)
//
// k: wavenumber
// r: radial distance (real)
// theta: zenith angle
// zRe, zIm: distributed source coordinates
// m: azimuthal mode
// Nrank: number of distributed sources
template<typename Real>
void svwf_distributed(int index, std::complex<Real> k, Real r, Real theta,
                      const std::vector<Real>& zRe, const std::vector<Real>& zIm,
                      int m, int Nrank,
                      std::vector<std::vector<std::complex<Real>>>& MV,
                      std::vector<std::vector<std::complex<Real>>>& NV);

// Complete vector spherical wave functions for ALL azimuthal modes m = 0..Mrank
// Port of SVWF.f90::MN_complete (lines 774-869)
//
// index: 1 = regular (Bessel j), 3 = radiating (Hankel h)
// z = k*r (complex wavenumber * radial distance)
// theta, phi: polar angles
// Mrank: maximum azimuthal order
// Nrank: maximum expansion order
// Nmax: total number of modes = Nrank + Mrank*(2*Nrank - Mrank + 1)
// plus: if true, m ordering is 0,+1,-1,+2,-2,...; if false, 0,-1,+1,-2,+2,...
// sym: if true, omit exp(i*m*phi); if false, include it
//
// Output: MV[3][Nmax], NV[3][Nmax] (3 = r,theta,phi components)
template<typename Real>
void svwf_complete(int index, std::complex<Real> z, Real theta, Real phi,
                   int Mrank, int Nrank, int Nmax, bool plus, bool sym,
                   std::vector<std::vector<std::complex<Real>>>& MV,
                   std::vector<std::vector<std::complex<Real>>>& NV);

} // namespace tmatrix

#endif
