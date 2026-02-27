#ifndef GUARD_tmatrix_ebcm_h
#define GUARD_tmatrix_ebcm_h

#include "tmatrix_types.hpp"
#include "tmatrix_geometry.hpp"
#include <complex>
#include <vector>

namespace tmatrix {

// Compute the T-matrix for an axisymmetric particle using the EBCM method.
//
// Arguments:
//   geom        - particle geometry (spheroid, cylinder, etc.)
//   k           - wavenumber in the medium (real, = 2*pi*sqrt(eps_m)/wavelength)
//   n_rel       - relative refractive index (complex, = sqrt(eps/eps_m))
//   lmax        - maximum angular momentum order
//   Nint        - number of integration (quadrature) points
//   use_ds      - if true, use distributed sources
//   complex_plane - if true, distribute sources in complex plane (oblate particles)
//   eps_z       - controls distributed source extent (fraction of particle extent)
//
// Returns:
//   T-matrix in MiePy convention: complex double array [2, rmax, 2, rmax]
//   where rmax = lmax*(lmax+2), first index [2] = electric(0)/magnetic(1)
//
// Internal precision controlled by Real template parameter.
// Output is always in double precision (downcast if Real = __float128).
template<typename Real>
std::vector<std::complex<double>> compute_axisymmetric_tmatrix(
    const AxialGeometry<Real>& geom,
    double k, std::complex<double> n_rel,
    int lmax, int Nint,
    bool use_ds, bool complex_plane, double eps_z,
    bool conducting = false);

} // namespace tmatrix

#endif
