#ifndef GUARD_tmatrix_ebcm_nonaxial_h
#define GUARD_tmatrix_ebcm_nonaxial_h

#include "tmatrix_types.hpp"
#include "tmatrix_ebcm.hpp"
#include "tmatrix_geometry.hpp"
#include <complex>
#include <vector>

namespace tmatrix {

// Compute the T-matrix for a non-axisymmetric particle using the EBCM method.
//
// Unlike the axisymmetric solver where each m-mode is independent, all m-modes
// couple into a single large matrix of size Nmax = Nrank + Mrank*(2*Nrank - Mrank + 1).
//
// Arguments:
//   geom        - non-axisymmetric particle geometry
//   k           - wavenumber in the medium (real)
//   n_rel       - relative refractive index (complex)
//   lmax        - maximum angular momentum order (= Nrank = Mrank)
//   Nint1, Nint2 - number of quadrature points in each surface direction
//   conducting  - if true, use PEC boundary conditions
//
// Returns:
//   TmatrixResult containing T-matrix in MiePy convention [2, rmax, 2, rmax]
//   and conditioning diagnostics.
template<typename Real>
TmatrixResult compute_nonaxisymmetric_tmatrix(
    const NonAxialGeometry<Real>& geom,
    double k, std::complex<double> n_rel,
    int lmax, int Nint1, int Nint2,
    bool conducting = false);

} // namespace tmatrix

#endif
