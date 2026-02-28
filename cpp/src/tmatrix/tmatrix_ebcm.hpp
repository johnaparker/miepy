#ifndef GUARD_tmatrix_ebcm_h
#define GUARD_tmatrix_ebcm_h

#include "tmatrix_types.hpp"
#include "tmatrix_geometry.hpp"
#include <complex>
#include <vector>

namespace tmatrix {

// Result of T-matrix computation, including conditioning diagnostics.
struct TmatrixResult {
    std::vector<std::complex<double>> T;  // T-matrix data [2, rmax, 2, rmax] flattened
    double worst_rcond;    // smallest reciprocal condition number across all m
    int worst_m;           // azimuthal mode with worst conditioning
};

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
//   TmatrixResult containing T-matrix data and conditioning diagnostics.
//   T-matrix in MiePy convention: complex double array [2, rmax, 2, rmax]
//   where rmax = lmax*(lmax+2), first index [2] = electric(0)/magnetic(1)
//
// Internal precision controlled by Real template parameter.
// Output is always in double precision (downcast if Real = __float128).
template<typename Real>
TmatrixResult compute_axisymmetric_tmatrix(
    const AxialGeometry<Real>& geom,
    double k, std::complex<double> n_rel,
    int lmax, int Nint,
    bool use_ds, bool complex_plane, double eps_z,
    bool conducting = false);

// Diagnostic: return Q31 and Q11 matrices for a single azimuthal mode m.
// Same code path as compute_axisymmetric_tmatrix, but returns raw Q matrices
// instead of the final T-matrix.
// Each matrix is flattened row-major [2*Nmax, 2*Nmax] where Nmax = Nrank-m+1 (or Nrank for m=0).
template<typename Real>
std::pair<std::vector<std::complex<double>>, std::vector<std::complex<double>>>
diagnostic_Q_matrices_m(
    const AxialGeometry<Real>& geom,
    double k, std::complex<double> n_rel,
    int m, int lmax, int Nint,
    bool use_ds, bool complex_plane, double eps_z,
    bool conducting = false);

} // namespace tmatrix

#endif
