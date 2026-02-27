#include "tmatrix/tmatrix_svwf.hpp"
#include "tmatrix/tmatrix_special.hpp"

namespace tmatrix {

// ============================================================================
// svwf_localized: Localized VSWFs for azimuthal mode m
// Port of SVWF.f90::MN (lines 11-74)
// ============================================================================
template<typename Real>
void svwf_localized(int index, std::complex<Real> z, Real theta,
                    int m, int Nrank, int Nmax,
                    std::vector<std::vector<std::complex<Real>>>& MV,
                    std::vector<std::vector<std::complex<Real>>>& NV) {
    using complex_t = std::complex<Real>;
    const complex_t zero(0, 0);
    const complex_t im_unit(0, 1);
    const Real MachEps = std::is_same<Real, double>::value ? 1e-15 : Real(1e-30);

    MV.assign(3, std::vector<complex_t>(Nmax, zero));
    NV.assign(3, std::vector<complex_t>(Nmax, zero));

    complex_t zc = z;
    if (tm_abs<Real>(zc) < MachEps)
        zc = complex_t(MachEps, MachEps);

    // Compute radial functions
    std::vector<complex_t> jh(Nrank + 1), jhd(Nrank + 1);
    if (index == 1) {
        bessel_j<Real>(zc, Nrank, jh, jhd);
    } else {
        bessel_h<Real>(zc, Nrank, jh, jhd);
    }

    // Compute angular functions
    std::vector<Real> Pnm, dPnm, pinm, taunm;
    legendre_normalized<Real>(theta, std::abs(m), Nrank, Pnm, dPnm, pinm, taunm);

    Real mr = Real(m);

    for (int k = 0; k < Nmax; k++) {
        int n;
        if (m == 0) {
            n = k + 1;
        } else {
            n = std::abs(m) + k;
        }

        Real nr = Real(n * (n + 1));
        Real nm = Real(1) / tm_sqrt(Real(2) * nr);

        complex_t fp(mr * pinm[n] * nm, 0);
        complex_t ft(taunm[n] * nm, 0);
        complex_t fl(nr * Pnm[n] * nm, 0);

        // M vector: (0, im*fp*jh, -ft*jh)
        MV[0][k] = zero;
        MV[1][k] = jh[n] * im_unit * fp;
        MV[2][k] = -jh[n] * ft;

        // N vector: (fl*jh/z, ft*jhd/z, im*fp*jhd/z)
        NV[0][k] = jh[n] * fl / zc;
        NV[1][k] = jhd[n] * ft / zc;
        NV[2][k] = jhd[n] * im_unit * fp / zc;
    }
}

// ============================================================================
// svwf_distributed: Distributed source VSWFs
// Port of SVWF.f90::MN_DS (lines 128-208)
// ============================================================================
template<typename Real>
void svwf_distributed(int index, std::complex<Real> k, Real r, Real theta,
                      const std::vector<Real>& zRe, const std::vector<Real>& zIm,
                      int m, int Nrank,
                      std::vector<std::vector<std::complex<Real>>>& MV,
                      std::vector<std::vector<std::complex<Real>>>& NV) {
    using complex_t = std::complex<Real>;
    const complex_t zero(0, 0);
    const complex_t im_unit(0, 1);
    const Real MachEps = std::is_same<Real, double>::value ? 1e-15 : Real(1e-30);

    MV.assign(3, std::vector<complex_t>(Nrank, zero));
    NV.assign(3, std::vector<complex_t>(Nrank, zero));

    // Determine the single order n used for all sources
    int n = (m == 0) ? 1 : std::abs(m);
    int Nbes = n + 1;

    Real mr = Real(m);
    Real nr = Real(n * (n + 1));
    Real nm = Real(1) / tm_sqrt(Real(2) * nr);

    Real sth = tm_sin(theta);
    Real cth = tm_cos(theta);
    Real ro = r * sth;    // cylindrical rho
    Real z_coord = r * cth;  // z coordinate

    std::vector<complex_t> jh(Nbes + 1), jhd(Nbes + 1);

    for (int p = 0; p < Nrank; p++) {
        Real dz = z_coord - zRe[p];
        complex_t dzc(dz, 0);
        complex_t roc(ro, 0);

        // Distance from source p to field point, in complex plane
        complex_t RR = tm_sqrt<Real>(roc * roc + (dzc - im_unit * complex_t(zIm[p], 0)) *
                                                  (dzc - im_unit * complex_t(zIm[p], 0)));
        if (tm_abs<Real>(RR) < MachEps)
            RR = complex_t(MachEps, MachEps);

        // Complex angles
        complex_t sint = roc / RR;
        complex_t cost = (dzc - im_unit * complex_t(zIm[p], 0)) / RR;
        complex_t argJ = k * RR;

        if (index == 1) {
            bessel_j<Real>(argJ, Nbes, jh, jhd);
        } else {
            bessel_h<Real>(argJ, Nbes, jh, jhd);
        }

        // Complex-argument Legendre functions
        complex_t Pmm = P_mm_complex<Real>(sint, cost, std::abs(m));
        complex_t pimm = pi_mm_complex<Real>(sint, std::abs(m));
        complex_t taumm = tau_mm_complex<Real>(sint, cost, std::abs(m));

        // Rotation from source-centered to global coordinates
        complex_t sinc = complex_t(sth, 0) * cost - complex_t(cth, 0) * sint;
        complex_t cosc = complex_t(cth, 0) * cost + complex_t(sth, 0) * sint;

        complex_t fp = im_unit * complex_t(mr, 0) * pimm * complex_t(nm, 0);
        complex_t ft = taumm * complex_t(nm, 0);
        complex_t fl = complex_t(nr, 0) * Pmm * complex_t(nm, 0);

        // M vector (rotated into global frame)
        complex_t factp = jh[n] * fp;
        MV[0][p] = factp * sinc;
        MV[1][p] = factp * cosc;
        MV[2][p] = -jh[n] * ft;

        // N vector (rotated into global frame)
        complex_t factl = jh[n] * fl;
        complex_t factt = jhd[n] * ft;
        NV[0][p] = (factl * cosc + factt * sinc) / argJ;
        NV[1][p] = (-factl * sinc + factt * cosc) / argJ;
        NV[2][p] = jhd[n] * fp / argJ;
    }
}

// Explicit instantiations
template void svwf_localized<double>(int, std::complex<double>, double, int, int, int,
    std::vector<std::vector<std::complex<double>>>&, std::vector<std::vector<std::complex<double>>>&);
template void svwf_distributed<double>(int, std::complex<double>, double, double,
    const std::vector<double>&, const std::vector<double>&, int, int,
    std::vector<std::vector<std::complex<double>>>&, std::vector<std::vector<std::complex<double>>>&);

#if MIEPY_HAS_QUAD
template void svwf_localized<__float128>(int, std::complex<__float128>, __float128, int, int, int,
    std::vector<std::vector<std::complex<__float128>>>&, std::vector<std::vector<std::complex<__float128>>>&);
template void svwf_distributed<__float128>(int, std::complex<__float128>, __float128, __float128,
    const std::vector<__float128>&, const std::vector<__float128>&, int, int,
    std::vector<std::vector<std::complex<__float128>>>&, std::vector<std::vector<std::complex<__float128>>>&);
#endif

} // namespace tmatrix
