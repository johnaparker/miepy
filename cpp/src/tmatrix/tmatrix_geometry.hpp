#ifndef GUARD_tmatrix_geometry_h
#define GUARD_tmatrix_geometry_h

#include "tmatrix_types.hpp"
#include "tmatrix_special.hpp"
#include <vector>
#include <stdexcept>

namespace tmatrix {

// Surface geometry point for axisymmetric particles
template<typename Real>
struct GeometryPoint {
    Real r;        // magnitude of position vector
    Real theta;    // zenith angle
    Real dA;       // surface element
    Real n_r;      // outward normal, r-component
    Real n_theta;  // outward normal, theta-component
};

// Base class for axisymmetric particle geometry
template<typename Real>
class AxialGeometry {
public:
    virtual ~AxialGeometry() = default;

    // Number of smooth curve segments (e.g. 1 for spheroid, 3 for cylinder)
    virtual int num_curves() const = 0;

    // Evaluate geometry at parameter value on curve segment
    virtual GeometryPoint<Real> evaluate(Real param, int curve) const = 0;

    // Whether the particle has mirror symmetry about z=0
    virtual bool is_mirror_symmetric() const = 0;

    // Generate distributed source coordinates
    // For prolate: sources on real z-axis, zIm = 0
    // For oblate: sources on imaginary z-axis, zRe = 0
    virtual void distributed_sources(int Nrank, bool complex_plane, Real eps_z,
                                     std::vector<Real>& zRe,
                                     std::vector<Real>& zIm) const = 0;

    // Get quadrature nodes and weights for all curves
    // Returns (paramG, weightsG, Nintparam) for each curve segment
    virtual void quadrature_points(int Nint, bool mirror,
                                   std::vector<std::vector<Real>>& paramG,
                                   std::vector<std::vector<Real>>& weightsG) const = 0;
};

// ============================================================================
// Spheroid: TypeGeom=1
// surf(1) = semi-axis along z (a), surf(2) = semi-axis perpendicular (b)
// Single curve, mirror symmetric
// ============================================================================
template<typename Real>
class SpheroidGeometry : public AxialGeometry<Real> {
public:
    SpheroidGeometry(Real axis_z, Real axis_xy)
        : a_(axis_z), b_(axis_xy) {}

    int num_curves() const override { return 1; }
    bool is_mirror_symmetric() const override { return true; }

    GeometryPoint<Real> evaluate(Real param, int /*curve*/) const override {
        GeometryPoint<Real> pt;
        Real e = a_ / b_;
        Real e2 = e * e;
        Real theta = param;
        Real cth = tm_cos(theta);
        Real sth = tm_sin(theta);

        Real tamp = cth * cth + sth * sth * e2;
        Real invtamp = Real(1) / tm_sqrt(tamp);
        Real r = a_ * invtamp;
        Real dr = -a_ * cth * sth * (e2 - Real(1)) * invtamp / tamp;

        Real norm = tm_sqrt(r * r + dr * dr);
        Real invnorm = Real(1) / norm;

        pt.r = r;
        pt.theta = theta;
        pt.dA = norm * r * sth;
        pt.n_r = r * invnorm;
        pt.n_theta = -dr * invnorm;
        return pt;
    }

    void distributed_sources(int Nrank, bool complex_plane, Real eps_z,
                             std::vector<Real>& zRe,
                             std::vector<Real>& zIm) const override {
        zRe.resize(Nrank);
        zIm.resize(Nrank);
        if (!complex_plane) {
            // Prolate: sources on real z-axis
            Real zmin = -eps_z * a_;
            Real zmax = eps_z * a_;
            for (int i = 0; i < Nrank; i++) {
                Real t = (Nrank > 1) ? Real(i) / Real(Nrank - 1) : Real(0.5);
                zRe[i] = zmin + t * (zmax - zmin);
                zIm[i] = Real(0);
            }
        } else {
            // Oblate: sources on imaginary z-axis
            Real zmin = -eps_z * b_;
            Real zmax = eps_z * b_;
            for (int i = 0; i < Nrank; i++) {
                Real t = (Nrank > 1) ? Real(i) / Real(Nrank - 1) : Real(0.5);
                zRe[i] = Real(0);
                zIm[i] = zmin + t * (zmax - zmin);
            }
        }
    }

    void quadrature_points(int Nint, bool mirror,
                           std::vector<std::vector<Real>>& paramG,
                           std::vector<std::vector<Real>>& weightsG) const override {
        paramG.resize(1);
        weightsG.resize(1);
        Real a = Real(0);
        Real b = mirror ? pi<Real>() / Real(2) : pi<Real>();
        gauss_legendre<Real>(a, b, Nint, weightsG[0], paramG[0]);
    }

    Real axis_z() const { return a_; }
    Real axis_xy() const { return b_; }

private:
    Real a_;  // semi-axis along z
    Real b_;  // semi-axis perpendicular
};

// ============================================================================
// Cylinder: TypeGeom=2
// surf(1) = half-length (a), surf(2) = radius (b)
// Three curves: top cap, side, bottom cap. Mirror symmetric.
// ============================================================================
template<typename Real>
class CylinderGeometry : public AxialGeometry<Real> {
public:
    CylinderGeometry(Real half_height, Real radius)
        : a_(half_height), b_(radius) {}

    int num_curves() const override { return 3; }
    bool is_mirror_symmetric() const override { return true; }

    GeometryPoint<Real> evaluate(Real param, int curve) const override {
        GeometryPoint<Real> pt;
        Real theta = param;
        Real cth = tm_cos(theta);
        Real sth = tm_sin(theta);
        Real r, dr;

        if (curve == 0) {
            // Top cap: r = a/cos(theta)
            Real invcos = Real(1) / cth;
            r = a_ * invcos;
            dr = a_ * sth * invcos * invcos;
        } else if (curve == 1) {
            // Side: r = b/sin(theta)
            Real invsin = Real(1) / sth;
            r = b_ * invsin;
            dr = -b_ * cth * invsin * invsin;
        } else {
            // Bottom cap: r = -a/cos(theta)
            Real invcos = Real(1) / cth;
            r = -a_ * invcos;
            dr = -a_ * sth * invcos * invcos;
        }

        Real norm = tm_sqrt(r * r + dr * dr);
        Real invnorm = Real(1) / norm;

        pt.r = r;
        pt.theta = theta;
        pt.dA = norm * r * sth;
        pt.n_r = r * invnorm;
        pt.n_theta = -dr * invnorm;
        return pt;
    }

    void distributed_sources(int Nrank, bool complex_plane, Real eps_z,
                             std::vector<Real>& zRe,
                             std::vector<Real>& zIm) const override {
        zRe.resize(Nrank);
        zIm.resize(Nrank);
        if (!complex_plane) {
            Real zmin = -eps_z * a_;
            Real zmax = eps_z * a_;
            for (int i = 0; i < Nrank; i++) {
                Real t = (Nrank > 1) ? Real(i) / Real(Nrank - 1) : Real(0.5);
                zRe[i] = zmin + t * (zmax - zmin);
                zIm[i] = Real(0);
            }
        } else {
            Real zmin = -eps_z * b_;
            Real zmax = eps_z * b_;
            for (int i = 0; i < Nrank; i++) {
                Real t = (Nrank > 1) ? Real(i) / Real(Nrank - 1) : Real(0.5);
                zRe[i] = Real(0);
                zIm[i] = zmin + t * (zmax - zmin);
            }
        }
    }

    void quadrature_points(int Nint, bool mirror,
                           std::vector<std::vector<Real>>& paramG,
                           std::vector<std::vector<Real>>& weightsG) const override {
        (void)mirror;  // Cylinder always uses full range (3 segments handle it)
        paramG.resize(3);
        weightsG.resize(3);

        Real theta0 = tm_atan(b_ / a_);

        // Distribute points proportionally
        int Nintb = int(b_ * Nint / Real(2) / (a_ + b_));
        int Ninta = Nint - 2 * Nintb;
        if (Nintb < 20) Nintb = 20;
        if (Ninta < 20) Ninta = 20;

        // Top cap: [0, theta0]
        gauss_legendre<Real>(Real(0), theta0, Nintb, weightsG[0], paramG[0]);
        // Side: [theta0, pi - theta0]
        gauss_legendre<Real>(theta0, pi<Real>() - theta0, Ninta, weightsG[1], paramG[1]);
        // Bottom cap: [pi - theta0, pi]
        gauss_legendre<Real>(pi<Real>() - theta0, pi<Real>(), Nintb, weightsG[2], paramG[2]);
    }

private:
    Real a_;  // half-height
    Real b_;  // radius
};

// ============================================================================
// Rounded Cylinder: TypeGeom=3
// surf(1) = half-length (a), surf(2) = radius including rounded part (b)
// Three curves with hemispherical caps. Mirror symmetric.
// ============================================================================
template<typename Real>
class RoundedCylinderGeometry : public AxialGeometry<Real> {
public:
    RoundedCylinderGeometry(Real half_height, Real radius)
        : a_(half_height), b_(radius) {
        if (half_height >= radius) {
            throw std::runtime_error("Rounded cylinder requires radius > half_height (oblate only)");
        }
    }

    int num_curves() const override { return 3; }
    bool is_mirror_symmetric() const override { return true; }

    GeometryPoint<Real> evaluate(Real param, int curve) const override {
        GeometryPoint<Real> pt;
        Real theta = param;
        Real cth = tm_cos(theta);
        Real sth = tm_sin(theta);
        Real r, dr;

        Real e = b_ - a_;
        Real e2 = e * e;

        if (curve == 0) {
            // Top cap
            Real invcos = Real(1) / cth;
            r = a_ * invcos;
            dr = a_ * sth * invcos * invcos;
        } else if (curve == 1) {
            // Rounded side
            Real tamp = tm_sqrt(a_ * a_ - e2 * cth * cth);
            r = e * sth + tamp;
            dr = e * cth + e2 * sth * cth / tamp;
        } else {
            // Bottom cap
            Real invcos = Real(1) / cth;
            r = -a_ * invcos;
            dr = -a_ * sth * invcos * invcos;
        }

        Real norm = tm_sqrt(r * r + dr * dr);
        Real invnorm = Real(1) / norm;

        pt.r = r;
        pt.theta = theta;
        pt.dA = norm * r * sth;
        pt.n_r = r * invnorm;
        pt.n_theta = -dr * invnorm;
        return pt;
    }

    void distributed_sources(int Nrank, bool complex_plane, Real eps_z,
                             std::vector<Real>& zRe,
                             std::vector<Real>& zIm) const override {
        zRe.resize(Nrank);
        zIm.resize(Nrank);
        if (!complex_plane) {
            Real zmin = -eps_z * a_;
            Real zmax = eps_z * a_;
            for (int i = 0; i < Nrank; i++) {
                Real t = (Nrank > 1) ? Real(i) / Real(Nrank - 1) : Real(0.5);
                zRe[i] = zmin + t * (zmax - zmin);
                zIm[i] = Real(0);
            }
        } else {
            Real zmin = -eps_z * b_;
            Real zmax = eps_z * b_;
            for (int i = 0; i < Nrank; i++) {
                Real t = (Nrank > 1) ? Real(i) / Real(Nrank - 1) : Real(0.5);
                zRe[i] = Real(0);
                zIm[i] = zmin + t * (zmax - zmin);
            }
        }
    }

    void quadrature_points(int Nint, bool mirror,
                           std::vector<std::vector<Real>>& paramG,
                           std::vector<std::vector<Real>>& weightsG) const override {
        (void)mirror;
        paramG.resize(3);
        weightsG.resize(3);

        Real e = b_ - a_;
        Real theta0 = tm_atan(e / a_);

        int Nintb = int(b_ * Nint / Real(2) / (a_ + b_));
        int Ninta = Nint - 2 * Nintb;
        if (Nintb < 20) Nintb = 20;
        if (Ninta < 20) Ninta = 20;

        gauss_legendre<Real>(Real(0), theta0, Nintb, weightsG[0], paramG[0]);
        gauss_legendre<Real>(theta0, pi<Real>() - theta0, Ninta, weightsG[1], paramG[1]);
        gauss_legendre<Real>(pi<Real>() - theta0, pi<Real>(), Nintb, weightsG[2], paramG[2]);
    }

private:
    Real a_;  // half-height
    Real b_;  // radius (including rounded part)
};

// ============================================================================
// 3D (non-axisymmetric) geometry support
// ============================================================================

// Surface geometry point for non-axisymmetric particles
template<typename Real>
struct GeometryPoint3D {
    Real r;                        // magnitude of position vector
    Real theta;                    // zenith angle
    Real phi;                      // azimuthal angle
    Real dA;                       // surface element
    Real n_r, n_theta, n_phi;     // outward normal in spherical coords
};

// Cartesian to spherical coordinate conversion
template<typename Real>
static void cartesian_to_spherical(Real x, Real y, Real z,
                                   Real& r, Real& theta, Real& phi) {
    r = tm_sqrt(x * x + y * y + z * z);
    Real rho = tm_sqrt(x * x + y * y);
    theta = tm_atan2(rho, z);
    phi = tm_atan2(y, x);
}

// Convert Cartesian normal to spherical normal components
template<typename Real>
static void normal_cart_to_sph(Real theta, Real phi,
                               Real nc_x, Real nc_y, Real nc_z,
                               Real& n_r, Real& n_theta, Real& n_phi) {
    Real st = tm_sin(theta), ct = tm_cos(theta);
    Real sp = tm_sin(phi), cp = tm_cos(phi);
    n_r     =  st * cp * nc_x + st * sp * nc_y + ct * nc_z;
    n_theta =  ct * cp * nc_x + ct * sp * nc_y - st * nc_z;
    n_phi   = -sp * nc_x + cp * nc_y;
}

// Base class for non-axisymmetric particle geometry
template<typename Real>
class NonAxialGeometry {
public:
    virtual ~NonAxialGeometry() = default;

    // Number of surface patches
    virtual int num_patches() const = 0;

    // Evaluate geometry at (param1, param2) on patch
    virtual GeometryPoint3D<Real> evaluate(Real param1, Real param2, int patch) const = 0;

    // Whether the particle has mirror symmetry about z=0
    virtual bool is_mirror_symmetric() const = 0;

    // Azimuthal N-fold symmetry (0 = none)
    virtual int azimuthal_symmetry() const = 0;

    // Generate 2D quadrature points for all patches
    virtual void quadrature_points(int Nint1, int Nint2, bool mirror,
                                   std::vector<std::vector<Real>>& param1G,
                                   std::vector<std::vector<Real>>& param2G,
                                   std::vector<std::vector<Real>>& weightsG) const = 0;
};

// ============================================================================
// Ellipsoid: 3D geometry
// Semi-axes a (x), b (y), c (z)
// Single patch parameterized by (u, v) on [0,pi] x [0,2pi]
// ============================================================================
template<typename Real>
class EllipsoidGeometry3D : public NonAxialGeometry<Real> {
public:
    EllipsoidGeometry3D(Real rx, Real ry, Real rz)
        : a_(rx), b_(ry), c_(rz) {}

    int num_patches() const override { return 1; }
    bool is_mirror_symmetric() const override { return true; }
    int azimuthal_symmetry() const override { return 0; }

    GeometryPoint3D<Real> evaluate(Real u, Real v, int /*patch*/) const override {
        GeometryPoint3D<Real> pt;

        Real su = tm_sin(u), cu = tm_cos(u);
        Real sv = tm_sin(v), cv = tm_cos(v);

        Real x = a_ * su * cv;
        Real y = b_ * su * sv;
        Real z = c_ * cu;

        cartesian_to_spherical(x, y, z, pt.r, pt.theta, pt.phi);

        // Tangent vectors
        Real xt =  a_ * cu * cv;
        Real yt =  b_ * cu * sv;
        Real zt = -c_ * su;
        Real xp = -a_ * su * sv;
        Real yp =  b_ * su * cv;
        Real zp =  Real(0);

        // First fundamental form
        Real E = xt * xt + yt * yt + zt * zt;
        Real G = xp * xp + yp * yp + zp * zp;
        Real F = xt * xp + yt * yp + zt * zp;
        pt.dA = tm_sqrt(E * G - F * F);

        // Normalized tangent vectors
        Real t1_mag = tm_sqrt(E);
        Real tau1x = xt / t1_mag;
        Real tau1y = yt / t1_mag;
        Real tau1z = zt / t1_mag;

        Real t2_mag = tm_sqrt(G);
        Real tau2x, tau2y, tau2z;
        if (t2_mag < Real(1e-30)) {
            tau2x = Real(0);
            tau2y = Real(1);
            tau2z = Real(0);
        } else {
            tau2x = xp / t2_mag;
            tau2y = yp / t2_mag;
            tau2z = Real(0);
        }

        // Cartesian normal = tau1 x tau2 (normalized)
        Real nc_x = tau1y * tau2z - tau1z * tau2y;
        Real nc_y = tau1z * tau2x - tau1x * tau2z;
        Real nc_z = tau1x * tau2y - tau1y * tau2x;
        Real nc_mag = tm_sqrt(nc_x * nc_x + nc_y * nc_y + nc_z * nc_z);
        nc_x /= nc_mag;
        nc_y /= nc_mag;
        nc_z /= nc_mag;

        // Convert normal to spherical coordinates
        normal_cart_to_sph(pt.theta, pt.phi, nc_x, nc_y, nc_z,
                          pt.n_r, pt.n_theta, pt.n_phi);

        return pt;
    }

    void quadrature_points(int Nint1, int Nint2, bool mirror,
                           std::vector<std::vector<Real>>& param1G,
                           std::vector<std::vector<Real>>& param2G,
                           std::vector<std::vector<Real>>& weightsG) const override {
        param1G.resize(1);
        param2G.resize(1);
        weightsG.resize(1);

        int N1 = mirror ? (Nint1 / 2 + 1) : Nint1;
        int N2 = Nint2;

        std::vector<Real> ut, wt, vp, wp;
        Real u_lo = Real(0);
        Real u_hi = mirror ? pi<Real>() / Real(2) : pi<Real>();
        gauss_legendre<Real>(u_lo, u_hi, N1, wt, ut);
        gauss_legendre<Real>(Real(0), Real(2) * pi<Real>(), N2, wp, vp);

        int Ntotal = N1 * N2;
        param1G[0].resize(Ntotal);
        param2G[0].resize(Ntotal);
        weightsG[0].resize(Ntotal);

        for (int j = 0; j < N1; j++) {
            for (int i = 0; i < N2; i++) {
                int idx = j * N2 + i;
                param1G[0][idx] = ut[j];   // u
                param2G[0][idx] = vp[i];   // v
                weightsG[0][idx] = wt[j] * wp[i];
            }
        }
    }

private:
    Real a_;  // semi-axis along x
    Real b_;  // semi-axis along y
    Real c_;  // semi-axis along z
};

// ============================================================================
// Regular N-hedral Prism: 3D geometry
// half-side a, half-height l, N sides
// Patches: lateral face (1), top cap (1), bottom cap (1, only if !mirror)
// ============================================================================
template<typename Real>
class RegularPrismGeometry3D : public NonAxialGeometry<Real> {
public:
    RegularPrismGeometry3D(int N, Real half_side, Real half_height)
        : N_(N), a_(half_side), l_(half_height) {}

    int num_patches() const override {
        // mirror: 2 patches (lateral + top); no mirror: 3 (+bottom)
        // Always return 3; the quadrature_points method handles the mirror case
        return 3;
    }

    bool is_mirror_symmetric() const override { return true; }
    int azimuthal_symmetry() const override { return N_; }

    GeometryPoint3D<Real> evaluate(Real param1, Real param2, int patch) const override {
        GeometryPoint3D<Real> pt;
        Real alpha = pi<Real>() / Real(N_);
        Real b = a_ / tm_tan(alpha);

        Real x, y, z, nc_x, nc_y, nc_z;

        if (patch == 0) {
            // Lateral face: param1 = z, param2 = side_coord
            x = param2;
            y = b;
            z = param1;
            pt.dA = Real(1);
            nc_x = Real(0);
            nc_y = Real(1);
            nc_z = Real(0);
        } else if (patch == 1) {
            // Top cap: param1 = phi, param2 = rho
            Real ro = param2;
            Real phi_loc = param1;
            x = ro * tm_sin(phi_loc);
            y = ro * tm_cos(phi_loc);
            z = l_;
            pt.dA = ro;
            nc_x = Real(0);
            nc_y = Real(0);
            nc_z = Real(1);
        } else {
            // Bottom cap: param1 = phi, param2 = rho
            Real ro = param2;
            Real phi_loc = param1;
            x =  ro * tm_sin(phi_loc);
            y =  ro * tm_cos(phi_loc);
            z = -l_;
            pt.dA = ro;
            nc_x =  Real(0);
            nc_y =  Real(0);
            nc_z = -Real(1);
        }

        cartesian_to_spherical(x, y, z, pt.r, pt.theta, pt.phi);
        normal_cart_to_sph(pt.theta, pt.phi, nc_x, nc_y, nc_z,
                          pt.n_r, pt.n_theta, pt.n_phi);

        return pt;
    }

    void quadrature_points(int Nint1, int Nint2, bool mirror,
                           std::vector<std::vector<Real>>& param1G,
                           std::vector<std::vector<Real>>& param2G,
                           std::vector<std::vector<Real>>& weightsG) const override {
        Real alpha = pi<Real>() / Real(N_);
        Real b = a_ / tm_tan(alpha);

        int Npatch = mirror ? 2 : 3;
        param1G.resize(Npatch);
        param2G.resize(Npatch);
        weightsG.resize(Npatch);

        // Lateral face (patch 0): param1=z, param2=side_coord
        int N1_lat = mirror ? (Nint1 / 2 + 1) : Nint1;
        int N2_lat = Nint2;

        std::vector<Real> zt_nodes, zt_weights, st_nodes, st_weights;
        if (mirror) {
            gauss_legendre<Real>(Real(0), l_, N1_lat, zt_weights, zt_nodes);
        } else {
            gauss_legendre<Real>(-l_, l_, N1_lat, zt_weights, zt_nodes);
        }
        gauss_legendre<Real>(-a_, a_, N2_lat, st_weights, st_nodes);

        int Nlat = N1_lat * N2_lat;
        param1G[0].resize(Nlat);
        param2G[0].resize(Nlat);
        weightsG[0].resize(Nlat);

        for (int i = 0; i < N1_lat; i++) {
            for (int j = 0; j < N2_lat; j++) {
                int idx = i * N2_lat + j;
                param1G[0][idx] = zt_nodes[i];   // z
                param2G[0][idx] = st_nodes[j];    // side_coord
                weightsG[0][idx] = zt_weights[i] * st_weights[j];
            }
        }

        // Cap: param1=phi, param2=rho
        // phi in [-alpha, alpha], rho in [0, b/cos(phi)] (phi-dependent!)
        int N2_cap = Nint2;
        int N3_cap = Nint2;  // same as N2 per Fortran convention

        std::vector<Real> phi_nodes, phi_weights;
        gauss_legendre<Real>(-alpha, alpha, N2_cap, phi_weights, phi_nodes);

        int Ncap = N2_cap * N3_cap;

        // Top cap (patch 1)
        param1G[1].resize(Ncap);
        param2G[1].resize(Ncap);
        weightsG[1].resize(Ncap);

        for (int i = 0; i < N2_cap; i++) {
            Real phi_val = phi_nodes[i];
            Real ro_max = b / tm_cos(phi_val);

            std::vector<Real> rho_nodes, rho_weights;
            gauss_legendre<Real>(Real(0), ro_max, N3_cap, rho_weights, rho_nodes);

            for (int j = 0; j < N3_cap; j++) {
                int idx = i * N3_cap + j;
                param1G[1][idx] = phi_nodes[i];   // phi
                param2G[1][idx] = rho_nodes[j];    // rho
                weightsG[1][idx] = phi_weights[i] * rho_weights[j];
            }
        }

        // Bottom cap (patch 2), only if !mirror
        if (!mirror) {
            param1G[2].resize(Ncap);
            param2G[2].resize(Ncap);
            weightsG[2].resize(Ncap);

            for (int i = 0; i < Ncap; i++) {
                param1G[2][i] = param1G[1][i];
                param2G[2][i] = param2G[1][i];
                weightsG[2][i] = weightsG[1][i];
            }
        }
    }

    Real half_side() const { return a_; }
    Real half_height() const { return l_; }
    int sides() const { return N_; }

private:
    int N_;    // number of sides
    Real a_;   // half-side length
    Real l_;   // half-height

    static Real tm_tan(Real x) { return tm_sin(x) / tm_cos(x); }
};

} // namespace tmatrix

#endif
