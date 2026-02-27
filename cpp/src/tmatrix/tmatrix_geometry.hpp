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

} // namespace tmatrix

#endif
