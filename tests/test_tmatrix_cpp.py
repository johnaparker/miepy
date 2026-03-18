"""Tests for the C++ EBCM T-matrix solver.

Verifies correctness of:
- Sphere degenerate case (spheroid with equal axes matches Mie theory)
- Cylinder geometry
- PEC (conducting) particles
- Non-degenerate spheroids
- Q31 conditioning diagnostics
"""

import warnings

import numpy as np
import pytest

import miepy

nm = 1e-9

Ag = miepy.materials.Ag()
metal = miepy.materials.metal()
eps_b = 1.5
medium = miepy.constant_material(index=eps_b**2)

radius = 75 * nm
wavelength = 600 * nm


class TestSphereDegenerate:
    """Spheroid with equal axes should match Mie theory T-matrix exactly."""

    @pytest.mark.parametrize("lmax", [1, 2, 3, 4])
    def test_dielectric_sphere(self, lmax):
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T_mie = miepy.tmatrix.tmatrix_sphere(radius, wavelength, eps, eps_m, lmax)
        T_ebcm = miepy.tmatrix.tmatrix_spheroid(
            radius, radius, wavelength, eps, eps_m, lmax, use_ds=False
        )

        assert np.allclose(T_mie, T_ebcm, rtol=0, atol=1e-10), (
            f"lmax={lmax}: max error = {np.max(np.abs(T_mie - T_ebcm))}"
        )

    def test_conducting_sphere(self):
        lmax = 2
        eps = metal.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T_mie = miepy.tmatrix.tmatrix_sphere(
            radius, wavelength, eps, eps_m, lmax, conducting=True
        )
        T_ebcm = miepy.tmatrix.tmatrix_spheroid(
            radius, radius, wavelength, eps, eps_m, lmax, conducting=True, use_ds=False
        )

        assert np.allclose(T_mie, T_ebcm, rtol=0, atol=5e-9), (
            f"max error = {np.max(np.abs(T_mie - T_ebcm))}"
        )


class TestTmatrixDiagonal:
    """For a sphere, T-matrix should be diagonal with only (n,m)=(n,m) entries."""

    def test_sphere_diagonal(self):
        lmax = 3
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_spheroid(radius, radius, wavelength, eps, eps_m, lmax, use_ds=False)

        # Off-diagonal in polarization should be zero for sphere
        assert np.allclose(T[0, :, 1, :], 0, atol=1e-12)
        assert np.allclose(T[1, :, 0, :], 0, atol=1e-12)

        # Within each polarization, should be diagonal
        rmax = lmax * (lmax + 2)
        for a in range(2):
            for i in range(rmax):
                for j in range(rmax):
                    if i != j:
                        assert abs(T[a, i, a, j]) < 1e-12, (
                            f"T[{a},{i},{a},{j}] = {T[a,i,a,j]} (should be zero)"
                        )


class TestNonDegenerateSpheroid:
    """Non-degenerate spheroids should produce physically reasonable T-matrices."""

    def test_prolate_spheroid_nonzero(self):
        """Prolate spheroid (axis_z > axis_xy) should have nonzero T-matrix."""
        lmax = 2
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_spheroid(
            radius, 2 * radius, wavelength, eps, eps_m, lmax
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_oblate_spheroid_nonzero(self):
        """Oblate spheroid (axis_xy > axis_z) should have nonzero T-matrix."""
        lmax = 2
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_spheroid(
            2 * radius, radius, wavelength, eps, eps_m, lmax
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_prolate_off_diagonal(self):
        """Prolate spheroid should have off-diagonal elements (not a sphere)."""
        lmax = 3
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_spheroid(
            radius, 2 * radius, wavelength, eps, eps_m, lmax
        )

        # For a non-spherical particle, some off-diagonal elements should be nonzero
        rmax = lmax * (lmax + 2)
        off_diag_norm = 0
        for a in range(2):
            for i in range(rmax):
                for j in range(rmax):
                    if i != j:
                        off_diag_norm += abs(T[a, i, a, j]) ** 2
        off_diag_norm = np.sqrt(off_diag_norm)

        assert off_diag_norm > 1e-10, "Non-spherical particle has no off-diagonal T-matrix elements"


class TestCylinder:
    """Tests for cylinder T-matrix computation."""

    def test_cylinder_nonzero(self):
        lmax = 2
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_cylinder(
            radius, 2 * radius, wavelength, eps, eps_m, lmax
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_cylinder_conducting(self):
        lmax = 2
        eps = metal.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_cylinder(
            radius, 2 * radius, wavelength, eps, eps_m, lmax, conducting=True
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"


class TestOpticalTheorem:
    """Extinction = scattering + absorption (optical theorem check via T-matrix)."""

    def test_passive_particle(self):
        """For a passive particle, absorption should be non-negative."""
        lmax = 3
        source = miepy.sources.plane_wave([1, 0])

        cluster = miepy.cluster(
            particles=miepy.spheroid([0, 0, 0], radius, 2 * radius, Ag),
            source=source,
            wavelength=wavelength,
            lmax=lmax,
            medium=medium,
        )

        C = cluster.cross_sections()
        # Allow tiny negative absorption from numerical noise (absolute tolerance)
        assert C.absorption >= -1e-12, f"Absorption is negative: {C.absorption}"
        assert abs(C.extinction - C.scattering - C.absorption) < 1e-12, (
            "Optical theorem violated"
        )


class TestDistributedSources:
    """Tests for distributed sources (DS) in EBCM T-matrix computation."""

    def test_ds_prolate_valid(self):
        """DS should produce a valid T-matrix for a 2:1 prolate spheroid."""
        lmax = 3
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_spheroid(
            radius, 2 * radius, wavelength, eps, eps_m, lmax, use_ds=True
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_ds_oblate_valid(self):
        """DS should produce a valid T-matrix for a 1:2 oblate spheroid."""
        lmax = 3
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_spheroid(
            2 * radius, radius, wavelength, eps, eps_m, lmax, use_ds=True
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_ds_vs_localized_moderate_aspect(self):
        """DS T-matrix should closely match localized T-matrix for a 2:1 prolate dielectric."""
        lmax = 3
        eps_diel = complex(4.0)
        eps_m = medium.eps(wavelength)

        T_loc = miepy.tmatrix.tmatrix_spheroid(
            radius, 2 * radius, wavelength, eps_diel, eps_m, lmax, use_ds=False
        )
        T_ds = miepy.tmatrix.tmatrix_spheroid(
            radius, 2 * radius, wavelength, eps_diel, eps_m, lmax, use_ds=True
        )

        frob_error = np.linalg.norm(T_ds - T_loc) / np.linalg.norm(T_loc)
        assert frob_error < 0.10, (
            f"DS vs localized Frobenius error = {frob_error:.4f} (target: < 10%)"
        )

    def test_ds_sphere_matches_mie(self):
        """Sphere with DS enabled should approximately match Mie theory.

        DS is designed for non-spherical particles; for spheres it introduces
        small numerical errors from the distributed source basis. The DS
        implementation shows ~1e-2 max error vs Mie for this case.
        """
        lmax = 3
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T_mie = miepy.tmatrix.tmatrix_sphere(radius, wavelength, eps, eps_m, lmax)
        T_ds = miepy.tmatrix.tmatrix_spheroid(
            radius, radius, wavelength, eps, eps_m, lmax, use_ds=True
        )

        max_error = np.max(np.abs(T_mie - T_ds))
        assert max_error < 0.1, (
            f"DS sphere vs Mie max error = {max_error:.6f} (target: < 0.1)"
        )

    def test_high_aspect_ratio_convergence(self):
        """A 10:1 prolate spheroid with DS should produce a valid T-matrix."""
        lmax = 3
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_spheroid(
            radius / 3, 10 * radius / 3, wavelength, eps, eps_m, lmax, use_ds=True
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

        # Optical theorem check: use in a cluster
        source = miepy.sources.plane_wave([1, 0])
        cluster = miepy.cluster(
            particles=miepy.spheroid([0, 0, 0], radius / 3, 10 * radius / 3, Ag),
            source=source,
            wavelength=wavelength,
            lmax=lmax,
            medium=medium,
        )

        C = cluster.cross_sections()
        assert C.absorption >= -1e-12, f"Absorption is negative: {C.absorption}"
        assert abs(C.extinction - C.scattering - C.absorption) < 1e-12, (
            "Optical theorem violated"
        )


has_quad = miepy.cpp.tmatrix.has_extended_precision


class TestExtendedPrecisionAvailability:
    """Verify the has_extended_precision attribute is exposed correctly."""

    def test_cpp_attribute_is_bool(self):
        assert isinstance(miepy.cpp.tmatrix.has_extended_precision, bool)

    def test_python_function_matches(self):
        assert miepy.has_extended_precision() == miepy.cpp.tmatrix.has_extended_precision


@pytest.mark.skipif(has_quad, reason="Platform has quad — test fallback on platforms without it")
class TestExtendedPrecisionFallback:
    """On platforms without __float128, extended_precision=True should raise."""

    def test_spheroid_raises(self):
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)
        with pytest.raises(RuntimeError, match="not available"):
            miepy.tmatrix.tmatrix_spheroid(
                radius, 2 * radius, wavelength, eps, eps_m, 2, extended_precision=True
            )

    def test_cylinder_raises(self):
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)
        with pytest.raises(RuntimeError, match="not available"):
            miepy.tmatrix.tmatrix_cylinder(
                radius, 2 * radius, wavelength, eps, eps_m, 2, extended_precision=True
            )


@pytest.mark.skipif(not has_quad, reason="__float128 not available")
class TestExtendedPrecision:
    """Verify quad produces valid T-matrices and matches double for easy cases."""

    def test_spheroid_valid(self):
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_spheroid(
            radius, 2 * radius, wavelength, eps, eps_m, 3, extended_precision=True
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_cylinder_valid(self):
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_cylinder(
            radius, 2 * radius, wavelength, eps, eps_m, 3, extended_precision=True
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_quad_matches_double_easy_case(self):
        """For a well-conditioned case, quad and double should agree closely."""
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T_double = miepy.tmatrix.tmatrix_spheroid(
            radius, 1.5 * radius, wavelength, eps, eps_m, 2, extended_precision=False
        )
        T_quad = miepy.tmatrix.tmatrix_spheroid(
            radius, 1.5 * radius, wavelength, eps, eps_m, 2, extended_precision=True
        )

        assert np.allclose(T_double, T_quad, rtol=1e-10), (
            f"max relative error = {np.max(np.abs(T_double - T_quad) / np.abs(T_quad).clip(1e-30))}"
        )


@pytest.mark.skipif(not has_quad, reason="__float128 not available")
class TestExtendedPrecisionChallengingCase:
    """High-aspect-ratio metallic spheroid — quad should satisfy optical theorem."""

    def test_high_aspect_ratio_metal_optical_theorem(self):
        """Quad should give non-negative absorption for a challenging case."""
        lmax = 4
        source = miepy.sources.plane_wave([1, 0])

        # 5:1 prolate metallic spheroid — a challenging case for double
        cluster = miepy.cluster(
            particles=miepy.spheroid(
                [0, 0, 0], radius / 2, 5 * radius / 2, Ag,
                extended_precision=True,
            ),
            source=source,
            wavelength=wavelength,
            lmax=lmax,
            medium=medium,
        )

        C = cluster.cross_sections()
        assert not np.isnan(C.extinction), "Extinction is NaN"
        assert not np.isnan(C.scattering), "Scattering is NaN"
        assert not np.isnan(C.absorption), "Absorption is NaN"
        assert C.absorption >= -1e-12, f"Absorption is negative: {C.absorption}"
        assert abs(C.extinction - C.scattering - C.absorption) < 1e-12, (
            "Optical theorem violated"
        )


@pytest.mark.skipif(not has_quad, reason="__float128 not available")
class TestExtendedPrecisionParticleAPI:
    """Verify extended_precision flows through particle constructors correctly."""

    def test_spheroid_default_false(self):
        p = miepy.spheroid([0, 0, 0], radius, 2 * radius, Ag)
        assert p.extended_precision is False

    def test_spheroid_accepts_true(self):
        p = miepy.spheroid([0, 0, 0], radius, 2 * radius, Ag, extended_precision=True)
        assert p.extended_precision is True

    def test_cylinder_default_false(self):
        p = miepy.cylinder([0, 0, 0], radius, 2 * radius, Ag)
        assert p.extended_precision is False

    def test_cylinder_accepts_true(self):
        p = miepy.cylinder([0, 0, 0], radius, 2 * radius, Ag, extended_precision=True)
        assert p.extended_precision is True

    def test_spheroid_dict_key_differs(self):
        p1 = miepy.spheroid([0, 0, 0], radius, 2 * radius, Ag, extended_precision=False)
        p2 = miepy.spheroid([0, 0, 0], radius, 2 * radius, Ag, extended_precision=True)
        assert p1._dict_key(wavelength) != p2._dict_key(wavelength)

    def test_cylinder_dict_key_differs(self):
        p1 = miepy.cylinder([0, 0, 0], radius, 2 * radius, Ag, extended_precision=False)
        p2 = miepy.cylinder([0, 0, 0], radius, 2 * radius, Ag, extended_precision=True)
        assert p1._dict_key(wavelength) != p2._dict_key(wavelength)

    def test_spheroid_cluster_computation(self):
        """Extended precision spheroid should work in a full cluster computation."""
        source = miepy.sources.plane_wave([1, 0])
        cluster = miepy.cluster(
            particles=miepy.spheroid(
                [0, 0, 0], radius, 2 * radius, Ag, extended_precision=True
            ),
            source=source,
            wavelength=wavelength,
            lmax=2,
            medium=medium,
        )

        C = cluster.cross_sections()
        assert C.scattering > 0
        assert C.extinction > 0


class TestConditioningWarning:
    """Tests for Q31 conditioning diagnostics in EBCM T-matrix."""

    def test_well_conditioned_no_warning(self):
        """A moderate dielectric sphere should not trigger a conditioning warning."""
        eps_diel = complex(4.0)
        eps_m = medium.eps(wavelength)

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            T = miepy.tmatrix.tmatrix_spheroid(
                radius, radius, wavelength, eps_diel, eps_m, 2, use_ds=False
            )
            assert not np.any(np.isnan(T))

    def test_ill_conditioned_warning(self):
        """A very high aspect ratio metallic spheroid at high lmax without DS
        should trigger a RuntimeWarning about Q31 conditioning."""
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        # 200:1 aspect ratio at lmax=6 without distributed sources
        axis_xy = radius / 10
        axis_z = axis_xy * 200
        with pytest.warns(RuntimeWarning, match="ill-conditioned"):
            miepy.tmatrix.tmatrix_spheroid(
                axis_z, axis_xy, wavelength, eps, eps_m, 6, use_ds=False
            )


# ============================================================================
# Non-axisymmetric T-matrix tests
# ============================================================================


class TestEllipsoidDegenerate:
    """Ellipsoid with rx=ry should match axisymmetric spheroid solver."""

    @pytest.mark.parametrize("lmax", [1, 2, 3])
    def test_ellipsoid_matches_spheroid(self, lmax):
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T_sph = miepy.tmatrix.tmatrix_spheroid(
            radius, 2 * radius, wavelength, eps, eps_m, lmax, use_ds=False
        )
        T_ell = miepy.tmatrix.tmatrix_ellipsoid(
            radius, radius, 2 * radius, wavelength, eps, eps_m, lmax,
            Nint1=100, Nint2=100,
        )

        assert np.allclose(T_sph, T_ell, rtol=0, atol=1e-10), (
            f"lmax={lmax}: max error = {np.max(np.abs(T_sph - T_ell))}"
        )

    def test_ellipsoid_matches_sphere(self):
        """Ellipsoid with rx=ry=rz should match Mie theory."""
        lmax = 2
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T_mie = miepy.tmatrix.tmatrix_sphere(radius, wavelength, eps, eps_m, lmax)
        T_ell = miepy.tmatrix.tmatrix_ellipsoid(
            radius, radius, radius, wavelength, eps, eps_m, lmax,
            Nint1=100, Nint2=100,
        )

        assert np.allclose(T_mie, T_ell, rtol=0, atol=1e-10), (
            f"max error = {np.max(np.abs(T_mie - T_ell))}"
        )


class TestEllipsoidNonDegenerate:
    """True 3-axis ellipsoid (rx != ry != rz) should produce valid results."""

    def test_ellipsoid_valid(self):
        lmax = 2
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_ellipsoid(
            60 * nm, 80 * nm, 100 * nm, wavelength, eps, eps_m, lmax,
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_ellipsoid_off_diagonal(self):
        """True ellipsoid should have cross-polarization elements."""
        lmax = 3
        eps = complex(4.0)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_ellipsoid(
            60 * nm, 80 * nm, 100 * nm, wavelength, eps, eps_m, lmax,
            Nint1=80, Nint2=80,
        )

        # Cross-polarization should be nonzero for a true ellipsoid
        cross_pol = np.max(np.abs(T[0, :, 1, :]))
        assert cross_pol > 1e-10, (
            f"Cross-polarization is too small: {cross_pol}"
        )


class TestRegularPrismBasic:
    """Basic validity tests for regular prism T-matrices."""

    def test_hexagonal_prism_valid(self):
        lmax = 2
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_regular_prism(
            6, 100 * nm, 150 * nm, wavelength, eps, eps_m, lmax,
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"

    def test_cube_valid(self):
        lmax = 2
        eps = Ag.eps(wavelength)
        eps_m = medium.eps(wavelength)

        T = miepy.tmatrix.tmatrix_square_prism(
            100 * nm, 100 * nm, wavelength, eps, eps_m, lmax,
        )

        assert not np.any(np.isnan(T)), "T-matrix contains NaN"
        assert not np.any(np.isinf(T)), "T-matrix contains Inf"
        assert np.max(np.abs(T)) > 0, "T-matrix is all zeros"


class TestNonaxialOpticalTheorem:
    """Optical theorem checks for non-axisymmetric particles."""

    def test_ellipsoid_optical_theorem(self):
        """Extinction = scattering + absorption for a dielectric ellipsoid."""
        lmax = 3
        source = miepy.sources.plane_wave([1, 0])

        cluster = miepy.cluster(
            particles=miepy.ellipsoid(
                [0, 0, 0], 60 * nm, 80 * nm, 100 * nm,
                miepy.constant_material(index=2.0),
            ),
            source=source,
            wavelength=wavelength,
            lmax=lmax,
            medium=medium,
        )

        C = cluster.cross_sections()
        assert C.absorption >= -1e-12, f"Absorption is negative: {C.absorption}"
        assert abs(C.extinction - C.scattering - C.absorption) < 1e-12, (
            "Optical theorem violated"
        )

    def test_cube_optical_theorem(self):
        """Optical theorem for a dielectric cube."""
        lmax = 2
        source = miepy.sources.plane_wave([1, 0])

        cluster = miepy.cluster(
            particles=miepy.cube(
                [0, 0, 0], 100 * nm,
                miepy.constant_material(index=2.0),
            ),
            source=source,
            wavelength=wavelength,
            lmax=lmax,
            medium=medium,
        )

        C = cluster.cross_sections()
        assert C.absorption >= -1e-12, f"Absorption is negative: {C.absorption}"
        assert abs(C.extinction - C.scattering - C.absorption) < 1e-12, (
            "Optical theorem violated"
        )
