"""Tests for the C++ EBCM T-matrix solver.

Verifies correctness of:
- Sphere degenerate case (spheroid with equal axes matches Mie theory)
- Cylinder geometry
- PEC (conducting) particles
- Non-degenerate spheroids
"""

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
        small numerical errors from the distributed source basis. Both C++ and
        Fortran DS implementations show ~1e-2 max error vs Mie for this case.
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
