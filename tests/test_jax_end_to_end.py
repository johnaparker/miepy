"""End-to-end tests comparing JAX backend against C++ backend for sphere_cluster.

Runs representative sphere_cluster computations with both backends and verifies
that all observables (cross-sections, forces, torques, expansion coefficients)
agree to within solver tolerance.
"""

import numpy as np
import pytest

jax = pytest.importorskip('jax')

import miepy

nm = 1e-9

Ag = miepy.materials.Ag()
radius = 75 * nm
source_rhc = miepy.sources.plane_wave.from_string(polarization="rhc")
source_x = miepy.sources.plane_wave.from_string(polarization="x")


def _solve_both_backends(positions, radius, material, source, wavelength, lmax,
                         interactions=True):
    """Solve the same problem with C++ and JAX backends, return both clusters."""
    # C++ backend (default)
    with miepy.backends.backend('cpp'):
        cpp_cluster = miepy.sphere_cluster(
            position=positions,
            radius=radius,
            material=material,
            source=source,
            wavelength=wavelength,
            lmax=lmax,
            interactions=interactions,
        )

    # JAX backend
    with miepy.backends.backend('jax'):
        jax_cluster = miepy.sphere_cluster(
            position=positions,
            radius=radius,
            material=material,
            source=source,
            wavelength=wavelength,
            lmax=lmax,
            interactions=interactions,
        )

    return cpp_cluster, jax_cluster


# ---------------------------------------------------------------------------
# Single particle (no interactions)
# ---------------------------------------------------------------------------

class TestSingleParticle:
    """Single sphere — no interaction matrix, tests source decomposition + Mie coefficients."""

    def test_cross_sections(self):
        cpp, jax_c = _solve_both_backends(
            [0, 0, 0], radius, Ag, source_x, 600*nm, lmax=2, interactions=False)
        C_cpp = cpp.cross_sections()
        C_jax = jax_c.cross_sections()
        for a, b in zip(C_cpp, C_jax):
            np.testing.assert_allclose(a, b, rtol=1e-10)

    def test_expansion_coefficients(self):
        cpp, jax_c = _solve_both_backends(
            [0, 0, 0], radius, Ag, source_x, 600*nm, lmax=2, interactions=False)
        np.testing.assert_allclose(jax_c.p_src, cpp.p_src, rtol=1e-10)
        np.testing.assert_allclose(jax_c.p_inc, cpp.p_inc, rtol=1e-10)
        np.testing.assert_allclose(jax_c.p_scat, cpp.p_scat, rtol=1e-10)


# ---------------------------------------------------------------------------
# Two particles — exact solver
# ---------------------------------------------------------------------------

class TestTwoParticlesExact:
    """Two-sphere dimer with exact (LU) solver."""

    @pytest.fixture(scope="class")
    def clusters(self):
        positions = [[100*nm, 0, 0], [-100*nm, 0, 0]]
        return _solve_both_backends(
            positions, radius, Ag, source_rhc, 600*nm, lmax=2)

    def test_p_inc(self, clusters):
        cpp, jax_c = clusters
        np.testing.assert_allclose(jax_c.p_inc, cpp.p_inc, rtol=1e-10)

    def test_p_scat(self, clusters):
        cpp, jax_c = clusters
        np.testing.assert_allclose(jax_c.p_scat, cpp.p_scat, rtol=1e-10)

    def test_cross_sections(self, clusters):
        cpp, jax_c = clusters
        C_cpp = cpp.cross_sections()
        C_jax = jax_c.cross_sections()
        for a, b in zip(C_cpp, C_jax):
            np.testing.assert_allclose(a, b, rtol=1e-10)

    def test_per_particle_cross_sections(self, clusters):
        cpp, jax_c = clusters
        for i in range(2):
            C_cpp = cpp.cross_sections_of_particle(i)
            C_jax = jax_c.cross_sections_of_particle(i)
            for a, b in zip(C_cpp, C_jax):
                np.testing.assert_allclose(a, b, rtol=1e-10)

    def test_force(self, clusters):
        cpp, jax_c = clusters
        for i in range(2):
            F_cpp = cpp.force_on_particle(i)
            F_jax = jax_c.force_on_particle(i)
            np.testing.assert_allclose(F_jax, F_cpp, rtol=1e-8, atol=1e-30)

    def test_torque(self, clusters):
        cpp, jax_c = clusters
        for i in range(2):
            T_cpp = cpp.torque_on_particle(i)
            T_jax = jax_c.torque_on_particle(i)
            np.testing.assert_allclose(T_jax, T_cpp, rtol=1e-8, atol=1e-30)


# ---------------------------------------------------------------------------
# Multi-particle — BiCGSTAB solver
# ---------------------------------------------------------------------------

class TestMultiParticleBiCGSTAB:
    """Five-sphere system with BiCGSTAB iterative solver."""

    @pytest.fixture(scope="class")
    def clusters(self):
        positions = [
            [0, 0, 0],
            [250*nm, 0, 0],
            [-250*nm, 0, 0],
            [0, 250*nm, 0],
            [0, -250*nm, 0],
        ]
        return _solve_both_backends(
            positions, radius, Ag, source_rhc, 600*nm, lmax=2)

    def test_p_inc(self, clusters):
        cpp, jax_c = clusters
        # BiCGSTAB tolerance allows some relative error; use atol for near-zero entries
        np.testing.assert_allclose(jax_c.p_inc, cpp.p_inc, rtol=1e-4, atol=1e-14)

    def test_cross_sections(self, clusters):
        cpp, jax_c = clusters
        C_cpp = cpp.cross_sections()
        C_jax = jax_c.cross_sections()
        for a, b in zip(C_cpp, C_jax):
            np.testing.assert_allclose(a, b, rtol=1e-4)

    def test_force(self, clusters):
        cpp, jax_c = clusters
        for i in range(5):
            F_cpp = cpp.force_on_particle(i)
            F_jax = jax_c.force_on_particle(i)
            np.testing.assert_allclose(F_jax, F_cpp, rtol=1e-3, atol=1e-30)


# ---------------------------------------------------------------------------
# Higher lmax
# ---------------------------------------------------------------------------

class TestHigherLmax:
    """Two particles at lmax=4 to test correctness at higher expansion orders."""

    @pytest.fixture(scope="class")
    def clusters(self):
        positions = [[125*nm, 0, 0], [-125*nm, 0, 0]]
        return _solve_both_backends(
            positions, radius, Ag, source_x, 600*nm, lmax=4)

    def test_p_inc(self, clusters):
        cpp, jax_c = clusters
        # Use atol for near-zero entries where relative error is meaningless
        np.testing.assert_allclose(jax_c.p_inc, cpp.p_inc, rtol=1e-9, atol=1e-14)

    def test_cross_sections(self, clusters):
        cpp, jax_c = clusters
        C_cpp = cpp.cross_sections()
        C_jax = jax_c.cross_sections()
        for a, b in zip(C_cpp, C_jax):
            np.testing.assert_allclose(a, b, rtol=1e-9)


# ---------------------------------------------------------------------------
# Displaced single particle (origin independence)
# ---------------------------------------------------------------------------

class TestDisplacedParticle:
    """A single sphere displaced from origin — cross-sections should be origin-independent."""

    def test_origin_independence(self):
        positions_origin = [0, 0, 0]
        positions_displaced = [40*nm, 50*nm, 60*nm]

        with miepy.backends.backend('jax'):
            c1 = miepy.sphere_cluster(
                position=positions_origin, radius=radius, material=Ag,
                source=source_x, wavelength=600*nm, lmax=2)
            C1 = c1.cross_sections()

            c2 = miepy.sphere_cluster(
                position=positions_displaced, radius=radius, material=Ag,
                source=source_x, wavelength=600*nm, lmax=2)
            C2 = c2.cross_sections()

        for a, b in zip(C1, C2):
            np.testing.assert_allclose(a, b, rtol=1e-12)


# ---------------------------------------------------------------------------
# Interactions off — N * single Mie
# ---------------------------------------------------------------------------

class TestInteractionsOff:
    """Two spheres with interactions disabled should give 2x single Mie cross-sections."""

    def test_extinction_equals_2x_single(self):
        sep = 200 * nm
        wavelength = 600 * nm

        with miepy.backends.backend('jax'):
            sol = miepy.sphere_cluster(
                position=[[sep/2, 0, 0], [-sep/2, 0, 0]],
                radius=radius, material=Ag, source=source_x,
                wavelength=wavelength, lmax=2, interactions=False)
            _, _, E_cluster = sol.cross_sections()

        sphere = miepy.single_mie_sphere(radius, Ag, np.array([wavelength]), lmax=2)
        _, _, E_single = sphere.cross_sections()

        np.testing.assert_allclose(E_cluster, 2 * E_single, rtol=1e-12)
