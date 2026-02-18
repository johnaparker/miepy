"""Tests for JAX backend special functions against C++ reference implementations."""

import numpy as np
import pytest

jax = pytest.importorskip('jax')
jnp = jax.numpy

from miepy.backends.jax import special
from miepy.backends.jax import mie as jax_mie
import miepy.cpp.special as cpp_special
import miepy.cpp.vsh_functions as cpp_vsh
from miepy.mie_single.mie_sphere import (
    mie_sphere_scattering_coefficients as np_mie_scat,
    mie_sphere_interior_coefficients as np_mie_int,
)


# ---------------------------------------------------------------------------
# Spherical Bessel j_n
# ---------------------------------------------------------------------------

class TestSphericalJn:
    """Test spherical_jn against C++ spherical_jn."""

    @pytest.mark.parametrize("n", [0, 1, 2, 5, 10])
    def test_real_z(self, n):
        z = np.linspace(0.5, 20.0, 50)
        z_complex = z.astype(complex)
        jax_vals = np.asarray(special.spherical_jn(n, jnp.array(z_complex)))
        cpp_vals = np.array([cpp_special.spherical_jn(n, complex(zi)) for zi in z_complex])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12)

    @pytest.mark.parametrize("n", [0, 1, 3, 7])
    def test_complex_z(self, n):
        # Use |z| > n to avoid forward-recursion instability (C++ has same issue)
        z = np.array([3.0 + 1.0j, 5.0 - 2.0j, 10.0 + 3.0j, 15.0 + 1.0j])
        jax_vals = np.asarray(special.spherical_jn(n, jnp.array(z)))
        cpp_vals = np.array([cpp_special.spherical_jn(n, complex(zi)) for zi in z])
        # Slightly relaxed tolerance for complex z due to FP operation ordering
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-11)

    @pytest.mark.parametrize("n", [1, 2, 5])
    def test_derivative(self, n):
        z = np.array([1.0 + 0.0j, 3.0 + 0.0j, 8.0 + 0.0j])
        jax_vals = np.asarray(special.spherical_jn(n, jnp.array(z), derivative=True))
        cpp_vals = np.array([cpp_special.spherical_jn(n, complex(zi), True) for zi in z])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12)


# ---------------------------------------------------------------------------
# Spherical Bessel y_n
# ---------------------------------------------------------------------------

class TestSphericalYn:
    @pytest.mark.parametrize("n", [0, 1, 2, 5, 10])
    def test_values(self, n):
        z = np.linspace(0.5, 20.0, 50)
        jax_vals = np.asarray(special.spherical_yn(n, jnp.array(z)))
        cpp_vals = np.array([cpp_special.spherical_yn(n, float(zi)) for zi in z])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12)

    @pytest.mark.parametrize("n", [1, 3, 7])
    def test_derivative(self, n):
        z = np.linspace(1.0, 15.0, 30)
        jax_vals = np.asarray(special.spherical_yn(n, jnp.array(z), derivative=True))
        cpp_vals = np.array([cpp_special.spherical_yn(n, float(zi), True) for zi in z])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12)


# ---------------------------------------------------------------------------
# Spherical Hankel h_n
# ---------------------------------------------------------------------------

class TestSphericalHn:
    @pytest.mark.parametrize("n", [0, 1, 2, 5, 10])
    def test_values(self, n):
        z = np.linspace(0.5, 20.0, 50)
        jax_vals = np.asarray(special.spherical_hn(n, jnp.array(z)))
        cpp_vals = np.array([cpp_special.spherical_hn(n, float(zi)) for zi in z])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12)

    @pytest.mark.parametrize("n", [1, 3, 7])
    def test_derivative(self, n):
        z = np.linspace(1.0, 15.0, 30)
        jax_vals = np.asarray(special.spherical_hn(n, jnp.array(z), derivative=True))
        cpp_vals = np.array([cpp_special.spherical_hn(n, float(zi), True) for zi in z])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12)


class TestSphericalHn2:
    @pytest.mark.parametrize("n", [0, 1, 5])
    def test_conjugate(self, n):
        z = np.linspace(0.5, 20.0, 50)
        jax_vals = np.asarray(special.spherical_hn_2(n, jnp.array(z)))
        cpp_vals = np.array([cpp_special.spherical_hn_2(n, float(zi)) for zi in z])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12)


class TestSphericalHnRecursion:
    @pytest.mark.parametrize("nmax", [1, 5, 10, 15])
    def test_all_orders(self, nmax):
        """Test recursion result against individual spherical_hn calls."""
        z_vals = [0.5, 2.0, 5.0, 15.0]
        for z in z_vals:
            jax_arr = np.asarray(special.spherical_hn_recursion(nmax, jnp.array(z)))
            ref = np.array([cpp_special.spherical_hn(n, z) for n in range(nmax + 1)])
            np.testing.assert_allclose(jax_arr, ref, rtol=1e-12,
                                       err_msg=f"nmax={nmax}, z={z}")


# ---------------------------------------------------------------------------
# Riccati-Bessel functions
# ---------------------------------------------------------------------------

class TestRiccati:
    @pytest.mark.parametrize("nmax", [1, 3, 6, 10])
    def test_riccati_1(self, nmax):
        """Compare riccati_1 against existing NumPy/scipy implementation.

        Forward recursion diverges from scipy for x < nmax; use x >> nmax.
        """
        from miepy.special_functions import riccati_1 as np_riccati_1
        x = np.array([float(nmax + 5), 15.0, 25.0])
        jax_vals = np.asarray(special.riccati_1(nmax, jnp.array(x)))
        np_vals = np_riccati_1(nmax, x)
        np.testing.assert_allclose(jax_vals, np_vals, rtol=1e-12)

    @pytest.mark.parametrize("nmax", [1, 3, 6, 10])
    def test_riccati_3(self, nmax):
        from miepy.special_functions import riccati_3 as np_riccati_3
        x = np.array([float(nmax + 5), 15.0, 25.0])
        jax_vals = np.asarray(special.riccati_3(nmax, jnp.array(x)))
        np_vals = np_riccati_3(nmax, x)
        np.testing.assert_allclose(jax_vals, np_vals, rtol=1e-12)

    @pytest.mark.parametrize("nmax", [1, 3, 6])
    def test_riccati_2(self, nmax):
        from miepy.special_functions import riccati_2 as np_riccati_2
        x = np.array([float(nmax + 5), 15.0, 25.0])
        jax_vals = np.asarray(special.riccati_2(nmax, jnp.array(x)))
        np_vals = np_riccati_2(nmax, x)
        np.testing.assert_allclose(jax_vals, np_vals, rtol=1e-12)


class TestRiccatiSingle:
    @pytest.mark.parametrize("n", [1, 3, 5, 10])
    def test_riccati_1_single(self, n):
        from miepy.special_functions import riccati_1_single as np_riccati_1s
        # Use x > n to keep forward recursion in stable regime
        x = np.array([float(n + 3), float(n + 10), 25.0, 40.0])
        jax_vals = np.asarray(special.riccati_1_single(n, jnp.array(x.astype(complex))))
        np_vals = np_riccati_1s(n, x)
        # JAX returns complex (input is complex); compare real parts for real input
        np.testing.assert_allclose(jax_vals.real, np_vals, rtol=1e-12)
        np.testing.assert_allclose(jax_vals.imag, 0, atol=1e-14)

    @pytest.mark.parametrize("n", [1, 3, 5])
    def test_riccati_3_single(self, n):
        from miepy.special_functions import riccati_3_single as np_riccati_3s
        x = np.array([float(n + 3), float(n + 10), 20.0])
        jax_vals = np.asarray(special.riccati_3_single(n, jnp.array(x)))
        np_vals = np_riccati_3s(n, x)
        np.testing.assert_allclose(jax_vals, np_vals, rtol=1e-12)


# ---------------------------------------------------------------------------
# Associated Legendre functions
# ---------------------------------------------------------------------------

class TestAssociatedLegendre:
    @pytest.mark.parametrize("n,m", [
        (0, 0), (1, 0), (1, 1), (1, -1),
        (2, 0), (2, 1), (2, 2), (2, -1), (2, -2),
        (3, 0), (3, 1), (3, 2), (3, 3), (3, -1), (3, -2), (3, -3),
        (5, 3), (5, -3), (7, 0), (7, 4),
    ])
    def test_values(self, n, m):
        z = np.array([-0.9, -0.5, 0.0, 0.3, 0.7, 0.99])
        jax_vals = np.asarray(special.associated_legendre(n, m, jnp.array(z)))
        cpp_vals = np.array([cpp_special.associated_legendre(n, m, float(zi)) for zi in z])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12, atol=1e-14)


class TestAssociatedLegendreRecursion:
    @pytest.mark.parametrize("nmax", [1, 3, 5, 8])
    def test_all_modes(self, nmax):
        """Test recursion output against individual associated_legendre calls."""
        z_vals = [-0.8, 0.0, 0.5, 0.95]
        for z in z_vals:
            jax_arr = np.asarray(special.associated_legendre_recursion(nmax, jnp.array(z)))
            # Build reference from individual calls
            total_size = nmax * (nmax + 2) + 1
            ref = np.zeros(total_size)
            for nn in range(0, nmax + 1):
                for m in range(-nn, nn + 1):
                    idx = nn * (nn + 2) - nn + m
                    ref[idx] = cpp_special.associated_legendre(nn, m, z)
            np.testing.assert_allclose(jax_arr, ref, rtol=1e-12, atol=1e-14,
                                       err_msg=f"nmax={nmax}, z={z}")


# ---------------------------------------------------------------------------
# Pi and Tau angular functions
# ---------------------------------------------------------------------------

class TestPiFunc:
    @pytest.mark.parametrize("n,m", [
        (1, 1), (1, -1), (1, 0),
        (2, 1), (2, -1), (2, 2), (2, 0),
        (3, 1), (3, -1), (3, 2), (3, 3),
        (5, 3), (5, -3),
    ])
    def test_general(self, n, m):
        theta = np.array([0.3, 0.7, 1.0, 1.5, 2.5])
        jax_vals = np.asarray(special.pi_func(n, m, jnp.array(theta)))
        cpp_vals = np.array([cpp_special.pi_func(n, m, float(t)) for t in theta])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-12, atol=1e-14)

    @pytest.mark.parametrize("n", [1, 2, 3, 5])
    def test_theta_zero(self, n):
        theta = jnp.array([0.0])
        for m in [-1, 0, 1]:
            jax_val = float(special.pi_func(n, m, theta)[0])
            cpp_val = cpp_special.pi_func(n, m, 0.0)
            np.testing.assert_allclose(jax_val, cpp_val, rtol=1e-12,
                                       err_msg=f"n={n}, m={m}, theta=0")

    @pytest.mark.parametrize("n", [1, 2, 3, 5])
    def test_theta_pi(self, n):
        theta = jnp.array([np.pi])
        for m in [-1, 0, 1]:
            jax_val = float(special.pi_func(n, m, theta)[0])
            cpp_val = cpp_special.pi_func(n, m, np.pi)
            np.testing.assert_allclose(jax_val, cpp_val, rtol=1e-12,
                                       err_msg=f"n={n}, m={m}, theta=pi")


class TestTauFunc:
    @pytest.mark.parametrize("n,m", [
        (1, 1), (1, -1), (1, 0),
        (2, 1), (2, -1), (2, 2), (2, 0),
        (3, 1), (3, -1), (3, 2), (3, 3),
        (5, 3), (5, -3),
    ])
    def test_general(self, n, m):
        theta = np.array([0.3, 0.7, 1.0, 1.5, 2.5])
        jax_vals = np.asarray(special.tau_func(n, m, jnp.array(theta)))
        cpp_vals = np.array([cpp_special.tau_func(n, m, float(t)) for t in theta])
        np.testing.assert_allclose(jax_vals, cpp_vals, rtol=1e-11, atol=1e-13)

    @pytest.mark.parametrize("n", [1, 2, 3, 5])
    def test_theta_zero(self, n):
        theta = jnp.array([0.0])
        for m in [-1, 0, 1]:
            jax_val = float(special.tau_func(n, m, theta)[0])
            cpp_val = cpp_special.tau_func(n, m, 0.0)
            np.testing.assert_allclose(jax_val, cpp_val, rtol=1e-12,
                                       err_msg=f"n={n}, m={m}, theta=0")

    @pytest.mark.parametrize("n", [1, 2, 3, 5])
    def test_theta_pi(self, n):
        theta = jnp.array([np.pi])
        for m in [-1, 0, 1]:
            jax_val = float(special.tau_func(n, m, theta)[0])
            cpp_val = cpp_special.tau_func(n, m, np.pi)
            np.testing.assert_allclose(jax_val, cpp_val, rtol=1e-12,
                                       err_msg=f"n={n}, m={m}, theta=pi")


# ---------------------------------------------------------------------------
# Emn normalization
# ---------------------------------------------------------------------------

class TestEmn:
    @pytest.mark.parametrize("n", [1, 2, 3, 5, 8])
    def test_values(self, n):
        for m in range(-n, n + 1):
            jax_val = special.Emn(m, n)
            cpp_val = cpp_vsh.Emn(m, n)
            np.testing.assert_allclose(jax_val, cpp_val, rtol=1e-14,
                                       err_msg=f"m={m}, n={n}")


# ---------------------------------------------------------------------------
# Wigner 3j symbols
# ---------------------------------------------------------------------------

class TestWigner3j:
    def test_known_values(self):
        """Test a few known Wigner 3j values against C++."""
        test_cases = [
            (1, 1, 0, 0, 0, 0),
            (1, 1, 2, 0, 0, 0),
            (2, 1, 3, 0, 0, 0),
            (2, 1, 3, 2, -1, -1),
            (1, 1, 2, 1, -1, 0),
            (3, 2, 1, 0, 0, 0),
            (3, 2, 3, 1, -2, 1),
            (5, 3, 4, 2, -1, -1),
        ]
        for j1, j2, j3, m1, m2, m3 in test_cases:
            jax_val = special.wigner_3j(j1, j2, j3, m1, m2, m3)
            cpp_val = cpp_special.wigner_3j(j1, j2, j3, m1, m2, m3)
            np.testing.assert_allclose(jax_val, cpp_val, rtol=1e-12,
                                       err_msg=f"({j1},{j2},{j3},{m1},{m2},{m3})")

    def test_selection_rules(self):
        """Wigner 3j should be zero when selection rules are violated."""
        assert special.wigner_3j(1, 1, 1, 0, 0, 1) == 0.0   # m1+m2+m3 != 0
        assert special.wigner_3j(1, 1, 5, 0, 0, 0) == 0.0   # triangle violated


class TestWigner3jBatch:
    @pytest.mark.parametrize("j2,j3,m1,m2,m3", [
        (1, 1, 0, 0, 0),
        (2, 1, 0, 0, 0),
        (3, 2, -1, 0, 1),
        (3, 2, 0, 0, 0),
        (5, 3, -2, 1, 1),
        (4, 4, 0, 0, 0),
    ])
    def test_batch_vs_individual(self, j2, j3, m1, m2, m3):
        """Test batch against individual wigner_3j calls via C++."""
        jax_jmin, jax_jmax, jax_vals = special.wigner_3j_batch(j2, j3, m1, m2, m3)

        if jax_jmax < jax_jmin:
            return  # empty result

        ref = np.array([
            cpp_special.wigner_3j(j1, j2, j3, m1, m2, m3)
            for j1 in range(jax_jmin, jax_jmax + 1)
        ])
        np.testing.assert_allclose(jax_vals, ref, rtol=1e-12,
                                   err_msg=f"j2={j2}, j3={j3}, m1={m1}, m2={m2}, m3={m3}")


# ---------------------------------------------------------------------------
# Gaunt coefficients
# ---------------------------------------------------------------------------

class TestGauntBatch:
    @pytest.mark.parametrize("m,n,u,v", [
        (0, 1, 0, 1),
        (1, 1, -1, 1),
        (1, 2, 0, 1),
        (1, 3, -1, 2),
        (0, 2, 0, 2),
        (2, 3, -1, 2),
        (1, 5, -1, 5),
        (0, 4, 0, 3),
        (-1, 2, 1, 3),
        (2, 4, -2, 4),
    ])
    def test_vs_individual(self, m, n, u, v):
        """Test batch Gaunt against individual a_func/b_func from C++."""
        jax_a, jax_b = special.gaunt_batch(m, n, u, v)

        qmax_A = min(n, v, (n + v - abs(m + u)) // 2)
        qmax_B = min(n, v, (n + v + 1 - abs(m + u)) // 2)

        assert len(jax_a) == qmax_A + 1
        assert len(jax_b) == qmax_B + 1

        # Compare a_vals
        ref_a = np.array([cpp_special.a_func(m, n, u, v, n + v - 2 * q) for q in range(qmax_A + 1)])
        if len(ref_a) > 0:
            np.testing.assert_allclose(jax_a, ref_a, rtol=1e-10, atol=1e-14,
                                       err_msg=f"a_vals: m={m},n={n},u={u},v={v}")

        # Compare b_vals (b_vals[0] = 0 always)
        ref_b = np.zeros(qmax_B + 1)
        for q in range(1, qmax_B + 1):
            p = n + v - 2 * q
            ref_b[q] = cpp_special.b_func(m, n, u, v, p)
        if len(ref_b) > 0:
            np.testing.assert_allclose(jax_b, ref_b, rtol=1e-10, atol=1e-14,
                                       err_msg=f"b_vals: m={m},n={n},u={u},v={v}")


# ---------------------------------------------------------------------------
# Mie scattering coefficients
# ---------------------------------------------------------------------------

class TestMieCoefficients:
    """Test JAX Mie coefficients against the NumPy reference."""

    def _make_params(self):
        """Common test parameters: gold sphere in vacuum."""
        radius = 75e-9
        wavelength = 600e-9
        eps = -11.04 + 0.78j   # approximate gold at 600 nm
        mu = 1.0
        eps_b = 1.0
        mu_b = 1.0
        k = 2 * np.pi / wavelength
        return radius, eps, mu, eps_b, mu_b, k

    @pytest.mark.parametrize("n", [1, 2, 3, 5])
    def test_scattering_coefficients(self, n):
        radius, eps, mu, eps_b, mu_b, k = self._make_params()
        jax_a, jax_b = jax_mie.mie_sphere_scattering_coefficients(
            radius, n, eps, mu, eps_b, mu_b, k)
        np_a, np_b = np_mie_scat(radius, n, eps, mu, eps_b, mu_b, k)
        # rtol=1e-8: for high n where x < n, forward recursion diverges from scipy
        np.testing.assert_allclose(complex(jax_a), complex(np_a), rtol=1e-8)
        np.testing.assert_allclose(complex(jax_b), complex(np_b), rtol=1e-8)

    @pytest.mark.parametrize("n", [1, 2, 3, 5])
    def test_interior_coefficients(self, n):
        radius, eps, mu, eps_b, mu_b, k = self._make_params()
        jax_c, jax_d = jax_mie.mie_sphere_interior_coefficients(
            radius, n, eps, mu, eps_b, mu_b, k)
        np_c, np_d = np_mie_int(radius, n, eps, mu, eps_b, mu_b, k)
        np.testing.assert_allclose(complex(jax_c), complex(np_c), rtol=1e-10)
        np.testing.assert_allclose(complex(jax_d), complex(np_d), rtol=1e-10)

    def test_conducting_sphere(self):
        radius, eps, mu, eps_b, mu_b, k = self._make_params()
        for n in [1, 2, 3]:
            jax_a, jax_b = jax_mie.mie_sphere_scattering_coefficients(
                radius, n, eps, mu, eps_b, mu_b, k, conducting=True)
            np_a, np_b = np_mie_scat(
                radius, n, eps, mu, eps_b, mu_b, k, conducting=True)
            np.testing.assert_allclose(complex(jax_a), complex(np_a), rtol=1e-10)
            np.testing.assert_allclose(complex(jax_b), complex(np_b), rtol=1e-10)

    def test_array_wavelengths(self):
        """Test with array-valued wavenumber."""
        radius = 75e-9
        wavelengths = np.array([400e-9, 500e-9, 600e-9, 700e-9])
        eps = np.array([-1.7 + 5.3j, -4.6 + 2.2j, -11.0 + 0.78j, -18.4 + 0.63j])
        mu = 1.0
        eps_b = 1.0
        mu_b = 1.0
        k = 2 * np.pi / wavelengths

        for n in [1, 2, 3]:
            jax_a, jax_b = jax_mie.mie_sphere_scattering_coefficients(
                radius, n, eps, mu, eps_b, mu_b, k)
            np_a, np_b = np_mie_scat(radius, n, eps, mu, eps_b, mu_b, k)
            np.testing.assert_allclose(np.asarray(jax_a), np_a, rtol=1e-10)
            np.testing.assert_allclose(np.asarray(jax_b), np_b, rtol=1e-10)
