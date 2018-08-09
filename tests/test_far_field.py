"""
tests for far-fields
"""
import numpy as np
import miepy

nm = 1e-9

dimer = miepy.sphere_cluster(position=[[-100*nm,0,0], [100*nm, 0, 0]],
                      radius=75*nm,
                      material=miepy.constant_material(3.6**2),
                      source=miepy.sources.plane_wave.from_string(polarization='y'),
                      wavelength=600*nm,
                      lmax=2)

theta = np.linspace(0, np.pi, 5)
phi = np.linspace(0, 2*np.pi, 5)
THETA, PHI = np.meshgrid(theta, phi)
R = 1e6
X, Y, Z = miepy.coordinates.sph_to_cart(R, THETA, PHI)
E_exact = dimer.E_field(R, THETA, PHI, source=False, interior=False, spherical=True)

def test_far_field_convergence():
    """far-field E and H field should agree with exact field in the large radius limit"""
    E_far = dimer.E_field(R, THETA, PHI, far=True, source=False, interior=False, spherical=True)
    np.testing.assert_allclose(E_exact, E_far, rtol=0, atol=1e-16)

def test_far_field_cluster_coefficient():
    """far-fields calculated from the cluster coefficients should be the same as the sum-over particle coefficients"""
    dimer.solve_cluster_coefficients(lmax=4)
    E_func = miepy.vsh.expand_E_far(dimer.p_cluster, dimer.material_data.k_b)
    E_far = E_func(R, THETA, PHI)

    np.testing.assert_allclose(E_exact, E_far, rtol=0, atol=1e-15)

def test_far_field_directly():
    """far-field function compared directly to total field function for n=2, m=-1"""
    x = 1e6

    n = 2
    m = -1
    Nfunc, Mfunc = miepy.vsh.VSH(n, m)

    rad, theta, phi = 1e6, 0.9, -0.6
    k = 1
    N = Nfunc(rad, theta, phi, k)
    M = Mfunc(rad, theta, phi, k)
    tau = miepy.vsh.special.tau_func(n, m, theta)
    pi = miepy.vsh.special.pi_func(n, m, theta)

    E_theta_1 = 1j*N[1]
    factor = np.exp(1j*k*rad)/(k*rad)
    E_theta_2 = 1j*factor*(-1j)**n*tau*np.exp(1j*m*phi)
    np.testing.assert_allclose(E_theta_1, E_theta_2, rtol=0, atol=1e-12)

    E_theta_1 = 1j*M[1]
    factor = np.exp(1j*k*rad)/(k*rad)
    E_theta_2 = 1j*factor*(-1j)**n*pi*np.exp(1j*m*phi)
    np.testing.assert_allclose(E_theta_1, E_theta_2, rtol=0, atol=1e-12)
