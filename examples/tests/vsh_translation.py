import numpy as np

import miepy

nm = 1e-9

v = 10
u = 0

n = 10
m = 0

ftype = "electric"

N, M = miepy.vsh.VSH(n, m)
if ftype == "magnetic":
    func = M
elif ftype == "electric":
    func = N

k = 2 * np.pi / (600 * nm)
r = 600 * nm

origin_1 = np.array([0, 0, 0])
THETA, PHI = miepy.coordinates.sphere_mesh(800)
E = func(r, THETA, PHI, k)
p = miepy.vsh.integral_project_fields_onto(E, r, k, ftype, n, m, spherical=True)

origin_2 = np.array([100 * nm, 0, 0])
origin_2 = miepy.coordinates.sph_to_cart(2 / k, 0.5, 0.5)
x, y, z = miepy.coordinates.sph_to_cart(r, THETA, PHI, origin=origin_2)
R_p, THETA_p, PHI_p = miepy.coordinates.cart_to_sph(x, y, z)
E_prime = func(R_p, THETA_p, PHI_p, k)
E_prime = miepy.coordinates.vec_sph_to_cart(E_prime, THETA_p, PHI_p)
E_prime = miepy.coordinates.vec_cart_to_sph(E_prime, THETA, PHI)
p_prime = miepy.vsh.integral_project_fields_onto(
    E_prime, r, k, ftype, v, u, spherical=True, mode=miepy.vsh.vsh_mode.incident
)

A = p_prime / p
print(A)

r_ij, theta_ij, phi_ij = miepy.coordinates.cart_to_sph(*(origin_2 - origin_1))
A = miepy.vsh.A_translation(m, n, u, v, r_ij, theta_ij, phi_ij, k)
# A = miepy.vsh.A_translation(-2, 6, -2, 10, r_ij, theta_ij, phi_ij, k)
# A = miepy.vsh.A_translation(0, 10, 0, 10, r_ij, theta_ij, phi_ij, k)
print(A)
