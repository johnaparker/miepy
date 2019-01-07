import numpy as np
import miepy
from math import factorial
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def get_tmatrix(rad, dS, eps, eps_m, wavelength, lmax):
    _, Ntheta, Nphi = dS.shape
    theta = np.linspace(0, np.pi, Ntheta)
    phi = np.linspace(0, 2*np.pi, Nphi)
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    dS = miepy.coordinates.vec_cart_to_sph(dS, THETA, PHI)

    rmax = lmax*(lmax+2)
    k1 = 2*np.pi/wavelength
    k2 = 2*np.pi*np.sqrt(eps)/wavelength

    np.conj1 = np.conj
    np.conj1 = lambda x: x
    np.conj2 = np.conj
    np.conj2 = lambda x: x
    def get_Q(p, q):
        modes = {1: miepy.vsh_mode.incident, 3: miepy.vsh_mode.outgoing}
        Q = np.zeros([2,rmax,2,rmax], dtype=complex)
        f = factorial

        for r1,n1,m1 in miepy.mode_indices(lmax):
            Nfunc, Mfunc = miepy.vsh.VSH(n1, m1, modes[p])
            Emn = miepy.vsh.Emn(m1, n1)
            M1_p_k2 = Emn*Mfunc(rad, THETA, PHI, k2)
            N1_p_k2 = Emn*Nfunc(rad, THETA, PHI, k2)
            M1_p_k1 = Emn*Mfunc(rad, THETA, PHI, k1)
            N1_p_k1 = Emn*Nfunc(rad, THETA, PHI, k1)

            for r2,n2,m2 in miepy.mode_indices(lmax):
                factor = (-1)**(m2-m1)
                factor = 1
                k3 = k1 if p == 1 else k1

                Nfunc, Mfunc = miepy.vsh.VSH(n2, m2, modes[q])
                Emn = miepy.vsh.Emn(m2, n2)
                M2_q_k2 = Emn*Mfunc(rad, THETA, PHI, k2)
                N2_q_k2 = Emn*Nfunc(rad, THETA, PHI, k2)
                M2_q_k1 = Emn*Mfunc(rad, THETA, PHI, k1)
                N2_q_k1 = Emn*Nfunc(rad, THETA, PHI, k1)

                integrand = np.sum(np.cross(dS, np.conj2(M2_q_k2), axis=0)*np.conj1(N1_p_k1), axis=0) \
                            + np.sqrt(eps/eps_m)*np.sum(np.cross(dS, np.conj2(N2_q_k2), axis=0)*np.conj1(M1_p_k1), axis=0)
                integrand *= 1j*k3**2/np.pi*factor
                Q[1,r1,1,r2] = miepy.vsh.misc.trapz_2d(theta, phi, integrand)

                integrand = np.sum(np.cross(dS, np.conj2(N2_q_k2), axis=0)*np.conj1(N1_p_k1), axis=0) \
                            + np.sqrt(eps/eps_m)*np.sum(np.cross(dS, np.conj2(M2_q_k2), axis=0)*np.conj1(M1_p_k1), axis=0)
                integrand *= 1j*k3**2/np.pi*factor
                Q[1,r1,0,r2] = miepy.vsh.misc.trapz_2d(theta, phi, integrand)

                integrand = np.sum(np.cross(dS, np.conj2(M2_q_k2), axis=0)*np.conj1(M1_p_k1), axis=0) \
                            + np.sqrt(eps/eps_m)*np.sum(np.cross(dS, np.conj2(N2_q_k2), axis=0)*np.conj1(N1_p_k1), axis=0)
                integrand *= 1j*k3**2/np.pi*factor
                Q[0,r1,1,r2] = miepy.vsh.misc.trapz_2d(theta, phi, integrand)

                integrand = np.sum(np.cross(dS, np.conj2(N2_q_k2), axis=0)*np.conj1(M1_p_k1), axis=0) \
                            + np.sqrt(eps/eps_m)*np.sum(np.cross(dS, np.conj2(M2_q_k2), axis=0)*np.conj1(N1_p_k1), axis=0)
                integrand *= 1j*k3**2/np.pi*factor
                Q[0,r1,0,r2] = miepy.vsh.misc.trapz_2d(theta, phi, integrand)

        return Q

    Q31 = get_Q(3, 1)
    Q11 = get_Q(1, 1)

    # T = -Q11*np.linalg.tensorinv(Q31)
    # Q31 = np.transpose(Q31, (2,3,0,1))
    T = -np.einsum('aibj,bjck->aick', Q11, np.linalg.tensorinv(Q31))
    # T = -Q11.reshape([2*rmax,2*rmax]) @ np.linalg.inv(Q31.reshape([2*rmax, 2*rmax]))
    # T = T.reshape([2,rmax,2,rmax])

    return T

nm = 1e-9

lmax = 4

wavelength = 600*nm
radius = 60*nm
eps = 4

material = miepy.constant_material(eps)
sphere = miepy.sphere([0,0,0], radius, material)
T = sphere.compute_tmatrix(lmax, wavelength, 1)
# print(T[0,:,0,0])

theta = np.linspace(0, np.pi, 50)
phi = np.linspace(0, 2*np.pi, 50)
THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
rhat, that, phat = miepy.coordinates.sph_basis_vectors(THETA, PHI)
dS = rhat*np.sin(THETA)*radius**2

# T = get_tmatrix(radius, dS, eps, 1, wavelength, lmax)
# print(T[0,:,0,0])

def ellipsoid_dS(a, b, c, theta, phi):
    rad = 1/np.sqrt(np.sin(theta)**2*np.cos(phi)**2/a**2 + np.sin(theta)**2*np.sin(phi)**2/b**2 + np.cos(theta)**2/c**2)
    dr_theta = -rad**3*(np.sin(theta)*np.cos(theta)*(np.cos(phi)**2/a**2 + np.sin(phi)**2/b**2 - 1/c**2))
    dr_phi = -rad**3*np.sin(phi)*np.cos(phi)*np.sin(theta)*(1/b**2 - 1/a**2)

    rhat, that, phat = miepy.coordinates.sph_basis_vectors(theta, phi)
    sigma = (rhat - 1/rad*dr_theta*that - 1/rad*dr_phi*phat)*rad**2*np.sin(theta)

    return rad, sigma

rad, dS = ellipsoid_dS(radius, radius, 1.3*radius, THETA, PHI)
T = get_tmatrix(rad, dS, eps, 1, wavelength, lmax)
print(T[0,:,0,0])

spheroid = miepy.spheroid([0,0,0], radius, 1.3*radius, material)
T = spheroid.compute_tmatrix(lmax, wavelength, 1)
print(T[0,:,0,0])


plt.show()
