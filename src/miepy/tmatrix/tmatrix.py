import numpy as np
import miepy

def get_tmatrix(rad, dS, eps, eps_m, wavelength, lmax):
    _, Ntheta, Nphi = dS.shape
    theta = np.linspace(0, np.pi, Ntheta)
    phi = np.linspace(0, 2*np.pi, Nphi)
    THETA, PHI = np.meshgrid(theta, phi, indexing='ij')
    dS = miepy.coordinates.vec_cart_to_sph(dS, THETA, PHI)

    rmax = lmax*(lmax+2)
    k1 = 2*np.pi/wavelength
    k2 = 2*np.pi*np.sqrt(eps)/wavelength

    class ingoing:
        M_k2  = np.empty((rmax,3) + THETA.shape, dtype=complex)
        N_k2  = np.empty((rmax,3) + THETA.shape, dtype=complex)
        M_k1  = np.empty((rmax,3) + THETA.shape, dtype=complex)
        N_k1  = np.empty((rmax,3) + THETA.shape, dtype=complex)    

    class outgoing:
        M_k2 = np.empty((rmax,3) + THETA.shape, dtype=complex)
        N_k2 = np.empty((rmax,3) + THETA.shape, dtype=complex)
        M_k1 = np.empty((rmax,3) + THETA.shape, dtype=complex)
        N_k1 = np.empty((rmax,3) + THETA.shape, dtype=complex)    

    for r,n,m in miepy.mode_indices(lmax):
        Emn = miepy.vsh.Emn(m, n)

        Nfunc, Mfunc = miepy.vsh.VSH(n, m, miepy.vsh_mode.incident)
        ingoing.M_k2[r] = Emn*Mfunc(rad, THETA, PHI, k2)
        ingoing.N_k2[r] = Emn*Nfunc(rad, THETA, PHI, k2)
        ingoing.M_k1[r] = Emn*Mfunc(rad, THETA, PHI, k1)
        ingoing.N_k1[r] = Emn*Nfunc(rad, THETA, PHI, k1)

        Nfunc, Mfunc = miepy.vsh.VSH(n, m, miepy.vsh_mode.outgoing)
        outgoing.M_k2[r] = Emn*Mfunc(rad, THETA, PHI, k2)
        outgoing.N_k2[r] = Emn*Nfunc(rad, THETA, PHI, k2)
        outgoing.M_k1[r] = Emn*Mfunc(rad, THETA, PHI, k1)
        outgoing.N_k1[r] = Emn*Nfunc(rad, THETA, PHI, k1)

    def get_Q(p, q):
        Q = np.zeros([2,rmax,2,rmax], dtype=complex)
        p_modes = ingoing if p == 1 else outgoing
        q_modes = ingoing if q == 1 else outgoing

        for r1,n1,m1 in miepy.mode_indices(lmax):
            M1_p_k2 = p_modes.M_k2[r1]
            N1_p_k2 = p_modes.N_k2[r1]
            M1_p_k1 = p_modes.M_k1[r1]
            N1_p_k1 = p_modes.N_k1[r1]

            for r2,n2,m2 in miepy.mode_indices(lmax):
                M2_q_k2 = q_modes.M_k2[r2]
                N2_q_k2 = q_modes.N_k2[r2]
                M2_q_k1 = q_modes.M_k1[r2]
                N2_q_k1 = q_modes.N_k1[r2]

                # factor = (-1)**(m2 - m1)
                factor = 1

                integrand = np.sum(np.cross(dS, M2_q_k2, axis=0)*N1_p_k1, axis=0) \
                            + np.sqrt(eps/eps_m)*np.sum(np.cross(dS, N2_q_k2, axis=0)*M1_p_k1, axis=0)
                integrand *= 1j*k1**2/np.pi*factor
                Q[1,r1,1,r2] = miepy.vsh.misc.trapz_2d(theta, phi, integrand)

                integrand = np.sum(np.cross(dS, N2_q_k2, axis=0)*N1_p_k1, axis=0) \
                            + np.sqrt(eps/eps_m)*np.sum(np.cross(dS, M2_q_k2, axis=0)*M1_p_k1, axis=0)
                integrand *= 1j*k1**2/np.pi*factor
                Q[1,r1,0,r2] = miepy.vsh.misc.trapz_2d(theta, phi, integrand)

                integrand = np.sum(np.cross(dS, M2_q_k2, axis=0)*M1_p_k1, axis=0) \
                            + np.sqrt(eps/eps_m)*np.sum(np.cross(dS, N2_q_k2, axis=0)*N1_p_k1, axis=0)
                integrand *= 1j*k1**2/np.pi*factor
                Q[0,r1,1,r2] = miepy.vsh.misc.trapz_2d(theta, phi, integrand)

                integrand = np.sum(np.cross(dS, N2_q_k2, axis=0)*M1_p_k1, axis=0) \
                            + np.sqrt(eps/eps_m)*np.sum(np.cross(dS, M2_q_k2, axis=0)*N1_p_k1, axis=0)
                integrand *= 1j*k1**2/np.pi*factor
                Q[0,r1,0,r2] = miepy.vsh.misc.trapz_2d(theta, phi, integrand)


        return Q

    Q31 = -get_Q(3, 1)
    Q11 = get_Q(1, 1)

    T = -np.einsum('aibj,bjck->aick', Q11, np.linalg.tensorinv(Q31))

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

theta = np.linspace(0, np.pi, 80)
phi = np.linspace(0, 2*np.pi, 80)
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

rad, dS = ellipsoid_dS(1*radius, 1.3*radius, 1.6*radius, THETA, PHI)
T = get_tmatrix(rad, dS, eps, 1, wavelength, lmax)
print(T[1,:,0,0])
T1 = np.copy(T)

spheroid = miepy.spheroid([0,0,0], radius, 1.3*radius, material)
spheroid = miepy.ellipsoid([0,0,0], radius, 1.3*radius, 1.6*radius, material)
T = spheroid.compute_tmatrix(lmax, wavelength, 1)
print(T[1,:,0,0])
T2 = np.copy(T)

import matplotlib.pyplot as plt
rmax = lmax*(lmax+2)
T1 = np.reshape(T1, [2*rmax, 2*rmax])
T2 = np.reshape(T2, [2*rmax, 2*rmax])

y = np.abs(T1-T2)
y /= np.max(y)
plt.pcolormesh(y, vmax=.1)
plt.gca().set_aspect('equal')
plt.colorbar()

plt.show()
