"""
Expressions for the force and torque given the expansion coefficients or the fields
"""

import numpy as np
from scipy import constants
import miepy
from my_pytools.my_numpy.integrate import simps_2d

def levi_civita():
    """return the levi-civita symbol"""

    eijk = np.zeros((3, 3, 3))
    eijk[0, 1, 2] = eijk[1, 2, 0] = eijk[2, 0, 1] = 1
    eijk[0, 2, 1] = eijk[2, 1, 0] = eijk[1, 0, 2] = -1
    return eijk

def force(p, q, p_inc, q_inc, k, E0, eps_b, mu_b, Lmax):
    """force from the expansion coefficients
    
       Arguments:
           p        p coefficients (scattered electric)
           q        q coefficients (scattered magnetic)
           p_inc    p_inc coefficients (incident electric)
           q_inc    q_inc coefficients (incident magnetic)
           k        wavenumber
           E0       field amplitude
           eps_b    background relative permitvitty
           mu_b     background relative permeability
           Lmax     maximum number of terms
    """
    Fxy = 0
    Fz = 0
    Axy = np.pi*E0**2/k**2*constants.epsilon_0*eps_b
    Az = -2*np.pi*E0**2/k**2*constants.epsilon_0*eps_b

    for n in range(1,Lmax+1):
        for m in range(-n,n+1):
            r = n**2 + n - 1 + m

            # Fxy, term 1/3
            if m != n:
                factor = Axy*np.sqrt((n+m+1)*(n-m))/(n*(n+1))
                r1 = n**2 + n - 1 + m + 1
                Fxy += factor*(2*p[r]*np.conj(q[r1]) \
                         - p[r]*np.conj(q_inc[r1]) \
                         - p_inc[r]*np.conj(q[r1]) \
                         + 2*q[r]*np.conj(p[r1]) \
                         - q[r]*np.conj(p_inc[r1]) \
                         - q_inc[r]*np.conj(p[r1]))

            # Fz, term 1/2
            factor = Az*m/(n*(n+1))
            Fz += factor*(2*p[r]*np.conj(q[r]) \
                    - p[r]*np.conj(q_inc[r]) \
                    - p_inc[r]*np.conj(q[r]))


            if n < Lmax:
                # Fxy, term 2/3
                factor = -Axy*np.sqrt((n+m+2)*(n+m+1)*n*(n+2)/((2*n+3)*(2*n+1)))/(n+1)
                r1 = (n+1)**2 + (n+1) - 1 + m + 1
                Fxy += factor*(2*p[r]*np.conj(p[r1]) \
                         - p[r]*np.conj(p_inc[r1]) \
                         - p_inc[r]*np.conj(p[r1]) \
                         + 2*q[r]*np.conj(q[r1]) \
                         - q[r]*np.conj(q_inc[r1]) \
                         - q_inc[r]*np.conj(q[r1]))

                # Fxy, term 3/3
                factor = Axy*np.sqrt((n-m+1)*(n-m+2)*n*(n+2)/((2*n+3)*(2*n+1)))/(n+1)
                r1 = (n+1)**2 + (n+1) - 1 + m - 1
                Fxy += factor*(2*p[r1]*np.conj(p[r]) \
                         - p[r1]*np.conj(p_inc[r]) \
                         - p_inc[r1]*np.conj(p[r]) \
                         + 2*q[r1]*np.conj(q[r]) \
                         - q[r1]*np.conj(q_inc[r]) \
                         - q_inc[r1]*np.conj(q[r]))

                # Fz, term 2/2
                factor = Az/(n+1)*np.sqrt((n-m+1)*(n+m+1)*n*(n+2)/(2*n+3)/(2*n+1))
                r1 = (n+1)**2 + (n+1) - 1 + m
                Fz += factor*(2*p[r1]*np.conj(p[r]) \
                        - p[r1]*np.conj(p_inc[r]) \
                        - p_inc[r1]*np.conj(p[r]) \
                        + 2*q[r1]*np.conj(q[r]) \
                        - q[r1]*np.conj(q_inc[r]) \
                        - q_inc[r1]*np.conj(q[r]))

    return np.array([np.real(Fxy), np.imag(Fxy), np.real(Fz)])

def torque(p, q, p_inc, q_inc, k, E0, eps_b, mu_b, Lmax):
    """torque from the expansion coefficients
    
       Arguments:
           p        p coefficients (scattered electric)
           q        q coefficients (scattered magnetic)
           p_inc    p_inc coefficients (incident electric)
           q_inc    q_inc coefficients (incident magnetic)
           k        wavenumber
           E0       field amplitude
           eps_b    background relative permitvitty
           mu_b     background relative permeability
           Lmax     maximum number of terms
    """
    T = np.zeros(3, dtype=float)
    A = -2*np.pi*E0**2/k**3*constants.epsilon_0*eps_b

    for n in range(1,Lmax+1):
        for m in range(-n,n+1):
            r = n**2 + n - 1 + m

            if m != n:
                # Tx
                factor = -A*np.sqrt((n-m)*(n+m+1))
                r1 = n**2 + n - 1 + m + 1
                T[0] += factor*np.real(p[r]*np.conj(p[r1]) \
                        + q[r]*np.conj(q[r1]) \
                        -0.5*(p[r1]*np.conj(p_inc[r]) \
                        + p[r]*np.conj(p_inc[r1]) \
                        + q[r1]*np.conj(q_inc[r]) \
                        + q[r]*np.conj(q_inc[r1])))

                # Ty
                T[1] += factor*np.imag(p[r]*np.conj(p[r1]) \
                        + q[r]*np.conj(q[r1]) \
                        +0.5*(p[r1]*np.conj(p_inc[r]) \
                        - p[r]*np.conj(p_inc[r1]) \
                        + q[r1]*np.conj(q_inc[r]) \
                        - q[r]*np.conj(q_inc[r1])))

            # Tz
            factor = A*m
            T[2] += factor* (np.abs(p[r])**2 + np.abs(q[r])**2 \
                    - np.real(p[r]*np.conj(p_inc[r]) \
                    + q[r]*np.conj(q_inc[r])))

    return T

def maxwell_stress_tensor(E, H, eps=1, mu=1):
    """Compute the Maxwell stress tensor
    
       Arguments:
           E[3,...]   electric field data
           H[3,...]   magnetic field data
           eps        medium permitvitty (default: 1)
           mu         medium permeability (default: 1)

       Returns T[3,3,...]
    """
    sigma = eps*np.einsum('i...,j...->ij...', E, np.conj(E)) \
            + mu*np.einsum('i...,j...->ij...', H, np.conj(H)) \
            - 0.5*np.einsum('ij,...->ij...', np.identity(3), eps*np.sum(np.abs(E)**2, axis=0)) \
            - 0.5*np.einsum('ij,...->ij...', np.identity(3), mu*np.sum(np.abs(H)**2, axis=0))
    sigma *= constants.epsilon_0/2

    return sigma

def force_and_torque_from_mst(E, H, radius, eps=1, mu=1):
    """Compute the force and torque on a particle using the Maxwell stress tensor

       Arguments:
           E[3,Ntheta,Nphi]     electric field values on the surface of a sphere
           H[3,Ntheta,Nphi]     magnetic field values on the surface of a sphere
           radius               radius of sphere 
           eps                  medium permitvitty (default: 1)
           mu                   medium permeability (default: 1)

       Returns (F[3], T[3])
    """
    sigma = maxwell_stress_tensor(E, H, eps, mu)
    levi = levi_civita()
    dA = radius**2

    Ntheta, Nphi = E.shape[1:]
    THETA, PHI = miepy.vsh.sphere_mesh(Ntheta)
    rhat,*_ = miepy.vsh.sph_basis_vectors(THETA, PHI)

    tau = np.linspace(-1, 1, Ntheta)
    phi = np.linspace(0, 2*np.pi, Nphi)

    integrand = np.einsum('ijxy,jxy->ixy', sigma, rhat)*dA
    F = np.array([simps_2d(tau, phi, integrand[x].real) for x in range(3)])

    integrand = np.einsum('imn,mxy,njxy,jxy->ixy', levi, radius*rhat, sigma, rhat)*dA
    T = np.array([simps_2d(tau, phi, integrand[x].real) for x in range(3)])

    return F,T

def _gmt_force_and_torque_from_mst(gmt, i, sampling=30):
    """FOR TESTING ONLY!
    Given GMT object and particle number i, return F,T using MST
    """
    radius = gmt.spheres.radius[i]
    X,Y,Z,THETA,PHI,tau,phi = miepy.vsh.cart_sphere_mesh(radius, gmt.spheres.position[i], sampling)

    E = gmt.E_field_from_particle(i, X, Y, Z)
    H = gmt.H_field_from_particle(i, X, Y, Z)

    F = np.zeros([3, gmt.Nfreq])
    T = np.zeros([3, gmt.Nfreq])

    for k in range(gmt.Nfreq):
        eps_b = gmt.material_data['eps_b'][k]
        mu_b = gmt.material_data['mu_b'][k]
        F[:,k], T[:,k] = force_and_torque_from_mst(E[:,k], H[:,k], radius, eps_b, mu_b)

    return F,T
