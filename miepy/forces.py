"""
Expressions for the force and torque given the expansion coefficients or the fields
"""

import numpy as np
from scipy import constants

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
    Axy = np.pi*E0**2/k**2
    Az = -2*np.pi*E0**2/k**2

    for n in range(1,Lmax+1):
        for m in range(-n,n+1):
            r = n**2 + n - 1 + m

            # Fxy, term 1/3
            if m != n:
                factor = Axy*np.sqrt((n+m+1)*(n-m))*(eps_b/mu_b)/(n+1)
                r1 = n**2 + n - 1 + m + 1
                Fxy += factor*(2*p[r]*np.conj(q[r1]) \
                         - p[r]*np.conj(q_inc[r1]) \
                         - p_inc[r]*np.conj(q[r1]) \
                         + 2*q[r]*np.conj(p[r1]) \
                         - q[r]*np.conj(p_inc[r1]) \
                         - q_inc[r]*np.conj(p[r1]))

            # Fz, term 1/2
            factor = Az*m*(eps_b/mu_b)/(n*(n+1))
            Fz += factor*(2*p[r]*np.conj(q[r]) \
                    - p[r]*np.conj(q_inc[r]) \
                    - p_inc[r]*np.conj(q[r]))


            if n < Lmax:
                # Fxy, term 2/3
                factor = -Axy*np.sqrt((n+m+2)*(n+m+1)*n*(n+2)/((2*n+3)*(2*n+1)))/(n+1)
                r1 = (n+1)**2 + (n+1) - 1 + m + 1
                Fxy += factor*(2*eps_b*p[r]*np.conj(p[r1]) \
                         - eps_b*p[r]*np.conj(p_inc[r1]) \
                         - eps_b*p_inc[r]*np.conj(p[r1]) \
                         + 2*eps_b/mu_b*q[r]*np.conj(q[r1]) \
                         - eps_b/mu_b*q[r]*np.conj(q_inc[r1]) \
                         - eps_b/mu_b*q_inc[r]*np.conj(q[r1]))

                # Fxy, term 3/3
                factor = Axy*np.sqrt((n-m+1)*(n-m+2)*n*(n+2)/((2*n+3)*(2*n+1)))/(n+1)
                r1 = (n+1)**2 + (n+1) - 1 + m - 1
                Fxy += factor*(2*eps_b*p[r1]*np.conj(p[r]) \
                         - eps_b*p[r1]*np.conj(p_inc[r]) \
                         - eps_b*p_inc[r1]*np.conj(p[r]) \
                         + 2*eps_b/mu_b*q[r1]*np.conj(q[r]) \
                         - eps_b/mu_b*q[r1]*np.conj(q_inc[r]) \
                         - eps_b/mu_b*q_inc[r1]*np.conj(q[r]))

                # Fz, term 2/2
                factor = Az/(n+1)*np.sqrt((n-m+1)*(n+m+1)*n*(n+2)/(2*n+3)/(2*n+1))
                r1 = (n+1)**2 + (n+1) - 1 + m
                Fz += factor*(2*eps_b*p[r1]*np.conj(p[r]) \
                        - eps_b*p[r1]*np.conj(p_inc[r]) \
                        - eps_b*p_inc[r1]*np.conj(p[r]) \
                        + eps_b/mu_b*2*q[r1]*np.conj(q[r]) \
                        - eps_b/mu_b*q[r1]*np.conj(q_inc[r]) \
                        - eps_b/mu_b*q_inc[r1]*np.conj(q[r]))

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
    A = -2*np.pi*E0**2/k**3*constants.epsilon_0

    for n in range(1,Lmax+1):
        for m in range(-n,n+1):
            r = n**2 + n - 1 + m

            if m != n:
                # Tx
                factor = -A*np.sqrt((n-m)*(n+m+1))
                r1 = n**2 + n - 1 + m + 1
                T[0] += factor*np.real(eps_b*p[r]*np.conj(p[r1]) \
                        + mu_b*q[r]*np.conj(q[r1]) \
                        -0.5*(eps_b*p[r1]*np.conj(p_inc[r]) \
                        + eps_b*p[r]*np.conj(p_inc[r1]) \
                        + mu_b*q[r1]*np.conj(q_inc[r]) \
                        + mu_b*q[r]*np.conj(q_inc[r1])))

                # Ty
                T[1] += factor*np.imag(eps_b*p[r]*np.conj(p[r1]) \
                        + mu_b*q[r]*np.conj(q[r1]) \
                        +0.5*(eps_b*p[r1]*np.conj(p_inc[r]) \
                        - eps_b*p[r]*np.conj(p_inc[r1]) \
                        + mu_b*q[r1]*np.conj(q_inc[r]) \
                        - mu_b*q[r]*np.conj(q_inc[r1])))

            # Tz
            factor = A*m/n
            T[2] += factor* (eps_b*np.abs(p[r])**2 + mu_b*np.abs(q[r])**2 \
                    - np.real(eps_b*p[r]*np.conj(p_inc[r]) \
                    + mu_b*q[r]*np.conj(q_inc[r])))

    return T
