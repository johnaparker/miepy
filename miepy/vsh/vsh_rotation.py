"""
VSH rotation functions
"""

import numpy as np
from math import factorial
import miepy
from miepy import vsh
from scipy import special
import sympy

n = 4
m = 3
mp = 1
beta = 0.4

def wigner_d(n, m, mp, beta):
    """recursion formulation (Mishchenko, Scattering, Absorption and Emission of Light by small Particles, p.365 (B.22 - B.24))"""
    wig_d = np.zeros(n + 1, dtype=complex)
    
    if beta < 0:
        aa = m
        bb = mp
        m = bb
        mp = aa
       
    if m == 0 and mp == 0:
        for nn in range(1, n + 1):
            wig_d[nn] = sympy.legendre(nn, np.cos(beta))          
    else:
        l_min = max(abs(m), abs(mp))
        wig_d[l_min - 1] = 0
        if mp >= m:
            zeta = 1
        else:
            zeta = (-1) ** (m - mp)
    
        wig_d[l_min] = (zeta * 2.0 ** (-l_min) * (factorial(2 * l_min) / (factorial(abs(m - mp)) 
                                                                          * factorial(abs(m + mp)))) ** 0.5
                        * (1 - np.cos(beta)) ** (abs(m - mp) / 2) 
                        * (1 + np.cos(beta)) ** (abs(m + mp) / 2 ))

        for ll in range(l_min, n):
            wig_d[ll + 1] = (((2 * ll + 1) * (ll * (ll + 1) * np.cos(beta) - m * mp) * wig_d[ll] 
                            - (ll + 1) * (ll ** 2 - m ** 2) ** 0.5 * (ll ** 2 - mp ** 2) ** 0.5 
                            * wig_d[ll - 1]) / (ll * ((ll + 1) ** 2 - m ** 2) ** 0.5 
                                                * ((ll + 1) ** 2 - mp ** 2) ** 0.5))

    return wig_d[n].real

    # f = factorial
    # factor = (f(n-m)*f(n+m)/f(n+mp)/f(n-mp))**0.5
    # factor = 1
    # d = factor*np.sin(beta/2)**a*np.cos(beta/2)**b \
        # * special.jacobi(n-m, a, b)(np.cos(beta))

def wigner_D(n, m, mp, alpha, beta, gamma):
    d = np.exp(1j*m*alpha)*wigner_d(n, m, mp, beta)*np.exp(1j*mp*gamma)
    return d

def vsh_rotation_matrix(n, alpha, beta, gamma):
    l = 2*n + 1
    R = np.zeros((l, l), dtype=complex)

    for i,m in enumerate(range(-n, n+1)):
        for j,mp in enumerate(range(-n, n+1)):
            R[i,j] = wigner_D(n, m, mp, alpha, beta, gamma) 

    return R

def rotate_coefficients(p_scat, alpha, beta, gamma):
    p_rot = np.empty_like(p_scat)
    rmax = p_scat.shape[1]
    lmax = miepy.vsh.rmax_to_lmax(rmax)


    for n in range(1,lmax+1):
        R = vsh_rotation_matrix(n, -alpha, -beta, -gamma)
        rmax = miepy.vsh.lmax_to_rmax(n)
        idx = np.s_[rmax-(2*n+1):rmax]
        p_rot[:,idx] = np.einsum('ab,pa->pb', R, p_scat[:,idx])
        p_rot[:,idx] = np.einsum('ab,pb->pa', R, p_scat[:,idx])

    return p_rot

# print(wigner_d(m, mp, n, beta))
# print(spherical_functions.wigner_d(n, m, mp, beta))
