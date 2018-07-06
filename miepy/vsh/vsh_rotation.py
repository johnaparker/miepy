"""
VSH rotation function
"""

import numpy as np
from math import factorial
import miepy
from miepy import vsh
from scipy import special
from smuthi import spherical_functions, field_expansion

n = 4
m = 3
mp = 1
beta = 0.4

def wigner_d(n, m, mp, beta):
    d = spherical_functions.wigner_d(n, m, mp, beta)

    # f = factorial
    # factor = (f(n-m)*f(n+m)/f(n+mp)/f(n-mp))**0.5

    # if m-mp < 0 or m+mp < 0:
        # return 0

    # d = factor*np.sin(beta/2)**(m-mp)*np.cos(beta/2)**(m+mp) \
        # * special.jacobi(n-m, m-mp, m+mp)(np.cos(beta))

    return d

def wigner_D(n, m, mp, alpha, beta, gamma):
    d = spherical_functions.wigner_D(n, m, mp, alpha, beta, gamma)
    d = np.exp(1j*m*alpha)*wigner_d(n, m, mp, beta)*np.exp(1j*mp*gamma)
    return d

def vsh_rotation_matrix(n, alpha, beta, gamma):
    l = 2*n + 1
    # return field_expansion.block_rotation_matrix_D_svwf(n, n, alpha, beta, gamma)[:l,:l]

    R = np.zeros((l, l), dtype=complex)

    for i,m in enumerate(range(-n, n+1)):
        for j,mp in enumerate(range(-n, n+1)):
            R[i,j] = wigner_D(n, m, mp, alpha, beta, gamma) 

    return R

def rotate_coefficients(p_scat, alpha, beta, gamma):
    beta *= -1
    alpha *= -1
    gamma *= -1
    p_rot = np.empty_like(p_scat)
    rmax = p_scat.shape[1]
    lmax = miepy.vsh.rmax_to_lmax(rmax)


    for n in range(1,lmax+1):
        R = vsh_rotation_matrix(n, alpha, beta, gamma)
        rmax = miepy.vsh.lmax_to_rmax(n)
        idx = np.s_[rmax-(2*n+1):rmax]
        p_rot[:,idx] = np.einsum('ab,pa->pb', R, p_scat[:,idx])
        p_rot[:,idx] = np.einsum('ab,pb->pa', R, p_scat[:,idx])

    return p_rot

# print(wigner_d(m, mp, n, beta))
# print(spherical_functions.wigner_d(n, m, mp, beta))
