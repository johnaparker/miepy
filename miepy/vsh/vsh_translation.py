"""
VSH translation function
"""

import numpy as np
from math import factorial
from miepy import vsh

def vsh_translation(m, n, u, v, r, theta, phi, k, mode):
    """VSH translation coefficients"""
    m *= -1
    f = factorial
    zn = vsh.get_zn(mode)

    factor = 0.5 * (-1.)**m * np.sqrt((2*v+1)*(2*n+1)*f(v-u)*f(n-m)
            /(v*(v+1)*n*(n+1)*f(v+u)*f(n+m))) * np.exp(1j*(u+m)*phi)

    qmax = min(n, v, (n+v - abs(m+u))//2)
    sum_term = 0
    for q in range(0, qmax+1):
        p = n + v - 2*q
        aq = vsh.special.a_func(m,n,u,v,p)
        A = 1j**p*(n*(n+1) + v*(v+1) - p*(p+1))*aq

        Pnm = vsh.special.associated_legendre(p,u+m)
        sum_term += A*zn(p, k*r)*Pnm(np.cos(theta))

    A_translation = factor*sum_term

    qmax = min(n, v, (n+v+1 - abs(m+u))//2)
    sum_term = 0
    for q in range(1, qmax+1):
        p = n + v - 2*q
        bq = vsh.special.b_func(m,n,u,v,p)
        A = 1j**(p+1)*(((p+1)**2 - (n-v)**2)*((n+v+1)**2 - (p+1)**2))**0.5*bq

        Pnm = vsh.special.associated_legendre(p+1,u+m)
        sum_term += A*zn(p+1, k*r)*Pnm(np.cos(theta))

    B_translation = -factor*sum_term

    return A_translation, B_translation
