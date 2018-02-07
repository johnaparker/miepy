import numpy as np
import miepy
from sympy.physics.wigner import wigner_3j, gaunt
from scipy import special
import math

# w = miepy.vsh.wigner_3j(1,1,1, 1,-1,0)
# print(w)
# w = wigner_3j(1,1,1,1,-1,0)
# print(w)
# from IPython import embed; embed()

def gaunt(m,n,u,v,p):
    """gaunt coefficient"""

    f = lambda n: special.gamma(n+1)
    numerator = f(n+m)*f(v+u)*f(p-m-u)
    denominator = f(n-m)*f(v-u)*f(p+m+u)
    factor = (-1)**(m+u)*(2*p+1)*(numerator/denominator)**0.5

    w1 = miepy.vsh.wigner_3j(n,v,p,0,0,0)
    w2 = miepy.vsh.wigner_3j(n,v,p,m,u,-m-u)

    return factor*w1*w2

# g = gaunt(1,2,0,1,1)
# print(g)
# g = miepy.vsh.a_func(1,2,0,1,1)
# print(g)
# from IPython import embed; embed()

def A_translation(m, n, u, v, r, theta, phi, k):
    f = lambda n: special.gamma(n+1)
    numerator = (2*v+1)*f(n-m)*f(v-u)
    denominator = 2*n*(n+1)*f(n+m)*f(v+u)

    factor = (-1)**m * numerator/denominator*np.exp(1j*(u+m)*phi)

    qmax = min(n, v, math.floor((n+v - abs(m+u))/2))
    sum_term = 0
    for q in range(0, qmax+1):
        p = n + v - 2*q
        aq = gaunt(m,n,u,v,p)
        A = 1j**p*(n*(n+1) + v*(v+1) - p*(p+1))*aq

        Pnm = miepy.vsh.associated_legendre(p,u+m)
        sum_term += A*miepy.vsh.spherical_hn(p, k*r)*Pnm(np.cos(theta))

    return factor*sum_term 

A = A_translation(0,10,0,10, 2, 0.5, 0.5, 1)
A = A_translation(2, 6, -2, 10, 2, 0.5, 0.5, 1)
A = A_translation(0, 1, 1, 1, 2, 0.5, 0.5, 1)
print(A)
A = miepy.vsh.A_translation(0, 1, 1, 1, 2, 0.5, 0.5, 1)
print(A)