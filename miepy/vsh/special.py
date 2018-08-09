"""
special functions required by vsh calculations
"""

import numpy as np
import sympy
from scipy import special
from functools import lru_cache
from math import factorial

from miepy.cpp.special import (spherical_hn, spherical_hn_2,
         wigner_3j, associated_legendre, pi_func, tau_func)

def riccati_1(n,z, derivative = False):
    jn = special.spherical_jn(n, z)

    if derivative:
        jn_p = special.spherical_jn(n, z, derivative=True)
        return z*jn_p + jn
    return z*jn

def riccati_2(n,z, derivative = False):
    yn = special.spherical_yn(n, z, derivative = derivative) 

    if derivative:
        yn_p = special.spherical_yn(n, z, derivative=True)
        return -z*yn_p - yn
    return -z*yn

def riccati_3(n,z, derivative = False):
    return riccati_2(n,z, derivative) - riccati_1(n,z, derivative)

@lru_cache(None)
def a_func(m,n,u,v,p):
    """gaunt coefficient"""

    f = factorial
    numerator = f(n+m)*f(v+u)*f(p-m-u)
    denominator = f(n-m)*f(v-u)*f(p+m+u)
    factor = (-1.)**(m+u)*(2*p+1)*(numerator/denominator)**0.5

    w1 = wigner_3j(n,v,p,0,0,0)
    w2 = wigner_3j(n,v,p,m,u,-m-u)

    return factor*w1*w2

@lru_cache(None)
def b_func(m,n,u,v,p):
    """b function"""

    f = factorial
    numerator = f(n+m)*f(v+u)*f(p-m-u+1)
    denominator = f(n-m)*f(v-u)*f(p+m+u+1)
    factor = (-1.)**(m+u)*(2*p+3)*(numerator/denominator)**0.5

    w1 = wigner_3j(n,v,p,0,0,0)
    w2 = wigner_3j(n,v,p+1,m,u,-m-u)

    return factor*w1*w2
