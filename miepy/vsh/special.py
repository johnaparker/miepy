"""
special functions required by vsh calculations
"""

import numpy as np
import sympy
from scipy import special
from functools import lru_cache
from math import factorial

def spherical_hn(n, z, derivative=False):
    """spherical hankel function of the first kind or its derivative

            n: int,array-like        order of the bessel function
            z: array[complex/float]  argument
            derivative: bool         If True, compute the derivative instead    """

    return special.spherical_jn(n,z,derivative) + 1j*special.spherical_yn(n,z,derivative)

def spherical_hn_2(n, z, derivative=False):
    """spherical hankel function of the second kind or its derivative

            n: int,array-like        order of the bessel function
            z: array[complex/float]  argument
            derivative: bool         If True, compute the derivative instead    """

    return special.spherical_jn(n,z,derivative) - 1j*special.spherical_yn(n,z,derivative)

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

@lru_cache(maxsize=None)
def associated_legendre(n,m, deriv=0):
    """associated legendre function of integer order and degree

            n: int         order
            m: int         degree
            deriv: int     derivative to take

        returns lpmv(x) function     """

    x = sympy.symbols('x')
    legfun_sym = (-1)**abs(m)*sympy.functions.special.polynomials.assoc_legendre(n,m,x)
    legfunc_sym_deriv = sympy.diff(legfun_sym, x, deriv)

    legfun_num_deriv = sympy.lambdify(x,legfunc_sym_deriv, modules='numpy')
    return legfun_num_deriv

@lru_cache(maxsize=None)
def pi_func(n, m):
    """pi special function that appears in the vector spherical harmonics (VSH)

            n: int         order
            m: int         degree

       returns pi(theta)                                                           """

    x = sympy.symbols('x')
    legfunc_sym = m*(-1)**abs(m)*sympy.functions.special.polynomials.assoc_legendre(n,m,x)
    legfunc_sym /= sympy.sqrt(1 - x**2)
    legfunc_sym = sympy.simplify(legfunc_sym)
    legfunc_sym = legfunc_sym.subs(x, sympy.cos(x))

    legfun_num = sympy.lambdify(x,legfunc_sym, modules='numpy')
    return legfun_num

@lru_cache(maxsize=None)
def tau_func(n, m):
    """pi special function that appears in the vector spherical harmonics (VSH)

            n: int         order
            m: int         degree

       returns tau(theta)                                                           """

    x = sympy.symbols('x')
    legfun_sym = (-1)**abs(m)*sympy.functions.special.polynomials.assoc_legendre(n,m,x)
    legfunc_sym_deriv = sympy.diff(legfun_sym, x, 1)
    legfunc_sym_deriv *= -1*sympy.sqrt(1 - x**2)
    legfunc_sym_deriv = sympy.expand(legfunc_sym_deriv)
    legfunc_sym_deriv = legfunc_sym_deriv.subs(x, sympy.cos(x))

    legfun_num_deriv = sympy.lambdify(x,legfunc_sym_deriv, modules='numpy')
    return legfun_num_deriv

def wigner_3j(j1, j2, j3, m1, m2, m3):
    """wigner3j coefficients"""
    kmin = max(0, j2-j3-m1, j1-j3+m2)
    kmax = max(j1+j2-j3, j1-m1, j2+m2)

    if m1 + m2 + m3 != 0:
        return 0

    j3min = max(abs(j1-j2), abs(m1+m2))
    if not j3min <= j3 <= j1 + j2:
        return 0

    f = factorial
    binomial = special.binom
    numerator = f(j1-m1)*f(j1+m1)*f(j2-m2)*f(j2+m2)*f(j3-m3)*f(j3+m3)
    denominator = f(j1+j2-j3)*f(j1-j2+j3)*f(-j1+j2+j3)*f(j1+j2+j3+1)
    factor = (-1)**(j1+j2+m3)*(numerator/denominator)**0.5

    sum_term = 0
    for k in range(kmin, kmax+1):
        b1 = binomial(j1+j2-j3, k)
        b2 = binomial(j1-j2+j3, j1-m1-k)
        b3 = binomial(-j1+j2+j3, j2+m2-k)
        sum_term += (-1)**k * b1 * b2 * b3
    
    return factor*sum_term

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
