"""
tests for the VSH translation functions
"""

import numpy as np
import miepy
from math import factorial

r_ji = 1
theta_ji = 0.7
phi_ji = 0.3

Lmax = 2
k = 0.8


for r,n,m in miepy.mode_indices(Lmax):
    for s,v,u in miepy.mode_indices(Lmax):

        ### original A and B
        A, B = miepy.vsh.vsh_translation(m, n, u, v, 
                r_ji, theta_ji, phi_ji, k, miepy.vsh.VSH_mode.outgoing)

        ### same mode, but reverse position
        A_calc, B_calc = miepy.vsh.vsh_translation(m, n, u, v, 
                r_ji, np.pi - theta_ji, np.pi + phi_ji, k, miepy.vsh.VSH_mode.outgoing)
        A_expect, B_expect = (-1)**(n+v)*A, (-1)**(n+v+1)*B

        # assert np.allclose(A_calc, A_expect)
        # assert np.allclose(B_calc, B_expect)

        ### swap (m,n) with (u,v) and negate m and u
        A_calc, B_calc = miepy.vsh.vsh_translation(-u, v, -m, n, 
                r_ji, theta_ji, phi_ji, k, miepy.vsh.VSH_mode.outgoing)
        A_expect, B_expect = (-1)**(m+u)*A, (-1)**(m+u+1)*B

        # assert np.allclose(A_calc, A_expect)
        # assert np.allclose(B_calc, B_expect)

        ### swap (m,n) with (u,v)
        A_calc, B_calc = miepy.vsh.vsh_translation(u, v, m, n, 
                r_ji, np.pi - theta_ji, np.pi + phi_ji, k, miepy.vsh.VSH_mode.outgoing)
        factor = (-1)**(n+v) #*factorial(n-m)*factorial(v+u)/(factorial(n+m)*factorial(v-u))
        A_expect, B_expect = factor*A, factor*B

        print(f'{(n,m)} -> {(v,u)}: A = {A_calc}, Ap = {A_expect}, B = {B_calc}, Bp = {B_expect}')
