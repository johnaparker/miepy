"""
performance test
"""

import numpy as np
import matplotlib.pyplot as plt
import miepy
from tqdm import tqdm

from timeit import default_timer as timer
import timeit


nm = 1e-9

Ag = miepy.materials.predefined.Ag()
radius = 75*nm
source = miepy.sources.rhc_polarized_plane_wave(amplitude=1)
separation = 250*nm

def force_func(N):
    spheres = miepy.spheres([[n*separation, 0, 0] for n in range(N)], radius, Ag)
    mie = miepy.gmt(spheres, source, 600*nm, 1)
    
    start = timer()
    mie.force()
    end = timer()

    return end - start

def flux_func(N):
    spheres = miepy.spheres([[n*separation, 0, 0] for n in range(N)], radius, Ag)
    mie = miepy.gmt(spheres, source, 600*nm, 1, origin=[separation/2,0,0])
    
    start = timer()
    mie.cross_sections()
    end = timer()

    return end - start

def interactions_func(N):
    spheres = miepy.spheres([[n*separation, 0, 0] for n in range(N)], radius, Ag)
    mie = miepy.gmt(spheres, source, 600*nm, 1)
    
    start = timer()
    mie._solve_interactions()
    end = timer()

    return end - start

def N_particle_test(func, Nmax = 20):
    print(func.__name__)
    for N in range(1, Nmax+1):
        time = func(N)
        print('  ', f'N = {N}: ', f'{time:.3f}s')

N_particle_test(interactions_func, 10)
N_particle_test(force_func, 10)
N_particle_test(flux_func, 10)
