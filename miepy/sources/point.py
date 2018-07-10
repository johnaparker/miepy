"""
point sources
"""

import numpy as np
import miepy
from miepy.sources.source_base import source

#TODO: implement for any polarization
class point_dipole(source):
    def __init__(self, location, polarization, amplitude=1, phase=0):
        super().__init__(amplitude, phase)

        self.location = location

        polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization = polarization
        self.polarization /= np.linalg.norm(polarization)
    
    #TODO: VSH should allow spherical=False or cartesian=True
    def E_field(self, x, y, z, k):
        N, _ = miepy.vsh.VSH(1,0)
        r, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.location)
        N_cart = miepy.coordinates.vec_sph_to_cart(N(r, theta, phi, k), theta, phi) 
        return self.amplitude*N_cart*np.exp(1j*self.phase)

    def H_field(self, x, y, z, k):
        _, M = miepy.vsh.VSH(1,0)
        r, theta, phi = miepy.coordinates.cart_to_sph(x, y , z, origin=self.location)
        M_cart = miepy.coordinates.vec_sph_to_cart(M(r, theta, phi, k), theta, phi) 
        return self.amplitude*M_cart*np.exp(1j*self.phase)
