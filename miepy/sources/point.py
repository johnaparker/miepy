"""
point sources
"""

import numpy as np
import miepy
from miepy.sources import source

#TODO: implement for any polarization
class point_dipole(source):
    def __init__(self, position, direction, amplitude=1, phase=0, mode='electric'):
        super().__init__(amplitude, phase)

        self.position = np.asarray(position)

        direction = np.asarray(direction, dtype=np.complex)
        self.direction = direction
        self.direction /= np.linalg.norm(direction)

        self.mode = mode

        self.weight = {}
        self.weight[1] = -1j*(self.direction[0] - 1j*self.direction[1])/2
        self.weight[0] = -1j*self.direction[2]
        self.weight[-1] = 1j*(self.direction[0] + 1j*self.direction[1])/2

    def structure(self, position, k, lmax, radius=None):
        rmax = miepy.vsh.lmax_to_rmax(lmax)
        p_src = np.zeros([2, rmax], dtype=complex)
        factor = self.amplitude*np.exp(1j*self.phase)

        dr = position - self.position
        rad, theta, phi = miepy.coordinates.cart_to_sph(*dr)

        for r,n,m in miepy.mode_indices(lmax):
            for u in [-1,0,1]:
                A, B = miepy.cpp.vsh_translation.vsh_translation(m, n, u, 1, rad, theta, phi, k, miepy.vsh_mode.outgoing)
                p_src[0,r] += -A*self.weight[u]
                p_src[1,r] += -B*self.weight[u]

        if self.mode == 'magnetic':
            p_src = p_src[::-1]

        return factor*p_src
    
    def E_field(self, x, y, z, k):
        r, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position)

        Esph = np.zeros((3,) + x.shape, dtype=complex)

        p_src = np.zeros([2,3], dtype=complex)
        p_src[0] = (self.weight[-1], self.weight[0], self.weight[1])

        if self.mode == 'magnetic':
            p_src = p_src[::-1]

        Efunc = miepy.vsh.expand_E(p_src, k, miepy.vsh_mode.outgoing)

        Esph = Efunc(r, theta, phi)
        Ecart = miepy.coordinates.vec_sph_to_cart(Esph, theta, phi) 

        return self.amplitude*Ecart*np.exp(1j*self.phase)

    def H_field(self, x, y, z, k):
        r, theta, phi = miepy.coordinates.cart_to_sph(x, y, z, origin=self.position)

        Esph = np.zeros((3,) + x.shape, dtype=complex)

        p_src = np.zeros([2,3], dtype=complex)
        p_src[0] = (self.weight[-1], self.weight[0], self.weight[1])

        if self.mode == 'magnetic':
            p_src = p_src[::-1]

        Efunc = miepy.vsh.expand_H(p_src, k, miepy.vsh_mode.outgoing, 1, 1)

        Esph = Efunc(r, theta, phi)
        Ecart = miepy.coordinates.vec_sph_to_cart(Esph, theta, phi) 

        return self.amplitude*Ecart*np.exp(1j*self.phase)
