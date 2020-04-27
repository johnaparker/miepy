"""
point sources
"""

import numpy as np
import miepy
from miepy.sources import source
from functools import partial

class point_dipole(source):
    def __init__(self, position, direction, amplitude=1, phase=0, mode='electric'):
        """
        Arguments:
            position     (x,y,z) position of the dipole
            direction    (Dx, Dy, Dz) direction the dipole points (can be complex)
            amplitude    amplitude of the dipole (default: 1)
            phase        additional phase factor (default: 0)
            mode         'electric' or 'magnetic' dipole (default: electric)
        """
        super().__init__(amplitude, phase, origin=position)

        self.position = np.asarray(position, dtype=float)

        direction = np.asarray(direction, dtype=complex)
        self.direction = direction
        self.direction /= np.linalg.norm(direction)

        self.mode = mode

        self.weight = {}
        self.weight[1] = -1j*(self.direction[0] - 1j*self.direction[1])/2
        self.weight[0] = -1j*self.direction[2]
        self.weight[-1] = 1j*(self.direction[0] + 1j*self.direction[1])/2

    def angular_spectrum(self, theta, phi, k):
        return self.E_angular(theta, phi, k)

    def E_field(self, x1, x2, x3, k, far=False, spherical=False):
        if far:
            Efunc = miepy.vsh.expand_E_far
        else:
            Efunc = partial(miepy.vsh.expand_E, mode=miepy.vsh_mode.outgoing)

        return self._field(Efunc, x1, x2, x3, k, spherical=spherical)

    def H_field(self, x1, x2, x3, k, far=False, spherical=False):
        if far:
            Hfunc = partial(miepy.vsh.expand_H_far, eps=1, mu=1)
        else:
            Hfunc = partial(miepy.vsh.expand_H, mode=miepy.vsh_mode.outgoing, eps=1, mu=1)

        return self._field(Hfunc, x1, x2, x3, k, spherical=spherical)

    #TODO: implement origin
    def E_angular(self, theta, phi, k, radius=None, origin=None):
        if radius is None:
            radius = 1e6*2*np.pi/k

        return self.E_field(radius, radius, theta, phi, far=True, spherical=False)[1:]

    #TODO: implement origin
    def H_angular(self, theta, phi, k, radius=None, origin=None):
        if radius is None:
            radius = 1e6*2*np.pi/k

        return self.H_field(radius, radius, theta, phi, far=True, spherical=False)[1:]

    def structure(self, position, k, lmax):
        position = np.asarray(position)
        Nparticles = len(position)

        rmax = miepy.vsh.lmax_to_rmax(lmax)
        p_src = np.zeros([Nparticles, 2, rmax], dtype=complex)
        factor = self.amplitude*np.exp(1j*self.phase)

        for i in range(Nparticles):
            dr = position[i] - self.position
            rad, theta, phi = miepy.coordinates.cart_to_sph(*dr)

            for r,n,m in miepy.mode_indices(lmax):
                for u in [-1,0,1]:
                    A, B = miepy.cpp.vsh_translation.vsh_translation(m, n, u, 1, rad, theta, phi, k, miepy.vsh_mode.outgoing)
                    p_src[i,0,r] += -A*self.weight[u]
                    p_src[i,1,r] += -B*self.weight[u]

        if self.mode == 'magnetic':
            p_src = p_src[:, ::-1]

        return factor*p_src

    def reflect(self, interface, medium, wavelength):
        raise NotImplementedError('Dipole interface reflections have not yet been implemented')

    def transmit(self, interface, medium, wavelength):
        raise NotImplementedError('Dipole interface transmission has not yet been implemented')
    
    def _field(self, Efunc, x1, x2, x3, k, spherical=False):
        if not spherical:
            rad, theta, phi = miepy.coordinates.cart_to_sph(x1, x2, x3, origin=self.position)
        else:
            rad, theta, phi = x1, x2, x3

        p_src = np.zeros([2,3], dtype=complex)
        p_src[0] = (self.weight[-1], self.weight[0], self.weight[1])

        if self.mode == 'magnetic':
            p_src = p_src[::-1]

        E = Efunc(p_src, k)(rad, theta, phi)
        if not spherical:
            E = miepy.coordinates.vec_sph_to_cart(E, theta, phi) 

        return E*self.amplitude*np.exp(1j*self.phase)
