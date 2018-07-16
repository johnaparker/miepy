import miepy
import numpy as np
from .particle_base import particle

class cylinder(particle):
    def __init__(self, position, radius, height, material, orientation=None, rounded=False, extended_precision=False, Nint=200,
            lmax=None):
        """A cylinder object

        Arguments:
            position[3]   x,y,z position of particle
            radius        cylinder radius
            height        cylinder height
            material      particle material (miepy.material object)
            orientation   particle orientation
        """
        super().__init__(position, orientation, material)
        self.radius = radius
        self.height = height
        self.extended_precision = extended_precision
        self.Nint = 200
        self.rounded = rounded
        self.lmax = lmax

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        if self.lmax is None:
            lmax_compute = lmax
        else:
            lmax_compute = self.lmax

        self.tmatrix = miepy.tmatrix.tmatrix_cylinder(self.radius, self.height, wavelength, 
                self.material.eps(wavelength), eps_m, lmax_compute, rounded=self.rounded, extended_precision=self.extended_precision,
                Nint=self.Nint)

        if self.lmax is not None:
            self.tmatrix = miepy.tmatrix.tmatrix_reduce_lmax(self.tmatrix, lmax)

        return self.tmatrix

    def enclosed_radius(self):
        return np.sqrt((self.height/2)**2 + self.radius**2)
