import miepy
import numpy as np
from .particle_base import particle

class cylinder(particle):
    def __init__(self, position, radius, height, material, orientation=None):
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

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        self.tmatrix = miepy.tmatrix.tmatrix_cylinder(self.radius, self.height, wavelength, 
                self.material.eps(wavelength), eps_m, lmax, rounded=False, extended_precision=False)
        return self.tmatrix

    def enclosed_radius(self):
        return np.sqrt((self.height/2)**2 + self.radius**2)
