import miepy
from .particle_base import particle

class cylinder(particle):
    def __init__(position, radius, height, material, orientation=None):
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

    def is_inside(pos):
        pass

    def compute_tmatrix(lmax, wavelength, eps_m, **kwargs):
        self.tmatrix = miepy.tmatrix.tmatrix_cylinder(self.radius, self.height, wavelength, 
                self.material.eps(wavelength), eps_m, lmax, rounded=False, extended_precision=False)
        return self.tmatrix
