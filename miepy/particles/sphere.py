import miepy
from .particle_base import particle

class sphere(particle):
    def __init__(position, radius, material):
        """A sphere object

        Arguments:
            position[3]   x,y,z position of particle
            radius        sphere radius
            material      particle material (miepy.material object)
        """
        super().__init__(position, None, material)
        self.radius = radius

    def is_inside(pos):
        pass

    def compute_tmatrix(lmax, wavelength, eps_m, **kwargs):
        self.tmatrix = miepy.tmatrix.tmatrix_sphere()
        return self.tmatrix
