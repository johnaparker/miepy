import miepy
from .particle_base import particle

class sphere(particle):
    def __init__(self, position, radius, material):
        """A sphere object

        Arguments:
            position[3]   x,y,z position of particle
            radius        sphere radius
            material      particle material (miepy.material object)
        """
        super().__init__(position, None, material)
        self.radius = radius

    def is_inside(self, pos):
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        self.tmatrix = miepy.tmatrix.tmatrix_sphere(self.radius, wavelength, 
                self.material.eps(wavelength), eps_m, lmax)
        self.tmatrix_fixed = self.tmatrix

        return self.tmatrix

    def enclosed_radius(self):
        return self.radius

    def _dict_key(self, wavelength):
        return (sphere, self.radius, self.material.eps(wavelength).item(), self.material.mu(wavelength).item())
