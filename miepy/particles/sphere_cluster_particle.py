import miepy
import numpy as np
from .particle_base import particle

class sphere_cluster_particle(particle):
    def __init__(self, position, radius, material, lmax, orientation=None):
        """A 'particle' that is a rigid collection of non-overlaping spheres

        Arguments:
            position[N,3]   x,y,z position of the spheres
            radius[N]       radii of the spheres
            material[N]     sphere materials
            lmax[N]         sphere lmax values
            orientation     orientation of the collection (relative to as created)
        """
        self.com = np.average(position, axis=0)
        self.Nparticles = len(position)

        self.p_position = np.atleast_2d(position)

        self.p_radii = np.empty(self.Nparticles, dtype=float)
        self.p_radii[...] = radius

        self.p_material = np.empty(self.Nparticles, dtype=object)
        self.p_material[...] = material

        self.p_lmax = np.empty(self.Nparticles, dtype=int)
        self.p_lmax[...] = lmax

        self.lmax_cluster = np.max(self.p_lmax)
        super().__init__(self.com, orientation, self.p_material[0])

    def __repr__(self):
        return f'''{self.__class__.__name__}:
    Nparticles = self.Nparticles
    position = {self.position}
    orientation = {self.orientation}
    '''

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        eps = np.empty(self.Nparticles, dtype=complex)
        for i in range(self.Nparticles):
            eps[i] = self.p_material[i].eps(wavelength)
        
        self.tmatrix_fixed = miepy.tmatrix.tmatrix_sphere_cluster(self.p_position, self.p_radii, self.p_lmax,
                self.lmax_cluster, wavelength, eps, eps_m, extended_precision=False)

        self.tmatrix_fixed = miepy.tmatrix.tmatrix_reduce_lmax(self.tmatrix_fixed, lmax)

        self._rotate_fixed_tmatrix()
        return self.tmatrix

    def enclosed_radius(self):
        return np.max(np.linalg.norm(self.p_position - self.position[:,np.newaxis], axis=1)) \
               + np.max(self.radii)
