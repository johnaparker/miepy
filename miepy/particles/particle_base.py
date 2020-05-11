import numpy as np
import miepy

#TODO: lmax per particle
#TODO: position and orientation should be properties
class particle:
    def __init__(self, position, orientation, material):
        """A particle consists of a position, orientation, material, and a lazily evaluated T-matrix

        Arguments:
            position[3]   x,y,z position of particle
            orientation   particle orientation
            material      particle material (miepy.material object)
        """

        self._position = np.asarray(position, dtype=float)

        if orientation is None:
            self._orientation = miepy.quaternion.one
        else:
            self._orientation = orientation

        self.material = material
        self.tmatrix_fixed = None
        self.tmatrix = None

        if self.material is not None and self.material.name == 'metal':
            self.conducting = True
        else:
            self.conducting = False

    @property
    def position(self):
        return self._position

    @position.setter
    def position(self, p):
        self.position[...] = p

    @property
    def orientation(self):
        return self._orientation

    @orientation.setter
    def orientation(self, n):
        self._orientation = n
        self._rotate_fixed_tmatrix()

    def is_inside(self, pos):
        """Return true if pos is inside the particle"""
        pass

    def enclosed_radius(self):
        """Return the radius of the smallest circumscribing sphere"""
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        """Compute the T-matrix of the particle
        
        Arguments:
            lmax         maximum number of multipoles
            wavelength   incident wavelength
            eps_m        permitvittiy of the medium
            kwargs       additional kwargs to pass
        """
        pass

    def _rotate_fixed_tmatrix(self):
        if self.tmatrix_fixed is not None:
            self.tmatrix = miepy.tmatrix.rotate_tmatrix(self.tmatrix_fixed, self.orientation)

    def _dict_key(self, wavelength):
        pass
