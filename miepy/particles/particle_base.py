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
        self.position = position
        self.orientation = orientation
        self.material = material
        self.tmatrix = None

    def is_inside(self, pos):
        """Return true if pos is inside the particle"""
        pass

    def enclosed_radius(self):
        """Return the radius of the smallest circumscribing sphere"""
        pass

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        """Compute the tmatrix of the particle
        
        Arguments:
            lmax         maximum number of multipoles
            wavelength   incident wavelength
            eps_m        permitvittiy of the medium
            kwargs       additional kwargs to pass
        """
        pass

    def translate(self, dr):
        """Translate the particle position by dr"""
        pass

    def rotate(self, dn):
        """Rotate the particle orientation by dn (and rotate the Tmatrix if it has been calculated)"""
        pass
