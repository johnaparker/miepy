import numpy as np

class material_struct:
    """struct to hold material data for N particles, and update it with given wavelength"""
    def __init__(self, materials, medium, wavelength=None):
        """Arguments:
               materials    list of materials
               medium       medium material
               wavelength   wavelength (default: None)
        """
        Nparticles = len(materials)

        self.materials = materials
        self.medium = medium
        self._wavelength = wavelength

        self.eps   = np.zeros(Nparticles, dtype=complex) 
        self.mu    = np.zeros(Nparticles, dtype=complex) 
        self.n     = np.zeros(Nparticles, dtype=complex) 
        self.eps_b = None
        self.mu_b  = None
        self.n_b   = None
        self.k_b   = None

        self.wavelength = wavelength

    @property
    def wavelength(self):
        """get the wavelength"""
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value):
        """set the wavelength to some value, changing all eps and mu data with it"""
        self._wavelength = value

        if value is not None:
            self.eps_b = self.medium.eps(self._wavelength)
            self.mu_b  = self.medium.mu(self._wavelength)
            self.n_b = (self.eps_b*self.mu_b)**0.5
            self.k_b = 2*np.pi*self.n_b/self._wavelength

            for i,material in enumerate(self.materials):
                self.eps[i] = material.eps(self._wavelength)
                self.mu[i]  = material.mu(self._wavelength)

            self.n[...]   = (self.eps*self.mu)**0.5
