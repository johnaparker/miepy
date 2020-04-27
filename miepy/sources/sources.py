"""
Abstract base classes for sources. Defines:

    source__________________________interface for all sources
    propagating_source______________sources that have a direction of propagation
    polarized_propagating_source____propagating_source with a global polarization state (TE, TM pair)
    combined_source_________________a suporposition of sources
"""

import numpy as np
from abc import ABCMeta, abstractmethod
import miepy

class source:
    """abstract base class for source objects"""
    __metaclass__ = ABCMeta

    def __init__(self, amplitude=1, phase=0, origin=None):
        """
        Arguments:
            amplitude   global amplitude factor (default 1)
            phase       global phase shift (default 0)
            origin      a reference point for the center of the source (default: [0,0,0])
        """
        self.amplitude = amplitude
        self.phase = phase

        if origin is None:
            self.origin = np.zeros(3, dtype=float)
        else:
            self.origin = np.asarray(origin, dtype=float)

    @abstractmethod
    def angular_spectrum(self, theta, phi, k):
        """Return the angular spectrum representation of the far-field source
        
        Arguments:
            theta    far-field theta angle
            phi      far-field phi angle
            k        wavenumber (in medium)
        """
        pass

    @abstractmethod
    def E_field(self, x1, x2, x3, k, far=False, spherical=False):
        """Compute the electric field of the source

        Arguments:
            x1        x (or r) position (array-like) 
            x2        y (or theta) position (array-like) 
            x3        z (or phi) position (array-like) 
            k         wavenumber (in medium)
            far       (optional) use the angular sprectrum for far-field calculations (bool, default=False)
            spherical (optional) input/output in spherical coordinates (bool, default=False)

        Returns: E[3,...]
        """
        pass

    @abstractmethod
    def H_field(self, x1, x2, x3, k, far=False, spherical=False):
        """Compute the magnetic field of the source

        Arguments:
            x1        x (or r) position (array-like) 
            x2        y (or theta) position (array-like) 
            x3        z (or phi) position (array-like) 
            k         wavenumber (in medium)
            far       (optional) use the angular sprectrum for far-field calculations (bool, default=False)
            spherical (optional) input/output in spherical coordinates (bool, default=False)

        Returns: E[3,...]
        """
        pass

    @abstractmethod
    def E_angular(self, theta, phi, k, radius=None, origin=None):
        """Compute the electric field in the far-field in spherical coordinates
             
        Arguments:
            theta    theta position (array-like) 
            phi      phi position (array-like) 
            k        wavenumber (in medium)
            radius   r position (default: large value)
            origin   origin around which to compute angular fields (default: self.origin)
        """
        pass

    @abstractmethod
    def H_angular(self, theta, phi, k, radius=None, origin=None):
        """Compute the magnetic field in the far-field in spherical coordinates
             
        Arguments:
            theta    theta position (array-like) 
            phi      phi position (array-like) 
            k        wavenumber (in medium)
            radius   r position (default: large value)
            origin   origin around which to compute angular fields (default: self.origin)
        """
        pass

    @abstractmethod
    def structure(self, position, k, lmax):
        """Compute the expansion coefficients of the source at a given position

        Arguments:
            position[N,3]   (x,y,z) position of the expansion origin for N points
            k               wavenumber (in medium)
            lmax            maximum expansion order
        """
        pass

    @abstractmethod
    def reflect(self, interface, medium, wavelength):
        """
        Create a source reflected by an interface
        """
        pass

    @abstractmethod
    def transmit(self, interface, medium, wavelength):
        """
        Create a source transmitted by an interface
        """
        pass

    def power_density(self, x1, x2, x3, k, far=False, spherical=False):
        """Compute the power density of the source

        Arguments:
            x1        x (or r) position (array-like) 
            x2        y (or theta) position (array-like) 
            x3        z (or phi) position (array-like) 
            far       (optional) use the angular sprectrum for far-field calculations (bool, default=False)
            spherical (optional) input/output in spherical coordinates (bool, default=False)

        Returns: E[3,...]
        """
        pass

    def H_angular_spectrum(self, theta, phi, k):
        H_inf = self.angular_spectrum(theta, phi, k)[::-1]
        H_inf[0] *= -1
        return H_inf

    def __add__(self, other):
        """Add two sources together"""
        if type(other) is combined_source:
            return combined_source(self, *other.sources)
        elif isinstance(other, source):
            return combined_source(self, other)
        else:
            raise ValueError('cannot add source with non-source of type {}'.format(type(other)))

class propagating_source(source):
    """abstract base class for propagating sources"""
    __metaclass__ = ABCMeta

    def __init__(self, amplitude=1, phase=0, origin=None, theta=0, phi=0, standing=False):
        source.__init__(self, amplitude=amplitude, phase=phase, origin=origin)

        self.theta = theta
        self.phi   = phi
        self.orientation = miepy.quaternion.from_spherical_coords(self.theta, self.phi)
        self.standing = standing

        ### TM and TE vectors
        self.k_hat, self.n_tm, self.n_te = miepy.coordinates.sph_basis_vectors(theta, phi)

class polarized_propagating_source(propagating_source):
    """abstract base class for polarized, propagating sources"""
    __metaclass__ = ABCMeta

    def __init__(self, polarization, amplitude=1, phase=0, origin=None, theta=0, phi=0, standing=False):
        propagating_source.__init__(self, amplitude=amplitude, phase=phase, origin=origin, theta=theta, phi=phi, standing=standing)

        self.polarization = np.asarray(polarization, dtype=np.complex)
        self.polarization /= np.linalg.norm(self.polarization)

    @abstractmethod
    def scalar_angular_spectrum(self, theta, phi, k):
        pass

    def angular_spectrum(self, theta, phi, k):
        U = self.scalar_angular_spectrum(theta, phi, k)

        Ex = U*self.polarization[0]
        Ey = U*self.polarization[1]
        Esph = np.array([-Ex*np.cos(phi) - Ey*np.sin(phi),
                         -Ex*np.sin(phi) + Ey*np.cos(phi)], dtype=complex)

        return Esph

class combined_source(source):
    """sources added together"""

    def __init__(self, *sources):
        self.sources = sources

    def __repr__(self):
        source_types = [type(s).__name__ for s in self.sources]
        return 'combined_source({})'.format(', '.join(source_types))

    def angular_spectrum(self, theta, phi, k):
        return sum((source.angular_spectrum(theta, phi, k) for source in self.sources))

    #TODO: alternatively, allow adding fields and performing structure only once
    def structure(self, position, k, lmax):
        return sum((source.structure(position, k, lmax) for source in self.sources))

    def E_field(self, x1, x2, x3, k, far=False, spherical=False):
        return sum((source.E_field(x1, x2, x3, k, far, spherical) for source in self.sources))

    def H_field(self, x1, x2, x3, k, far=False, spherical=False):
        return sum((source.H_field(x1, x2, x3, k, far, spherical) for source in self.sources))

    def E_angular(self, theta, phi, k, radius=None, origin=None):
        return sum((source.E_angular(theta, phi, k, radius, origin) for source in self.sources))

    def H_angular(self, theta, phi, k, radius=None, origin=None):
        return sum((source.H_angular(theta, phi, k, radius, origin) for source in self.sources))

    def reflect(self, interface, medium, wavelength):
        return combined_source(*[src.reflect(interface, medium, wavelength) for src in self.sources])

    def transmit(self, interface, medium, wavelength):
        return combined_source(*[src.transmit(interface, medium, wavelength) for src in self.sources])

    def __add__(self, other):
        if type(other) is combined_source:
            return combined_source(*self.sources, *other.sources)
        elif isinstance(other, source):
            return combined_source(*self.sources, other)
        else:
            raise ValueError('cannot add source with non-source of type {}'.format(type(other)))
