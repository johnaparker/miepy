"""
Source abstract base class
"""

import numpy as np
from abc import ABCMeta, abstractmethod
import miepy

#TODO: in addition to an amplitude, all sources have a phase, polarization (maybe), and direction of propagation (maybe)
class source:
    """source interface base class"""
    __metaclass__ = ABCMeta

    def __init__(self, amplitude):
        self.amplitude = amplitude

    @abstractmethod
    def E_field(self, x, y, z, k): pass

    @abstractmethod
    def H_field(self, x, y, z, k): pass

    @abstractmethod
    def structure(self, position, k, Lmax, radius): pass

    @abstractmethod
    def is_paraxial(self, k): pass

    def spherical_ingoing(self, theta, phi, k):
        raise NotImplementedError('source has not defined a spherical ingoing function')

    def __add__(self, other):
        return combined_source(self, other)

class combined_source(source):
    """sources added together"""

    def __init__(self, *sources):
        self.sources = sources
        self.amplitude = sources[0].amplitude

    def E_field(self, x, y, z, k):
        return sum(map(lambda source: source.E_field(x, y, z, k), self.sources))

    def H_field(self, x, y, z, k):
        return sum(map(lambda source: source.H_field(x, y, z, k), self.sources))

    #TODO: alternatively, allow adding fields and performing structure only once
    def structure(self, position, k, Lmax, size):
        return sum((source.structure(position, k, Lmax, size) for source in self.sources))

    def spherical_ingoing(self, theta, phi, k):
        return sum((source.spherical_ingoing(theta, phi, k) for source in self.sources))

    def is_paraxial(self, k):
        return all((source.is_paraxial(k) for source in sources))
