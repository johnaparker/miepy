"""
Source abstract base class
"""

import numpy as np
from abc import ABCMeta, abstractmethod
import miepy

class source:
    """source interface base class"""
    __metaclass__ = ABCMeta

    def __init__(self, amplitude):
        self.amplitude = amplitude

    @abstractmethod
    def E_field(self, x, y, z, k): pass

    @abstractmethod
    def H_field(self, x, y, z, k): pass

    def structure(self, position, k, Lmax):
        p_src = miepy.sources.decomposition.point_matching(self,
                      position, radius=1e-9, k=k, Lmax=Lmax, sampling=6)
        return p_src

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
