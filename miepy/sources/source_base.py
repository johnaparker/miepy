"""
Source abstract base class
"""

import numpy as np
from abc import ABCMeta, abstractmethod

class source:
    """source interface base class"""
    __metaclass__ = ABCMeta

    def __init__(self, amplitude):
        self.amplitude = amplitude

    @abstractmethod
    def E(self, r, k): pass

    @abstractmethod
    def H(self, r, k): pass

    @abstractmethod
    def structure_of_mode(self, n, m, r, k): pass

    def structure(self, position, k, Nmax):
        p = np.zeros([Nmax, 2*Nmax+1], dtype=complex)
        q = np.zeros([Nmax, 2*Nmax+1], dtype=complex)
        for n in range(1, Nmax+1):
            for m in range(-n,n+1):
                p[n-1,m+n], q[n-1,m+n] = self.structure_of_mode(n, m, position, k)

        return p,q

    def __add__(self, other):
        return combined_source(self, other)

class combined_source(source):
    """sources added together"""

    def __init__(self, *sources):
        self.sources = sources
        self.amplitude = sources[0].amplitude

    def E(self, r, k):
        return sum(map(lambda source: source.E(r,k), self.sources))

    def H(self, r, k):
        return sum(map(lambda source: source.H(r,k), self.sources))
    
    def structure_of_mode(self, n, m, r, k):
        structure = map(lambda source: source.structure_of_mode(n, m, r, k), self.sources)
        p_all = sum([p for p,q in structure])
        q_all = sum([q for p,q in structure])

        return p_all, q_all
