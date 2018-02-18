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
