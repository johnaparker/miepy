"""
Utility functions
"""

import numpy as np

def atleast(array, dim, length, dtype=None):
    """Given an n-dimensional array, return either an n or n+1 dimensional repeated array
        array      input array (or scalar)
        dim        dimension of the returned array (must be n or n+1)
        dtype      array datatype (default None)
    """

    ret = np.asarray(np.atleast_1d(array), dtype=dtype)
    if (dim not in (ret.ndim,ret.ndim+1)):
        raise ValueError('dim = {0} is invalid with input dim = {1}'.format(dim, ret.ndim))

    if len(ret.shape) == dim and ret.shape[0] != length:
        ret = np.repeat(ret, length)
    elif len(ret.shape) != dim:
        ret = np.array([ret]*length)

    return ret
