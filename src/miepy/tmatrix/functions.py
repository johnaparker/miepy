import numpy as np

import miepy


def tmatrix_reduce_lmax(tmatrix, lmax):
    """Reduce the dimensions of a tmatrix to lmax.

    Arguments:
        tmatrix[2,rmax,2,rmax]    the tmatrix object
        lmax    the new value of lmax

    Returns:
        A new tmatrix[2,rmax',2,rmax']
    """
    rmax_original = tmatrix.shape[1]
    lmax_original = miepy.vsh.rmax_to_lmax(rmax_original)

    if lmax > lmax_original:
        raise ValueError(f"cannot reduce tmatrix to lmax={lmax} since it is already smaller than this")

    rmax = miepy.vsh.lmax_to_rmax(lmax)

    return tmatrix[:, :rmax, :, :rmax]


def rotate_tmatrix(tmatrix, quat):
    """Rotate a T-matrix.

    Arguments:
        tmatrix    tmatrix to be rotated
        quat       quaternion representing the rotation

    Returns:
        The rotated T-matrix
    """
    rmax = tmatrix.shape[1]
    lmax = miepy.vsh.rmax_to_lmax(rmax)

    R = np.zeros([rmax, rmax], dtype=complex)

    for _i, n, _m in miepy.mode_indices(lmax):
        r = miepy.vsh.lmax_to_rmax(n)
        idx = np.s_[r - (2 * n + 1) : r]
        R[idx, idx] = miepy.vsh.vsh_rotation_matrix(n, quat)

    tmatrix_rot = np.einsum("ms,uw,asbw->ambu", R, np.conjugate(R), tmatrix)

    return tmatrix_rot
