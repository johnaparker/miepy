"""
VSH rotation functions
"""

import numpy as np
import miepy

def wigner_D(n, m, mp, quat):
    """Wigner-D function
    
    Arguments:
        n       multipole order
        m       multipole orientation (from)
        mp      multipole orientation (to)
        quat    quaternion representing the rotation
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from spherical_functions import Wigner_D_element

    return Wigner_D_element(quat, n, mp, m)

def vsh_rotation_matrix(n, quat):
    """Rotation matrix for a given multipole order

    Arguments:
        n       multipole order
        quat    quaternion representing the rotation

    Returns:
        Rotation matrix R[2n+1,2n+1], such that p' = R*p
    """
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        from spherical_functions import Wigner_D_matrices

    l = 2*n + 1
    R = Wigner_D_matrices(quat, n, n).reshape((l,l))

    m = np.arange(-n, n+1)
    R *= np.power(-1.0, np.subtract.outer(m, m))

    return np.conj(R)

def rotate_expansion_coefficients(p_exp, quat):
    """Rotate a set of expansion coefficients to a new reference frame

    Arguments:
        p_exp[2,rmax]   expansion coefficients
        quat            quaternion representing the rotation

    Returns:
        The rotated expansion coefficients, p_rot[2,rmax]
    """
    p_rot = np.empty_like(p_exp)
    rmax = p_exp.shape[-1]
    lmax = miepy.vsh.rmax_to_lmax(rmax)

    for n in range(1,lmax+1):
        R = vsh_rotation_matrix(n, quat)
        rmax = miepy.vsh.lmax_to_rmax(n)
        idx = np.s_[rmax-(2*n+1):rmax]

        if p_rot.ndim == 3:
            p_rot[...,idx] = np.einsum('ab,Npb->Npa', R, p_exp[...,idx])
        else:
            p_rot[:,idx] = np.einsum('ab,pb->pa', R, p_exp[:,idx])

    return p_rot
