import miepy
import numpy as np
import quaternion

def tmatrix_reduce_lmax(tmatrix, lmax):
    """Reduce the dimensions of a tmatrix to lmax
    
    Arguments:
        tmatrix[2,rmax,2,rmax]    the tmatrix object
        lmax    the new value of lmax

    Returns:
        A new tmatrix[2,rmax',2,rmax']
    """

    rmax_original = tmatrix.shape[1]
    lmax_original = miepy.vsh.rmax_to_lmax(rmax_original)

    if lmax > lmax_original:
        raise ValueError(f'cannot reduce tmatrix to lmax={lmax} since it is already smaller than this') 

    rmax = miepy.vsh.lmax_to_rmax(lmax)

    return tmatrix[:, :rmax, :, :rmax]

def rotate_tmatrix(tmatrix, quat):
    """Rotate a T-matrix
    
    Arguments:
        tmatrix    tmatrix to be rotated
        quat       quaternion representing the rotation

    Returns:
        The rotated T-matrix
    """
    alpha, beta, gamma = -quaternion.as_euler_angles(quat)

    rmax = tmatrix.shape[1]
    lmax = miepy.vsh.rmax_to_lmax(rmax)

    R = np.zeros([rmax, rmax], dtype=complex)

    for i,n,m in miepy.mode_indices(lmax):
        #TODO: remove this for loop, fill in R block diagonal style
        for j,v,u in miepy.mode_indices(lmax):
            if n == v:
                R[i,j] = miepy.vsh.vsh_rotation.wigner_D(n, m, u, alpha, beta, gamma) 

    #TODO: which should be correct?
    # tmatrix_rot = np.einsum('sm,wu,asbw->ambu', R, np.conjugate(R), tmatrix)
    tmatrix_rot = np.einsum('sm,wu,ambu->asbw', R, np.conjugate(R), tmatrix)

    return tmatrix_rot
