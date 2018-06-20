import miepy

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
