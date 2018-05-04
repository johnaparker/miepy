"""
Functions related to the (m,n) vsh mode indices
"""

def rmax_to_Lmax(rmax):
    """obtain Lmax from rmax"""
    Lmax = int(-1 + (1+rmax)**0.5)
    return Lmax

def Lmax_to_rmax(Lmax):
    """obtain rmax from Lmax"""
    rmax = Lmax*(Lmax + 2)
    return rmax

def mode_indices(Lmax):
    """generate (n,m) index pairs up to n=Lmax"""
    n = 1
    m = -1

    for n in range(1,Lmax+1):
        for m in range(-n,n+1):
            yield n, m
