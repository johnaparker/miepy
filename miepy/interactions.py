"""
functions for building interaction matrices and solving them
"""

import numpy as np
import miepy
from scipy.sparse.linalg import bicgstab
from miepy.cpp.vsh_translation import vsh_translation_numpy as vsh_translation

def solve_linear_system(tmatrix, p_src, method):
    """Solve the linear system p_inc = p_src - tmatrix*p_inc
        size
       Arguments:
           tmatrix[N,2,rmax,N,2,rmax]   particle aggregate tmatrix
           p_src[N,2,rmax]   source scattering coefficients
           method    solver method ('exact', 'bicgstab')
    """
    Nparticles = tmatrix.shape[0]
    rmax = p_src.shape[-1]
    size = Nparticles*2*rmax

    return miepy.cpp.interactions.solve_linear_system(tmatrix.reshape(size, size), p_src.reshape(-1), method=method).reshape([Nparticles,2,rmax])

def interactions_precomputation(positions, k, lmax):
    """Get the relative r,theta,phi positions of the particles and precomputed zn function

       Arguments:
           positions[N,3]      particles positions
           k      medium wavenumber
           lmax   maximum number of multipoles
    
       Returns:
           r[idx], theta[idx], phi[idx], zn[nmax,idx] (idx enumerates i,j and j > i)
    """
    Nparticles = positions.shape[0]
    size = int(Nparticles*(Nparticles-1)/2)
    r_ji     = np.zeros(size, dtype=float)
    theta_ji = np.zeros(size, dtype=float)
    phi_ji   = np.zeros(size, dtype=float)

    idx = 0
    for i in range(Nparticles):
        for j in range(i+1, Nparticles):
            pi = positions[i]
            pj = positions[j]
            dji = pi -  pj
            r_ji[idx] = np.linalg.norm(dji)
            theta_ji[idx] = np.arccos(dji[2]/r_ji[idx])
            phi_ji[idx] = np.arctan2(dji[1], dji[0])

            idx +=1

    nmax = 2*lmax + 1
    zn_values = np.zeros((nmax, size), dtype=complex)
    for n in range(nmax):
        zn_values[n] = miepy.vsh.special.spherical_hn(n, k*r_ji)

    return r_ji, theta_ji, phi_ji, zn_values

def sphere_aggregate_tmatrix(positions, mie, k):
    """Obtain the particle-centered aggregate T-matrix for a cluster of spheres
       Returns T[N,2,rmax,N,2,rmax]
    
       Arguments:
           positions[N,3]      particles positions
           mie[N,2,lmax]       mie scattering coefficients
           k                   medium wavenumber
    """
    Nparticles = positions.shape[0]
    lmax = mie.shape[-1]
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    return miepy.cpp.interactions.sphere_aggregate_tmatrix(positions, mie.reshape([Nparticles,-1]), k).reshape([Nparticles,2,rmax,Nparticles,2,rmax])


def sphere_aggregate_tmatrix_periodic(positions, mie, k, symmetry, k_hat):
    """Obtain the particle-centered aggregate T-matrix for a cluster of spheres with a given periodic symmetry
       Returns T[N,2,rmax,N,2,rmax]
    
       Arguments:
           positions[N,3]      particles positions
           mie[N,2,lmax]       mie scattering coefficients
           k                   medium wavenumber
           symmetry            type of symmetry
           k_hat               unit k-vector
    """

    Nparticles = positions.shape[0]
    lmax = mie.shape[-1]
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    # agg_tmatrix = np.zeros(shape=(Nparticles, 2, rmax, Nparticles, 2, rmax), dtype=complex)
    agg_tmatrix = sphere_aggregate_tmatrix(positions, mie, k)

    xpos, ypos, zpos = symmetry.generate(5000) 
    cells = np.array([xpos,ypos,zpos]).T


    for i in range(Nparticles):
        for j in range(Nparticles):
            dr = positions[i] - (cells + positions[j])
            rad, theta, phi = miepy.coordinates.cart_to_sph(dr[:,0], dr[:,1], dr[:,2])
            for r,n,m in miepy.mode_indices(lmax):
                for s,v,u in miepy.mode_indices(lmax):
                    phase_factor = np.exp(-1j*k*(k_hat[0]*xpos + k_hat[1]*ypos + k_hat[2]*zpos))
                    A_transfer, B_transfer = vsh_translation(m, n, u, v, 
                            rad, theta, phi, k, miepy.vsh_mode.outgoing)
                    for a in range(2):
                        for b in range(2):
                            val = (A_transfer, B_transfer)[(a+b)%2]*phase_factor
                            agg_tmatrix[i,a,r,j,b,s] += np.sum(val)*mie[j,b,v-1]

    return agg_tmatrix

#TODO this function is more general than above and can be used for both cases (change only the einsum)
def particle_aggregate_tmatrix(positions, tmatrix, k):
    """Obtain the particle-centered aggregate T-matrix for a cluster of particles
       Returns T[N,2,rmax,N,2,rmax]
    
       Arguments:
           positions[N,3]      particles positions
           tmatrix[N,2,rmax,2,rmax]   single particle T-matrices
           k                   medium wavenumber
    """

    Nparticles = positions.shape[0]
    rmax = tmatrix.shape[-1]
    lmax = miepy.vsh.rmax_to_lmax(rmax)

    return miepy.cpp.interactions.particle_aggregate_tmatrix(positions, tmatrix.reshape([Nparticles,-1]), k).reshape([Nparticles,2,rmax,Nparticles,2,rmax])

def reflection_matrix_nia(positions, mie, k, reflected, z):
    """Obtain the particle-centered aggregate T-matrix for a cluster of spheres
       Returns T[N,2,rmax,N,2,rmax]
    
       Arguments:
           positions[N,3]      particles positions
           mie[N,2,lmax]       mie scattering coefficients
           k                   medium wavenumber
           z                   z coordinate of interface
    """
    Nparticles = positions.shape[0]
    lmax = mie.shape[-1]
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    R = miepy.cpp.interactions.reflection_matrix_nia(positions, mie.reshape([Nparticles,-1]), k, reflected, z).reshape([Nparticles,2,rmax,Nparticles,2,rmax])
    return R

