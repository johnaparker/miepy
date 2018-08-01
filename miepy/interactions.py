"""
functions for building interaction matrices and solving them
"""

import numpy as np
import miepy
from scipy.sparse.linalg import bicgstab

def solve_linear_system(tmatrix, p_src, method):
    """Solve the linear system p_inc = p_src - tmatrix*p_inc
        
       Arguments:
           tmatrix[N,2,rmax,N,2,rmax]   particle aggregate tmatrix
           p_src[N,2,rmax]   source scattering coefficients
           method    solver method ('exact', 'bicgstab')
    """
    interaction_matrix = np.copy(tmatrix)
    np.einsum('airair->air', interaction_matrix)[...] += 1

    if method == 'exact':
        return np.linalg.tensorsolve(interaction_matrix, p_src)

    elif method == 'bicgstab':
        N = p_src.shape[0]
        rmax = p_src.shape[-1]

        sol = bicgstab(interaction_matrix.reshape((2*N*rmax, -1)), p_src.reshape((-1,)))[0]
        return sol.reshape((N,2,rmax))

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

def sphere_aggregate_tmatrix(positions, a, k):
    """Obtain the particle-centered aggregate T-matrix for a cluster of spheres
       Returns T[N,2,rmax,N,2,rmax]
    
       Arguments:
           positions[N,3]      particles positions
           a[N,2,lmax]         mie scattering coefficients
           k                   medium wavenumber
    """
    Nparticles = positions.shape[0]
    lmax = a.shape[-1]
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    agg_tmatrix = np.zeros(shape=(Nparticles, 2, rmax, Nparticles, 2, rmax), dtype=complex)

    if Nparticles == 1:
        return agg_tmatrix
    
    r_ji, theta_ji, phi_ji, zn_values = interactions_precomputation(positions, k, lmax)

    for r,n,m in miepy.mode_indices(lmax):
        for s,v,u in miepy.mode_indices(lmax):
            # if s - 2*u < r: continue
            A_transfer, B_transfer = miepy.vsh.vsh_translation(m, n, u, v, 
                    r_ji, theta_ji, phi_ji, k, miepy.vsh_mode.outgoing, zn_values=zn_values)

            upper_idx = np.triu_indices(Nparticles, 1)
            lower_idx = upper_idx[::-1]

            agg_tmatrix[:,0,r,:,0,s][upper_idx] = A_transfer
            agg_tmatrix[:,0,r,:,1,s][upper_idx] = B_transfer
            agg_tmatrix[:,1,r,:,0,s][upper_idx] = B_transfer
            agg_tmatrix[:,1,r,:,1,s][upper_idx] = A_transfer

            agg_tmatrix[:,0,r,:,0,s][lower_idx] = (-1)**(n+v)  *A_transfer
            agg_tmatrix[:,0,r,:,1,s][lower_idx] = (-1)**(n+v+1)*B_transfer
            agg_tmatrix[:,1,r,:,0,s][lower_idx] = (-1)**(n+v+1)*B_transfer
            agg_tmatrix[:,1,r,:,1,s][lower_idx] = (-1)**(n+v)  *A_transfer

            # agg_tmatrix[:,0,s-2*u,:,0,r-2*m][upper_idx] = (-1)**(m+u)  *A_transfer
            # agg_tmatrix[:,0,s-2*u,:,1,r-2*m][upper_idx] = (-1)**(m+u+1)*B_transfer
            # agg_tmatrix[:,1,s-2*u,:,0,r-2*m][upper_idx] = (-1)**(m+u+1)*B_transfer
            # agg_tmatrix[:,1,s-2*u,:,1,r-2*m][upper_idx] = (-1)**(m+u)  *A_transfer

            # agg_tmatrix[:,0,s-2*u,:,0,r-2*m][lower_idx] = (-1)**(m+u+n+v)*A_transfer
            # agg_tmatrix[:,0,s-2*u,:,1,r-2*m][lower_idx] = (-1)**(m+u+n+v)*B_transfer
            # agg_tmatrix[:,1,s-2*u,:,0,r-2*m][lower_idx] = (-1)**(m+u+n+v)*B_transfer
            # agg_tmatrix[:,1,s-2*u,:,1,r-2*m][lower_idx] = (-1)**(m+u+n+v)*A_transfer

            agg_tmatrix[:,:,r,:,:,s] *= a[:,:,v-1]
            # agg_tmatrix[:,:,s-2*u,:,:,r-2*m] *= a[:,:,n-1]

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
    agg_tmatrix = np.zeros(shape=(Nparticles, 2, rmax, Nparticles, 2, rmax), dtype=complex)

    if Nparticles == 1:
        return agg_tmatrix
    
    r_ji, theta_ji, phi_ji, zn_values = interactions_precomputation(positions, k, lmax)

    for r,n,m in miepy.mode_indices(lmax):
        for s,v,u in miepy.mode_indices(lmax):
            # if s - 2*u < r: continue
            A_transfer, B_transfer = miepy.vsh.vsh_translation(m, n, u, v, 
                    r_ji, theta_ji, phi_ji, k, miepy.vsh_mode.outgoing, zn_values=zn_values)

            upper_idx = np.triu_indices(Nparticles, 1)
            lower_idx = upper_idx[::-1]

            agg_tmatrix[:,0,r,:,0,s][upper_idx] = A_transfer
            agg_tmatrix[:,0,r,:,1,s][upper_idx] = B_transfer
            agg_tmatrix[:,1,r,:,0,s][upper_idx] = B_transfer
            agg_tmatrix[:,1,r,:,1,s][upper_idx] = A_transfer

            agg_tmatrix[:,0,r,:,0,s][lower_idx] = (-1)**(n+v)  *A_transfer
            agg_tmatrix[:,0,r,:,1,s][lower_idx] = (-1)**(n+v+1)*B_transfer
            agg_tmatrix[:,1,r,:,0,s][lower_idx] = (-1)**(n+v+1)*B_transfer
            agg_tmatrix[:,1,r,:,1,s][lower_idx] = (-1)**(n+v)  *A_transfer

            # agg_tmatrix[:,0,s-2*u,:,0,r-2*m][upper_idx] = (-1)**(m+u)  *A_transfer
            # agg_tmatrix[:,0,s-2*u,:,1,r-2*m][upper_idx] = (-1)**(m+u+1)*B_transfer
            # agg_tmatrix[:,1,s-2*u,:,0,r-2*m][upper_idx] = (-1)**(m+u+1)*B_transfer
            # agg_tmatrix[:,1,s-2*u,:,1,r-2*m][upper_idx] = (-1)**(m+u)  *A_transfer

            # agg_tmatrix[:,0,s-2*u,:,0,r-2*m][lower_idx] = (-1)**(m+u+n+v)*A_transfer
            # agg_tmatrix[:,0,s-2*u,:,1,r-2*m][lower_idx] = (-1)**(m+u+n+v)*B_transfer
            # agg_tmatrix[:,1,s-2*u,:,0,r-2*m][lower_idx] = (-1)**(m+u+n+v)*B_transfer
            # agg_tmatrix[:,1,s-2*u,:,1,r-2*m][lower_idx] = (-1)**(m+u+n+v)*A_transfer

    agg_tmatrix = np.einsum('iabjcd,jcdef->iabjef', agg_tmatrix, tmatrix)
    return agg_tmatrix
