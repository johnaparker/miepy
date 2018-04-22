"""
functions for building interaction matrices and solving them
"""

import numpy as np
import miepy

#TODO vectorize for loops. Avoid transpose of position->pass x,y,z to source instead...?
def sphere_cluster_t_matrix(positions, a, k):
    """Obtain the T-matrix for a cluster of spheres.
       Returns T[2,N,rmax,2,N,rmax]
    
       Arguments:
           positions[N,3]      particles positions
           a[N,2,Lmax]         mie scattering coefficients
           k                   medium wavenumber
    """
    Nparticles = positions.shape[0]
    Lmax = a.shape[-1]
    rmax = miepy.vsh.Lmax_to_rmax(Lmax)
    identity = np.zeros(shape = (Nparticles, 2, rmax, Nparticles, 2, rmax), dtype=complex)
    np.einsum('airair->air', identity)[...] = 1
    
    interaction_matrix = np.zeros(shape = (Nparticles, 2, rmax, Nparticles, 2, rmax), dtype=complex)

    for i in range(Nparticles):
        for j in range(Nparticles):
            if i == j: continue

            pi = positions[i]
            pj = positions[j]
            dji = pi -  pj
            r_ji = np.linalg.norm(dji)
            theta_ji = np.arccos(dji[2]/r_ji)
            phi_ji = np.arctan2(dji[1], dji[0])

            for r,(n,m) in enumerate(miepy.vsh.mode_indices(Lmax)):
                for s,(v,u) in enumerate(miepy.vsh.mode_indices(Lmax)):
                    A_transfer = miepy.vsh.A_translation(m,n,u,v,r_ji,theta_ji,phi_ji,
                            k, miepy.vsh.VSH_mode.outgoing)
                    B_transfer = miepy.vsh.B_translation(m,n,u,v,r_ji,theta_ji,phi_ji,
                            k, miepy.vsh.VSH_mode.outgoing)

                    interaction_matrix[i,0,r,j,0,s] = A_transfer*a[j,0,v-1]
                    interaction_matrix[i,0,r,j,1,s] = B_transfer*a[j,1,v-1]
                    interaction_matrix[i,1,r,j,0,s] = B_transfer*a[j,0,v-1]
                    interaction_matrix[i,1,r,j,1,s] = A_transfer*a[j,1,v-1]

    t_matrix = identity + interaction_matrix
    return t_matrix

def solve_sphere_cluster(positions, a, p_src, k):
    """Solve interactions of a collection of spheres.
       Returns [p_inc, q_inc]
    
       Arguments:
           positions[N,3]      particles positions
           a[N,2,Lmax]         mie scattering coefficients
           p_src[N,2,rmax]     source scattering coefficients
           k                   medium wavenumber
    """
    A = sphere_cluster_t_matrix(positions, a,  k)
    sol = np.linalg.tensorsolve(A, p_src)

    return sol
