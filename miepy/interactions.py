"""
functions for building interaction matrices and solving them
"""

import numpy as np
import miepy

#TODO vectorize for loops. Avoid transpose of position->pass x,y,z to source instead...?
def sphere_cluster_t_matrix(positions, a, b, k):
    """Obtain the T-matrix for a cluster of spheres.
       Returns T[2,N,rmax,2,N,rmax]
    
       Arguments:
           positions[N,3]      particles positions
           a[N,Lmax]           an mie coefficients
           b[N,Lmax]           bn mie coefficients
           k                   medium wavenumber
    """
    Nparticles = positions.shape[0]
    Lmax = a.shape[1]
    rmax = Lmax*(Lmax + 2)
    identity = np.zeros(shape = (2, Nparticles, rmax, 2, Nparticles, rmax), dtype=complex)
    np.einsum('airair->air', identity)[...] = 1
    
    interaction_matrix = np.zeros(shape = (2, Nparticles, rmax, 2, Nparticles, rmax), dtype=complex)

    for i in range(Nparticles):
        for j in range(Nparticles):
            if i == j: continue

            pi = positions[i]
            pj = positions[j]
            dji = pi -  pj
            r_ji = np.linalg.norm(dji)
            theta_ji = np.arccos(dji[2]/r_ji)
            phi_ji = np.arctan2(dji[1], dji[0])

            for n in range(1, Lmax+1):
                for m in range(-n, n+1):
                    r = n**2 + n - 1 + m
                    for v in range(1, Lmax+1):
                        for u in range(-v, v+1):
                            s = v**2 + v - 1 + u

                            A_transfer = miepy.vsh.A_translation(m,n,u,v,r_ji,theta_ji,phi_ji,
                                    k, miepy.vsh.VSH_mode.outgoing)
                            B_transfer = miepy.vsh.B_translation(m,n,u,v,r_ji,theta_ji,phi_ji,
                                    k, miepy.vsh.VSH_mode.outgoing)

                            interaction_matrix[0,i,r,0,j,s] = A_transfer*a[j,v-1]
                            interaction_matrix[0,i,r,1,j,s] = B_transfer*b[j,v-1]
                            interaction_matrix[1,i,r,0,j,s] = B_transfer*a[j,v-1]
                            interaction_matrix[1,i,r,1,j,s] = A_transfer*b[j,v-1]

    t_matrix = identity + interaction_matrix
    return t_matrix

def solve_sphere_cluster(positions, a, b, p_src, q_src, k):
    """Solve interactions of a collection of spheres.
       Returns [p_inc, q_inc]
    
       Arguments:
           positions[N,3]      particles positions
           a[N,Lmax]           an mie coefficients
           b[N,Lmax]           bn mie coefficients
           p_src[N,rmax]       p coefficients of the source
           q_src[N,rmax]       q coefficients of the source
           k                   medium wavenumber
    """

    A = sphere_cluster_t_matrix(positions, a, b, k)
    b = np.array([p_src, q_src])
    sol = np.linalg.tensorsolve(A, b)

    return sol

def cluster_coefficients():
    pass
