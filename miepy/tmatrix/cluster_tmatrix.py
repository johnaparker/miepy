import miepy
import numpy as np
from miepy.cpp.vsh_translation import vsh_translation_numpy as vsh_translation

def cluster_tmatrix(cluster):
    if isinstance(cluster, miepy.cluster):
        T = miepy.interactions.particle_aggregate_tmatrix(cluster.position, cluster.tmatrix,
                                  cluster.material_data.k_b)
    elif isinstance(cluster, miepy.sphere_cluster):
        T = miepy.interactions.sphere_aggregate_tmatrix(cluster.position, cluster.mie_scat,
                                  cluster.material_data.k_b)

    lmax = cluster.lmax
    rmax = miepy.vsh.lmax_to_rmax(lmax)
    N = cluster.Nparticles
    A1 = np.empty([N,2,rmax,2,rmax], dtype=complex)
    A2 = np.empty([N,2,rmax,2,rmax], dtype=complex)
    k = cluster.material_data.k_b

    dr = cluster.position
    rad, theta, phi = miepy.coordinates.cart_to_sph(dr[:,0], dr[:,1], dr[:,2])
    for r,n,m in miepy.mode_indices(lmax):
        for s,v,u in miepy.mode_indices(lmax):
            Emn = miepy.vsh.Emn(m, n)
            Euv = miepy.vsh.Emn(u, v)
            A_transfer, B_transfer = vsh_translation(m, n, u, v,
                    rad, theta, phi, k, miepy.vsh_mode.incident)

            A1[:,0,r,0,s] = 1/(Emn/Euv)*A_transfer
            A1[:,1,r,1,s] = 1/(Emn/Euv)*A_transfer
            A1[:,0,r,1,s] = 1/(Emn/Euv)*B_transfer
            A1[:,1,r,0,s] = 1/(Emn/Euv)*B_transfer

    for i in range(N):
        for j in range(N):
            for r,n,m in miepy.mode_indices(lmax):
                for s,v,u in miepy.mode_indices(lmax):
                    for a in range(2):
                        for b in range(2):
                            T[i,a,r,j,b,s] *= cluster.mie_scat[i,a,n-1]
                            # T[i,a,r,j,b,s] *= cluster.mie_scat[j,b,v-1]

    dr = -cluster.position
    rad, theta, phi = miepy.coordinates.cart_to_sph(dr[:,0], dr[:,1], dr[:,2])
    for r,n,m in miepy.mode_indices(lmax):
        for s,v,u in miepy.mode_indices(lmax):
            Emn = miepy.vsh.Emn(m, n)
            Euv = miepy.vsh.Emn(u, v)
            A_transfer, B_transfer = vsh_translation(m, n, u, v,
                    rad, theta, phi, k, miepy.vsh_mode.incident)

            A2[:,0,r,0,s] = 1/(Emn/Euv)*A_transfer
            A2[:,1,r,1,s] = 1/(Emn/Euv)*A_transfer
            A2[:,0,r,1,s] = 1/(Emn/Euv)*B_transfer
            A2[:,1,r,0,s] = 1/(Emn/Euv)*B_transfer

    np.einsum('abcabc->abc', T)[...] += 1
    T = np.linalg.tensorinv(T, 3)
    # T_cluster = np.einsum('Xabcd,XcdYef,Yefgh->abgh', A, T, A)
    # T_cluster = np.einsum('Xabcd,XcdYef,Yefgh->abgh', A1, T, A2)
    T_cluster = np.einsum('Xabcd,XcdYef,Yefgh->abgh', A1, T, A2)

    return T_cluster

class custom_particle(miepy.particles.particle):
    def __init__(self, position, T, orientation=None):
        super().__init__(position, orientation, miepy.materials.Ag())
        self.T = T

    def compute_tmatrix(self, lmax, wavelength, eps_m, **kwargs):
        self.tmatrix_fixed = miepy.tmatrix.tmatrix_reduce_lmax(self.T, lmax)
        self._rotate_fixed_tmatrix()

        return self.tmatrix

    def enclosed_radius(self):
        return 150*nm

    def _dict_key(self, wavelength):
        return (custom_particle,)

if __name__ == '__main__':
    nm = 1e-9
    Ag = miepy.materials.Ag()
    source = miepy.sources.plane_wave([1,1])
    L = 155*nm

    cluster = miepy.sphere_cluster(position=[[-L/2, 0,0], [L/2,0,0]],
                                   material=Ag,
                                   radius=75*nm,
                                   source=source,
                                   lmax=14,
                                   wavelength=800*nm)
    print(cluster.cross_sections())

    pos = cluster.position
    radii = cluster.radius
    lmax = np.empty_like(radii, dtype=int)
    lmax[...] = cluster.lmax
    lmax_cluster = 4
    wavelength = cluster.wavelength
    eps = cluster.material_data.eps
    eps_m = cluster.material_data.eps_b
    # T = miepy.tmatrix.tmatrix_sphere_cluster(pos, radii, lmax, lmax_cluster, wavelength, eps=eps, eps_m=eps_m)
    # T = cluster_tmatrix(cluster)
    # particle = custom_particle([0,0,0], T)
    particle = miepy.sphere_cluster_particle(cluster.position, cluster.radius, Ag, lmax=14)
    cluster = miepy.cluster(particles=particle,
                            source=source,
                            lmax=4,
                            wavelength=800*nm)
    print(cluster.cross_sections())
