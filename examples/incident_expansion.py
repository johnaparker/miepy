"""Decompose an x-polarized plane wave into the VSHs"""
import numpy as np
import matplotlib.pyplot as plt
import miepy

### source definition
source = miepy.sources.x_polarized_plane_wave()
k = 2*np.pi/1
Nmax = 10

### grid plane
x = np.linspace(-.3,.3,30)
y = np.linspace(-.3,.3,30)
z = 0 
X,Y,Z = np.meshgrid(x,y,z, indexing='ij')
coords = np.array([X,Y,Z])
E = source.E(coords, k).squeeze()

### exact fields
def exact():
    I = np.sum(np.abs(E)**2, axis=0)

    fig,ax = plt.subplots()
    im = ax.pcolormesh(X.squeeze(),Y.squeeze(),I, vmin=0.9, vmax=1.1)
    plt.colorbar(im)
    ax.quiver(X.squeeze(),Y.squeeze(),E[0].real, E[1].real)

### analytic expansion
def analytic():
    expanded_E = np.zeros_like(E, dtype=complex)

    R, THETA, PHI = miepy.vsh.cart_to_sph(X,Y,Z)
    rhat,that,phat = miepy.vsh.sph_basis_vectors(THETA,PHI)

    for n in range(1,Nmax+1):
        for m in range(-n, n+1):
            Nfunc,Mfunc = miepy.vsh.VSH(n,m, mode=miepy.vsh.VSH_mode.incident)

            Emn = miepy.vsh.Emn(m,n,source.amplitude)
            p,q = source.structure(n,m,[0,0,0],k)

            # phase = np.exp(1j*k*Z)
            # if m == 1:
                # p,q = (phase*1/2, phase*1/2)
            # elif m == -1:
                # p,q = (-phase*1/(2*n*(n+1)), phase*1/(2*n*(n+1)))
            # else:
                # p,q = (np.zeros_like(coords[0]), np.zeros_like(coords[0]))

            N = Nfunc(R,THETA,PHI,k)
            M = Mfunc(R,THETA,PHI,k)

            N = rhat*N[0] + that*N[1] + phat*N[2]
            M = rhat*M[0] + that*M[1] + phat*M[2]

            expanded_E += -1j*Emn*(p*N + q*M).squeeze()

    fig,ax = plt.subplots()
    I = np.sum(np.abs(expanded_E)**2, axis=0)
    im = ax.pcolormesh(X.squeeze(),Y.squeeze(),I)
    plt.colorbar(im)
    ax.quiver(X.squeeze(),Y.squeeze(),expanded_E[0].real, expanded_E[1].real)

### numerical decomposition and expansion
p,q = miepy.vsh.decompose_source(source, k, Nmax, sampling=60)
f =  miepy.vsh.expand(p, q, k, miepy.vsh.VSH_mode.incident)
expanded_E = f(X,Y,Z).squeeze()

fig,ax = plt.subplots()
I = np.sum(np.abs(expanded_E)**2, axis=0)
im = ax.pcolormesh(X.squeeze(),Y.squeeze(),I)
plt.colorbar(im)
ax.quiver(X.squeeze(),Y.squeeze(),expanded_E[0].real, expanded_E[1].real)


plt.show()
