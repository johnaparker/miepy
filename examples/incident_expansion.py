import numpy as np
import matplotlib.pyplot as plt
import miepy
import my_pytools.my_numpy.special as special

source = miepy.sources.x_polarized_plane_wave()

x = np.linspace(-.3,.3,30)
y = np.linspace(-.3,.3,30)
z = 1.5
X,Y,Z = np.meshgrid(x,y,z, indexing='ij')
k = 2*np.pi/1
coords = np.array([X,Y,Z])

E = source.E(coords, k).squeeze()
I = np.sum(np.abs(E)**2, axis=0)

fig,ax = plt.subplots()
im = ax.pcolormesh(X.squeeze(),Y.squeeze(),I, vmin=0.9, vmax=1.1)
plt.colorbar(im)
ax.quiver(X.squeeze(),Y.squeeze(),E[0].real, E[1].real)

expanded_E = np.zeros_like(E, dtype=complex)
Nmax = 10

R = (X**2 + Y**2 + Z**2)**0.5
THETA = np.arccos(Z/R)
PHI = np.arctan2(Y,X)

rhat = np.array([np.sin(THETA)*np.cos(PHI), np.sin(THETA)*np.sin(PHI), np.cos(THETA)])
that = np.array([np.cos(THETA)*np.cos(PHI), np.cos(THETA)*np.sin(PHI), -np.sin(THETA)])
phat = np.array([-np.sin(PHI), np.cos(PHI), np.zeros_like(THETA)])

for n in range(1,Nmax+1):
    for m in range(-n, n+1):
        Nfunc,Mfunc = special.VSH(n,m, mode=special.VSH_mode.incident)

        Emn = special.Emn(m,n,source.amplitude)
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
# im = ax.pcolormesh(X.squeeze(),Y.squeeze(),I, vmin=0.9, vmax=1.1)
im = ax.pcolormesh(X.squeeze(),Y.squeeze(),I, vmin=0.9, vmax=1.1)
plt.colorbar(im)
ax.quiver(X.squeeze(),Y.squeeze(),expanded_E[0].real, expanded_E[1].real)

plt.show()
