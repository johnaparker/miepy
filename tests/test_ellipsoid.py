import miepy
import numpy as np

nm = 1e-9
rx = ry = 40*nm
rz = 70*nm
material = miepy.materials.Ag()
wavelength = 700*nm
lmax = 4

e = miepy.spheroid([0,0,0], rx, rz, material)
T1 = e.compute_tmatrix(lmax, wavelength, 1.0)
idx = np.abs(T1) > 1e-4
print(T1[idx])
e = miepy.ellipsoid([0,0,0], rx, ry, rz, material)
T2 = e.compute_tmatrix(lmax, wavelength, 1.0)
idx = np.abs(T2) > 1e-4
print(T2[idx])
