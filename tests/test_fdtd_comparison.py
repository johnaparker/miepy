"""
Comparison to FDTD for the following simulation:

    gold nanoparticle dimer
    radius: 75 nm
    separation: 400 nm
    source: rhc plane wave
    medium: vacuum
    wavelengths: 400 - 1000 nm

Compares cross-sections and forces
"""

import numpy as np
import miepy

nm = 1e-9

### cross-sections from FDTD
scat = np.array([1.07751721e-14, 1.26679525e-14, 1.47235693e-14, 1.69336991e-14,
 1.93312588e-14, 2.20001815e-14, 2.50714233e-14, 2.87037141e-14,
 3.30512937e-14, 3.82243390e-14, 4.42501808e-14, 5.10440839e-14,
 5.83971075e-14, 6.59855440e-14, 7.34022363e-14, 8.02055770e-14,
 8.59781467e-14, 9.03846091e-14, 9.32182172e-14, 9.44272008e-14,
 9.41160899e-14, 9.25219015e-14, 8.99700925e-14, 8.68191764e-14,
 8.34050560e-14, 7.99959190e-14, 7.67659970e-14, 7.37921204e-14,
 7.10717753e-14, 6.85564654e-14, 6.61907578e-14, 6.39463002e-14,
 6.18416674e-14, 5.99428258e-14, 5.83443961e-14, 5.71374210e-14,
 5.63735369e-14, 5.60371068e-14, 5.60353184e-14, 5.62116381e-14])

absorb = np.array([-1.29930581e-16,  4.26142796e-16,  1.29811152e-15,  2.29521847e-15, 
  3.21331503e-15,  3.90907817e-15,  4.36107031e-15,  4.70134889e-15, 
  5.20877786e-15,  6.26448011e-15,  8.27876128e-15,  1.16053238e-14, 
  1.64613307e-14,  2.28704352e-14,  3.06407777e-14,  3.93824443e-14, 
  4.85607018e-14,  5.75742483e-14,  6.58431557e-14,  7.28899615e-14, 
  7.83995365e-14,  8.22482748e-14,  8.44996373e-14,  8.53697211e-14, 
  8.51719553e-14,  8.42532409e-14,  8.29343604e-14,  8.14653854e-14, 
  8.00027303e-14,  7.86094162e-14,  7.72751338e-14,  7.59488477e-14, 
  7.45746768e-14,  7.31219302e-14,  7.16023180e-14,  7.00709520e-14, 
  6.86119790e-14,  6.73136192e-14,  6.62401429e-14,  6.54092670e-14])

extinct = scat + absorb

### forces from FDTD
Fx = np.array([-4.07243862e-27, -3.91900859e-27, -3.61153493e-27, -3.12542938e-27,
 -2.42888016e-27, -1.47864254e-27, -2.23234359e-28,  1.39133044e-27,
  3.42151870e-27,  5.93599328e-27,  9.02917753e-27,  1.28237179e-26,
  1.74459514e-26,  2.29688894e-26,  2.93354587e-26,  3.62916914e-26,
  4.33651252e-26,  4.99124358e-26,  5.52350930e-26,  5.87330386e-26,
  6.00473138e-26,  5.91428340e-26,  5.63030681e-26,  5.20407523e-26,
  4.69589369e-26,  4.16117613e-26,  3.64087144e-26,  3.15837372e-26,
  2.72226358e-26,  2.33216744e-26,  1.98449691e-26,  1.67583103e-26,
  1.40349324e-26,  1.16443521e-26,  9.54127173e-27,  7.66647560e-27,
  5.96052308e-27,  4.38136450e-27,  2.91438605e-27,  1.56853445e-27])

Fy = np.array([ 3.51470526e-28,  3.44618449e-28,  3.47563103e-28,  3.94093803e-28,
  4.86803364e-28,  5.97580415e-28,  6.86839430e-28,  7.26996647e-28,
  7.14346542e-28,  6.62211672e-28,  5.80974422e-28,  4.59312742e-28,
  2.60480585e-28, -6.16644205e-29, -5.34355420e-28, -1.14079652e-27,
 -1.81291296e-27, -2.44762161e-27, -2.94010930e-27, -3.21899139e-27,
 -3.26762340e-27, -3.12318603e-27, -2.85615547e-27, -2.54147888e-27,
 -2.23496721e-27, -1.96370123e-27, -1.73093348e-27, -1.52889845e-27,
 -1.35065863e-27, -1.19503691e-27, -1.06425495e-27, -9.58420137e-28,
 -8.71979149e-28, -7.94728284e-28, -7.16179896e-28, -6.29768478e-28,
 -5.33921003e-28, -4.29657555e-28, -3.16968808e-28, -1.92834908e-28])

Fz = np.array([ 2.37269916e-26, 2.81872570e-26, 3.31290757e-26, 3.86162382e-26,
  4.48058406e-26, 5.19166367e-26, 6.01677554e-26, 6.97313326e-26,
  8.07393588e-26, 9.33536687e-26, 1.07863410e-25, 1.24741446e-25,
  1.44591854e-25, 1.67961108e-25, 1.95051053e-25, 2.25432340e-25,
  2.57881407e-25, 2.90435624e-25, 3.20687215e-25, 3.46246252e-25,
  3.65234952e-25, 3.76659539e-25, 3.80549425e-25, 3.77837694e-25,
  3.70045934e-25, 3.58894003e-25, 3.45960925e-25, 3.32480526e-25,
  3.19289117e-25, 3.06884205e-25, 2.95526515e-25, 2.85327587e-25,
  2.76298431e-25, 2.68368699e-25, 2.61401900e-25, 2.55225434e-25,
  2.49673969e-25, 2.44627455e-25, 2.40023333e-25, 2.35836991e-25])

### gold eps data in FDTD
eps = np.array([ -1.00220232+5.44762485j,  -1.35782477+5.24761921j,
                 -1.78271609+5.04627076j,  -2.27027164+4.84897369j,
                 -2.81446263+4.65913201j,  -3.4100077 +4.47873681j,
                 -4.05239683+4.30880799j,  -4.73783996+4.14971756j,
                 -5.46318344+4.00141902j,  -6.22581829+3.86360643j,
                 -7.02359246+3.73582256j,  -7.8547325 +3.6175312j ,
                 -8.71777635+3.50816443j,  -9.61151695+3.40715291j,
                -10.53495565+3.31394448j, -11.48726407+3.22801512j,
                -12.4677531 +3.14887463j, -13.47584767+3.07606916j,
                -14.51106626+3.0091814j , -15.57300441+2.94782964j,
                -16.66132118+2.8916659j , -17.77572823+2.84037382j,
                -18.91598088+2.79366621j, -20.08187082+2.75128268j,
                -21.2732201 +2.71298731j, -22.48987626+2.67856644j,
                -23.73170824+2.64782657j, -24.99860304+2.62059251j,
                -26.29046293+2.5967056j , -27.60720313+2.57602214j,
                -28.94874993+2.55841196j, -30.31503899+2.54375714j,
                -31.70601405+2.53195081j, -33.12162576+2.52289615j,
                -34.56183069+2.51650543j, -36.02659057+2.51269917j,
                -37.51587153+2.51140537j, -39.02964358+2.51255884j,
                -40.56788003+2.51610058j, -42.13055709+2.52197726j,])

wavelengths = np.linspace(400*nm, 1000*nm, 40)
frequency = np.linspace(1, 1/0.4, 40)

separation = 400*nm
radius = 75*nm
Au = miepy.data_material(wavelengths, eps)
source = miepy.sources.plane_wave.from_string(polarization='rhc')
gmtF = np.zeros((3,) + wavelengths.shape)
gmtC = np.zeros((3,) + wavelengths.shape)

for i,wavelength in enumerate(wavelengths):
    sol = miepy.sphere_cluster(position=[[-separation/2,0,0],[separation/2,0,0]],
                               radius=radius,
                               material=Au,
                               source=source,
                               wavelength=wavelength,
                               lmax=2)

    gmtF[:,i] = sol.force_on_particle(1) 
    gmtC[:,i] = sol.cross_sections()

def test_fdtd_cross_sections(plot=False):
    """compare cross-sections to fdtd"""
    if not plot:
        L2 = np.linalg.norm(scat - gmtC[0])/scat.shape[0]
        assert np.all(L2 < 8e-15)

        L2 = np.linalg.norm(absorb - gmtC[1])/absorb.shape[0]
        assert np.all(L2 < 1.1e-14)

        L2 = np.linalg.norm(extinct - gmtC[2])/extinct.shape[0]
        assert np.all(L2 < 1.9e-14)
    else:
        fig, ax = plt.subplots()
        ax.plot(1000/frequency, scat, 'o', color='C0', label='scattering (FDTD)')
        ax.plot(1000/frequency, absorb, 'o', color='C1', label='absorbption (FDTD)')
        ax.plot(1000/frequency, extinct, 'o', color='C2', label='extinction (FDTD)')

        ax.plot(wavelengths/nm, gmtC[0], color='C0', label='scattering (FDTD)')
        ax.plot(wavelengths/nm, gmtC[1], color='C1', label='absorbption (FDTD)')
        ax.plot(wavelengths/nm, gmtC[2], color='C2', label='extinction (FDTD)')

        ax.legend()
        ax.set(ylabel='cross-section', xlabel='wavelength (nm)')
        fig.suptitle(test_fdtd_cross_sections.__name__, weight='bold')

def test_fdtd_forces(plot=False):
    """compare forces to fdtd"""
    if not plot:
        L2 = np.linalg.norm(Fx - gmtF[0])/Fx.shape[0]
        assert np.all(L2 < 5e-27)

        L2 = np.linalg.norm(Fy - gmtF[1])/Fy.shape[0]
        assert np.all(L2 < 3.4e-28)

        L2 = np.linalg.norm(Fz - gmtF[2])/Fz.shape[0]
        assert np.all(L2 < 3.7e-26)
    else:
        fig, axes = plt.subplots(nrows=2, figsize=(7,6), sharex=True,
                      gridspec_kw=dict(height_ratios=[2,1], hspace=0.05))

        for ax in axes:
            ax.plot(1000/frequency, Fx, 'o', color='C0', label='Fx (FDTD)')
            ax.plot(1000/frequency, Fy, 'o', color='C1', label='Fy (FDTD)')
            ax.plot(1000/frequency, Fz, 'o', color='C2', label='Fz (FDTD)')

        for ax in axes:
            ax.axhline(0, linestyle='--', color='black')
            ax.plot(wavelengths/nm, gmtF[0], color='C0', label='Fx (GMT)')
            ax.plot(wavelengths/nm, gmtF[1], color='C1', label='Fy (GMT)')
            ax.plot(wavelengths/nm, gmtF[2], color='C2', label='Fz (GMT)')

        axes[0].legend()
        axes[0].set(ylabel='force')
        axes[1].set(xlabel='wavelength (nm)', ylabel='force', ylim=[-0.035e-25,0.01e-25])
        fig.suptitle(test_fdtd_forces.__name__, weight='bold')

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    test_fdtd_cross_sections(plot=True)
    test_fdtd_forces(plot=True)
    plt.show()
