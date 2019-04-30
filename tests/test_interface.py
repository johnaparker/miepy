"""
Tests related to miepy.interface
"""

import numpy as np
import miepy

nm = 1e-9
wavelength = 600*nm
k = 2*np.pi/wavelength
radius = 75*nm
medium = miepy.materials.water()
material = miepy.materials.Ag()

width = 200*nm
polarization = [1,0]

zpos = 400*nm


def test_interface_z_translation():
    """
    Moving the source and particle is identical to moving the interface (cross-section comparison)
    """
    source = miepy.sources.gaussian_beam(width=width, polarization=polarization, center=[0,0,-zpos])
    interface = miepy.interface(miepy.materials.constant_material(index=1.7))
    cluster = miepy.sphere_cluster(position=[0,0,-zpos],
                                   radius=radius,
                                   material=material,
                                   medium=medium,
                                   lmax=2,
                                   source=source,
                                   interface=interface,
                                   wavelength=wavelength)
    C1 = np.array(cluster.cross_sections())


    source = miepy.sources.gaussian_beam(width=width, polarization=polarization)
    interface = miepy.interface(miepy.materials.constant_material(index=1.7), z=zpos)
    cluster = miepy.sphere_cluster(position=[0,0,0],
                                   radius=radius,
                                   material=material,
                                   medium=medium,
                                   lmax=2,
                                   source=source,
                                   interface=interface,
                                   wavelength=wavelength)
    C2 = np.array(cluster.cross_sections())

    assert np.allclose(C1, C2, atol=0, rtol=0)

def test_index_matched_interface():
    """
    An interface that is index-matched with the medium is identical to not having an interface (cross-section comparison)
    """
    source = miepy.sources.gaussian_beam(width=width, polarization=polarization)
    interface = miepy.interface(medium, z=zpos)
    cluster = miepy.sphere_cluster(position=[0,0,0],
                                   radius=radius,
                                   material=material,
                                   medium=medium,
                                   lmax=2,
                                   source=source,
                                   interface=interface,
                                   wavelength=wavelength)
    C1 = np.array(cluster.cross_sections())


    cluster = miepy.sphere_cluster(position=[0,0,0],
                                   radius=radius,
                                   material=material,
                                   medium=medium,
                                   lmax=2,
                                   source=source,
                                   wavelength=wavelength)
    C2 = np.array(cluster.cross_sections())

    assert np.allclose(C1, C2, atol=0, rtol=1e-15)
