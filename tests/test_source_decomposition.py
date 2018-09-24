import miepy
import numpy as np

def test_plane_wave_point_matching():
    """point matching a plane wave agrees with analytic results"""

    k = 1
    wav = 2*np.pi/k
    lmax = 3

    source = miepy.sources.plane_wave([1,0])
    p_src = source.structure([0,0,0], k, lmax)

    p_src_numeric = miepy.vsh.decomposition.near_field_point_matching(source, [0,0,0], wav/4, k, lmax, sampling=10)

    print(p_src)
    print(p_src_numeric)
     
    L2 = np.linalg.norm(p_src - p_src_numeric)/p_src.size
    avg = np.average(np.abs(p_src) + np.abs(p_src_numeric))/2
    print(L2/avg)

    assert L2 < 7e-3*avg

def test_paraxial_gaussian_beam_approx_plane_wave():
    """a wide gaussian beam is approximately a plane wave at the center"""

    k = 1
    wav = 2*np.pi/k
    lmax = 3

    source = miepy.sources.plane_wave([1,0])
    p_plane = source.structure([0,0,0], k, lmax)

    source = miepy.sources.gaussian_beam(10*wav, [1,0], amplitude=1)
    p_gauss = miepy.vsh.decomposition.near_field_point_matching(source, [0,0,0], wav/4, k, lmax, sampling=10)

    print(p_plane)
    print(p_gauss)
     
    L2 = np.linalg.norm(p_plane - p_gauss)/p_plane.size
    avg = np.average(np.abs(p_plane) + np.abs(p_gauss))/2
    print(L2/avg)

    assert L2 < 6e-3*avg
