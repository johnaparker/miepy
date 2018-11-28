import miepy
from collections import namedtuple
import quaternion
import numpy as np

def get_texture(material):
    """Return a (texture, color) tuple for a given material"""
    import vpython
    paint = namedtuple('paint', ['texture', 'color'])

    if material.name == 'Ag':
        return paint(vpython.textures.metal, vpython.vec(*[.9]*3))
    if material.name == 'Au':
        return paint(vpython.textures.metal, vpython.vec(1, .6, .3))

    return paint(None, vpython.vec(*[1]*3))

def draw_particle(particle, paint, **kwargs):
    import vpython
    vec = vpython.vec

    if type(particle) is miepy.sphere:
        return vpython.sphere(pos=vec(*particle.position),
                              radius=particle.radius,
                              texture=paint.texture,
                              color=paint.color, **kwargs)

    else:
        q = particle.orientation
        theta, phi = quaternion.as_spherical_coords(q)
        axis = vec(np.sin(theta)*np.cos(phi),
                   np.sin(theta)*np.sin(phi),
                   np.cos(theta))

        if type(particle) is miepy.cylinder:
            return vpython.cylinder(pos=vec(*particle.position) - particle.height/2*axis,
                                    radius=particle.radius,
                                    length=particle.height,
                                    texture=paint.texture,
                                    axis=axis,
                                    color=paint.color, **kwargs)

        elif type(particle) is miepy.spheroid:
            return vpython.ellipsoid(pos=vec(*particle.position),
                                    width=2*particle.axis_xy,
                                    height=2*particle.axis_xy,
                                    length=2*particle.axis_z,
                                    texture=paint.texture,
                                    axis=axis,
                                    color=paint.color, **kwargs)

#TODO: animation: pan camera around
#TODO: scale bar: auto-detect position and size
#TODO: center: COM + FOV. Options: 'origin', 'auto', np.array
def visualize(cluster, animation=False, scale=None, center='auto'):
    """
    Create a 3D visualization of a particle cluster using VPython
    
    Arguments:
        cluster     sphere or particle cluster
    """
    import vpython
    vec = vpython.vec
    scene = vpython.canvas(width=1000, height=800, background=vec(1,1,1))

    if type(cluster) == miepy.sphere_cluster:
        for i in range(cluster.Nparticles):
            position = cluster.position[i]
            radius = cluster.radius[i]
            paint = get_texture(cluster.material[i])

            sphere = vpython.sphere(pos=vec(*position), radius=radius, texture=paint.texture, color=paint.color)

    elif type(cluster) == miepy.cluster:
        for i in range(cluster.Nparticles):
            particle = cluster.particles[i]
            paint = get_texture(particle.material)

            obj = draw_particle(particle, paint)
