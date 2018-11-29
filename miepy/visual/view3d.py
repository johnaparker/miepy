import miepy
from collections import namedtuple
import numpy as np
from itertools import count
import matplotlib.colors as cm

def get_texture(material):
    """Return a (texture, color) tuple for a given material"""
    import vpython
    paint = namedtuple('paint', ['texture', 'color', 'shininess'])

    if material.name == 'Ag':
        return paint(vpython.textures.metal, vpython.vec(*[.9]*3), shininess=0.15)
    if material.name == 'Au':
        return paint(vpython.textures.metal, vpython.vec(255/255, 210/255, 20/255), shininess=0.10)

    default_color = vpython.vec(*cm.to_rgb('C0'))
    return paint(None, default_color, shininess=1)

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
        theta, phi = miepy.quaternion.as_spherical_coords(q)
        axis = vec(np.sin(theta)*np.cos(phi),
                   np.sin(theta)*np.sin(phi),
                   np.cos(theta))

        if type(particle) is miepy.cylinder:
            return vpython.cylinder(pos=vec(*particle.position) - particle.height/2*axis,
                                    radius=particle.radius,
                                    length=particle.height,
                                    texture=paint.texture,
                                    axis=axis,
                                    color=paint.color,
                                    shininess=paint.shininess, **kwargs)

        elif type(particle) is miepy.spheroid:
            return vpython.ellipsoid(pos=vec(*particle.position),
                                    width=2*particle.axis_xy,
                                    height=2*particle.axis_xy,
                                    length=2*particle.axis_z,
                                    texture=paint.texture,
                                    axis=axis,
                                    color=paint.color,
                                    shininess=paint.shininess, **kwargs)

#TODO: animation: pan camera around (widgets for rotation speed, transparency, pause, play, stop)
#TODO: scale bar: auto-detect position and size
#TODO: center: COM + FOV. Options: 'origin', 'auto', np.array
def visualize(cluster, animation=False, transparent=False, scale=None, center='auto'):
    """
    Create a 3D visualization of a particle cluster using VPython
    
    Arguments:
        cluster      sphere or particle cluster
        animation    (bool) an animation visualization
        transparent  (bool) transparent particles
        scale        (str) if set, display a scale bar
        origin       origin of camera view
    """
    import vpython
    vec = vpython.vec
    # scene = vpython.canvas(width=750, height=600, background=vec(1,1,1))
    scene = vpython.canvas(width=500, height=400, background=vec(1,1,1))

    if type(cluster) == miepy.sphere_cluster:
        for i in range(cluster.Nparticles):
            position = cluster.position[i]
            radius = cluster.radius[i]
            paint = get_texture(cluster.material[i])

            sphere = vpython.sphere(pos=vec(*position), radius=radius, texture=paint.texture, color=paint.color, shininess=paint.shininess)

    elif type(cluster) == miepy.cluster:
        for i in range(cluster.Nparticles):
            particle = cluster.particles[i]
            paint = get_texture(particle.material)

            obj = draw_particle(particle, paint)

    if transparent:
        for obj in scene.objects:
            obj.opacity = .5

    if animation:
        T = 300
        # cluster = vpython.compound(scene.objects)
        # cluster.rotate(angle=2*np.pi/T, origin=vec(0,0,0), axis=vec(0,1,0))

        for t in count():
            vpython.rate(30)

            # phi = 2*np.pi*t/T
            # scene.forward = vec(-np.sin(phi), 0, -np.cos(phi))
            
            for obj in scene.objects:
                obj.rotate(angle=2*np.pi/T, origin=vec(0,0,0), axis=vec(0,1,0))
