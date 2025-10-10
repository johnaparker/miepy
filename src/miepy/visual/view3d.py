import miepy
from collections import namedtuple
import numpy as np
from itertools import count
import matplotlib.colors as cm

nm = 1e-9

def get_paint(material):
    """Return a (texture, color, shininess) tuple for a given material"""
    import vpython
    paint = namedtuple('paint', ['texture', 'color', 'shininess'])

    if material.name == 'Ag':
        return paint(vpython.textures.metal, vpython.vec(*[.9]*3), shininess=0.15)
    if material.name == 'Au':
        return paint(vpython.textures.metal, vpython.vec(255/255, 210/255, 0/255), shininess=0.07)

    default_color = vpython.vec(*cm.to_rgb('C0'))
    return paint(None, default_color, shininess=1)

def draw_particle(particle, paint, **kwargs):
    """Draw a (non-spherical) particle with a given paint"""
    import vpython
    vec = vpython.vec

    if type(particle) is miepy.sphere:
        return vpython.sphere(pos=vec(*particle.position)/nm,
                              radius=particle.radius/nm,
                              texture=paint.texture,
                              color=paint.color, **kwargs)

    else:
        q = particle.orientation
        theta, phi = miepy.quaternion.as_spherical_coords(q)
        axis = vec(np.sin(theta)*np.cos(phi),
                   np.sin(theta)*np.sin(phi),
                   np.cos(theta))
        rot_vec = miepy.quaternion.as_rotation_vector(q)
        angle = np.linalg.norm(rot_vec)
        if angle != 0:
            axis_rot = rot_vec/angle
        else:
            axis_rot = [0,0,1]

        if type(particle) is miepy.cylinder:
            return vpython.cylinder(pos=(vec(*particle.position) - particle.height/2*axis)/nm,
                                    radius=particle.radius/nm,
                                    length=particle.height/nm,
                                    texture=paint.texture,
                                    axis=axis,
                                    color=paint.color,
                                    shininess=paint.shininess, **kwargs)

        elif type(particle) is miepy.spheroid:
            return vpython.ellipsoid(pos=vec(*particle.position)/nm,
                                    width=2*particle.axis_xy/nm,
                                    height=2*particle.axis_xy/nm,
                                    length=2*particle.axis_z/nm,
                                    texture=paint.texture,
                                    axis=axis,
                                    color=paint.color,
                                    shininess=paint.shininess, **kwargs)

        elif type(particle) is miepy.ellipsoid:
            obj = vpython.ellipsoid(pos=vec(*particle.position)/nm,
                                    width=2*particle.rx/nm,
                                    height=2*particle.ry/nm,
                                    length=2*particle.rz/nm,
                                    texture=paint.texture,
                                    color=paint.color,
                                    axis=vec(0,0,1),
                                    shininess=paint.shininess, **kwargs)

            obj.rotate(angle=angle, axis=vec(*axis_rot))
            return obj

        elif type(particle) is miepy.regular_prism:
            rotate = -(2*np.pi - (particle.N-2)*np.pi/(particle.N))/2
            shape = vpython.shapes.ngon(np=particle.N, length=particle.width/nm, rotate=rotate)
            path = [vec(0,0,0), vec(0,0,particle.height/nm)]

            obj = vpython.extrusion(pos=vec(*particle.position)/nm,
                                    path=path,
                                    shape=shape,
                                    texture=paint.texture,
                                    color=paint.color,
                                    shininess=paint.shininess, **kwargs)

            obj.rotate(angle=angle, axis=vec(*axis_rot))
            return obj

#TODO: scale bar: auto-detect position and size
#TODO: center: COM + FOV. Options: 'origin', 'auto', np.array
def visualize(cluster, animation=False, transparent=False, scale=None, origin='auto'):
    """
    Create a 3D visualization of a particle cluster using VPython
    
    Arguments:
        cluster      sphere or particle cluster
        animation    (bool) an animation visualization
        transparent  (bool) transparent particles
        scale        (str) if set, display a scale bar
        origin       origin of camera view ('auto' for centroid, 'cluster' for cluster origin, or np.array)
    """
    import vpython
    vec = vpython.vec
    # scene = vpython.canvas(width=750, height=600, background=vec(1,1,1))

    if origin == 'auto':
        origin = vec(*np.average(cluster.position, axis=0))/nm
    elif origin == 'cluster':
        origin = vec(*cluster.origin)/nm
    else:
        origin = vec(*origin)/nm

    scene = vpython.canvas(width=500, height=325, background=vec(1,1,1), center=origin)

    if type(cluster) == miepy.sphere_cluster:
        for i in range(cluster.Nparticles):
            position = cluster.position[i]/nm
            radius = cluster.radius[i]/nm
            paint = get_paint(cluster.material[i])

            sphere = vpython.sphere(pos=vec(*position), radius=radius, texture=paint.texture, color=paint.color, shininess=paint.shininess)

    elif type(cluster) == miepy.cluster:
        for i in range(cluster.Nparticles):
            particle = cluster.particles[i]
            paint = get_paint(particle.material)

            obj = draw_particle(particle, paint)

    if transparent:
        for obj in scene.objects:
            obj.opacity = .5
            
    if cluster.interface is not None:
        scene.autoscale = False
        length = height = 1e-3/nm
        width = 1e-3/nm
        z = cluster.interface.z
        box = vpython.box(pos=vec(0,0,-z-width/2), width=width, height=height, length=length,
        color=vec(0.6,0.6,0.6), shininess=0)

    if animation:
        T = 300
        # cluster = vpython.compound(scene.objects)
        # cluster.rotate(angle=2*np.pi/T, origin=vec(0,0,0), axis=vec(0,1,0))

        scene.append_to_caption('\n')

        def setspeed(s):
            nonlocal T
            wt.text = 'f = {:1.2f} Hz'.format(s.value)
            T = 30/s.value
            
        sl = vpython.slider(min=0.01, max=.3, value=.1, length=180, bind=setspeed, right=15, top=7)
        wt = vpython.wtext(text='f = {:1.2f} Hz'.format(sl.value))

        def transparency(b):
            if b.checked:
                for obj in scene.objects:
                    obj.opacity = .5
            else:
                for obj in scene.objects:
                    obj.opacity = 1

        scene.append_to_caption('    ')
        vpython.checkbox(bind=transparency, text='Transparent')

        running = True
        def Run(b):
            nonlocal running
            running = not running
            if running:
                b.text = " ▌▌"
                b.color = vec(1,0,0)
            else:
                b.text = "  ►  "
                b.color = vec(0,.6,0)

        scene.append_to_caption('       ')
        vpython.button(text=" ▌▌", bind=Run, color=vec(1,0,0))


        for t in count():
            vpython.rate(30)

            # phi = 2*np.pi*t/T
            # scene.forward = vec(-np.sin(phi), 0, -np.cos(phi))
            
            if running:
                for obj in scene.objects:
                    obj.rotate(angle=2*np.pi/T, origin=origin, axis=vec(0,1,0))
