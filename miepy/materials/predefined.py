from miepy.materials.load import load_material

def Ag(author='Johnson'):
    """Return silver material from MiePy data"""
    return load_material("Ag", author)
    # return load_material(miepy.__path__[0] + "/materials/ag.npy")

def Au(author='Johnson'):
    """Return gold material from MiePy data"""
    return load_material("Au", author)