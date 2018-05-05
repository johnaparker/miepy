from miepy.materials.load import load_material
from miepy.materials.create import constant_material

### mediums
def vacuum():
    return constant_material(1)

def water():
    return constant_material(1.33**2)

### metals
def Ag(author='Johnson'):
    """Return silver material from MiePy data"""
    return load_material("Ag", author)

def Au(author='Johnson'):
    """Return gold material from MiePy data"""
    return load_material("Au", author)

def Al(author='Rakic'):
    """Return aluminum material from MiePy data"""
    return load_material("Al", author)
