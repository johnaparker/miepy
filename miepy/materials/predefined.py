import miepy
from miepy.materials.load import load_material

### mediums
def vacuum():
    return miepy.constant_material(1)

def water():
    return miepy.constant_material(1.33**2)

### metals
def Ag(author='Johnson'):
    """Return silver material from MiePy data"""
    return load_material("Ag", author)

def Au(author='Johnson'):
    """Return gold material from MiePy data"""
    return load_material("Au", author)
