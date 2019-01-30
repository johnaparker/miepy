from miepy.materials.load import load_material
from miepy.materials.create import constant_material
import numpy as np

### mediums
def vacuum():
    return constant_material(1, name='vacuum')

def air():
    return constant_material(1, name='air')

def water():
    return constant_material(1.33**2, name='water')

### metals
def metal():
    return constant_material(1 + 1j*1e40, name='metal')

def Ag(author='Johnson'):
    """Return silver material from MiePy data"""
    return load_material("Ag", author)

def Au(author='Johnson'):
    """Return gold material from MiePy data"""
    return load_material("Au", author)

def Al(author='Rakic'):
    """Return aluminum material from MiePy data"""
    return load_material("Al", author)

def Ni(author='Johnson'):
    """Return Nickel material from MiePy data"""
    return load_material("Ni", author)

def Cu(author='Johnson'):
    """Return Copper material from MiePy data"""
    return load_material("Cu", author)

def Co(author='Johnson'):
    """Return Cobalt material from MiePy data"""
    return load_material("Co", author)

def Pt(author='Rakic-LD'):
    """Return Platinum material from MiePy data"""
    return load_material("Pt", author)

def TiO2(author='Siefke'):
    """Return TiO2 material from MiePy data"""
    return load_material("TiO2", author)

def silica(author='Gao'):
    """Return silica material from MiePy data"""
    return load_material("SiO2", author)
