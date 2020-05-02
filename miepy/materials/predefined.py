from miepy.material_functions.load import load_material
from miepy.material_functions.create import constant_material

### mediums
class materials:
    @staticmethod
    def vacuum():
        return constant_material(index=1, name='vacuum')

    @staticmethod
    def air():
        return constant_material(index=1, name='air')

    @staticmethod
    def water():
        return constant_material(index=1.33, name='water')

    ### metals
    @staticmethod
    def metal():
        return constant_material(eps=complex(1,1e40), name='metal')

    @staticmethod
    def Ag(author='Johnson'):
        """Return silver material from MiePy data"""
        return load_material("Ag", author)

    @staticmethod
    def Au(author='Johnson'):
        """Return gold material from MiePy data"""
        return load_material("Au", author)

    @staticmethod
    def Al(author='Rakic'):
        """Return aluminum material from MiePy data"""
        return load_material("Al", author)

    @staticmethod
    def Ni(author='Johnson'):
        """Return Nickel material from MiePy data"""
        return load_material("Ni", author)

    @staticmethod
    def Cu(author='Johnson'):
        """Return Copper material from MiePy data"""
        return load_material("Cu", author)

    @staticmethod
    def Co(author='Johnson'):
        """Return Cobalt material from MiePy data"""
        return load_material("Co", author)

    @staticmethod
    def Pt(author='Rakic-LD'):
        """Return Platinum material from MiePy data"""
        return load_material("Pt", author)

    @staticmethod
    def TiO2(author='Siefke'):
        """Return TiO2 material from MiePy data"""
        return load_material("TiO2", author)

    @staticmethod
    def silica(author='Gao'):
        """Return silica material from MiePy data"""
        return load_material("SiO2", author)

    @staticmethod
    def Si(author='Green-2008'):
        """Return silicon material from MiePy data"""
        return load_material("Si", author)
