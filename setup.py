import os
import subprocess
from setuptools import setup
from setuptools.command.build_ext import build_ext as _build_ext

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

class build_ext(_build_ext):
    def run(self):
        if not os.path.isdir('miepy/bin'):
            subprocess.call(['mkdir', 'miepy/bin'])

        if not os.path.exists('miepy/bin/tmatrix'):
            subprocess.call(['make clean'], cwd='./miepy/tmatrix/nfmds', shell=True)
            subprocess.call(['make precision=double'], cwd='./miepy/tmatrix/nfmds', shell=True)
            subprocess.call(['make clean'], cwd='./miepy/tmatrix/nfmds', shell=True)
            subprocess.call(['mv', './miepy/tmatrix/nfmds/tmatrix', 'miepy/bin'])

        if not os.path.exists('miepy/bin/tmatrix_extended'):
            subprocess.call(['make clean'], cwd='./miepy/tmatrix/nfmds', shell=True)
            subprocess.call(['make precision=quad'], cwd='./miepy/tmatrix/nfmds', shell=True)
            subprocess.call(['make clean'], cwd='./miepy/tmatrix/nfmds', shell=True)
            subprocess.call(['mv', './miepy/tmatrix/nfmds/tmatrix_extended', 'miepy/bin'])

        if not os.path.isdir('./miepy/materials/database'):
            command = ["./download_materials.sh"]
            print('Downloading material database')
            if subprocess.call(protoc_command) != 0:
                sys.exit(-1)
            install.run(self)
        else:
            print('Material database already downloaded')

setup(
    cmdclass={
        'build_ext': build_ext,
    },
    name = "miepy",
    version = "0.2",
    author = "John Parker",
    author_email = "japarker@uchicago.com",
    description = ("Python module to solve Maxwell's equations for a cluster of spheres using the generalized multiparicle Mie theory"),
    license = "MIT",
    keywords = "electrodynamics mie scattering",
    url = "",
    packages=['miepy'],
    long_description=read('README.md'),
    install_requires=['numpy', 'scipy', 'matplotlib'],
    include_package_data = True,
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
    ],
)
