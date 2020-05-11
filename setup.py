import os
import re
import sys
import platform
import subprocess

from setuptools import setup, find_packages, Extension
from distutils.command.build import build
from setuptools.command.build_ext import build_ext
from setuptools import Command
from distutils.version import LooseVersion


NAME = 'miepy'
DESCRIPTION = "Solve Maxwell's equations for a cluster of particles using the generalized multiparticle Mie theory (GMMT)"
URL = ''
EMAIL = 'japarker@uchicago.edu'
AUTHOR = 'John Parker'
KEYWORDS = 'electrodynamics mie scattering'
# REQUIRES_PYTHON = '>=3.6.0'
VERSION = '0.5.0'
LICENSE = 'GPLv3'

REQUIRED = [
    'numpy', 
    'scipy',
    'matplotlib',
    'tqdm',
    'sympy',
    'pandas',
    'pyyaml',
    'numpy_quaternion',
    'spherical_functions',
    'vpython',
]


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


def unzip_material_database():
    import zipfile
    path = 'miepy/materials/database.zip'
    with zipfile.ZipFile(path, 'r') as zip_ref:
        zip_ref.extractall('miepy/materials')


def build_nfmds(build_direc, lib_dir):
    import pathlib
    if platform.system() == "Windows":
        build_direc = build_direc.replace('\\', r'/')
        lib_dir = lib_dir.replace('\\', r'/')
        
    src_dir = 'miepy/tmatrix/nfmds'
    obj_dir = '../../../{direc}/nfmds'.format(direc=build_direc)
    exe_dir = '../../../{direc}/miepy/bin'.format(direc=lib_dir)
    exe_dir = '../../bin'

    exe_dir_root = 'miepy/bin'
    obj_dir_root = '{build_direc}/nfmds'.format(build_direc=build_direc)
    pathlib.Path(exe_dir_root).mkdir(exist_ok=True) 
    pathlib.Path(obj_dir_root).mkdir(exist_ok=True) 

    command = ['make', 'objdir={obj_dir}'.format(obj_dir=obj_dir),
           'exedir={exe_dir}'.format(exe_dir=exe_dir)]
    subprocess.check_call(' '.join(command), cwd=src_dir, shell=True)


class builder(build):
    def run(self):
        if not os.path.isdir('miepy/materials/database'):
            unzip_material_database()

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        build_nfmds(self.build_temp, self.build_lib)

        super().run()


class builder_ext(build_ext):
    def run(self):
        try:
            out = subprocess.check_output(['cmake', '--version'])
        except OSError:
            raise RuntimeError("CMake must be installed to build the following extensions: " +
                               ", ".join(e.name for e in self.extensions))

        if platform.system() == "Windows":
            cmake_version = LooseVersion(re.search(r'version\s*([\d.]+)', out.decode()).group(1))
            if cmake_version < '3.1.0':
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=' + extdir,
                      '-DPYTHON_EXECUTABLE=' + sys.executable]

        cfg = 'Debug' if self.debug else 'Release'
        build_args = ['--config', cfg]

        if platform.system() == "Windows":
            cmake_args += ['-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}'.format(cfg.upper(), extdir)]
            if sys.maxsize > 2**32:
                cmake_args += ['-A', 'x64']
            build_args += ['--', '/m']
        else:
            cmake_args += ['-DCMAKE_BUILD_TYPE=' + cfg]
            build_args += ['--', '-j2']

        env = os.environ.copy()
        env['CXXFLAGS'] = '{} -DVERSION_INFO=\\"{}\\"'.format(env.get('CXXFLAGS', ''),
                                                              self.distribution.get_version())
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(['cmake', ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env)
        subprocess.check_call(['cmake', '--build', '.'] + build_args, cwd=self.build_temp)

setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    description=DESCRIPTION,
    license=LICENSE,
    keywords=KEYWORDS,
    url=URL,
    packages=find_packages(),
    long_description=read('README.md'),
    long_description_content_type='text/markdown',
    install_requires=REQUIRED,
    include_package_data = True,
    ext_modules=[CMakeExtension('miepy/cpp', './cpp')],
    cmdclass={
        'build': builder,
        'build_ext': builder_ext,
    },
    classifiers=[
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python',
        'Programming Language :: C++',
        'Programming Language :: Fortran',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: MacOS',
        'Operating System :: Microsoft :: Windows',
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
    ],
    zip_safe=False,
)
