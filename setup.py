import os
import re
import sys
import platform
import subprocess

from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext
from setuptools import Command
from distutils.version import LooseVersion

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=''):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)

class nfmds(Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        print('building nfmds...')

class builder(build_ext):
    def run(self):
        ### Fortran compile
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

        ### Material database
        if not os.path.isdir('./miepy/materials/database'):
            command = ["./download_materials.sh"]
            print('Downloading material database')
            if subprocess.call(protoc_command) != 0:
                sys.exit(-1)
            install.run(self)
        else:
            print('Material database already downloaded')

        ### CMake compile
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

        subprocess.call(['cp', f'{extdir}/cpp.cpython-37m-x86_64-linux-gnu.so', 'miepy'])

setup(
    name = "miepy",
    version = "0.3",
    author = "John Parker",
    author_email = "japarker@uchicago.com",
    description = ("Python module to solve Maxwell's equations for a cluster of particles using the generalized multiparicle Mie theory (GMMT)"),
    license = "MIT",
    keywords = "electrodynamics mie scattering",
    url = "",
    packages=find_packages(),
    long_description=read('README.md'),
    install_requires=['numpy', 
                      'scipy',
                      'matplotlib',
                      'tqdm',
                      'sympy',
                      'pandas',
                      'pyyaml',
                      'numpy_quaternion',
                      'spherical_functions'],
    include_package_data = True,
    ext_modules=[CMakeExtension('miepy/cpp', './cpp')],
    cmdclass={
        'build_ext': builder,
        'build_nfmds': nfmds,
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
    ],
)
