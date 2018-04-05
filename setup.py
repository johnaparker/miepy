import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "miepy",
    version = "0.1",
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
