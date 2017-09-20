import os
from setuptools import setup

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "miepy",
    version = "0.1.0",
    author = "John Parker",
    author_email = "japarker@uchicago.com",
    description = ("Python module to calcuate scattering coefficients of a plane wave incident on a sphere or core-shell structure using Mie theory"),
    license = "MIT",
    keywords = "mie scattering bohren huffman core-shell",
    url = "http://packages.python.org/an_example_pypi_project",
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
