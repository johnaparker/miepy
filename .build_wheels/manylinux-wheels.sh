#!/bin/bash

# pull docker image
#       docker pull quay.io/pypa/manylinux1_x86_64
# run docker environment interactively (from MiePy root directory)
#       docker run -it -v `pwd`:/io quay.io/pypa/manylinux1_x86_64
# run this script in docker (from MiePy root directory)
#       docker run --rm -v `pwd`:/io quay.io/pypa/manylinux1_x86_64 /io/.build_wheels/manylinux-wheels.sh
set -e -x

# install cmake
/opt/python/cp37-cp37m/bin/pip install cmake
#alias cmake='/opt/python/cp37-cp37m/bin/cmake'
export PATH=/opt/python/cp37-cp37m/bin/:$PATH

# unzip  libraries to build
cd /root
cp /io/.build_wheels/centos_build.zip .
unzip centos_build.zip

# configure GCC paths
export PATH=/root/gcc-8.2/bin:$PATH
export CC=/root/gcc-8.2/bin/gcc
export CXX=/root/gcc-8.2/bin/g++
export LD_LIBRARY_PATH=/root/gcc-8.2/lib64:$LD_LIBRARY_PATH

# install GSL
cd gsl-2.5
make install

# install Eigen
cd ../eigen3
mkdir build
cd build
cmake ..
make install

# Compile wheels
export PATH_BASE=$PATH
for VERSION in {cp35-cp35m,cp36-cp36m,cp37-cp37m}; do
    export PYBIN="/opt/python/${VERSION}/bin/"
    export PATH="${PYBIN}:${PATH_BASE}"
    export CPATH="/opt/python/${VERSION}/include/python${VERSION:2:1}.${VERSION:3:1}m:${CPATH}"

    "${PYBIN}/pip" install numpy
    "${PYBIN}/pip" install cmake
    #"${PYBIN}/pip" install -r /io/requirements.txt
    "${PYBIN}/pip" wheel /io/ -w /io/.build_wheels/wheelhouse/
    
    # Bundle external shared libraries into the wheel
    for whl in /io/.build_wheels/wheelhouse/miepy*.whl; do
        if [[ "$whl" = *"${VERSION}"* ]] && [[ "$whl" != *"manylinux1"* ]] ; then
            auditwheel repair "$whl" -w /io/.build_wheels/wheelhouse/
        fi
    done

done

# Bundle external shared libraries into the wheels
#for whl in /io/.build_wheels/wheelhouse/miepy*.whl; do
    #auditwheel repair "$whl" -w /io/.build_wheels/wheelhouse/
#done

# Install packages and test
#for PYBIN in /opt/python/{cp35-cp35m,cp36-cp36m,cp37-cp37m}/bin/; do
    #"${PYBIN}/pip" install miepy --no-index -f /io/.build_wheels/wheelhouse
    #(cd "$HOME"; "${PYBIN}/pytest" /io/tests)
#done
