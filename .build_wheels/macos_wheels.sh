#!/bin/bash

export MACOSX_DEPLOYMENT_TARGET=10.9
export CPATH=/usr/local/include/:$CPATH

export CC=/usr/local/bin/gcc-8
export CXX=/usr/local/bin/g++-8

for VERSION in {3.6,3.7,3.8}; do
    export PATH=/Library/Frameworks/Python.framework/Versions/${VERSION}/bin:$PATH
    pip3 wheel ./ -w .build_wheels/temp
    delocate-wheel -w .build_wheels/wheelhouse .build_wheels/temp/*.whl
    rm .build_wheels/temp/*
done

rmdir .build_wheels/temp

