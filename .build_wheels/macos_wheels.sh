#!/bin/bash

export CC=/usr/local/bin/gcc-8
export CXX=/usr/local/bin/g++-8

for VERSION in {3.5,3.6,3.7}; do
    export PATH=/Library/Frameworks/Python.framework/Versions/${VERSION}/bin:$PATH
    pip3 wheel ./ -w .build_wheels/temp
    delocate-wheel -w .build_wheels/wheelhouse .build_wheels/temp/*.whl
    rm .build_wheels/temp/*
done

rmdir .build_wheels/temp

