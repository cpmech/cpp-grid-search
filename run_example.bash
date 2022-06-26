#!/bin/bash

mkdir -p build
cd build
cmake ../
make
#./example
./example_triangles

