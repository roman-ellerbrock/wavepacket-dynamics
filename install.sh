#!/bin/sh
# initialize submodules
git submodule update --init --recursive

# Create build directory and run CMakeLists.txt
mkdir -p build
cd build
cmake ..

# build the executable WavepacketDynamics and return to original directory
make
cd ..
