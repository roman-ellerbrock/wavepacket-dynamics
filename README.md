# Standard Wavepacket Dynamics

This program simulates the quantum dynamics of low-dimensional systems by
expanding the wavefunction in a cartesian product basis. 

## Getting Started

Start by cloning the repo into your code directory. Make sure you check out
the branch 'tutorial'. 

0. Make sure you have all dependencies installed. On FIRE run: ml 
1. Create a build directory, e.g. named 'build' via: mkdir build
2. Change directory to the directory via: cd build
3. Run cmake on CMakeLists.txt via: cmake ..
4. Run make via: make

If everything worked out, you should have a binary named 'WavepacketDynamics' in
your build directory.

Dependencies:
a.) QuTree from https://github.com/roman-ellerbrock/QuTree
b.) yaml-cpp
c.) CMake >=3.2,
d.) C++ compilers
e.) A plotting software (e.g. gnuplot, matplotlib, ...)

## How to work with this tutorial

Different parts of the program are missing. The goal of the tutorial is to
complete the code parts that are missing and to run and analyze several systems.
Before jumping into solving the exercise, familiarize yourself with the codebase.

## Note on performance

This program is written for didactic purposes and readability of code
is preferred over performance. The code can with moderate effort be 
translated into efficient code, mostly by avoiding unnecessary copying of
objects.
