# Standard Wavepacket Dynamics

This program simulates the quantum dynamics of low-dimensional systems by
expanding the wavefunction in a cartesian product basis. 

## Getting Started

Start by cloning the repo into your code directory. Make sure you check out
the correct branch (see slides) before proceeding.



0. Make sure you have the dependencies installed.
1. Run ./install.sh

If everything worked out, you should have a binary named 'WavepacketDynamics' in
your build directory. Otherwise, try fixing the issues or ask me (romanellerbrock@gmail.com).

In order to re-build the project, simply run
1. cd build
2. make

FAQ:
- If compiling QuTree crashes, try checking out branch *feature-distribute* in external/QuTree. This branch is designed to have less requirements.

Dependencies:
- CMake >=3.16,
- C++ compilers
- A plotting software (e.g. gnuplot, matplotlib, ...)

Git submodules (project will automatically install these):
- QuTree (https://github.com/roman-ellerbrock/QuTree)
- yaml-cpp

## How to run the code
The code can be executed via *./WavepacketDynamics ${input}*. 
The binary is found in your build directory after compiling the code.
Try running the code on examples/basis_eigenstates.yaml.

## How to work with this tutorial

Different parts of the program are missing. The goal of the tutorial is to
complete the code parts that are missing and to run and analyze several systems.
Before jumping into solving the exercise, familiarize yourself with the codebase.

## Note on performance

This program is written for educational purposes and readability of code
is preferred over performance. The code can with moderate effort be 
translated into efficient code, mostly by avoiding unnecessary copying of
objects. The external QuTree
