#!/bin/sh
# Creates a .gif of the wavefunction plots in ../tmp created by the SIL in Lanczos.h
# Note: make sure to run ./clean.sh if you have previously propagated a wavefunction
# for more steps and you now run a simulation for less steps. Otherwise, your GIF
# will have frames from the last simulation at the end.
rm Psi.*.dat
rm frame.*.png
cp ../tmp/Psi.*.dat .
gnuplot frame.plt
convert -delay 10 *.png out.gif
rm Psi.*.dat
rm frame.*.png
