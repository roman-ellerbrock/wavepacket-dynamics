#!/bin/sh
rm Psi.*.dat
rm frame.*.png
cp ../tmp/Psi.*.dat .
gnuplot frame_1D.plt
convert -delay 10 *.png out.gif

