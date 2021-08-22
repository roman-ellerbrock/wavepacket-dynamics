#!/bin/sh
rm Psi.*.dat
rm frame.*.png
cp ../Psi.*.dat .
gnuplot frame.plt
convert -delay 8 *.png out.gif

