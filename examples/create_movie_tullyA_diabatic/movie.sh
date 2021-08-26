#!/bin/sh
rm Psi.*.dat
rm frame.*.png
cp ../tmp/Psi.*.dat .
gnuplot frame.plt
convert -delay 10 *.png out.gif

