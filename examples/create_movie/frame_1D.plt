reset
set terminal pngcairo

# for HO & FFT models
set xr [-20:20]
set yr [0:0.15]

set xl "x"
set yl "probability density"

f(x) = 0.1*exp(-x*x)
do for [t=0:120] {
    outfile = sprintf('frame.%03d.png',t)
    set output outfile
    q = t/10
    plot 'Psi.'.t.'.dat' title '1D Eckhard Barrier t = '.q w l, f(x) title 'V(x)'
}
