reset
set terminal pngcairo

set xr [-8:8]
set yr [-8:8]
set zr [0:0.14]


do for [t=0:149] {
    outfile = sprintf('frame.%03d.png',t)
    set output outfile
    q = t/5
    splot 'Psi.'.t.'.dat' title '2D Harmonic Osciallator t = '.q w l
}
