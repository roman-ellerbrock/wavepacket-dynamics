reset
set terminal pngcairo

# for HO & FFT models
#set xr [-8:8]
set yr [0:0.25]

do for [t=0:149] {
    outfile = sprintf('frame.%03d.png',t)
    set output outfile
    q = t/10
    plot 'Psi.'.t.'.dat' title '1D Harmonic Osciallator t = '.q w l
}
