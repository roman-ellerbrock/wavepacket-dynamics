reset
set terminal pngcairo

# for HO & FFT models
#set xr [-8:8]
#set yr [-8:8]
#set zr [0:0.14]

# for tully plotting
set view 90,0,1
set xr [-12:20]

do for [t=0:149] {
    outfile = sprintf('frame.%03d.png',t)
    set output outfile
    q = t/5
    splot 'Psi.'.t.'.dat' title '2D Harmonic Osciallator t = '.q w l
}
