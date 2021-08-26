reset
set terminal pngcairo

# for HO & FFT models
#set xr [-8:8]
#set yr [-8:8]
#set zr [0:0.14]

# for tully plotting
set view 90,0,1
#set xr [-450:900]

do for [t=0:149] {
    outfile = sprintf('frame.%03d.png',t)
    set output outfile
    out = 50
    q = t*out
    # splot 'Psi.'.t.'.dat' title '2D Harmonic Osciallator t = '.q w l
    splot 'Psi.'.t.'.dat' title 'Tully Model A, adiabatic g.s., t = '.q w l #, '../tully_adia1.dat' w l
}
