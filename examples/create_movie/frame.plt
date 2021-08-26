reset
set terminal pngcairo

# for HO & FFT models
#set xr [-8:8]
#set yr [-8:8]
#set zr [0:0.14]

# see only one side
#set view 90,0,1

do for [t=0:199] {
    outfile = sprintf('frame.%03d.png',t)
    set output outfile
    q = t / 10 
    splot 'Psi.'.t.'.dat' title '2D Harmonic Osciallator t = '.q w l
}
