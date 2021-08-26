reset
set terminal pngcairo

# for HO & FFT models
#set xr [-8:8]
#set yr [-8:8]
#set zr [0:0.14]

# for tully plotting
#set view 90,0,1
set xr [-1000:1000]
set zr [0:]

do for [t=0:149] {
    outfile = sprintf('frame.%03d.png',t)
    set output outfile
    out = 50
    q = t*out
    splot 'Psi.'.t.'.dat' title 'Tully Model A, diabatic g.s., t = '.q 
}
