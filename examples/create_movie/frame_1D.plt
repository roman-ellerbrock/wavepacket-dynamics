reset
set terminal pngcairo

# for HO & FFT models
set xr [-450:450]
#set yr [0:0.15]

set xl "x"
set yl "probability density"

do for [t=0:120] {
    outfile = sprintf('frame.%03d.png',t)
    set output outfile
    q = t/10
    plot 'Psi.'.t.'.dat' title 'Tully Model A, adiabatic g.s., t = '.q w l, '../tully_adia1.dat' w l
}
