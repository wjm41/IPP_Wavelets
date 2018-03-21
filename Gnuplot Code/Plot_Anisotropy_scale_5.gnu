Title="Velocity anisotropy against normalised timescale t_{j} for b=0.5"
set terminal postscript eps enhanced color size 25,15 font 'Helvetica,90' 
File_name="beta0p5_Anisotropy_Against_Scale"
set output File_name.'.eps'
set xlabel "t_{j}" 
set ylabel "Anisotropy"
#set grid
set zeroaxis
set logscale
set format y "10^{%T}"
set format x "10^{%T}"
set ytics scale 8, 3
set xtics scale 4
set yrange[1E-2:1E1]
set title "Plot of " . Title
E_VPAR_2= 'beta0p5_eps0p4/electrons/Moments/moments_vpar_t2.dat'
E_VPERP_2= 'beta0p5_eps0p4/electrons/Moments/moments_vperp_t2.dat'

I_VPAR_2= 'beta0p5_eps0p4/ions/Moments/moments_vpar_t2.dat'
I_VPERP_2= 'beta0p5_eps0p4/ions/Moments/moments_vperp_t2.dat'

plot '<paste '.E_VPERP_2.' '.E_VPAR_2 u 1:((($4/($2**(2)))-1)/(($10/($8**(2)))-1)) every ::3::20 with linespoints ps 10 lw 12 title "Electron Anisotropy",'<paste '.I_VPERP_2.' '.I_VPAR_2 u 1:((($4/($2**(2)))-1)/(($10/($8**(2)))-1)) every ::3::20 with linespoints ps 10 lw 12 title "Ion Anisotropy", 1 with lines dt 2 lt rgb "red" lw 12

#pause -1