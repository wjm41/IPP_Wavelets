Title="b = 0.5 for ".ARG1
set terminal postscript eps enhanced color size 25,15 font 'Helvetica,90' 
File_name="beta0p5_".ARG1."_Flatness_Against_Scale"
set output File_name.'.eps'
set xlabel "t_{j}" 
set ylabel "Flatness"
#set grid
set zeroaxis
set logscale 
set format y "10^{%T}"
set format x "10^{%T}"
set ytics scale 8, 3
set xtics scale 4
set title "Plot of " . Title
plot 'beta0p5_eps0p4/'.ARG1.'/Moments/moments_vpar_t2.dat' u 1:($4/($2**(2))) every ::1::20 with linespoints lw 12 ps 15 title "Parallel", 'beta0p5_eps0p4/'.ARG1.'/Moments/moments_vperp_t2.dat' u 1:($4/($2**(2))) every ::1::20 with linespoints lw 12 ps 15 title "Perpendicular", 3 with lines dt 2 lw 15 lt rgb "black" title "Gaussian Flatness = 3", 4.2 with lines dt 3 lw 15 lt rgb "red" title "Logistic Flatness = 4.2", 5 with lines dt 4 lw 15 lt rgb "blue" title "Hyperbolic Secant Flatness = 5"
#pause -1
