Name='Acceleration PDFs for b=0.5 at various scales'
set terminal postscript eps enhanced color size 25,15 font 'Helvetica,60'
File_name="beta0p5_Acc_species_comparison_mag_perp2_j".ARG1."_j".ARG2
set output File_name.'.eps'
set multiplot layout 1,2 title Name
#
set xlabel "No. of standard deviations" 
set ylabel "PDF"
set zeroaxis
set logscale y
set ytics scale 8, 3
set xtics scale 4
set format y "10^{%T}"
set xrange[-8:8]
set yrange[1E-6:2]

Title='Magnitude'
set title "Plot of " . Title

E_J1=ARG1+6
E_J2=ARG2+6

plot 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.E_J1.'_mag_t2.dat' w lines lw 12 lt rgb "green" title "Electrons j=".E_J1, 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.E_J2.'_mag_t2.dat' w lines lw 12 lt rgb "web-blue" title "Electrons j=".E_J2,'beta0p5_eps0p4/ions/PDFs/pdf_Acc_j'.ARG1.'_mag_t2.dat' w lines lw 12 lt rgb "red" title "Ions j=".ARG1, 'beta0p5_eps0p4/ions/PDFs/pdf_Acc_j'.ARG2.'_mag_t2.dat' w lines lw 12 lt rgb "black" title "Ions j=".ARG2
unset key
#
set xlabel "No. of standard deviations" 
set ylabel "PDF"
set zeroaxis
set logscale y
set format y "10^{%T}"
set xrange[-8:8]
set yrange[1E-6:2]

Title='Perpendicular (in direction of vxb)'
set title "Plot of " . Title

plot 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.E_J1.'_perp2_t2.dat' w lines lw 12 lt rgb "green" title "Electrons j=".E_J1, 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.E_J2.'_perp2_t2.dat' w lines lw 12 lt rgb "web-blue" title "Electrons j=".E_J2,'beta0p5_eps0p4/ions/PDFs/pdf_Acc_j'.ARG1.'_perp2_t2.dat' w lines lw 12 lt rgb "red" title "Ions j=".ARG1, 'beta0p5_eps0p4/ions/PDFs/pdf_Acc_j'.ARG2.'_perp2_t2.dat' w lines lw 12 lt rgb "black" title "Ions j=".ARG2
unset multiplot
#


#pause -1