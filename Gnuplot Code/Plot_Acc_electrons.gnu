Name='PDFs of acceleration decomposition for ELECTRONS quarterways through the simulation at scale 2^{-'.ARG1.'}'
set terminal postscript eps enhanced color size 25,15 font 'Helvetica,60'
File_name="Acc_electrons_comparison_j".ARG1
set output File_name.'.eps'
set multiplot layout 2,2 title Name
#
set xlabel "No. of standard deviations" 
set ylabel "PDF"
set zeroaxis
set logscale y
set format y "10^{%T}"
set xrange[-8:8]
set yrange[1E-6:2]

Title='Magnitude PDFs'
set title "Plot of " . Title

plot 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_mag_t1.dat' w lines lw 12 lt rgb "magenta" title "b=0.2 T1", 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_mag_t2.dat' w lines lw 12 lt rgb "green" title "b=0.2 T2", 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_mag_t3.dat' w lines lw 12 lt rgb "web-blue" title "b=0.2 T3",'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_mag_t1.dat' w lines lw 12 lt rgb "red" title "b=0.5 T1", 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_mag_t2.dat' w lines lw 12 lt rgb "orange" title "b=0.5 T2", 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_mag_t3.dat' w lines lw 12 lt rgb "black" title "b=0.5 T3"
unset key
#
set xlabel "No. of standard deviations" 
set ylabel "PDF"
set zeroaxis
set logscale y
set format y "10^{%T}"
set xrange[-8:8]
set yrange[1E-6:2]

Title='Perpendicular (in direction of vperp) PDFs'
set title "Plot of " . Title

plot 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp1_t1.dat' w lines lw 12 lt rgb "magenta" title "b=0.2 T1", 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp1_t2.dat' w lines lw 12 lt rgb "green" title "b=0.2 T2", 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp1_t3.dat' w lines lw 12 lt rgb "web-blue" title "b=0.2 T3",'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp1_t1.dat' w lines lw 12 lt rgb "red" title "b=0.5 T1", 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp1_t2.dat' w lines lw 12 lt rgb "orange" title "b=0.5 T2", 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp1_t3.dat' w lines lw 12 lt rgb "black" title "b=0.5 T3"
unset key
#
set xlabel "No. of standard deviations" 
set ylabel "PDF"
set zeroaxis
set logscale y
set format y "10^{%T}"
set xrange[-8:8]
set yrange[1E-6:2]

Title='Parallel PDFs'
set title "Plot of " . Title

plot 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_par_t1.dat' w lines lw 12 lt rgb "magenta" title "b=0.2 T1", 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_par_t2.dat' w lines lw 12 lt rgb "green" title "b=0.2 T2", 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_par_t3.dat' w lines lw 12 lt rgb "web-blue" title "b=0.2 T3",'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_par_t1.dat' w lines lw 12 lt rgb "red" title "b=0.5 T1", 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_par_t2.dat' w lines lw 12 lt rgb "orange" title "b=0.5 T2", 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_par_t3.dat' w lines lw 12 lt rgb "black" title "b=0.5 T3"
unset key
#
set xlabel "No. of standard deviations" 
set ylabel "PDF"
set zeroaxis
set logscale y
set format y "10^{%T}"
set xrange[-8:8]
set yrange[1E-6:2]

Title='Perpendicular (in direction of vxb) PDFs'
set title "Plot of " . Title

plot 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp2_t1.dat' w lines lw 12 lt rgb "magenta" title "b=0.2 T1", 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp2_t2.dat' w lines lw 12 lt rgb "green" title "b=0.2 T2", 'beta0p2_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp2_t3.dat' w lines lw 12 lt rgb "web-blue" title "b=0.2 T3",'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp2_t1.dat' w lines lw 12 lt rgb "red" title "b=0.5 T1", 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp2_t2.dat' w lines lw 12 lt rgb "orange" title "b=0.5 T2", 'beta0p5_eps0p4/electrons/PDFs/pdf_Acc_j'.ARG1.'_perp2_t3.dat' w lines lw 12 lt rgb "black" title "b=0.5 T3"
unset multiplot
#

#pause -1