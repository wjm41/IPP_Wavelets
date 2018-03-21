Name='Ion acceleration PDFs at various scales'
set terminal postscript eps enhanced color size 25,15 font 'Helvetica,60'
File_name="Acc_ions_par_perp1_j".ARG1."_j".ARG2
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

Title='Parallel PDFs'
set title "Plot of " . Title

plot 'beta0p2_eps0p4/ions/PDFs/pdf_Acc_j'.ARG1.'_par_t2.dat' w lines lw 12 lt rgb "green" title "b=0.2 j=".ARG1, 'beta0p2_eps0p4/ions/PDFs/pdf_Acc_j'.ARG2.'_par_t2.dat' w lines lw 12 lt rgb "web-blue" title "b=0.2 j=".ARG2,'beta0p5_eps0p4/ions/PDFs/pdf_Acc_j'.ARG1.'_par_t2.dat' w lines lw 12 lt rgb "red" title "b=0.5 j=".ARG1, 'beta0p5_eps0p4/ions/PDFs/pdf_Acc_j'.ARG2.'_par_t2.dat' w lines lw 12 lt rgb "black" title "b=0.5 j=".ARG2
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

plot 'beta0p2_eps0p4/ions/PDFs/pdf_Acc_j'.ARG1.'_perp1_t2.dat' w lines lw 12 lt rgb "green" title "b=0.2 j=".ARG1, 'beta0p2_eps0p4/ions/PDFs/pdf_Acc_j'.ARG2.'_perp1_t2.dat' w lines lw 12 lt rgb "web-blue" title "b=0.2 j=".ARG2,'beta0p5_eps0p4/ions/PDFs/pdf_Acc_j'.ARG1.'_perp1_t2.dat' w lines lw 12 lt rgb "red" title "b=0.5 j=".ARG1, 'beta0p5_eps0p4/ions/PDFs/pdf_Acc_j'.ARG2.'_perp1_t2.dat' w lines lw 12 lt rgb "black" title "b=0.5 j=".ARG2
unset multiplot
#


#pause -1