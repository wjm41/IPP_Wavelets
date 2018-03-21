Name='PDFs of acceleration magnitude for '.ARG1
set terminal postscript eps enhanced color size 25,15 font 'Helvetica,90'
File_name="beta0p2_Acc_multiple_j_".ARG1
set output File_name.'.eps'
set title Name
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

plot 'beta0p2_eps0p4/'.ARG1.'/PDFs/pdf_Acc_j2_mag_t2.dat' w lines lw 12 lt rgb "green" title ARG1." j=2", 'beta0p2_eps0p4/'.ARG1.'/PDFs/pdf_Acc_j4_mag_t2.dat' w lines lw 12 lt rgb "web-blue" title ARG1." j=4",'beta0p2_eps0p4/'.ARG1.'/PDFs/pdf_Acc_j6_mag_t2.dat' w lines lw 12 lt rgb "red" title ARG1." j=6", 'beta0p2_eps0p4/'.ARG1.'/PDFs/pdf_Acc_j8_mag_t2.dat' w lines lw 12 lt rgb "black" title ARG1." j=8"

#pause -1