#!/bin/bash  
echo "This is a shell script for outputting all the Acc PDFs"

gnuplot -c "Plot_Acc_par_perp1_2.gnu" 1 5
gnuplot -c "Plot_Acc_par_perp1_5.gnu" 1 5
gnuplot -c "Plot_Acc_mag_perp2_2.gnu" 1 5
gnuplot -c "Plot_Acc_mag_perp2_5.gnu" 1 5

gnuplot -c "Plot_Acc_par_perp1_electrons.gnu" 1 5
gnuplot -c "Plot_Acc_par_perp1_ions.gnu" 1 5
gnuplot -c "Plot_Acc_mag_perp2_electrons.gnu" 1 5
gnuplot -c "Plot_Acc_mag_perp2_ions.gnu" 1 5
