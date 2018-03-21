#!/bin/bash  
echo "This is a shell script for outputting all the power PDFs and skewnesses"

gnuplot -c "Plot_Power_diff_j_2.gnu" 1 5
gnuplot -c "Plot_Power_diff_j_5.gnu" 1 5

gnuplot -c "Plot_Power_diff_j_electrons.gnu" 1 6
gnuplot -c "Plot_Power_diff_j_electrons.gnu" 8 11
gnuplot -c "Plot_Power_diff_j_ions.gnu" 1 8
