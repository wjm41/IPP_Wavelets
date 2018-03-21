#!/bin/bash  
echo "This is a shell script for outputting all the Flatnesses as well as the anisotropy"

gnuplot -c "Plot_Flatness_scale_2.gnu" ions
gnuplot -c "Plot_Flatness_scale_2.gnu" electrons
gnuplot -c "Plot_Flatness_scale_5.gnu" ions
gnuplot -c "Plot_Flatness_scale_5.gnu" electrons

gnuplot -c "Plot_Anisotropy_scale_2.gnu" 
gnuplot -c "Plot_Anisotropy_scale_5.gnu" 