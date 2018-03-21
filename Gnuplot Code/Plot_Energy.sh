#!/bin/bash  
echo "This is a shell script for outputting all the energy PDFs and Skewnesses"
for i in {1..3}
do
gnuplot -c "Plot_Xu_PDF_2.gnu" electrons $i
gnuplot -c "Plot_Xu_PDF_2.gnu" ions $i
gnuplot -c "Plot_Xu_PDF_5.gnu" electrons $i
gnuplot -c "Plot_Xu_PDF_5.gnu" ions $i
gnuplot -c "Plot_Xu_Skewness_2.gnu" electrons $i
gnuplot -c "Plot_Xu_Skewness_2.gnu" ions $i
gnuplot -c "Plot_Xu_Skewness_5.gnu" electrons $i
gnuplot -c "Plot_Xu_Skewness_5.gnu" ions $i
done

