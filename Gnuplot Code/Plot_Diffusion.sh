#!/bin/bash  
echo "This is a shell script for outputting all the velocity diffusion traces"
gnuplot -c "Plot_Diffusion.gnu" beta0p2_eps0p4 v raw raw
gnuplot -c "Plot_Diffusion.gnu" beta0p5_eps0p4 v raw raw
for i in {1..5}
do
gnuplot -c "Plot_Diffusion.gnu" beta0p2_eps0p4 v j$(($i+6)) j$i
gnuplot -c "Plot_Diffusion.gnu" beta0p5_eps0p4 v j$(($i+6)) j$i
done

