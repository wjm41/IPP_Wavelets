set terminal postscript eps enhanced color size 25,15 font 'Helvetica,60'
File_name=ARG1."_Diffusion_".ARG2."_".ARG4
set output File_name.'.eps'
set xlabel "Time" 
set ylabel "RMS"
#set grid
set key top right
Title='RMS values of the Diffusion parallel and perpendicular to the magnetic field as a function of time.'
set title "Plot of " . Title
plot ARG1.'/electrons/Diffusion/'.ARG2.'_par_rms_'.ARG3.'.dat' with lines lw 12 title "Electrons Parallel", ARG1.'/electrons/Diffusion/'.ARG2.'_perp_rms_'.ARG3.'.dat' with lines lw 12 title "Electrons Perpendicular", ARG1.'/ions/Diffusion/'.ARG2.'_par_rms_'.ARG4.'.dat' with lines lw 12 title "Ions Parallel", ARG1.'/ions/Diffusion/'.ARG2.'_perp_rms_'.ARG4.'.dat' with lines lw 12 title "Ions Perpendicular"
#pause -1
