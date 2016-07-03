# Gnuplot script for plotting TKE results
# run it in the directory of tke.dat with
# the option -persist

#Uncomment this for saving the plot as svg
#set terminal svg size 755.9,453.5 fname 'Verdana' fsize 14
#set output 'tke.svg'


set multiplot
set key autotitle columnhead
set xrange [0:1]
set xlabel "y/H"
plot for [i=2:7] 'tke.dat' u 1:i w l lw 3, 'tke.dat' u 1:($2+$3-$4+$5+$6+$7) w l lw 3 title "RESIDUAL"
set origin 0.25,0.25
set size 0.5,0.7
set xrange [0:0.1]
set xlabel "y/H"
plot for [i=2:7] 'tke.dat' u 1:i w l lw 3, 'tke.dat' u 1:($2+$3-$4+$5+$6+$7) w l lw 3 title "RESIDUAL"
