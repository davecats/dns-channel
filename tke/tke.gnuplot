# Gnuplot script for plotting TKE results
# run it in the directory of tke.dat with
# the option -persist

set multiplot
set key autotitle columnhead
set xrange [0:1]
set xlabel "y/H"
plot for [i=2:6] 'tke.dat' u 1:i w l lw 3, 'tke.dat' u 1:($2-$3+$4+$5+$6) w l lw 3 title "RESIDUAL"
set origin 0.25,0.25
set size 0.5,0.7
set xrange [0:0.1]
set xlabel "y/H"
plot for [i=2:6] 'tke.dat' u 1:i w l lw 3, 'tke.dat' u 1:($2-$3+$4+$5+$6) w l lw 3 title "RESIDUAL"
