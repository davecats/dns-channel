# Gnuplot script for plotting MKE results
# run it in the directory of mke.dat with the
# option -persist

set multiplot
set key autotitle columnhead
set xrange [0:1]
set xlabel "y/H"
plot for [i=2:8] 'mke.dat' u 1:i w l lw 3,  'mke.dat' u 1:($2+$5+$6-$7+$3+$4-$8) w l lw 3 title "RESIDUAL"
set origin 0.25,0.35
set size 0.5,0.6
set xrange [0:0.1]
set xlabel "y/H"
plot for [i=2:8] 'mke.dat' u 1:i w l lw 3,  'mke.dat' u 1:($2+$5+$6-$7+$3+$4-$8) w l lw 3 title "RESIDUAL"
