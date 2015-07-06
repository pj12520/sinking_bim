#Gnuplot script to plot system configurations and sphere trajectory

set terminal postscript

set output traj_plot

unset key

plot traj_data u 2:4 w lines

set output config_plot

set xrange[0:15]
set yrange[-10:5]
set size square

set title time

plot interf_dat  u 2:3 w lines, sphere_dat u 2:3 w lines