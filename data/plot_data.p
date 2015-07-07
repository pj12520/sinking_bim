#Gnuplot script to plot system configurations and sphere trajectory

set terminal postscript

set output traj_plot

unset key

set xlabel "Dimensionless Time"
set ylabel "Dimensionless Position"

plot traj_data u 2:3 w lines

