set terminal postscript

unset key

set output config_plot

set xrange[-15:15]
set yrange[-15:15]
set size square

set title time

plot interf_dat u 2:3 w lines, sphere_dat u 2:3 w lines, interf_dat u (-$2):3 w lines, sphere_dat u (-$2):3 w lines