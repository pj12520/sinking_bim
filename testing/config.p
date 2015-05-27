# Gnuplot script to plot profiles showing the initial configuration of the system as produced from testing

# Plot to show the initial configuration

set xrange[0:15]
set yrange[-1:5]

unset key

set terminal postscript enhanced color
set output 'config.eps'

plot 'sphere_config.dat' u 2:3 w lines, 'interf_config.dat' u 2:3 w lines

#Plot to show r(theta) and z(theta)

set xrange[0:3.142]
set yrange[0:4]

set xlabel 'Theta'

set key below

set terminal postscript enhanced 
set output 'sphere.eps'

plot 'sphere_config.dat' u 1:2 w lines title 'r', 'sphere_config.dat' u 1:3 w lines title 'z'

#Plot to show r(arc) and z(arc)

set xrange[0:15]
set yrange[-1:15]

set xlabel 'Arc'

set key below

set terminal postscript enhanced 
set output 'interf.eps'

plot 'interf_config.dat' u 1:2 w lines title 'r', 'interf_config.dat' u 1:3 w lines title 'z'