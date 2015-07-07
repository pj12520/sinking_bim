#Gnuplot script to plot regime diagram for tailing and film thinning regimes

   set terminal postscript enhanced color

   set xlabel 'Modified Density Ratio'
   set ylabel 'Bond Number'

   set key below
   
   set title 'Tailing Vs. Film Drainage regime diagram

   set logscale x
   set logscale y

   set xrange[0.1:10000]
   set yrange[0.0001:10000]

   set output 'regime.eps'

   plot 'regime.dat' u 1:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'regime.dat' u 1:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'
