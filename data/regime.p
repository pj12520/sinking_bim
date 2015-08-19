#Gnuplot script to plot regime diagram for tailing and film thinning regimes

   set terminal postscript enhanced color

   set logscale x
   set logscale y

   set key below
   
   set xlabel 'Modified Density Ratio'
   set ylabel 'Bond Number'

   set xrange[0.1:10000]
   set yrange[0.0001:10000]

   set title 'Viscosity Ratio = 0.001'

   set output 'const_viscos/viscos0.001.eps'

   plot 'const_viscos/viscos0.001.dat' u 1:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_viscos/viscos0.001.dat' u 1:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'


   set title 'Viscosity Ratio = 0.01'

   set output 'const_viscos/viscos0.01.eps'

   plot 'const_viscos/viscos0.01.dat' u 1:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_viscos/viscos0.01.dat' u 1:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'


   set title 'Viscosity Ratio = 0.1'

   set output 'const_viscos/viscos0.1.eps'

   plot 'const_viscos/viscos0.1.dat' u 1:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_viscos/viscos0.1.dat' u 1:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'

   
   set title 'Viscosity Ratio = 1'

   set output 'const_viscos/viscos1.eps'

   plot 'const_viscos/viscos1.dat' u 1:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_viscos/viscos1.dat' u 1:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'


   set title 'Viscosity Ratio = 10'

   set output 'const_viscos/viscos10.eps'

   plot 'const_viscos/viscos10.dat' u 1:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_viscos/viscos10.dat' u 1:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'


   set title 'Viscosity Ratio = 100'

   set output 'const_viscos/viscos100.eps'

   plot 'const_viscos/viscos100.dat' u 1:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_viscos/viscos100.dat' u 1:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'


   set title 'Viscosity Ratio = 1000'

   set output 'const_viscos/viscos1000.eps'

   plot 'const_viscos/viscos1000.dat' u 1:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_viscos/viscos1000.dat' u 1:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'




   set xlabel 'Modified Density Ratio'
   set ylabel 'Viscosity Ratio'

   set xrange[0.1:10000]
   set yrange[0.0001:10000]

   set title 'Bond Number = 0.001'

   set output 'const_bond/bond0.001.eps'

   plot 'const_bond/bond0.001.dat' u 1:(stringcolumn(4) eq 'film'? $3:1/0) title 'Film Drainage', 'const_bond/bond0.001.dat' u 1:(stringcolumn(4) eq 'tail'? $3:1/0) title 'Tailing'


   set title 'Bond Number = 0.01'

   set output 'const_bond/bond0.01.eps'

   plot 'const_bond/bond0.01.dat' u 1:(stringcolumn(4) eq 'film'? $3:1/0) title 'Film Drainage', 'const_bond/bond0.01.dat' u 1:(stringcolumn(4) eq 'tail'? $3:1/0) title 'Tailing'
   

   set title 'Bond Number = 0.1'

   set output 'const_bond/bond0.1.eps'

   plot 'const_bond/bond0.1.dat' u 1:(stringcolumn(4) eq 'film'? $3:1/0) title 'Film Drainage', 'const_bond/bond0.1.dat' u 1:(stringcolumn(4) eq 'tail'? $3:1/0) title 'Tailing'

   
   set title 'Bond Number = 1'

   set output 'const_bond/bond1.eps'

   plot 'const_bond/bond1.dat' u 1:(stringcolumn(4) eq 'film'? $3:1/0) title 'Film Drainage', 'const_bond/bond1.dat' u 1:(stringcolumn(4) eq 'tail'? $3:1/0) title 'Tailing'


   set title 'Bond Number = 10'

   set output 'const_bond/bond10.eps'

   plot 'const_bond/bond10.dat' u 1:(stringcolumn(4) eq 'film'? $3:1/0) title 'Film Drainage', 'const_bond/bond10.dat' u 1:(stringcolumn(4) eq 'tail'? $3:1/0) title 'Tailing'


   set title 'Bond Number = 100'

   set output 'const_bond/bond100.eps'

   plot 'const_bond/bond100.dat' u 1:(stringcolumn(4) eq 'film'? $3:1/0) title 'Film Drainage', 'const_bond/bond100.dat' u 1:(stringcolumn(4) eq 'tail'? $3:1/0) title 'Tailing'


   set title 'Bond Number = 1000'

   set output 'const_bond/bond1000.eps'

   plot 'const_bond/bond1000.dat' u 1:(stringcolumn(4) eq 'film'? $3:1/0) title 'Film Drainage', 'const_bond/bond1000.dat' u 1:(stringcolumn(4) eq 'tail'? $3:1/0) title 'Tailing'


   



   set xlabel 'Viscosity Ratio'
   set ylabel 'Bond Number'

   set xrange[0.0001:10000]
   set yrange[0.0001:10000]


   set title 'MDR = 1'

   set output 'const_mdr/mdr1.eps'

   plot 'const_mdr/mdr1.dat' u 3:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_mdr/mdr1.dat' u 3:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'


   set title 'MDR = 10'

   set output 'const_mdr/mdr10.eps'

   plot 'const_mdr/mdr10.dat' u 3:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_mdr/mdr10.dat' u 3:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'


   set title 'MDR = 100'

   set output 'const_mdr/mdr100.eps'

   plot 'const_mdr/mdr100.dat' u 3:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_mdr/mdr100.dat' u 3:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'


   set title 'MDR = 1000'

   set output 'const_mdr/mdr1000.eps'

   plot 'const_mdr/mdr1000.dat' u 3:(stringcolumn(4) eq 'film'? $2:1/0) title 'Film Drainage', 'const_mdr/mdr1000.dat' u 3:(stringcolumn(4) eq 'tail'? $2:1/0) title 'Tailing'




   set zlabel 'MDR'

   set zrange[0.1:10000]
   set logscale z

   set title '3D regime diagram'

   set view 16,286

   set output 'regime.eps'

   splot 'regime.dat' u 3:2:(stringcolumn(4) eq 'film'? $1:1/0) title 'Film Drainage', 'regime.dat' u 3:2:(stringcolumn(4) eq 'tail'? $1:1/0) title 'Tailing'