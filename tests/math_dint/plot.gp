#!/usr/bin/gnuplot

set title "int 1/(1+(kx)^2) from -1 to 1  (= 2/k atan(k))"
set xlabel "k"
set ylabel "log(rel.err.)"

plot "test_dint.dat" using 1:(log(abs($2))) with linespoints pt 6 title "G2",\
     "test_dint_gk.dat" using 1:(log(abs($2))) with linespoints pt 6 title "G7+K13",\
     "test_dint_gk.dat" using 1:(log(abs($3))) with linespoints pt 6 title "G7+K13 est err",\
     "test_dint_gka.dat" using 1:(log(abs($2))) with linespoints pt 6 title "G7+K13 adaptive"
pause -1