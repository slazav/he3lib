#!/usr/bin/gnuplot

plot "test_dint.dat" using 1:(log($2)) with linespoints pt 6 title "G2",\
     "test_dint_gk.dat" using 1:(log($2)) with linespoints pt 6 title "G7+K13",\
     "test_dint_gk.dat" using 1:(log($3)) with linespoints pt 6 title "G7+K13 est err",\
     "test_dint_gka.dat" using 1:(log($2)) with linespoints pt 6 title "G7+K13 adaptive"
pause -1