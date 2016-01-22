#!/usr/bin/gnuplot

r0=0.7
sc=4

v2=2*pi*r0/sc
v3=4*pi*r0**2/sc

plot\
 "test2d.dat" using 1:(100*($2-v2)/v2) with linespoints pt 6\
 title "2D relative error, %",\
 "test3d.dat" using 1:(100*($2-v3)/v3) with linespoints pt 6\
 title "3D relative error, %"

pause -1