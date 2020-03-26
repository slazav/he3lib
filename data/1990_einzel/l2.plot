#!/usr/bin/gnuplot

f(x) = a*x**3 + b*x**2 + c*x + d

a = 5e-6
b = -4e-4
c = 9.5e-3
d = 0.68

#fit f(x) "l2.txt" via a,b,c,d

set nokey
plot "l2.txt" with points pt 7, f(x)


print a,b,c,d


pause -1