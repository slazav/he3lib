#!/usr/bin/gnuplot

set nokey

y0 = 1.42613

f(x) = A*x**2 + B*x + 1
fit f(x) "wcp_coeff.txt" using ($1-y0):2 via A,B
plot "wcp_coeff.txt" using ($1-y0):2 pt 7 title "a", f(x)
print A,B,1
pause -1

f(x) = A*x**2 + B*x
fit f(x) "wcp_coeff.txt" using ($1-y0):3 via A,B
plot "wcp_coeff.txt" using ($1-y0):3 pt 7 title "b", f(x)
print A,B,0
pause -1

f(x) = A*x**2 + B*x + C
fit f(x) "wcp_coeff.txt" using ($1-y0):4 via A,B,C
plot "wcp_coeff.txt" using ($1-y0):4 pt 7 title "c", f(x)
print A,B,C
pause -1

fit f(x) "wcp_coeff.txt" using ($1-y0):5 via A,B,C
plot "wcp_coeff.txt" using ($1-y0):5 pt 7 title "d", f(x)
print A,B,C
pause -1
