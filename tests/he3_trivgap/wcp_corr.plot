#!/usr/bin/gnuplot

f(x) = a + b*exp(-x*c-x**2*d)
fit f(x) "wcp_corr.txt" using (1-$1):($7**2) via a,b,c,d

plot\
 "wcp_corr.txt" using (1-$1):($3**2) title "1.43",\
 "wcp_corr.txt" using (1-$1):($4**2) title "1.6",\
 "wcp_corr.txt" using (1-$1):($5**2) title "1.8",\
 "wcp_corr.txt" using (1-$1):($6**2) title "2.0",\
 "wcp_corr.txt" using (1-$1):($7**2) title "2.2",\
 f(x)

print a,b,c,d
pause -1

fit f(x) "wcp_corr.txt" using (1-$1):($6**2) via a,b,c,d
replot
print a,b,c,d
pause -1

fit f(x) "wcp_corr.txt" using (1-$1):($5**2) via a,b,c,d
replot
print a,b,c,d
pause -1

fit f(x) "wcp_corr.txt" using (1-$1):($4**2) via a,b,c,d
replot
print a,b,c,d
pause -1
