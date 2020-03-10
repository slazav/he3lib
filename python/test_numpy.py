#!/usr/bin/python3

import he3lib
import numpy
import matplotlib.pyplot as plt

# phase diagram, tc, tab
p = numpy.arange(0,45,0.1) # [bar]
tc   = he3lib.he3_tc(p)
tab  = he3lib.he3_tab(p)

plt.plot(tc,p, 'b-')
plt.plot(tab,p, 'b-')

# phase diagram, P_melt
t  = numpy.arange(0,3,0.01) # [mK]
pmelt = he3lib.he3_pmelt(t/1000)
plt.plot(t,pmelt, 'b-')

# tc in finite magnetic field
H = 1000; # [G]
tc_h = he3lib.he3_b2tab(p,H);
plt.plot(tc_h,p, 'r-')

plt.xlim(0,3)
plt.ylim(0,40)
plt.savefig('test.png', format='png', dpi='figure')
