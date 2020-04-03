#!/bin/sh -x

# it's a simplies script to build a command line program. For those who
# do not want to adjust the huge and strange Makefile.

./make_inc
gcc -c he3.c -o he3.o
gfortran -fno-range-check -O -I. -std=legacy functions/*.f he3.o -o he3
gfortran -fno-range-check -O -I. -std=legacy --shared -fPIC functions/*.f -o he3.dll

## matlab -- not tested.
## -rpath should contain folder where .so/.dll file is located (if is is not
##  a standard system path for libraries)

#matlab64 -r 'mex -output he3lib.mex -L. LDFLAGS="-Wl,-rpath=." -lhe3 he3lib_mex.c'

