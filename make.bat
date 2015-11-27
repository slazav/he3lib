#!/bin/sh -x

# it's a simplies script to build a command line program. For those who
# do not want to adjust the huge and strange Makefile.

gcc -c he3.c -o he3.o
gfortran -fno-range-check -O -std=legacy *.f he3.o -o he3
