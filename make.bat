#!/bin/sh -x
gcc -c he3.c -o he3.o
gfortran -fno-range-check -O -std=legacy *.f external/poly/*.f he3.o -o he3
