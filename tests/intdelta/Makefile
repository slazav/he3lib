# tests for functions from he3_math.f
# math_dint2d_delta
# math_dint3d_delta

# gfortran parameters
FC=gfortran
FFLAGS= -Werror -Wconversion\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O\
  -std=legacy -g

all: plot

test2d: test2d.f ../he3_math.f
test3d: test3d.f ../he3_math.f

%.dat: %
	./$< > $@

plot: plot.gp test2d.dat test3d.dat
	./plot.gp