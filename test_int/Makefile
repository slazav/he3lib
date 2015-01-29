# gfortran parameters
FC=gfortran
FFLAGS= -Werror -Wconversion\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O\
  -std=legacy -g

all: plot

test_dint: test_dint.f ../he3_math.f
test_dint_gk: test_dint_gk.f ../he3_math.f
test_dint_gka: test_dint_gka.f ../he3_math.f

%.dat: %
	./$< > $@

plot: plot.gp test_dint.dat test_dint_gk.dat test_dint_gka.dat
	./plot.gp