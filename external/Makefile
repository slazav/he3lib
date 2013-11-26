FFLAGS= -Wconversion -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -g\
  -std=legacy

all: libint.a libtn.a  libbtn.a libpde.a libpoly.a

FC=gfortran

# integration functions from http://www.netlib.org/slatec/index.html
INT_SRC=dqags dqagse dqelg dqk21 dqpsrt d1mach xermsg\
     j4save xerprn xersve xerhlt xercnt fdump\
     xgetua i1mach
INT_SRC:=$(patsubst %,int/%.o,$(INT_SRC))

# minimization
TN_SRC=blas tn
TN_SRC:=$(patsubst %,tn/%.o,$(TN_SRC))

BTN_SRC=btnlib chklib basics d1mach
BTN_SRC:=$(patsubst %,btn/%.o,$(BTN_SRC))

# pdecol
PDE_SRC=pde_dp
PDE_SRC:=$(patsubst %,pde/%.o,$(PDE_SRC))

POLY_SRC=E02AEE E02CBE M01AGE P01AAE X02AAE X04AAE
POLY_SRC:=$(patsubst %,poly/%.o,$(POLY_SRC))


libint.a: $(INT_SRC)
	ar rs $@ $+

libtn.a: $(TN_SRC)
	ar rs $@ $+

libbtn.a: $(BTN_SRC)
	ar rs $@ $+

libpde.a: $(PDE_SRC)
	ar rs $@ $+

libpoly.a: $(POLY_SRC)
	ar rs $@ $+


clean:
	rm -f *.a *.so */*.o
