
#FFLAGS= -Werror -Wconversion\
#  -Wintrinsic-shadow -Wline-truncation\
#  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
#  -Wline-truncation -Wintrinsics-std -Wsurprising -Wno-tabs -Wunderflow\
#  -Wno-unused-parameter -Wno-align-commons -fno-range-check

# gfortran parameters
FFLAGS= -Werror -Wconversion\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O\
  -std=legacy

# important flags:
# -std=legacy -- to allow blas.f compilation
# -fno-range-check -- to allow NaN values

LIBNAME=libhe3
all: external $(LIBNAME).a $(LIBNAME).so he3.f90h

FC=gfortran

###################################

# he3 constants and functions (see src/)
LIBOBJS=he3_const he3_phase he3_fermi he3_normal\
        he3_math he3_gap he3_dipole\
        he3_transp_n he3_transp_b\
        he3_text he3_other\
        he3_rota

# additional fitting functions used in libhe3
ADDOBJS=E02AEE E02CBE M01AGE P01AAE X02AAE X04AAE
#        dgesv dgetrs dlaswp dtrsm lsame xerbla

# h-file for f90 is created from he3.fh
he3.f90h: he3.fh
	echo "! This file is created automatically from $<" > $@
	sed -e 's/!F90_ONLY!//' $< >> $@

# Legget equations
LEGG_EQ_OBJS=he3b_legg_rot1d

OBJS=\
  $(patsubst %,%.o,$(LIBOBJS))\
  $(patsubst %,legg_eq/%.o,$(LEGG_EQ_OBJS))\
  external/poly/*.o\
  external/int/*.o\
  external/tn/*.o

$(LIBNAME).a: $(OBJS)
	ar rs $@ $+

$(LIBNAME).so: $(OBJS)
	$(FC) --shared -fPIC -o $@ $+

clean:
	rm -f *.a *.so *.o libs/*.o legg_eq/*.o he3.f90h
	make -C matlab clean
	make -C doc clean

test1: test.f libhe3.a
	$(FC) -g -fno-range-check $+ -o $@ ../external/libint.a 