
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
  -Wno-unused-parameter -fPIC -fno-range-check

LIBNAME=libhe3
all: $(LIBNAME).a $(LIBNAME).so

FC=gfortran

###################################

# he3 constants and functions (see src/)
LIBOBJS=he3_const he3_phase he3_fermi he3_gap\
        he3_transp_n he3_transp_b\
        he3_flegg he3_d_exp\
        he3_swvel he3_tau_lt

# additional fitting functions used in libhe3
ADDOBJS=E02AEE E02CBE M01AGE P01AAE X02AAE X04AAE
#        dgesv dgetrs dlaswp dtrsm lsame xerbla



# Legget equations
LEGG_EQ_OBJS=he3b_legg_rot1d

OBJS=\
  $(patsubst %,%.o,$(LIBOBJS))\
  $(patsubst %,libs/%.o,$(ADDOBJS))\
  $(patsubst %,legg_eq/%.o,$(LEGG_EQ_OBJS))

$(LIBNAME).a: $(OBJS)
	ar rs $@ $+

$(LIBNAME).so: $(OBJS)
	$(FC) --shared -fPIC -o $@ $+

clean:
	rm -f *.a *.so *.o libs/*.o legg_eq/*.o
