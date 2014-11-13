
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

all: build_headers\
     build_library\
     build_cmdline

#     build_doc
#     build_octave

FC=gfortran

###################################

# header files for C/F77/F90 are created from he3.def
build_headers: he3.f90h he3.fh he3.h he3tab.h
he3.f90h he3.fh he3.h he3tab.h: he3.def make_inc
	./make_inc

LIBNAME=libhe3
build_library: $(LIBNAME).a $(LIBNAME).so

# he3 constants and functions (see src/)
LIBOBJS=he3_const he3_phase he3_fermi he3_normal\
        he3_math he3_gap he3_dipole he3_grad he3_text\
        he3_transp_n he3_transp_b he3_other\
        he3_rota
OBJS= $(patsubst %,%.o,$(LIBOBJS))

$(LIBNAME).a: $(OBJS)
	ar rs $@ $+

$(LIBNAME).so: $(OBJS)
	$(FC) --shared -fPIC -o $@ $+

# cmdline program
build_cmdline: he3
he3.o: he3.c he3tab.h
	$(CC) -c he3.c -o he3.o
he3: he3.o libhe3.a
	$(FC) $+ -o $@

###################################
build_octave: build_library
	make -C matlab octave

build_doc: build_octave
	make -C doc
	make -C doc_tex

###################################
clean:
	rm -f *.a *.so *.o libs/*.o he3.f90h he3.fh he3.h he3
#	make -C matlab clean
#	make -C doc clean
