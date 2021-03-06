# fortran compiler parameters
FFLAGS= -Werror -Wconversion -I.\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O\
  -std=legacy

# important flags:
# -std=legacy -- to allow blas.f compilation
# -fno-range-check -- to allow NaN values

# fortran compiler
#FC=gfortran

all: headers\
     library\
     cmdline
#    doc
#    octave
#    octave-mex
#    matlab

.PHONY: octave matlab octave-mex matlab64 headers library cmdline doc


###################################

# he3 constants and functions (see src/)
LIBOBJS=he3_const he3_phase he3_fermi he3_normal\
        he3_math he3_gap he3_dipole he3_grad he3_text\
        he3_transp_n he3_transp_b he3_b2 he3_other he3_polar\
        he3_rota he3_bspec he3_a he4 he3_wire he34
OBJS= $(patsubst %,functions/%.o,$(LIBOBJS))
SRCS= $(patsubst %,functions/%.f,$(LIBOBJS))

$(OBJS): he3.fh

###################################

# header files for C/F77/F90 are created from source files
he3.f90h he3.fh he3.h he3tab.h: $(SRCS) make_inc
	./make_inc

LIBNAME=libhe3
library: $(LIBNAME).a $(LIBNAME).so

$(LIBNAME).a: $(OBJS)
	ar rs $@ $+

$(LIBNAME).so: $(OBJS)
	$(FC) --shared -fPIC -o $@ $+

# cmdline program
cmdline: he3
he3.o: he3.c he3tab.h
	$(CC) -c he3.c -o he3.o
he3: he3.o $(LIBNAME).a
	$(FC) $+ -o $@

###################################
## install rules

bindir  ?= /usr/bin
libdir  ?= /usr/lib64
datadir ?= /usr/share
includedir  ?= /usr/include
octdir = ${datadir}/octave/site/he3lib

install_headers: he3.f90h he3.fh he3.h
	mkdir -p ${includedir}
	install -m0644 $+  ${includedir}

install_library: $(LIBNAME).so
	mkdir -p ${libdir}
	install -m0644 $+ ${libdir}

install_cmdline: he3
	mkdir -p ${bindir}
	install -m0755 $+ ${bindir}

install_octave: octave install_library
	mkdir -p ${octdir}/packinfo
	install -m0644 *.oct ${octdir}
	install -m0644 m/*.m ${octdir}
	install -m0644 m/DESCRIPTION ${octdir}/packinfo

###################################

octave:     library he3lib.oct
octave-mex: library he3lib.mex
matlab:     library he3lib.mexglx
matlab64:   library he3lib.mexa64

he3lib.oct: he3lib_oct.cc he3.h he3tab.h
	mkoctfile $< -lhe3 -L. -W -std=c++11 -s -v -o $@

he3lib.mex: he3lib_mex.c he3.h he3tab.h
	mkoctfile $< -mex -lhe3 -L. -W -s -v -o $@

he3lib.mexglx: he3lib_mex.c he3.h he3tab.h
	matlab -nojvm -nosplash -r "mex -output $@ libhe3.so $<"

he3lib.mexa64: he3lib_mex.c he3.h he3tab.h
	matlab64 -nojvm -nosplash -r "mex -output $@ libhe3.so $<"

doc: octave
	make -C doc_misc
	make -C doc_tex

###################################
clean:
	rm -f *.a *.so *.o *.oct m/*.m he3.f90h he3.fh he3.h he3 he3tab.h functions/*.o
