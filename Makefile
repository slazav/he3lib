
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

all: headers\
     library\
     cmdline
#    doc
#    octave
#    matlab

.PHONY: octave matlab matlab64 headers library cmdline doc

#FC=gfortran

###################################

# header files for C/F77/F90 are created from he3.def
headers: he3.f90h he3.fh he3.h he3tab.h
he3.f90h he3.fh he3.h he3tab.h octave/he3.tab: he3.def make_inc
	./make_inc

LIBNAME=libhe3
library: $(LIBNAME).a $(LIBNAME).so

# he3 constants and functions (see src/)
LIBOBJS=he3_const he3_phase he3_fermi he3_normal\
        he3_math he3_gap he3_dipole he3_grad he3_text\
        he3_transp_n he3_transp_b he3_b2 he3_other he3_polar\
        he3_rota he3_bspec
OBJS= $(patsubst %,%.o,$(LIBOBJS))

$(LIBNAME).a: $(OBJS)
	ar rs $@ $+

$(LIBNAME).so: $(OBJS)
	$(FC) --shared -fPIC -o $@ $+

# cmdline program
cmdline: he3
he3.o: he3.c he3tab.h
	$(CC) -c he3.c -o he3.o
he3: he3.o libhe3.a
	$(FC) $+ -o $@

###################################
## install rules

bindir  ?= /usr/bin
libdir  ?= /usr/lib
datadir ?= /usr/share
includedir  ?= /usr/include
octdir = ${datadir}/octave/packages/he3lib

install: install_headers install_library install_cmdline

install_headers: he3.f90h he3.fh he3.h
	mkdir -p ${includedir}
	install -m0644 $+  ${includedir}

install_library: $(LIBNAME).a $(LIBNAME).so
	mkdir -p ${libdir}
	install -m0644 $+ ${libdir}

install_cmdline: he3
	mkdir -p ${bindir}
	install -m0755 $+ ${bindir}

install_octave: octave
	mkdir -p ${octdir}/packinfo
	install -m0644 octave/*.oct ${octdir}
	install -m0644 octave/DESCRIPTION ${octdir}/packinfo

###################################
octave: library octave/he3.tab
	make -C octave

matlab: library
	make -C matlab matlab

matlab64: library
	make -C matlab matlab64

doc: octave
	make -C doc
	make -C doc_tex

###################################
clean:
	rm -f *.a *.so *.o libs/*.o he3.f90h he3.fh he3.h he3
#	make -C matlab clean
#	make -C doc clean
