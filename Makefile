# fortran compiler parameters
FFLAGS= -Werror -Wconversion -I.\
  -Wline-truncation\
  -Waliasing  -Wampersand -Warray-bounds -Wcharacter-truncation\
  -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow\
  -Wno-unused-parameter -fPIC -fno-range-check -O\
  -std=legacy -frecursive -g

# important flags:
# -std=legacy -- to allow blas.f compilation
# -fno-range-check -- to allow NaN values
# -frecursive -- to allocate all arrays on the stack

# fortran compiler
#FC=gfortran

all: headers\
     library\
     cmdline
#    doc
#    octave
#    octave-pkg
#    octave-mex
#    matlab
#    matlab64

.PHONY: headers library cmdline doc\
        octave matlab octave-mex matlab64 octave-pkg\
        install_headers install_library install_cmdline\
        install_octave_local install_octave_global


###################################

# he3 constants and functions (see functions/)
LIBOBJS=he3_const he3_phase he3_fermi he3_normal\
        he3_math he3_gap he3_dipole he3_grad he3_text\
        he3_transp_n he3_transp_b he3_b2 he3_other he3_polar\
        he3_rota he3_bspec he3_a he4 he3_wire_orig he34
OBJS= $(patsubst %,functions/%.o,$(LIBOBJS))
SRCS= $(patsubst %,functions/%.f,$(LIBOBJS))

$(OBJS): he3.fh

###################################

# header files for C/F77/F90 are created from source files
# Header files for C/F77/F90 are created from he3.def by make_inc script.
# This script also creates matlab/octave wrappers in m/
headers: he3.f90h he3.fh he3.h he3tab.h
he3.f90h he3.fh he3.h he3tab.h: $(SRCS) make_inc
	./make_inc

# library
LIBNAME=libhe3
library: $(LIBNAME).a $(LIBNAME).so
$(LIBNAME).a: $(OBJS)
	ar rs $@ $+
$(LIBNAME).so: $(OBJS)
	$(FC) --shared -fPIC -o $@ $+

# cmdline program
cmdline: he3
he3.o: he3.c he3tab.h
he3: he3.o $(LIBNAME).a
	$(FC) $+ -o $@

# octave
octave: headers he3lib.oct
he3lib.oct: he3lib_oct.cc libhe3.a he3.h he3tab.h
	mkoctfile $< libhe3.a -W -std=c++11 -s -v -o $@

# octave-mex
octave-mex: headers he3lib.mex
he3lib.mex: he3lib_mex.c libhe3.a he3.h he3tab.h
	mkoctfile $< -mex libhe3.a -W -s -v -o $@

# matlab
matlab: headers he3lib.mexglx
he3lib.mexglx: he3lib_mex.c libhe3.a he3.h he3tab.h
	matlab -nojvm -nosplash -r "mex -output $@ libhe3.a $<"

# matlab64
matlab64: headers he3lib.mexa64
he3lib.mexa64: he3lib_mex.c libhe3.a he3.h he3tab.h
	matlab64 -nojvm -nosplash -r "mex -output $@ libhe3.a $<"

# octave package (see https://octave.org/doc/v6.4.0/Creating-Packages.html)
octave-pkg:octave-he3lib.tgz
octave-he3lib.tgz: octave headers
	tar -cvzh -f $@  m/*.m he3lib.oct DESCRIPTION LICENSE\
	 --xform 's|LICENSE|COPYING|'\
	 --xform 's|^m/|inst/|'\
	 --xform 's|^\(he3lib.oct\)|inst/\1|'\
	 --xform 's|^|pkg/|'

doc: octave
	make -C doc_misc
	make -C doc_tex

###################################
## install rules

bindir  ?= /usr/bin
libdir  ?= /usr/lib64
datadir ?= /usr/share
includedir  ?= /usr/include
octdir = ${libdir}/octave/site/

install_headers: he3.f90h he3.fh he3.h
	mkdir -p ${includedir}
	install -m0644 $+  ${includedir}

install_library: $(LIBNAME).so
	mkdir -p ${libdir}
	install -m0644 $+ ${libdir}

install_cmdline: he3
	mkdir -p ${bindir}
	install -m0755 $+ ${bindir}

install_octave_local: octave-he3lib.tgz
	octave-cli --eval "pkg install -local $<"

install_octave_global: octave-he3lib.tgz
	octave-cli --eval "pkg install -global $<"

###################################

clean:
	rm -f *.a *.so *.o *.dll *.oct m/*.m functions/*.o
	rm -f he3.f90h he3.fh he3.h he3 he3tab.h *.tgz
