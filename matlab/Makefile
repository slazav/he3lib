LIBNAME=libhe3
MEX=/opt/matlab-R2011b/bin/mex

FUNC_I1 = he3_pmelt he3_tc he3_tab
FUNC_I0 = he3_pa he3_pabn he3_tabn he3_amass

MEXS_I0=$(patsubst %,%.mexglx,$(FUNC_I0))
MEXS_I1=$(patsubst %,%.mexglx,$(FUNC_I1))

all: $(LIBNAME).so

$(LIBNAME).so:../$(LIBNAME).so
	ln -s $< $@

#$(MEXS_I0): mexfunc.c $(LIBNAME).so
#	name=$(subst .mexglx,,$@);\
#	$(MEX) -glnx86 -DFUNC=$${name}_ -DNARGIN=0\
#	   mexfunc.c ./$(LIBNAME).so -output $$name
#
#$(MEXS_I1): mexfunc.c ../$(LIBNAME).so
#	name=$(subst .mexglx,,$@);\
#	$(MEX) -glnx86 -DFUNC=$${name}_ -DNARGIN=1\
#	   mexfunc.c ./$(LIBNAME).so -output $$name
#
######

#%.mexa64: %.o mexfunc.c ../$(LIBNAME).so
#	$(MEX) -v -glnxa64 -DFUNC=$*_\
#	   mexfunc.c $(LIBNAME).so -output $*

clean:
	rm -f *.a *.so *.o *.mexglx
