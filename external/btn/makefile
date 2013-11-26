all: btn.a main1 main2 main3 main4

GOBJS =	btnlib.o chklib.o basics.o d1mach.o
FC = gfortran

.SUFFIXES: .f .o

.f.o:
	gfortran -c $<

main1:  $(GOBJS) main1.f
	$(FC) $+ -o $@
main2:  $(GOBJS) main2.f
	$(FC) $+ -o $@
main3:  $(GOBJS) main3.f
	$(FC) $+ -o $@
main4:  $(GOBJS) main4.f
	$(FC) $+ -o $@

btn.a: $(GOBJS)
	ar rs $@ $+

clean:
	rm -f -- *.o *.a main? fort.10 get_d1mach