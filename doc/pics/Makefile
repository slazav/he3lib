all: phase1.eps phase2.eps\
     gap1.eps\
     transp_tau.eps transp_diff.eps\
     other_lf.eps\
     yosida_int.eps

%.eps: %.m
	octave -q $<
	touch st

clean:
	rm -f st *.tex *.eps