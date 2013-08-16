all: phase1.eps phase2.eps gap1.eps gap_sus.eps\
     gap_lf.eps yosida_int.eps

%.eps: %.m
	octave -q $<
	touch st

clean:
	rm -f st *.tex *.eps