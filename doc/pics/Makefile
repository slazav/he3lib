all: phase.eps gap.eps yosida_int.eps

phase1.eps phase2.eps: phase.m
	octave -q $<
	touch st

%.eps: %.m
	octave -q $<
	touch st

clean:
	rm -f st *.tex *.eps