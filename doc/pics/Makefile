all: phase.eps gap.eps yosida_int.eps

%.eps: %.m
	octave -q $<
	touch st

clean:
	rm -f st *.tex *.eps