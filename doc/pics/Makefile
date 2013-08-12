all: fermi.tex phase.tex phase.eps gap.eps yosida_int.eps


%.tex: %.m
	octave -q $<
	touch st

%.eps: %.m
	octave -q $<
	touch st

clean:
	rm -f st *.tex *.eps