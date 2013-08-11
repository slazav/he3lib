all: fermi.tex phase.tex phase.eps


%.tex: %.m
	octave -q $<
	touch st

%.eps: %.m
	octave -q $<
	touch st

clean:
	rm -f st *.tex *.eps