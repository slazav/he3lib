all: phase1.eps phase2.eps\
     gap1.eps\
     nu_b.eps text_pars.eps text_pars2.eps text_gr.eps\
     yosida_int.eps\
     transp.eps sdiff.eps\
     qball.eps

%.eps: %.m
	octave-cli -q $<
	touch st

clean:
	rm -f st *.tex *.eps