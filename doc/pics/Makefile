all: phase1.eps phase2.eps\
     gap1.eps\
     transp_tau.eps transp_diff.eps\
     nu_b.eps text_pars.eps text_pars2.eps text_gr.eps\
     yosida_int.eps\
     sdiff_per.eps sdiff_par.eps

%.eps: %.m
	octave -q $<
	touch st

clean:
	rm -f st *.tex *.eps