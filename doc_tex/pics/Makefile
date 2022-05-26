all: phase1.png gap1.png\
     nu_b.png text_pars.png text_pars2.png text_gr.png\
     yosida_int.png\
     transp.png sdiff.png\
     qball.png

%.png: %.m
	octave-cli -q $<
	touch st

clean:
	rm -f st *.tex *.png