all: he3lib.pdf yosida.pdf spin_diff.pdf

he3lib.pdf: ptable.tex
ptable.tex: ptable.m
	octave -q $<

%.pdf: %.ps
	ps2pdf $<
%.ps: %.dvi
	dvips $<
%.dvi: %.tex
	latex $<

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.dvi *.eps
	rm -f ptable.tex
