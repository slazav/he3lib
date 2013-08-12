all:  pics/st he3lib.pdf yosida.pdf spin_diff.pdf

pics/st:
	makefile -C pics

%.pdf: %.ps
	ps2pdf $<
%.ps: %.dvi
	dvips $<
%.dvi: %.tex
	latex $<

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.dvi *.eps
